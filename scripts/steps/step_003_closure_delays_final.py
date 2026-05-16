#!/usr/bin/env python3
"""
STEP 003: Closure Delay Measurement — Final Production Version

This script computes closure delays from secondary spectra and tests for
synchronization holonomy as predicted by the Temporal Equivalence Principle.

Theory:
-------
For a triplet of scintillation arclets forming a closed loop with paths i, j, k,
the closure residual is:

   C_ijk = tau_hat_ij + tau_hat_jk + tau_hat_ki

where tau_hat_ab is the independently measured delay of the a,b cross-term
extracted via sub-pixel peak fitting in the secondary spectrum. In a purely
additive screen model, tau_hat_ab = tau_b - tau_a within measurement noise,
implying C_ijk = 0. The closure observable therefore tests the failure of this
additive representation.

Under standard physics (additive delays): C_ijk ≡ 0 within noise
Under TEP (non-additive): C_ijk ≠ 0 systematically

Implementation Notes:
---------------------
1. AXIS CONVENTION: Secondary spectrum S has shape (n_fD, n_tau)
   S[fD_idx, tau_idx] where fD_idx corresponds to fD_mHz[fD_idx]

2. ALIASING CHECK: Cross-terms at (tauⱼ-tauᵢ, fDⱼ-fDᵢ) must be within bounds

3. phase_delta_ns: per-triplet narrowband equivalent delay (ns) from phase_closure_rad
   and epoch centre frequency in MHz; see PHASE_DELTA_NS_FORMULA / PHASE_DELTA_NS_DEFINITION
   on closure summaries. Primary holonomy magnitudes remain H_* from geometric closure.
"""

import argparse
import functools
import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.binary_reflex_kinematics import reflex_binary_transverse_velocity_kms
from scripts.utils.config import (
    C_LIGHT_KM_S,
    J0437_A1_LC,
    J0437_DEC_RAD,
    J0437_ECC,
    J0437_INC_DEG,
    J0437_KOM_DEG,
    J0437_OM_DEG,
    J0437_PB_DAYS,
    J0437_PSI_DEG,
    J0437_RA_RAD,
    J0437_T0_MJD,
)
from scripts.utils.json_numpy import NpEncoder
from scripts.utils.logger import TEPLogger, print_status, set_step_logger
from scripts.utils.parallel_workers import configure_blas_thread_env, worker_count

# Import astropy at module level to avoid overhead/instability in process pool
from astropy.time import Time
from astropy.coordinates import get_body_barycentric_posvel

# =============================================================================
# CONFIGURABLE PARAMETERS WITH PHYSICAL JUSTIFICATION
# =============================================================================

# Parabolic interpolation: minimum denominator to avoid numerical instability
# Value chosen based on machine epsilon for float32 (~1e-7) squared
PARABOLIC_DENOM_TOLERANCE = 1e-10

# Search radius for sub-pixel peak finding (pixels)
# Based on typical arclet FWHM of 6-10 pixels in secondary spectrum
SEARCH_RADIUS_PIXELS = 8

# Margin from spectrum edges to avoid aliasing artifacts (microseconds / mHz)
# Physically corresponds to ~1% of typical spectrum bounds
EDGE_MARGIN = 0.01  # 0.01 us and 0.01 mHz respectively

# Maximum arclet deviation tolerance (in units of grid spacing)
# 10x grid spacing allows sub-pixel peaks while rejecting false positives
MAX_DEVIATION_FACTOR = 10

# Maximum arclets to consider for closure combinations
# 12 arclets yields 220 closure triplets (C(12,3) = 220)
# This provides sufficient statistics while limiting computation
MAX_ARCLETS_FOR_TRIPLETS = 12

# Maximum triplets to return per epoch
# Limits output size while preserving highest-SNR measurements
MAX_TRIPLETS_PER_EPOCH = 20

# NOTE: A legacy hard cap at 0.05 μs (50 ns) on |geometric_delta_us| was removed.
# That gate truncated the tail of the Stokes-weighted closure distribution and could
# discard physically large closures; outlier control is left to cross-term validity + SNR.

# Machine-readable semantics for per-triplet phase_delta_ns (see ClosureResult).
# tau_ns = psi / (2*pi*f_Hz) with f_Hz = nu_ref_mhz*1e6  =>  tau_ns = psi*1e3/(2*pi*nu_ref_mhz).
# psi comes from processed secondary cross-terms; it is not guaranteed to equal 2*pi*f*tau_geom.
# Primary holonomy magnitudes are H_* fields from geometric_delta_us (microseconds -> ns).
PHASE_DELTA_NS_DEFINITION = (
    "phase_delta_ns = psi * 1e3 / (2*pi*nu_ref_mhz), psi=phase_closure_rad (rad), "
    "nu_ref_mhz=epoch band centre (MHz); narrowband identity tau_ns=psi/(2*pi*f_Hz). "
    "Illustrative only: psi need not match geometric closure delay. "
    "Primary scalar observables are H_magnitude_ns / H_trim_magnitude_ns from geometric closure."
)
PHASE_DELTA_NS_FORMULA = "tau_ns = psi * 1e3 / (2*pi*nu_ref_mhz)"

SEC_DIR = PROJECT_ROOT / "data" / "secondary" / "j0437"
SEC_DIR_J1603 = PROJECT_ROOT / "data" / "secondary" / "j1603"
SEC_DIR_JIAMUSI = PROJECT_ROOT / "data" / "secondary" / "jiamusi"
SEC_DIR_MEERKAT = PROJECT_ROOT / "data" / "secondary" / "meerkat"
RESULTS_DIR = PROJECT_ROOT / "results"

# J0437-4715 Kinematic Parameters (distance/PM from Reardon et al. 2021; binary from 2024 PPTA-DR3 in config.py)
J0437_DIST_PC = 156.3
J0437_PM_RA = 121.439  # mas/yr
J0437_PM_DEC = -71.438  # mas/yr
J0437_S_SCREEN = 0.6  # Relative screen distance (D_s / D_p)

# J1603-7202 Kinematic Parameters (from Walker et al. 2022)
J1603_DIST_PC = 250.0
J1603_PM_RA = -6.6  # mas/yr
J1603_PM_DEC = -25.9  # mas/yr
J1603_PB_DAYS = 0.0  # Non-binary
J1603_T0_MJD = 0.0
J1603_K_KMS = 0.0
J1603_S_SCREEN = 0.5


@dataclass
class CrossTermMeasurement:
    """Result of measuring a cross-term in the secondary spectrum."""

    tau_meas: float
    fD_meas: float
    tau_pred: float
    fD_pred: float
    peak_value: float
    snr: float
    valid: bool
    tau_offset_frac: float


@dataclass
class ClosureResult:
    """Complete closure measurement for a triplet of arclets.

    delta_us: Closure delay in microseconds
    sigma_us: Uncertainty estimate from error propagation (sigma_δ = √3 * sigma_tau)
    snr: Triplet SNR from cross-term measurements (independent of closure delay)
    closure_snr: Closure SNR = |delta|/sigma (for statistical analysis only, NOT for filtering)
    geometric_delta_us: Geometrically aligned closure delay (delta * geom_sign)
    geom_sign: Signed 2D geometric orientation of the loop (+1 CCW, -1 CW)
    phase_closure_rad: Closure of cross-term phases, wrapped to [-pi, pi).
    phase_delta_ns: Narrowband equivalent delay (ns) if psi = 2*pi*f_Hz*tau with
        f_Hz = nu_ref_mhz*1e6; i.e. tau_ns = psi*1e3/(2*pi*nu_ref_mhz). Not the same
        observable as geometric_delta_us; see PHASE_DELTA_NS_DEFINITION on summaries.
    """

    tau_01: CrossTermMeasurement
    tau_12: CrossTermMeasurement
    tau_02: CrossTermMeasurement
    triplet_indices: List[int]
    arclet_snrs: List[float]
    delta_us: float
    sigma_us: float
    snr: float
    closure_snr: float
    geometric_delta_us: float
    geom_sign: float
    saa_sign: float
    phase_closure_rad: float
    phase_delta_ns: float


def parabolic_peak_interp(y_minus: float, y_peak: float, y_plus: float) -> float:
    """Parabolic interpolation for sub-pixel peak position.
    
    Returns NaN if denominator is near zero (flat parabola), indicating
    unreliable sub-pixel position. Caller must check with np.isnan().
    """
    denom = y_minus - 2 * y_peak + y_plus
    if abs(denom) < PARABOLIC_DENOM_TOLERANCE:
        return float('nan')  # Invalid measurement - flat parabola
    offset = (y_minus - y_plus) / (2 * denom)
    return float(np.clip(offset, -1.0, 1.0))


def measure_cross_term_subpixel(
    S: np.ndarray,
    tau_us: np.ndarray,
    fD_mHz: np.ndarray,
    arclet_i: np.ndarray,
    arclet_j: np.ndarray,
    tau_bounds: tuple,
    fD_bounds: tuple,
    margin: float = EDGE_MARGIN,
    verbose: bool = False,
) -> CrossTermMeasurement:
    """Measure cross-term position with sub-pixel precision."""
    tau_i, fD_i, _ = arclet_i
    tau_j, fD_j, _ = arclet_j

    tau_pred = tau_j - tau_i
    fD_pred = fD_j - fD_i

    tau_min, tau_max = tau_bounds
    fD_min, fD_max = fD_bounds

    in_bounds = (
        tau_min + margin <= tau_pred <= tau_max - margin
        and fD_min + margin <= fD_pred <= fD_max - margin
    )

    if not in_bounds:
        return CrossTermMeasurement(
            tau_meas=tau_pred,
            fD_meas=fD_pred,
            tau_pred=tau_pred,
            fD_pred=fD_pred,
            peak_value=0.0,
            snr=0.0,
            valid=False,
            tau_offset_frac=0.0,
        )

    tau_idx = np.argmin(np.abs(tau_us - tau_pred))
    fD_idx = np.argmin(np.abs(fD_mHz - fD_pred))

    # Use configurable search radius based on typical arclet FWHM
    search_radius = SEARCH_RADIUS_PIXELS

    tau_start = max(1, tau_idx - search_radius)
    tau_end = min(len(tau_us) - 1, tau_idx + search_radius + 1)
    fD_start = max(1, fD_idx - search_radius // 2)
    fD_end = min(len(fD_mHz) - 1, fD_idx + search_radius // 2 + 1)

    if tau_end <= tau_start or fD_end <= fD_start:
        return CrossTermMeasurement(
            tau_meas=tau_pred,
            fD_meas=fD_pred,
            tau_pred=tau_pred,
            fD_pred=fD_pred,
            peak_value=0.0,
            snr=0.0,
            valid=False,
            tau_offset_frac=0.0,
        )

    sub_S = S[fD_start:fD_end, tau_start:tau_end].copy()

    if sub_S.size == 0 or sub_S.max() < 1e-30:
        return CrossTermMeasurement(
            tau_meas=tau_pred,
            fD_meas=fD_pred,
            tau_pred=tau_pred,
            fD_pred=fD_pred,
            peak_value=0.0,
            snr=0.0,
            valid=False,
            tau_offset_frac=0.0,
        )

    peak_local = np.unravel_index(np.argmax(sub_S), sub_S.shape)
    peak_fD_local = peak_local[0]
    peak_tau_local = peak_local[1]
    peak_val = sub_S[peak_local]

    tau_offset_frac = 0.0
    parabola_valid = True
    if 1 <= peak_tau_local < sub_S.shape[1] - 1:
        y_minus = sub_S[peak_fD_local, peak_tau_local - 1]
        y_center = sub_S[peak_fD_local, peak_tau_local]
        y_plus = sub_S[peak_fD_local, peak_tau_local + 1]
        tau_offset_frac = parabolic_peak_interp(y_minus, y_center, y_plus)
        parabola_valid = not np.isnan(tau_offset_frac)

    peak_tau_idx_float = tau_start + peak_tau_local + tau_offset_frac
    peak_fD_idx = fD_start + peak_fD_local

    tau_idx_int = int(np.clip(np.floor(peak_tau_idx_float), 0, len(tau_us) - 2))
    tau_frac = peak_tau_idx_float - tau_idx_int
    tau_meas = tau_us[tau_idx_int] * (1 - tau_frac) + tau_us[tau_idx_int + 1] * tau_frac

    fD_meas = fD_mHz[peak_fD_idx]

    tau_deviation = abs(tau_meas - tau_pred)
    fD_deviation = abs(fD_meas - fD_pred)

    d_tau = np.median(np.abs(np.diff(tau_us)))
    d_fD = np.median(np.abs(np.diff(fD_mHz)))
    # Use configurable tolerance factor (grid spacing multiple)
    max_deviation_tau = MAX_DEVIATION_FACTOR * d_tau
    max_deviation_fD = MAX_DEVIATION_FACTOR * d_fD

    # Measurement is invalid if parabola was flat (NaN) OR deviation too large
    valid = parabola_valid and (tau_deviation <= max_deviation_tau) and (fD_deviation <= max_deviation_fD)

    # 1. GLOBAL SNR check
    global_bg_level = np.median(S)
    global_snr = peak_val / global_bg_level if global_bg_level > 0 else 0.0

    # 2. LOCAL SNR check (Noise Floor Rejection)
    # Check 11x11 region around peak for diffuse noise dominance
    win = 5
    loc_f_start = max(0, peak_fD_idx - win)
    loc_f_end = min(len(fD_mHz), peak_fD_idx + win + 1)
    loc_t_start = max(0, int(np.floor(peak_tau_idx_float)) - win)
    loc_t_end = min(len(tau_us), int(np.floor(peak_tau_idx_float)) + win + 1)
    
    local_patch = S[loc_f_start:loc_f_end, loc_t_start:loc_t_end].astype(float, copy=True)
    if local_patch.size > 9:
        # Median of patch excluding central 3×3 around the peak (avoids biasing local_bg high)
        pi_f = peak_fD_idx - loc_f_start
        pi_t = int(np.floor(peak_tau_idx_float)) - loc_t_start
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                fi, tj = pi_f + di, pi_t + dj
                if 0 <= fi < local_patch.shape[0] and 0 <= tj < local_patch.shape[1]:
                    local_patch[fi, tj] = np.nan
        local_bg = float(np.nanmedian(local_patch))
        if not np.isfinite(local_bg) or local_bg <= 0:
            local_snr = global_snr
        else:
            local_snr = peak_val / local_bg
    else:
        local_snr = global_snr

    # Requirement: Global SNR >= 5.0 AND Local SNR >= 3.0
    # This ensures the peak is not just a noise fluctuation on a diffuse background
    snr_valid = (global_snr >= 5.0) and (local_snr >= 3.0)

    # Combined validation
    valid = parabola_valid and (tau_deviation <= max_deviation_tau) and (fD_deviation <= max_deviation_fD) and snr_valid

    if verbose:
        print_status(f"Sub-pixel measurement: tau_pred={tau_pred:.4f}, fD_pred={fD_pred:.2f}", "DATA")
        print_status(f"  Refined peak: tau_meas={tau_meas:.4f} us (offset={tau_offset_frac:.3f}), SNR={global_snr:.2f}", "CALC")
        if not valid:
            if not parabola_valid:
                print_status(f"  INVALID: Flat parabola - unreliable sub-pixel position", "WARNING")
            else:
                print_status(f"  INVALID: Dev tau={tau_deviation:.3f} > {max_deviation_tau:.3f} or fD={fD_deviation:.2f} > {max_deviation_fD:.2f}", "WARNING")

    return CrossTermMeasurement(
        tau_meas=float(tau_meas),
        fD_meas=float(fD_meas),
        tau_pred=float(tau_pred),
        fD_pred=float(fD_pred),
        peak_value=float(peak_val),
        snr=float(global_snr),
        valid=valid,
        tau_offset_frac=float(tau_offset_frac),
    )


# Pulsar kinematic and ISM parameter database
PULSAR_PARAMS = {
    "J0437-4715": {
        "dist": J0437_DIST_PC,
        "pm_ra": J0437_PM_RA,
        "pm_dec": J0437_PM_DEC,
        "pb": J0437_PB_DAYS,
        "t0": J0437_T0_MJD,
        "a1": J0437_A1_LC,
        "s": J0437_S_SCREEN,
        "psi": float(J0437_PSI_DEG),
    },
    "J1603-7202": {
        "dist": J1603_DIST_PC,
        "pm_ra": J1603_PM_RA,
        "pm_dec": J1603_PM_DEC,
        "pb": J1603_PB_DAYS,
        "t0": J1603_T0_MJD,
        "a1": 0.0,
        "s": J1603_S_SCREEN,
        "psi": 105.0,  # From Walker et al. 2022
    },
    "B0329+54": {
        "dist": 1000.0,
        "pm_ra": 7.10,
        "pm_dec": -11.75,
        "pb": 0.0,
        "t0": 0.0,
        "a1": 0.0,
        "s": 0.99,
        "psi": 110.0,
    },
    "B0355+54": {
        "dist": 1000.0,
        "pm_ra": 9.17,
        "pm_dec": 0.70,
        "pb": 0.0,
        "t0": 0.0,
        "a1": 0.0,
        "s": 0.62,
        "psi": 95.0,
    },
    "B0540+23": {
        "dist": 1600.0,
        "pm_ra": 2.50,
        "pm_dec": -21.80,
        "pb": 0.0,
        "t0": 0.0,
        "a1": 0.0,
        "s": 0.99,
        "psi": 120.0,
    },
    "B0740-28": {
        "dist": 2000.0,
        "pm_ra": -2.44,
        "pm_dec": -0.09,
        "pb": 0.0,
        "t0": 0.0,
        "a1": 0.0,
        "s": 0.5,
        "psi": 0.0,
    },
    "B1508+55": {
        "dist": 2100.0,
        "pm_ra": -73.70,
        "pm_dec": -62.70,
        "pb": 0.0,
        "t0": 0.0,
        "a1": 0.0,
        "s": 0.34,
        "psi": 0.0,
    },
    "B1933+16": {
        "dist": 3700.0,
        "pm_ra": -2.00,
        "pm_dec": -0.10,
        "pb": 0.0,
        "t0": 0.0,
        "a1": 0.0,
        "s": 0.94,
        "psi": 0.0,
    },
    "B2154+40": {
        "dist": 2900.0,
        "pm_ra": 14.60,
        "pm_dec": -2.60,
        "pb": 0.0,
        "t0": 0.0,
        "a1": 0.0,
        "s": 0.75,
        "psi": 0.0,
    },
    "B2310+42": {
        "dist": 1060.0,
        "pm_ra": -3.00,
        "pm_dec": -6.00,
        "pb": 0.0,
        "t0": 0.0,
        "a1": 0.0,
        "s": 0.99,
        "psi": 0.0,
    },
    "B2324+60": {
        "dist": 2700.0,
        "pm_ra": 0.0,
        "pm_dec": 0.0,
        "pb": 0.0,
        "t0": 0.0,
        "a1": 0.0,
        "s": 0.5,
        "psi": 0.0,
    },
    "B2351+61": {
        "dist": 2400.0,
        "pm_ra": -0.19,
        "pm_dec": -0.01,
        "pb": 0.0,
        "t0": 0.0,
        "a1": 0.0,
        "s": 0.5,
        "psi": 0.0,
    },
    "J0908-1739": {
        "dist": 400.0,
        "pm_ra": 0.0,
        "pm_dec": 0.0,
        "pb": 0.0,
        "t0": 0.0,
        "a1": 0.0,
        "s": 0.5,
        "psi": 0.0,
    },
    "J0922-0638": {
        "dist": 1000.0,
        "pm_ra": 0.0,
        "pm_dec": 0.0,
        "pb": 0.0,
        "t0": 0.0,
        "a1": 0.0,
        "s": 0.5,
        "psi": 0.0,
    },
    "J1731-4744": {
        "dist": 400.0,
        "pm_ra": 0.0,
        "pm_dec": 0.0,
        "pb": 0.0,
        "t0": 0.0,
        "a1": 0.0,
        "s": 0.5,
        "psi": 0.0,
    },
}


@functools.lru_cache(maxsize=1024)
def calculate_velocity_vector(
    mjd: float, pulsar_name: str = "J0437-4715", verbose: bool = False
) -> np.ndarray:
    """Calculate the relative transverse velocity vector for a specific pulsar.
    
    Uses caching to avoid redundant astropy calculations for repeated MJDs.
    """
    # Strict parameter lookup: raise KeyError if pulsar not in database
    # This prevents "silent fallbacks" and ensures no synthetic results.
    if pulsar_name not in PULSAR_PARAMS:
        raise KeyError(f"Pulsar '{pulsar_name}' not found in kinematic database. Cannot perform authentic TEP alignment.")
    
    p = PULSAR_PARAMS[pulsar_name]

    # 1. Proper Motion (converted to km/s)
    v_pm_x = 4.74 * p["pm_ra"] * p["dist"] / 1000.0
    v_pm_y = 4.74 * p["pm_dec"] * p["dist"] / 1000.0
    v_pm = np.array([v_pm_x, v_pm_y])

    t = Time(mjd, format='mjd')
    _, v_earth_bary = get_body_barycentric_posvel('earth', t)
    # Convert AU/day to km/s
    v_earth_km_s = v_earth_bary.xyz.value * 149597870.7 / 86400.0
    # Project to transverse plane (approximate for J0437-4715 coordinates)
    # Using local RA/Dec projection:
    v_earth = np.array([v_earth_km_s[0], v_earth_km_s[1]])

    # 3. Pulsar binary reflex velocity in the same tangent basis as proper motion (km/s)
    v_orbit = np.array([0.0, 0.0])
    if p["pb"] > 0:
        if pulsar_name == "J0437-4715":
            v_orbit = reflex_binary_transverse_velocity_kms(
                float(mjd),
                float(p["pb"]),
                float(p["t0"]),
                float(p["a1"]),
                float(J0437_ECC),
                float(J0437_OM_DEG),
                float(J0437_INC_DEG),
                float(J0437_KOM_DEG),
                float(C_LIGHT_KM_S),
                float(J0437_RA_RAD),
                float(J0437_DEC_RAD),
            )
        else:
            raise NotImplementedError(
                f"Binary sky-plane reflex velocity is only wired for J0437-4715; "
                f"received '{pulsar_name}'. Add orbit geometry to PULSAR_PARAMS and "
                f"extend binary_reflex_kinematics if this pulsar is required."
            )

    # Effective Scintillation Velocity
    weight_pulsar = (1.0 - p["s"]) / p["s"]
    v_eff = weight_pulsar * (v_pm + v_orbit) + v_earth

    if verbose:
        print_status(f"Kinematic alignment for {pulsar_name} (MJD {mjd:.2f}):", "PHYSICS")
        print_status(f"  v_pm: {v_pm}, v_earth: {v_earth}, v_orbit: {v_orbit}", "PHYSICS")
        print_status(f"  v_eff (transverse): {v_eff} km/s, project_axis_psi: {p['psi']}°", "PHYSICS")

    return v_eff


def circular_mean_and_rbar(angles: np.ndarray, weights: Optional[np.ndarray] = None) -> tuple:
    """Circular mean and mean resultant length.

    Returns
    -------
    psi_mean : float
        Circular mean in [-pi, pi).
    r_bar : float
        Mean resultant length (0 <= r_bar <= 1).
    """
    if weights is None:
        weights = np.ones_like(angles)
    weights = np.asarray(weights)
    z = np.sum(weights * np.exp(1j * angles))
    w_sum = np.sum(weights)
    r_bar = float(np.abs(z) / w_sum) if w_sum > 0 else 0.0
    psi_mean = float(np.angle(z / w_sum)) if w_sum > 0 else 0.0
    return psi_mean, r_bar


def rayleigh_test(angles: np.ndarray, weights: Optional[np.ndarray] = None) -> tuple:
    """Rayleigh test for uniformity of circular data.

    Tests H0: angles are uniformly distributed (no preferred direction).
    For large n, 2n*R_bar^2 ~ chi-squared(2).

    Returns
    -------
    z_stat : float
        Rayleigh Z = 2n * R_bar^2.
    p_value : float
        Two-tailed p-value from chi-squared(2).
    """
    n = len(angles)
    if n < 3:
        return 0.0, 1.0
    _, r_bar = circular_mean_and_rbar(angles, weights)
    # Effective sample size with weights
    if weights is not None:
        w = np.asarray(weights)
        n_eff = float(np.sum(w)**2 / np.sum(w**2))
    else:
        n_eff = float(n)
    z_stat = float(2.0 * n_eff * r_bar**2)
    p_value = float(stats.chi2.sf(z_stat, 2))
    return z_stat, p_value


def v_test_circular(angles: np.ndarray, mu0: float = 0.0, weights: Optional[np.ndarray] = None) -> tuple:
    """V-test for a specified mean direction of circular data.

    Tests H0: angles are uniformly distributed against H1: concentration
    around the specified direction mu0. For large n, V ~ N(0,1) under H0.

    Returns
    -------
    v_stat : float
        V = sqrt(2n) * R_bar * cos(psi_mean - mu0).
    p_value : float
        Two-tailed p-value.
    """
    n = len(angles)
    if n < 3:
        return 0.0, 1.0
    psi_mean, r_bar = circular_mean_and_rbar(angles, weights)
    if weights is not None:
        w = np.asarray(weights)
        n_eff = float(np.sum(w)**2 / np.sum(w**2))
    else:
        n_eff = float(n)
    v_stat = float(np.sqrt(2.0 * n_eff) * r_bar * np.cos(psi_mean - mu0))
    p_value = float(2.0 * stats.norm.sf(abs(v_stat)))
    return v_stat, p_value


def circular_bootstrap_ci(
    angles: np.ndarray,
    weights: Optional[np.ndarray] = None,
    n_boot: int = 10000,
    seed: int = 42,
) -> tuple:
    """Bootstrap confidence interval for circular mean.

    Resamples angles with replacement and computes the circular mean each
    time. Returns a standard error and a 95% confidence interval that is
    robust to the -pi/pi branch cut.

    Returns
    -------
    se : float
        Bootstrap standard error (radians).
    ci_95 : tuple
        (lower, upper) 95% confidence interval in radians.
    """
    n = len(angles)
    if n < 3:
        return float("inf"), (float("-inf"), float("inf"))
    if weights is None:
        weights = np.ones(n)
    weights = np.asarray(weights)
    psi_mean, _ = circular_mean_and_rbar(angles, weights)
    rng = np.random.RandomState(seed)
    boot_means = []
    for _ in range(n_boot):
        idx = rng.choice(n, size=n, replace=True)
        boot_psi, _ = circular_mean_and_rbar(angles[idx], weights[idx])
        boot_means.append(boot_psi)
    boot_means = np.array(boot_means)
    # Center at observed mean to unwrap branch cut
    centered = (boot_means - psi_mean + np.pi) % (2.0 * np.pi) - np.pi
    se = float(np.std(centered, ddof=1))
    ci_low = float((psi_mean + np.percentile(centered, 2.5) + np.pi) % (2.0 * np.pi) - np.pi)
    ci_high = float((psi_mean + np.percentile(centered, 97.5) + np.pi) % (2.0 * np.pi) - np.pi)
    return se, (ci_low, ci_high)


def compute_closure(
    S: np.ndarray,
    S_complex: np.ndarray,  # New: Complex secondary spectrum
    tau_us: np.ndarray,
    fD_mHz: np.ndarray,
    arclet_i: np.ndarray,
    arclet_j: np.ndarray,
    arclet_k: np.ndarray,
    tau_bounds: tuple,
    fD_bounds: tuple,
    v_eff: np.ndarray,
    mjd: float = 0.0,
    nu_ref_mhz: float = None,  # Must be provided from epoch metadata
    pulsar_name: str = "J0437-4715",
) -> Optional[ClosureResult]:
    """Compute closure delay and phase closure for a triplet of arclets.
    
    Note: nu_ref_mhz is required and must come from epoch metadata.
    Hardcoding would compromise chromatic vs achromatic discrimination tests.
    """
    if nu_ref_mhz is None:
        raise ValueError("nu_ref_mhz is required - must be provided from epoch metadata (e.g., frequency_mhz from catalog). "
                        "Hardcoding 1400 MHz would compromise chromatic vs achromatic discrimination tests.")
    
    # Sort arclets by tau to ensure consistent orientation and measurement
    # regardless of input index order (from SNR-based sorting)
    # This eliminates "Index Bias" from the discovery statistics.
    sorted_arclets = sorted([arclet_i, arclet_j, arclet_k], key=lambda x: x[0])
    a0, a1, a2 = sorted_arclets

    ct_01 = measure_cross_term_subpixel(
        S, tau_us, fD_mHz, a0, a1, tau_bounds, fD_bounds, verbose=False
    )
    ct_12 = measure_cross_term_subpixel(
        S, tau_us, fD_mHz, a1, a2, tau_bounds, fD_bounds, verbose=False
    )
    ct_02 = measure_cross_term_subpixel(
        S, tau_us, fD_mHz, a0, a2, tau_bounds, fD_bounds, verbose=False
    )

    if not (ct_01.valid and ct_12.valid and ct_02.valid):
        return None

    # PHASE CLOSURE EXTRACTION (The true TEP signal)
    # Extract complex phase at the measured peak positions
    def get_phase_at_peak(ct):
        # Find closest indices for the complex grid
        f_idx = np.argmin(np.abs(fD_mHz - ct.fD_meas))
        t_idx = np.argmin(np.abs(tau_us - ct.tau_meas))

        # Average phase over a small 3x3 region around the peak to reduce noise
        f_start = max(0, f_idx - 1)
        f_end = min(len(fD_mHz), f_idx + 2)
        t_start = max(0, t_idx - 1)
        t_end = min(len(tau_us), t_idx + 2)

        region = S_complex[f_start:f_end, t_start:t_end]
        # Compute phase of the complex average
        avg_complex = np.mean(region)
        return np.angle(avg_complex)

    phi_01 = get_phase_at_peak(ct_01)
    phi_12 = get_phase_at_peak(ct_12)
    phi_02 = get_phase_at_peak(ct_02)

    # Wrap closure phase to [-pi, pi]
    psi = (phi_01 + phi_12 - phi_02 + np.pi) % (2 * np.pi) - np.pi

    # ALIGNMENT (Kinematic + Geometric)
    dx1 = a1[0] - a0[0]
    dy1 = a1[1] - a0[1]
    dx2 = a2[0] - a0[0]
    dy2 = a2[1] - a0[1]
    cross_prod = dx1 * dy2 - dy1 * dx2
    # Use > 0 instead of >= 0 to avoid bias (>= gives 52.5% positive for random data)
    geom_sign = 1.0 if cross_prod > 0 else -1.0

    if pulsar_name not in PULSAR_PARAMS:
        raise KeyError(
            f"Pulsar '{pulsar_name}' not found in kinematic database. "
            "Cannot perform authentic TEP alignment."
        )
    pulsar_info = PULSAR_PARAMS[pulsar_name]
    rad_psi = np.radians(pulsar_info["psi"])
    v_proj = v_eff[0] * np.cos(rad_psi) + v_eff[1] * np.sin(rad_psi)
    saa_sign = 1.0 if v_proj >= 0 else -1.0

    # Primary Measurement: Raw Delay Closure
    delta_us_raw = ct_01.tau_meas + ct_12.tau_meas - ct_02.tau_meas

    # Stokes-theorem alignment
    # Velocity normalization: 50 km/s is a characteristic ISS velocity scale
    # This represents typical interstellar scintillation velocities for pulsars,
    # providing a dimensionless weighting factor that scales the geometric alignment
    # with the projected pulsar velocity relative to the scattering screen
    v_weight = v_proj / 50.0  # characteristic velocity scale (typical ISS velocity)
    geometric_delta = float(delta_us_raw * geom_sign * v_weight)
    delta_us_val = float(delta_us_raw)

    # Secondary Measurement: narrowband phase-to-time equivalent (ns) at band centre.
    # SI: tau_s = psi / (2*pi*f_Hz), f_Hz = nu_ref_mhz * 1e6  =>  tau_ns = psi * 1e3 / (2*pi*nu_ref_mhz).
    phase_delta_ns = (psi / (2 * np.pi * nu_ref_mhz)) * 1e3

    d_tau = np.median(np.abs(np.diff(tau_us)))
    sigma = d_tau * np.sqrt(3)
    triplet_snr = np.mean([ct_01.snr, ct_12.snr, ct_02.snr])
    closure_snr = abs(delta_us_raw) / sigma if sigma > 0 else 0.0

    return ClosureResult(
        delta_us=delta_us_val,
        sigma_us=float(sigma),
        snr=float(triplet_snr),
        closure_snr=float(closure_snr),
        geometric_delta_us=geometric_delta,
        geom_sign=float(geom_sign),
        saa_sign=float(saa_sign),
        phase_closure_rad=float(psi),
        phase_delta_ns=float(phase_delta_ns),
        tau_01=ct_01,
        tau_12=ct_12,
        tau_02=ct_02,
        triplet_indices=[0, 1, 2], # Not used for logic, just for dataclass
        arclet_snrs=[float(a0[2]), float(a1[2]), float(a2[2])],
    )


def find_best_triplets(
    S: np.ndarray,
    S_complex: np.ndarray,
    tau_us: np.ndarray,
    fD_mHz: np.ndarray,
    arclets: np.ndarray,
    v_eff: np.ndarray,
    eta1: float = 0.0,
    eta2: float = 0.0,
    mjd: float = 0.0,
    max_triplets: int = MAX_TRIPLETS_PER_EPOCH,
    use_snr_threshold: bool = True,
    min_snr: float = 5.0,  # 5-sigma SNR threshold (particle physics standard)
    nu_ref_mhz: float = None,  # Must be provided from epoch metadata
    pulsar_name: str = "J0437-4715",
) -> List[ClosureResult]:
    """
    Find measurable triplets, with PRIORITY for cross-screen loops.
    
    Note: nu_ref_mhz is required and must come from epoch metadata.
    Hardcoding would compromise chromatic vs achromatic discrimination tests.
    """
    if nu_ref_mhz is None:
        raise ValueError("nu_ref_mhz is required - must be provided from epoch metadata (e.g., frequency_mhz from catalog). "
                        "Hardcoding 1400 MHz would compromise chromatic vs achromatic discrimination tests.")
    
    if arclets is None:
        return []
    
    tau_bounds = (float(tau_us.min()), float(tau_us.max()))
    fD_bounds = (float(fD_mHz.min()), float(fD_mHz.max()))

    n = len(arclets)
    if n < 3:
        return []

    # Assign arclets to screens based on curvature proximity
    # arclet = [tau, fD, snr]
    arc1_indices = []
    arc2_indices = []

    for idx, a in enumerate(arclets):
        tau, fD, _ = a
        if fD == 0:
            continue

        # Calculate residuals for both screens
        res1 = abs(tau - eta1 * fD**2)
        res2 = abs(tau - eta2 * fD**2) if eta2 > 0 else 1e9

        # Assign to closest arc if within tolerance (e.g. 1us)
        if res1 < res2 and res1 < 1.0:
            arc1_indices.append(idx)
        elif res2 < res1 and res2 < 1.0:
            arc2_indices.append(idx)
        # NO FALLBACK: Ambiguous arclets are NOT assigned
        # This prevents bias from arbitrary assignment to arc1

    # Prioritize triplets that span both screens (Cross-Screen Holonomy)
    cross_triplets = []
    # Combination 1: Two from Arc 1, One from Arc 2
    if len(arc1_indices) >= 2 and len(arc2_indices) >= 1:
        for i, j in combinations(arc1_indices[:6], 2):
            for k in arc2_indices[:6]:
                cross_triplets.append((i, j, k))

    # Combination 2: One from Arc 1, Two from Arc 2
    if len(arc1_indices) >= 1 and len(arc2_indices) >= 2:
        for i in arc1_indices[:6]:
            for j, k in combinations(arc2_indices[:6], 2):
                cross_triplets.append((i, j, k))

    print_status(f"Triplet selection: Found {len(arc1_indices)} arclets on Arc 1, {len(arc2_indices)} on Arc 2", "DATA")
    print_status(f"  Prioritized cross-screen triplets: {len(cross_triplets)}", "DATA")

    # NO FALLBACK: If no cross-screen triplets, epoch is skipped
    # Cross-screen triplets are required for authentic TEP measurement
    # Using same-screen triplets would measure standard ISS, not TEP
    if not cross_triplets:
        print_status(f"  No cross-screen triplets found - skipping epoch (requires cross-screen for TEP measurement)", "WARNING")
        return None

    results = []
    for i, j, k in cross_triplets:
        a0, a1, a2 = arclets[i], arclets[j], arclets[k]
        result = compute_closure(
            S,
            S_complex,
            tau_us,
            fD_mHz,
            a0,
            a1,
            a2,
            tau_bounds,
            fD_bounds,
            v_eff,
            mjd=mjd,
            nu_ref_mhz=nu_ref_mhz,
            pulsar_name=pulsar_name,
        )
        if result is not None:
            # Transfer triplet indices
            result.triplet_indices = [int(i), int(j), int(k)]
            results.append(result)

    if not results:
        return []

    # Apply selection criterion (cap triplets per epoch to limit correlated resampling)
    if use_snr_threshold:
        filtered = [r for r in results if r.snr >= min_snr]
        filtered.sort(key=lambda x: -x.snr)
        return filtered[:max_triplets]
    else:
        results.sort(key=lambda x: -x.snr)
        return results[:max_triplets]


def process_epoch(
    sec_path: Path,
    max_arclets: int = MAX_ARCLETS_FOR_TRIPLETS,
    use_snr_threshold: bool = True,
    min_snr: float = 5.0,  # 5-sigma SNR threshold (particle physics standard)
    pulsar_name: str = "J0437-4715",
) -> Optional[Dict[str, Any]]:
    """Process a single epoch and return closure measurements.
    
    Uses memory-mapped arrays to reduce memory usage when loading large secondary spectra.
    """
    try:
        # Use memory-mapped arrays to avoid loading full data into memory
        data = np.load(sec_path, allow_pickle=False, mmap_mode='r')
        
        # Copy all needed data before closing the memory-mapped file
        S = data["secondary"].copy()
        tau_us = data["tau_us"].copy()
        fD_mHz = data["fD_mHz"].copy()
        
        # Check if arclets key exists before accessing
        if "arclets" not in data:
            data.close()
            print_status(f"  [FAIL] {sec_path.stem}: Missing arclets data in secondary spectrum", "ERROR")
            return None
        
        arclets = data["arclets"].copy()
        mjd = float(data["mjd_start"])
        eta1 = data.get("eta_screen1")
        eta2 = data.get("eta_screen2")
        freq_mhz = data.get("frequency_mhz")
        
        # secondary_complex is only saved for files with >= 3 arclets
        # Load it conditionally - will be None if not present (checked later)
        S_complex = data.get("secondary_complex")
        if S_complex is not None:
            S_complex = S_complex.copy()
        
        # Close the memory-mapped file to release file handle
        data.close()
        
        if eta1 is None:
            print_status(f"  [FAIL] {sec_path.stem}: Missing screen curvature (eta_screen1). Cannot perform structural alignment.", "ERROR")
            return None
    except Exception as e:
        print_status(f"  [FAIL] {sec_path.stem}: Failed to load secondary spectrum: {e}", "ERROR")
        return None

    # Additional safety check - handle case where arclets might not be defined
    try:
        if arclets is None or len(arclets) < 3:
            return None
    except NameError:
        print_status(f"  [FAIL] {sec_path.stem}: arclets not defined after loading", "ERROR")
        return None
    
    # Now that the analysis has confirmed sufficient arclets, require complex data
    if S_complex is None:
        print_status(f"  [FAIL] {sec_path.stem}: Missing secondary_complex data (needed for phase closure)", "ERROR")
        return None

    # Additional safety checks for variables that might be None if exception occurred
    if mjd is None or S is None or tau_us is None or fD_mHz is None:
        print_status(f"  [FAIL] {sec_path.stem}: Required data missing after loading", "ERROR")
        return None

    v_eff = calculate_velocity_vector(mjd, pulsar_name=pulsar_name, verbose=True)

    results = find_best_triplets(
        S,
        S_complex,
        tau_us,
        fD_mHz,
        arclets,
        v_eff,
        eta1=eta1,
        eta2=eta2,
        mjd=mjd,
        max_triplets=MAX_TRIPLETS_PER_EPOCH,
        use_snr_threshold=use_snr_threshold,
        min_snr=min_snr,
        nu_ref_mhz=freq_mhz,  # Pass frequency from epoch metadata (may be None for older data)
        pulsar_name=pulsar_name,
    )

    # Store n_arclets before cleanup
    n_arclets = len(arclets)
    
    # Explicit memory cleanup for large arrays (removed arclets deletion to prevent UnboundLocalError)
    del S, S_complex, tau_us, fD_mHz

    if not results:
        return None

    def format_closure_result(r: ClosureResult) -> Dict[str, Any]:
        return {
            "delta_us": r.delta_us,
            "sigma_us": r.sigma_us,
            "snr": r.snr,  # Cross-term SNR (independent of closure delay, used for filtering)
            "closure_snr": r.closure_snr,  # Circular closure SNR (for statistical analysis only)
            "triplet": r.triplet_indices,
            "arclet_snrs": r.arclet_snrs,
            "tau_01": r.tau_01.tau_meas,
            "tau_12": r.tau_12.tau_meas,
            "tau_02": r.tau_02.tau_meas,
            "tau_01_pred": r.tau_01.tau_pred,
            "tau_12_pred": r.tau_12.tau_pred,
            "tau_02_pred": r.tau_02.tau_pred,
            "fD_01": r.tau_01.fD_meas,
            "fD_12": r.tau_12.fD_meas,
            "fD_02": r.tau_02.fD_meas,
            "cross_term_snrs": [
                r.tau_01.snr,
                r.tau_12.snr,
                r.tau_02.snr,
            ],  # Individual cross-term SNRs
            "geometric_delta_us": float(r.geometric_delta_us),
            "geom_sign": float(r.geom_sign),
            "saa_sign": float(r.saa_sign),
            "phase_closure_rad": float(r.phase_closure_rad),
            "phase_delta_ns": float(r.phase_delta_ns),
        }

    return {
        "epoch": sec_path.stem,
        "mjd": mjd,
        "n_arclets": n_arclets,
        "n_triplets": len(results),
        "triplets": [format_closure_result(r) for r in results],
        "primary": format_closure_result(results[0]),
    }


def step_main(logger=None, verbose=True):
    """Pipeline entry point for Step 003."""
    if logger:
        set_step_logger(logger)
    return main()


def process_pulsar(
    sec_dir: Path,
    pulsar_name: str,
    max_arclets: int,
    use_snr_threshold: bool,
    min_snr: float,
    workers: int = None,
    mode: str = "full",
):
    """Process closure delays for a single pulsar with parallel processing."""
    print_status(f"\nProcessing {pulsar_name}...", "TITLE")
    print_status("=" * 60, "TITLE")

    # Extract pulsar prefix for filtering
    pulsar_prefix = pulsar_name.split("+")[0].split("-")[0]

    # Get all secondary files
    if mode == "sb0":
        sec_files = sorted(sec_dir.glob("*_sb0_secondary.npz"))
    elif mode == "sb1":
        sec_files = sorted(sec_dir.glob("*_sb1_secondary.npz"))
    else:
        # Default to full band files (excluding sub-bands)
        sec_files = sorted([f for f in sec_dir.glob("*_secondary.npz") if "_sb" not in f.name])

    # Filter secondary files by pulsar prefix for multi-pulsar directories
    # (Jiamusi and MeerKAT both contain multiple pulsars in one directory)
    if "jiamusi" in str(sec_dir).lower() or "meerkat" in str(sec_dir).lower():
        sec_files = [f for f in sec_files if f.stem.startswith(pulsar_prefix)]

    if not sec_files:
        print_status(f"No secondary spectra found for {pulsar_name}", "WARNING")
        return None, None

    if workers is None:
        workers = worker_count(role="cpu_bound", reserve=2)

    print_status(
        f"Processing {len(sec_files)} epochs with {workers} parallel workers...", "INFO"
    )

    all_results = []
    all_deltas = []
    viable = 0

    # Serial or Parallel epoch processing
    if workers == 1:
        print_status("  Running in serial mode for debugging...", "INFO")
        skipped_count = 0
        for sf in sec_files:
            try:
                result = process_epoch(
                    sf,
                    max_arclets,
                    use_snr_threshold,
                    min_snr,
                    pulsar_name=pulsar_name,
                )
                if result is None:
                    skipped_count += 1
                    continue
                viable += 1
                all_results.append(result)
                delta = result["primary"]["geometric_delta_us"]
                all_deltas.append(delta)
                print_status(
                    f"  [OK] {result['epoch'][:35]:35s} Delta={delta:+.4f} us", "INFO"
                )
            except Exception as e:
                import traceback
                print_status(f"  [FAIL] {sf.stem[:40]:40s} Error: {e}", "ERROR")
                traceback.print_exc()
        
        if skipped_count > 0:
            print_status(f"  Skipped {skipped_count} epochs with no measurable triplets", "INFO")
    else:
        # Parallel epoch processing using ProcessPoolExecutor
        skipped_count = 0
        with ProcessPoolExecutor(max_workers=workers) as executor:
            # Submit all tasks
            future_to_file = {
                executor.submit(
                    process_epoch,
                    sf,
                    max_arclets,
                    use_snr_threshold,
                    min_snr,
                    pulsar_name=pulsar_name,
                ): sf
                for sf in sec_files
            }

            # Collect results as they complete
            for future in as_completed(future_to_file):
                sf = future_to_file[future]
                try:
                    result = future.result()
                    if result is None:
                        skipped_count += 1
                        continue

                    viable += 1
                    all_results.append(result)

                    delta = result["primary"]["geometric_delta_us"]
                    all_deltas.append(delta)

                    print_status(
                        f"  [OK] {result['epoch'][:35]:35s} Delta={delta:+.4f} us", "INFO"
                    )
                except Exception as e:
                    import traceback
                    print_status(f"  [FAIL] {sf.stem[:40]:40s} Error: {e}", "ERROR")
                    traceback.print_exc()
        
        if skipped_count > 0:
            print_status(f"  Skipped {skipped_count} epochs with no measurable triplets", "INFO")

    if not all_deltas:
        # Only print warning for primary target pulsars, not control pulsars
        if pulsar_name in ["j0437", "J0437-4715"]:
            print_status(f"No valid closure measurements for {pulsar_name}", "WARNING")
        return None, None

    # --- CONSERVATIVE STATISTICAL FRAMEWORK: Independence Correction ---
    # The true number of independent degrees of freedom is limited by the number of 
    # observational epochs (viable) and the number of independent arclets, not the 
    # total number of triplets (n_total_triplets), which share correlated arclets.
    
    epoch_psi_means = []
    epoch_abs_delta_means = []
    epoch_delta_means = []
    epoch_weights = []
    epoch_geom_signs = []
    
    n_total_triplets = sum(len(r["triplets"]) for r in all_results)
    
    for r in all_results:
        triplets = r["triplets"]
        if not triplets:
            continue

        # Skip epochs with very low triplet counts (< 5) to avoid unreliable statistics
        # This creates a discrepancy: n_epochs counts all viable epochs (>=3 triplets)
        # while n_independent_samples counts only epochs with >=5 triplets used in aggregation
        if len(triplets) < 5:
            continue

        # 1. Phase Closure (primary): SNR^2-weighted circular mean within epoch
        ep_psi = np.array([t["phase_closure_rad"] for t in triplets])
        ep_delta = np.array([t["geometric_delta_us"] for t in triplets])
        ep_abs_delta = np.abs(ep_delta)
        ep_snr = np.array([float(t["snr"]) for t in triplets], dtype=float)
        w_psi = np.square(np.maximum(ep_snr, 1e-6))
        ep_psi_mean, _ = circular_mean_and_rbar(ep_psi, w_psi)
        
        # 2. Mean for delays
        ep_delta_mean = np.mean(ep_delta)
        ep_abs_delta_mean = np.mean(ep_abs_delta)
        
        # 3. Calculate internal uncertainty (SEM of triplets within epoch)
        # Note: This still uses N_triplets internally, but we aggregate epochs conservatively next
        ep_delta_std = np.std(ep_delta, ddof=1) if len(ep_delta) > 1 else 1e-3
        ep_sem = ep_delta_std / np.sqrt(len(ep_delta))
        
        # Use inverse-variance weighting for global aggregation
        # Ensure we don't divide by zero; 0.001 us (1 ns) is a floor for triplet precision
        weight = 1.0 / (ep_sem**2 + (0.001)**2)
        
        # Track majority geometric orientation for monopole/bipole decomposition
        ep_geom_signs = np.array([t.get("geom_sign", 0) for t in triplets])
        majority_geom = float(np.sign(np.sum(ep_geom_signs)))
        
        epoch_psi_means.append(ep_psi_mean)
        epoch_abs_delta_means.append(ep_abs_delta_mean)
        epoch_delta_means.append(ep_delta_mean)
        epoch_weights.append(weight)
        epoch_geom_signs.append(majority_geom)

    epoch_psi_means = np.array(epoch_psi_means)
    epoch_abs_delta_means = np.array(epoch_abs_delta_means)
    epoch_delta_means = np.array(epoch_delta_means)
    epoch_weights = np.array(epoch_weights)
    epoch_geom_signs = np.array(epoch_geom_signs)
    
    n_independent = len(epoch_weights)
    
    # --- Global Weighted Aggregation ---
    sum_w = np.sum(epoch_weights)
    
    # 1. Phase Closure (Primary Metric)
    # Vector weighted mean using circular statistics
    psi_mean, psi_rbar = circular_mean_and_rbar(epoch_psi_means, epoch_weights)
    # Numerical safety: floating error can push R_bar marginally outside [0, 1]
    psi_rbar = float(np.clip(psi_rbar, 0.0, 1.0))

    # Circular standard deviation and analytical standard error
    # sigma_0 = sqrt(-2 * ln(R_bar)) for concentrated wrapped-normal-like data
    if 1e-12 < psi_rbar < 1.0:
        psi_circ_std = float(np.sqrt(-2.0 * np.log(psi_rbar)))
        psi_circ_se = float(psi_circ_std / np.sqrt(n_independent))
    elif psi_rbar >= 1.0:
        psi_circ_std = 0.0
        psi_circ_se = 0.0
    else:
        psi_circ_std = float(np.pi)
        psi_circ_se = float(np.pi)

    # Legacy t-statistic (kept for backward compatibility with downstream scripts)
    psi_std_err = np.std(epoch_psi_means, ddof=1) / np.sqrt(n_independent)
    psi_t = psi_mean / psi_std_err if psi_std_err > 0 else 0.0

    # Rayleigh test for uniformity (H0: no preferred direction)
    rayleigh_z, rayleigh_p = rayleigh_test(epoch_psi_means, epoch_weights)

    # V-test against mu0 = 0 (GR prediction)
    v_stat, v_p = v_test_circular(epoch_psi_means, mu0=0.0, weights=epoch_weights)

    # Epoch-level bootstrap for circular mean
    psi_boot_se, psi_boot_ci = circular_bootstrap_ci(
        epoch_psi_means, epoch_weights, n_boot=10000
    )
    
    # --- Monopole / Bipole Vector Decomposition ---
    # Separate epochs by majority geometric orientation and compute vector means
    pos_mask = epoch_geom_signs > 0
    neg_mask = epoch_geom_signs < 0
    
    if np.sum(pos_mask) > 0:
        pos_z = np.mean(np.exp(1j * epoch_psi_means[pos_mask]))
        pos_angle = float(np.angle(pos_z))
        pos_rbar = float(np.abs(pos_z))
        pos_n = int(np.sum(pos_mask))
    else:
        pos_angle, pos_rbar, pos_n = 0.0, 0.0, 0
    
    if np.sum(neg_mask) > 0:
        neg_z = np.mean(np.exp(1j * epoch_psi_means[neg_mask]))
        neg_angle = float(np.angle(neg_z))
        neg_rbar = float(np.abs(neg_z))
        neg_n = int(np.sum(neg_mask))
    else:
        neg_angle, neg_rbar, neg_n = 0.0, 0.0, 0
    
    # Monopole = midpoint, Bipole = half-difference
    # ψ = ψ₀ + ψ₁ * geom_sign  =>  ψ₀ = (pos + neg)/2,  ψ₁ = (pos - neg)/2
    psi_monopole = float(((pos_angle + neg_angle) / 2.0 + np.pi) % (2 * np.pi) - np.pi)
    psi_bipole = float(((pos_angle - neg_angle) / 2.0 + np.pi) % (2 * np.pi) - np.pi)
    
    # Angular separation (smallest arc)
    angular_sep = float(abs(pos_angle - neg_angle))
    if angular_sep > np.pi:
        angular_sep = 2 * np.pi - angular_sep
    
    bipole_ratio = float(abs(psi_bipole) / abs(psi_monopole)) if abs(psi_monopole) > 0.1 else float('inf')
    
    # 2. Holonomy Magnitude |H| (Unsigned)
    H_mean_us = np.sum(epoch_weights * epoch_abs_delta_means) / sum_w if sum_w > 0 else 0.0
    # SEM consistent with inverse-variance epoch weights: Var(sum w_i x_i / sum w_i) = 1/sum w_i
    H_sem_us = float(np.sqrt(1.0 / sum_w)) if sum_w > 0 else 0.0
    H_sem_between_epoch_us = (
        float(np.std(epoch_abs_delta_means, ddof=1) / np.sqrt(n_independent))
        if n_independent > 1
        else 0.0
    )
    H_t = H_mean_us / H_sem_us if H_sem_us > 0 else 0.0
    # Gaussian two-sided p under independent epoch-likelihood approximation (meta-analysis convention)
    H_p = float(2.0 * stats.norm.sf(abs(H_t))) if sum_w > 0 else 1.0

    # 3. Signed Mean (Diagnostic)
    H_signed_mean_us = np.sum(epoch_weights * epoch_delta_means) / sum_w if sum_w > 0 else 0.0
    H_signed_sem_us = H_sem_us
    H_signed_t = H_signed_mean_us / H_signed_sem_us if H_signed_sem_us > 0 else 0.0
    
    # 4. Robust Trimmed Magnitude
    # We trim the epoch means themselves to be truly robust to epoch-level artifacts
    trim_frac = 0.10
    sorted_ep_abs = np.sort(epoch_abs_delta_means)
    n_trim = int(np.floor(trim_frac * n_independent))
    trimmed_ep = sorted_ep_abs[n_trim : n_independent - n_trim]
    if len(trimmed_ep) > 1:
        H_trim_mean = float(np.mean(trimmed_ep))
        H_trim_sem = float(np.std(trimmed_ep, ddof=1) / np.sqrt(len(trimmed_ep)))
        H_trim_t = H_trim_mean / H_trim_sem if H_trim_sem > 0 else 0.0
    else:
        H_trim_mean = H_mean_us
        H_trim_sem = H_sem_us
        H_trim_t = H_t

    # --- Kinematic Parameters ---
    # Compute transverse velocity from proper motion and distance
    p = PULSAR_PARAMS.get(pulsar_name, PULSAR_PARAMS["J0437-4715"])
    v_pm_ra_kms = 4.74 * p["pm_ra"] * p["dist"] / 1000.0
    v_pm_dec_kms = 4.74 * p["pm_dec"] * p["dist"] / 1000.0
    v_transverse_kms = float(np.sqrt(v_pm_ra_kms**2 + v_pm_dec_kms**2))
    dist_pc = float(p["dist"])
    
    # --- Variance Decomposition ---
    # Wrapped-normal variance decomposition: sigma_circ^2 = sigma_signal^2 + k*(D/v)
    # k is the noise-variance scaling coefficient from D/v sampling
    sigma_circ_sq = psi_circ_std**2
    dv_ratio = dist_pc / v_transverse_kms if v_transverse_kms > 0 else 0.0
    
    # For single-pulsar fit, we cannot uniquely determine sigma_signal and k
    # Store D/v ratio for cross-pulsar comparison
    # (Full decomposition requires two pulsars with known distances/velocities)
    
    # --- Noise Floor & Excess ---
    # E[|X|] = sigma * sqrt(2/pi) for zero-mean Gaussian scatter. Per epoch, sigma is
    # a robust MAD scale of signed geometric_delta_us about the epoch median (less
    # inflated than std when holonomy varies across triplets).
    epoch_noise_floors_us = []
    for r in all_results:
        triplets = r["triplets"]
        # Same filter as the main aggregation loop (>=5 triplets for reliable stats)
        if len(triplets) < 5:
            continue
        ep_deltas = np.array([t["geometric_delta_us"] for t in triplets], dtype=float)
        med = float(np.median(ep_deltas))
        mad = float(np.median(np.abs(ep_deltas - med)))
        sigma_us = 1.4826 * mad
        if (not np.isfinite(sigma_us)) or sigma_us <= 0.0:
            sigma_us = float(np.std(ep_deltas, ddof=1))
        epoch_noise_floors_us.append(sigma_us * np.sqrt(2.0 / np.pi))

    epoch_noise_floors_us = np.array(epoch_noise_floors_us)
    if len(epoch_noise_floors_us) == 0 or np.any(np.isnan(epoch_noise_floors_us)):
        print_status(
            f"[WARN] Insufficient epoch statistics for noise-floor estimation "
            f"for pulsar '{pulsar_name}' (no epoch with >=5 triplets). "
            f"Setting noise floor to NaN and excess to 0.0.",
            "WARNING",
        )
        H_noise_bias_ns = float("nan")
        H_excess_ns = 0.0
        H_excess_t = 0.0
    else:
        # Weighted mean noise floor (us -> ns), using same weights as H_mean_us
        H_noise_bias_us = float(np.sum(epoch_weights * epoch_noise_floors_us) / sum_w)
        H_noise_bias_ns = float(H_noise_bias_us * 1e3)

        # Noise-subtracted excess: linear subtraction for transparency.
        # Under pure noise, H_mean_ns ≈ H_noise_bias_ns and excess → 0.
        # A positive excess indicates signal above the Rice floor.
        H_excess_ns = float(max(0.0, H_mean_us * 1e3 - H_noise_bias_ns))
        H_excess_t = H_excess_ns / (H_sem_us * 1e3) if H_excess_ns > 0 and H_sem_us > 0 else 0.0

    # Detection status based on proper circular statistics
    # The legacy linear t-statistic (psi_t) is inappropriate for circular data
    # and yields false positives for highly dispersed distributions (e.g. J1603).
    # We use Rayleigh/V-test significance with bootstrap CI confirmation.
    ci_excludes_zero = bool(not (psi_boot_ci[0] <= 0 <= psi_boot_ci[1])) if np.isfinite(psi_boot_ci[0]) and np.isfinite(psi_boot_ci[1]) else False
    detected_3sigma = bool((rayleigh_p < 0.003 or v_p < 0.003) and ci_excludes_zero)
    detected_5sigma = bool((rayleigh_p < 1e-6 or v_p < 1e-6) and ci_excludes_zero)

    # Unweighted mean for transparency (weighted mean uses inverse-variance weighting)
    H_unweighted_mean_us = float(np.mean(epoch_abs_delta_means))
    
    summary = {
        "phase_delta_ns_formula": PHASE_DELTA_NS_FORMULA,
        "phase_delta_ns_definition": PHASE_DELTA_NS_DEFINITION,
        "pulsar": pulsar_name,
        "n_epochs": viable,
        "n_total_triplets": n_total_triplets,
        "n_independent_samples": n_independent,
        "n_excluded_low_triplet_count": viable - n_independent,
        "mean_geometric_closure_us": float(H_mean_us),
        "std_geometric_closure_us": float(np.std(epoch_delta_means)),
        "sem_geometric_closure_us": float(H_sem_us),
        "H_magnitude_ns": float(H_mean_us * 1e3),
        "H_unweighted_mean_ns": float(H_unweighted_mean_us * 1e3),
        "H_sem_ns": float(H_sem_us * 1e3),
        "H_sem_between_epoch_unweighted_ns": float(H_sem_between_epoch_us * 1e3),
        "H_t_statistic": float(H_t),
        "H_p_value": float(H_p),
        "H_noise_bias_ns": H_noise_bias_ns,
        "H_excess_ns": H_excess_ns,
        "H_excess_t_statistic": float(H_excess_t),
        "H_signed_mean_ns": float(H_signed_mean_us * 1e3),
        "H_signed_sem_ns": float(H_signed_sem_us * 1e3),
        "H_signed_t_statistic": float(H_signed_t),
        "H_trim_magnitude_ns": float(H_trim_mean * 1e3),
        "H_trim_sem_ns": float(H_trim_sem * 1e3),
        "H_trim_t_statistic": float(H_trim_t),
        "H_trim_fraction": float(trim_frac),
        "H_noise_floor_method": "mad_median_folded_normal_E_abs",
        "phase_closure_epoch_mean_weighting": "triplet_snr_squared_circular_mean",
        "phase_closure_mean_rad": float(psi_mean),
        "phase_closure_t_statistic": float(psi_t),
        "phase_closure_rbar": float(psi_rbar),
        "phase_closure_circ_std_rad": float(psi_circ_std),
        "phase_closure_circ_se_rad": float(psi_circ_se),
        "phase_closure_rayleigh_z": float(rayleigh_z),
        "phase_closure_rayleigh_p": float(rayleigh_p),
        "phase_closure_v_stat": float(v_stat),
        "phase_closure_v_p": float(v_p),
        "phase_closure_bootstrap_se_rad": float(psi_boot_se) if np.isfinite(psi_boot_se) else None,
        "phase_closure_bootstrap_ci_95_lower_rad": float(psi_boot_ci[0]) if np.isfinite(psi_boot_ci[0]) else None,
        "phase_closure_bootstrap_ci_95_upper_rad": float(psi_boot_ci[1]) if np.isfinite(psi_boot_ci[1]) else None,
        "phase_closure_circ_var_rad2": float(sigma_circ_sq),
        "phase_closure_dv_ratio_pc_per_kms": float(dv_ratio),
        "phase_closure_monopole_rad": float(psi_monopole),
        "phase_closure_bipole_rad": float(psi_bipole),
        "phase_closure_bipole_ratio": float(bipole_ratio) if np.isfinite(bipole_ratio) else None,
        "phase_closure_angular_sep_deg": float(np.degrees(angular_sep)),
        "phase_closure_pos_orientation_n": pos_n,
        "phase_closure_pos_orientation_psi_rad": float(pos_angle),
        "phase_closure_pos_orientation_rbar": float(pos_rbar),
        "phase_closure_neg_orientation_n": neg_n,
        "phase_closure_neg_orientation_psi_rad": float(neg_angle),
        "phase_closure_neg_orientation_rbar": float(neg_rbar),
        "pulsar_v_transverse_kms": v_transverse_kms,
        "pulsar_dist_pc": dist_pc,
        "detected_3sigma": detected_3sigma,
        "detected_5sigma": detected_5sigma,
        "method": "Conservative Independence Framework (Epoch-level aggregation)",
        "note": "Phase Closure uses circular statistics: circular mean, circular SE, Rayleigh/V-tests, and epoch-level bootstrap CI. Legacy t-statistic retained for backward compatibility. n_independent_samples excludes epochs with <5 triplets (unreliable statistics). Within each epoch, triplet ψ uses SNR^2-weighted circular mean. H_magnitude_ns is inverse-variance weighted across epochs; H_sem_ns matches that weighting (1/sqrt(sum w_i)); H_sem_between_epoch_unweighted_ns is the simple between-epoch SEM of mean(|H|) for comparison. H_t and H_p use the weighted SEM with a Gaussian two-sided tail (not Student t on n_epochs). Triplet cap applies equally with or without SNR filtering. H_noise_bias_ns uses per-epoch MAD scale of signed geometric delays about the epoch median, then E|H|_floor = sigma*sqrt(2/pi). Binary reflex velocity uses Thiele–Innes projection (Reardon et al. 2024 PPTA-DR3 e, omega, i, Omega, T0, x).",
    }


    # Scientific output
    print_status("\n" + "=" * 70, "TITLE")
    print_status(f"GEOMETRIC CLOSURE DELAY ANALYSIS: {pulsar_name}", "TITLE")
    print_status("=" * 70, "TITLE")
    print_status(
        f"  Epochs analyzed:          {viable}/{len(sec_files)} viable", "INFO"
    )
    print_status(f"  Total triplets:         {n_total_triplets:,}", "INFO")
    print_status(f"  ", "INFO")
    print_status(f"  ── TEP Magnitude |H| (Mean of |geometric_delta|) ---", "INFO")
    print_status(
        f"  Total |H| Mean:    {H_mean_us * 1e3:.3f} +/- {H_sem_us * 1e3:.3f} ns (mean/weighted-SEM ratio = {H_t:.1f} sigma)", "INFO"
    )
    print_status(
        f"  Unweighted |H| Mean: {H_unweighted_mean_us * 1e3:.3f} ns", "INFO"
    )
    # Rice noise floor: expected |H| under zero-mean Gaussian noise (folded-normal)
    detection_label = "SUCCESS" if detected_3sigma else "INFO"
    print_status(f"  Noise floor E[|H|]:      {H_noise_bias_ns:.3f} ns (from epoch-level sigma)", "INFO")
    print_status(f"  Noise-subtracted excess:  {H_excess_ns:.3f} ns (t = {H_excess_t:.1f} sigma)", detection_label)
    print_status(f"  ", "INFO")
    print_status(f"  -- Signed Mean (Diagnostic: shows bipolar cancellation)", "INFO")
    print_status(
        f"  Signed mean:              {H_signed_mean_us * 1e3:+.3f} +/- {H_signed_sem_us * 1e3:.3f} ns",
        "INFO",
    )
    print_status(f"  Signed t-statistic:       {H_signed_t:.2f} sigma", "INFO")
    print_status(f"  ", "INFO")
    print_status(f"  -- Robust Estimator (10% trimmed mean of |H|) ------", "INFO")
    print_status(
        f"  |H|_trimmed:              {H_trim_mean * 1e3:.3f} +/- {H_trim_sem * 1e3:.3f} ns",
        "INFO",
    )
    print_status(f"  Trimmed t-statistic:      {H_trim_t:.2f} sigma", "INFO")
    print_status(f"  ", "INFO")
    print_status(f"  -- Phase Closure (Circular Statistics) -------------", "INFO")
    print_status(
        f"  Circular mean psi:        {psi_mean:+.4f} rad", "INFO"
    )
    print_status(
        f"  Circular SE:              {psi_circ_se:.4f} rad (R_bar = {psi_rbar:.3f})", "INFO"
    )
    print_status(
        f"  Bootstrap SE:             {psi_boot_se:.4f} rad", "INFO"
    )
    print_status(
        f"  95% Bootstrap CI:         [{psi_boot_ci[0]:+.4f}, {psi_boot_ci[1]:+.4f}] rad", "INFO"
    )
    print_status(
        f"  Rayleigh Z:               {rayleigh_z:.2f} (p = {rayleigh_p:.2e})", "INFO"
    )
    print_status(
        f"  V-test (mu0=0):           {v_stat:+.2f} (p = {v_p:.2e})", "INFO"
    )
    print_status(f"  Legacy t-statistic:       {psi_t:.2f} sigma", "INFO")

    if detected_5sigma:
        print_status(
            f"\n  High-significance TEP detection (Phase Closure): Rayleigh p = {rayleigh_p:.2e}, V-test p = {v_p:.2e}",
            "SUCCESS",
        )
    elif detected_3sigma:
        print_status(
            f"\n  Moderate-significance TEP detection (Phase Closure): Rayleigh p = {rayleigh_p:.2e}, V-test p = {v_p:.2e}; further confirmation warranted",
            "WARNING",
        )
    else:
        print_status(
            f"\n  Note: No significant phase closure coherence (Rayleigh p = {rayleigh_p:.2e}, V-test p = {v_p:.2e}) (Delay Excess: {H_excess_t:.1f} sigma)",
            "INFO",
        )

    print_status(f"\n" + "=" * 70, "TITLE")

    return all_results, summary


def main(
    max_arclets: int = MAX_ARCLETS_FOR_TRIPLETS,
    use_snr_threshold: bool = True,
    min_snr: float = 5.0,
    mode: str = "full",
    pulsars: list = ["all"]
):
    # Logger is set by run_pipeline.py via set_step_logger()
    # Do not create a new logger here to avoid overriding the pipeline's logger

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    configure_blas_thread_env()
    workers = worker_count(role="cpu_bound", reserve=2)

    print_status("=" * 60, "TITLE")
    print_status("TEP-J0437 Step 003: Closure Delay Extraction", "TITLE")
    print_status("=" * 60, "TITLE")
    print_status(
        f"Method: Sub-pixel cross-term measurement Parallel: {workers} workers",
        "INFO",
    )
    print_status(
        f"Parameters: max_arclets={max_arclets}, use_snr_threshold={use_snr_threshold}, min_snr={min_snr}",
        "INFO",
    )

    # Process J0437 in different modes (full, sb0, sb1)
    modes = ["full", "sb0", "sb1"]
    j0437_summary_full = None  # Save full-mode summary for cross-pulsar comparison
    for mode in modes:
        print_status(f"\nProcessing J0437-4715 (Mode: {mode})...", "INFO")
        j0437_results, j0437_summary = process_pulsar(
            SEC_DIR,
            "J0437-4715",
            max_arclets,
            use_snr_threshold,
            min_snr,
            workers=workers,
            mode=mode,
        )
        if j0437_results is not None:
            suffix = f"_{mode}" if mode != "full" else ""
            print_status(f"J0437 {mode}: Extracted {len(j0437_results)} epoch results", "SUCCESS")
            
            # Save results
            with open(RESULTS_DIR / f"step_003_closure_final_per_epoch_j0437{suffix}.json", "w") as f:
                json.dump(j0437_results, f, indent=2, cls=NpEncoder)
            with open(RESULTS_DIR / f"step_003_closure_final_summary_j0437{suffix}.json", "w") as f:
                json.dump(j0437_summary, f, indent=2, cls=NpEncoder)
                
            if mode == "full":
                # Legacy compatibility
                with open(RESULTS_DIR / "step_003_closure_final_per_epoch.json", "w") as f:
                    json.dump(j0437_results, f, indent=2, cls=NpEncoder)
                with open(RESULTS_DIR / "step_003_closure_final_summary.json", "w") as f:
                    json.dump(j0437_summary, f, indent=2, cls=NpEncoder)
                j0437_summary_full = j0437_summary
        else:
            print_status(f"J0437 {mode}: No valid results extracted", "WARNING")

    # Process J1603
    print_status("\nProcessing J1603-7202...", "INFO")
    j1603_results, j1603_summary = process_pulsar(
        SEC_DIR_J1603,
        "J1603-7202",
        max_arclets,
        use_snr_threshold,
        min_snr,
        workers=workers,
    )
    if j1603_results is not None:
        print_status(f"J1603: Extracted {len(j1603_results)} epoch results", "SUCCESS")
        print_status(
            f"J1603 Summary: |H| = {j1603_summary['H_magnitude_ns']:.3f} +/- {j1603_summary['H_sem_ns']:.3f} ns, t = {j1603_summary['H_t_statistic']:.2f}sigma",
            "INFO",
        )
        # Write with step_003_ prefix (new naming convention)
        with open(
            RESULTS_DIR / "step_003_closure_final_per_epoch_j1603.json", "w"
        ) as f:
            json.dump(j1603_results, f, indent=2, cls=NpEncoder)
        with open(RESULTS_DIR / "step_003_closure_final_summary_j1603.json", "w") as f:
            json.dump(j1603_summary, f, indent=2, cls=NpEncoder)
    else:
        print_status("J1603: No valid results extracted", "WARNING")

    # Process Jiamusi pulsars (B1933+16 has most epochs)
    if SEC_DIR_JIAMUSI.exists():
        print_status("\n" + "=" * 70, "TITLE")
        print_status("Jiamusi Pulsar Analysis", "TITLE")
        print_status("=" * 70, "TITLE")

        # Process B1933+16 (8 epochs - best candidate)
        print_status("\nProcessing B1933+16...", "INFO")
        b1933_results, b1933_summary = process_pulsar(
            SEC_DIR_JIAMUSI,
            "B1933+16",
            max_arclets,
            use_snr_threshold,
            3.0,
            workers=workers,
        )
        if b1933_results is not None:
            # Write with step_003_ prefix (new naming convention)
            with open(
                RESULTS_DIR / "step_003_closure_final_per_epoch_B1933.json", "w"
            ) as f:
                json.dump(b1933_results, f, indent=2, cls=NpEncoder)
            with open(
                RESULTS_DIR / "step_003_closure_final_summary_B1933.json", "w"
            ) as f:
                json.dump(b1933_summary, f, indent=2, cls=NpEncoder)
            print_status(
                f"B1933+16: {b1933_summary.get('n_epochs', 0)} viable epochs", "INFO"
            )

        # Process B0355+54 (5 epochs)
        print_status("\nProcessing B0355+54...", "INFO")
        b0355_results, b0355_summary = process_pulsar(
            SEC_DIR_JIAMUSI,
            "B0355+54",
            max_arclets,
            use_snr_threshold,
            3.0,
            workers=workers,
        )
        if b0355_results is not None:
            # Write with step_003_ prefix (new naming convention)
            with open(
                RESULTS_DIR / "step_003_closure_final_per_epoch_B0355.json", "w"
            ) as f:
                json.dump(b0355_results, f, indent=2, cls=NpEncoder)
            with open(
                RESULTS_DIR / "step_003_closure_final_summary_B0355.json", "w"
            ) as f:
                json.dump(b0355_summary, f, indent=2, cls=NpEncoder)
            print_status(
                f"B0355+54: {b0355_summary.get('n_epochs', 0)} viable epochs", "INFO"
            )

        # Process remaining Jiamusi pulsars
        jiamusi_pulsars = [
            "B0329+54",
            "B0540+23",
            "B0740-28",
            "B1508+55",
            "B2154+40",
            "B2310+42",
            "B2324+60",
            "B2351+61",
        ]

        print_status(
            f"\nProcessing {len(jiamusi_pulsars)} additional Jiamusi pulsars...", "INFO"
        )
        for pulsar_name in jiamusi_pulsars:
            pulsar_prefix = pulsar_name.split("+")[0].split("-")[0]
            print_status(f"  Processing {pulsar_name}...", "INFO")
            results, summary = process_pulsar(
                SEC_DIR_JIAMUSI,
                pulsar_name,
                max_arclets,
                use_snr_threshold,
                3.0,
                workers=workers,
            )
            if results is not None:
                # Write with step_003_ prefix (new naming convention)
                with open(
                    RESULTS_DIR
                    / f"step_003_closure_final_per_epoch_{pulsar_prefix}.json",
                    "w",
                ) as f:
                    json.dump(results, f, indent=2, cls=NpEncoder)
                with open(
                    RESULTS_DIR
                    / f"step_003_closure_final_summary_{pulsar_prefix}.json",
                    "w",
                ) as f:
                    json.dump(summary, f, indent=2, cls=NpEncoder)
                print_status(
                    f"{pulsar_name}: {summary.get('n_epochs', 0)} viable epochs", "INFO"
                )
            else:
                # Only print warning for primary target pulsars, not control pulsars
                if pulsar_name in ["J0437-4715", "J1603-7202"]:
                    print_status(f"{pulsar_name}: No valid results", "WARNING")

        # Process MeerKAT pulsars (L-band, different frequency from Parkes)
        if SEC_DIR_MEERKAT.exists():
            print_status("\n" + "=" * 70, "TITLE")
            print_status("MeerKAT Pulsar Analysis", "TITLE")
            print_status("=" * 70, "TITLE")

            meerkat_pulsars = [
                "J0908-1739",
                "J0922-0638",
                "J1731-4744",
            ]
            for pulsar_name in meerkat_pulsars:
                pulsar_prefix = pulsar_name.split("+")[0].split("-")[0]
                print_status(f"  Processing {pulsar_name}...", "INFO")
                results, summary = process_pulsar(
                    SEC_DIR_MEERKAT,
                    pulsar_name,
                    max_arclets,
                    use_snr_threshold,
                    2.0,
                    workers=workers,
                )
                if results is not None:
                    with open(
                        RESULTS_DIR
                        / f"step_003_closure_final_per_epoch_{pulsar_prefix}.json",
                        "w",
                    ) as f:
                        json.dump(results, f, indent=2, cls=NpEncoder)
                    with open(
                        RESULTS_DIR
                        / f"step_003_closure_final_summary_{pulsar_prefix}.json",
                        "w",
                    ) as f:
                        json.dump(summary, f, indent=2, cls=NpEncoder)
                    print_status(
                        f"{pulsar_name}: {summary.get('n_epochs', 0)} viable epochs", "INFO"
                    )
                else:
                    print_status(f"{pulsar_name}: No valid results", "INFO")

        print_status("\n" + "=" * 70, "TITLE")
        print_status("STEP 003 COMPLETED SUCCESSFULLY", "SUCCESS")
        print_status("=" * 70, "TITLE")
    else:
        print_status("Jiamusi data directory not found", "WARNING")

    # Combined summary for comparison
    if j0437_summary_full is not None and j1603_summary is not None:
        # Use full-mode J0437 summary for cross-pulsar comparison
        j0437_summary = j0437_summary_full
        
        # Calculate significance of difference in EXCESS holonomy
        h_diff = j0437_summary["H_excess_ns"] - j1603_summary["H_excess_ns"]
        h_diff_se = np.sqrt(j0437_summary["H_sem_ns"]**2 + j1603_summary["H_sem_ns"]**2)
        H_excess_diff_t = abs(h_diff) / h_diff_se

        print_status("\n" + "=" * 70, "TITLE")
        print_status("CONTROL PULSAR COMPARISON: J0437-4715 vs J1603-7202", "TITLE")
        print_status("=" * 70, "TITLE")
        print_status(
            f"  J0437-4715 Excess: {j0437_summary['H_excess_t_statistic']:.2f}sigma, H_ex = {j0437_summary['H_excess_ns']:.3f} ns",
            "INFO",
        )
        print_status(
            f"  J1603-7202 Excess: {j1603_summary['H_excess_t_statistic']:.2f}sigma, H_ex = {j1603_summary['H_excess_ns']:.3f} ns",
            "INFO",
        )
        print_status(f"  Difference Significance: {H_excess_diff_t:.2f}sigma", "INFO")
        print_status(f"  ", "INFO")
        if (
            H_excess_diff_t > 3.0
        ):
            print_status(
                f"  Pipeline specificity confirmed: Signal difference significant ({H_excess_diff_t:.1f}sigma)",
                "SUCCESS",
            )
        else:
            print_status(
                f"  Pipeline specificity inconclusive: Both pulsars show similar noise-bias levels",
                "WARNING",
            )
        print_status("=" * 70, "TITLE")

        # --- Cross-Pulsar Variance Decomposition ---
        # Model: sigma_circ^2 = sigma_signal^2 + k * (D/v)
        # Two pulsars allow unique determination of both parameters
        dv1 = j0437_summary["phase_closure_dv_ratio_pc_per_kms"]
        dv2 = j1603_summary["phase_closure_dv_ratio_pc_per_kms"]
        var1 = j0437_summary["phase_closure_circ_var_rad2"]
        var2 = j1603_summary["phase_closure_circ_var_rad2"]
        
        # Linear fit: var = a + b * dv  =>  a = sigma_signal^2, b = k
        b = (var2 - var1) / (dv2 - dv1)
        a = var1 - b * dv1
        
        sigma_signal = float(np.sqrt(max(0, a)))
        k_noise = float(b)
        
        # Model-predicted values
        var1_pred = a + b * dv1
        var2_pred = a + b * dv2
        resid1 = var1 - var1_pred
        resid2 = var2 - var2_pred
        
        # Variance explained by model
        var_mean = (var1 + var2) / 2.0
        ss_res = resid1**2 + resid2**2
        ss_tot = (var1 - var_mean)**2 + (var2 - var_mean)**2
        r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
        
        print_status("\n" + "=" * 70, "TITLE")
        print_status("CROSS-PULSAR VARIANCE DECOMPOSITION", "TITLE")
        print_status("=" * 70, "TITLE")
        print_status(f"  Model: sigma_circ^2 = sigma_signal^2 + k * (D/v)", "INFO")
        print_status(f"  J0437: D/v = {dv1:.2f}, sigma_circ^2 = {var1:.3f} rad^2", "INFO")
        print_status(f"  J1603: D/v = {dv2:.2f}, sigma_circ^2 = {var2:.3f} rad^2", "INFO")
        print_status(f"  ", "INFO")
        print_status(f"  Fit: sigma_signal = {sigma_signal:.3f} rad ({np.degrees(sigma_signal):.1f}°)", "INFO")
        print_status(f"  Fit: k = {k_noise:.3f} rad^2 / (pc / km/s)", "INFO")
        print_status(f"  R-squared = {r_squared:.3f}", "INFO")
        print_status(f"  ", "INFO")
        print_status(
            f"  Intrinsic holonomy coherence: {sigma_signal:.3f} rad ({np.degrees(sigma_signal):.1f}°)",
            "INFO",
        )
        if sigma_signal < 0.15:
            print_status(
                f"  Coherence < 0.15 rad: Phases are nearly perfectly coherent (intrinsic dispersion negligible)",
                "SUCCESS",
            )
        print_status("=" * 70, "TITLE")
        
        # Store cross-pulsar results in both summaries
        cross_pulsar_result = {
            "model": "sigma_circ^2 = sigma_signal^2 + k * (D/v)",
            "sigma_signal_rad": sigma_signal,
            "sigma_signal_deg": float(np.degrees(sigma_signal)),
            "k_noise_rad2_per_dv": k_noise,
            "r_squared": r_squared,
            "J0437_dv": dv1,
            "J0437_var_observed": var1,
            "J0437_var_predicted": float(var1_pred),
            "J1603_dv": dv2,
            "J1603_var_observed": var2,
            "J1603_var_predicted": float(var2_pred),
            "n_pulsars": 2,
            "note": "Two-pulsar linear fit determines intrinsic signal coherence and noise scaling coefficient. With only n=2 pulsars, the fit has zero degrees of freedom and R^2=1.0 is mathematically trivial, not evidence of model validity. Additional pulsars are required for empirical validation of the variance decomposition model.",
        }
        
        j0437_summary["cross_pulsar_variance_decomposition"] = cross_pulsar_result
        j1603_summary["cross_pulsar_variance_decomposition"] = cross_pulsar_result
        
        # Rewrite updated summaries
        with open(RESULTS_DIR / "step_003_closure_final_summary.json", "w") as f:
            json.dump(j0437_summary, f, indent=2, cls=NpEncoder)
        with open(RESULTS_DIR / "step_003_closure_final_summary_j1603.json", "w") as f:
            json.dump(j1603_summary, f, indent=2, cls=NpEncoder)

    print_status("Analysis complete. Results saved to results/", "SUCCESS")
    print_status("=" * 70, "TITLE")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute closure delays from secondary spectra"
    )
    parser.add_argument(
        "--max-arclets",
        type=int,
        default=MAX_ARCLETS_FOR_TRIPLETS,
        help="Maximum arclets to consider",
    )
    parser.add_argument(
        "--no-snr-threshold",
        action="store_true",
        help="Disable SNR threshold (use original max triplets method)",
    )
    parser.add_argument(
        "--min-snr",
        type=float,
        default=5.0,
        help="Minimum SNR threshold for triplet selection (5-sigma, particle physics standard)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Number of parallel workers (default: auto-detect)",
    )
    parser.add_argument(
        "--mode",
        type=str,
        default="full",
        choices=["full", "sb0", "sb1"],
        help="Processing mode (full band or sub-bands)",
    )
    parser.add_argument(
        "--pulsar",
        type=str,
        default="all",
        help="Specific pulsar to process (or 'all')",
    )
    args = parser.parse_args()

    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger("step_003_final", str(log_dir / "step_003_final.log"))
    set_step_logger(logger)

    main(
        max_arclets=args.max_arclets,
        use_snr_threshold=not args.no_snr_threshold,
        min_snr=args.min_snr,
        mode=args.mode,
        pulsars=[args.pulsar]
    )
