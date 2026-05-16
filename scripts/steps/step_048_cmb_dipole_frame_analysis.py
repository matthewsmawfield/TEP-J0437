#!/usr/bin/env python3
"""
================================================================================
STEP 048: CMB DIPOLE FRAME ANALYSIS
================================================================================

Purpose: Test whether Stokes-weighted closure statistics are stable when a
fixed CMB kinematic dipole bulk vector is added to the Step 003 effective
velocity, and whether phase-only (unweighted) Phase Closure is invariant as
required.

Step 003 uses solar-system barycentric Earth velocity (ICRS) combined with
pulsar proper motion and binary terms in a 2D equatorial (vx, vy) reduction.
This step does not perform a full Lorentz-frame transformation of the
scattering problem; it applies the Planck 2018 dipole vector in the same
2D convention as Step 003 for sensitivity and null-structure tests.

The Solar System barycentre moves at ~369.8 km/s toward Galactic
(l ~ 264.02 deg, b ~ 48.25 deg) relative to the CMB rest frame (Planck 2018).
That bulk speed is large compared to Earth's orbital speed (~30 km/s) and
comparable to or larger than some pulsars' transverse speeds, so slow
pulsars can show large changes in v_proj and Stokes weights.

Methodology:
------------
1. Build the CMB dipole velocity in ICRS Cartesian coordinates (Planck 2018).
2. For each epoch, v_eff_ssb = Step 003 effective velocity; sensitivity model
   v_eff_cmb_xy = v_eff_ssb_xy + v_cmb_xy (same 2D plane as Step 003).
3. Recompute geometric_delta_us = delta_us * geom_sign * (v_proj/50 km/s).
4. Recompute epoch-aggregated |H|, signed means, weighted and unweighted psi.
5. Enhanced tests: Earth velocity projected onto the *unit* CMB dipole in 3D
   ICRS for annual-modulation covariates; wrong-direction controls; joint
   linear models; multipole check; random-direction |H| specificity.

Triplet phase_closure_rad is unchanged by velocity rescaling; unweighted
circular summaries of epoch-mean psi from those phases are invariant under
bulk-vector swaps. Weighted psi can move because epoch inverse-variance
weights depend on Stokes-weighted delay scatter.
================================================================================
"""

import argparse
import json
import math
import sys
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
    J0437_RA_RAD,
    J0437_T0_MJD,
)
from scripts.utils.json_numpy import NpEncoder
from scripts.utils.logger import print_status, set_step_logger

from astropy.coordinates import Galactic, ICRS, UnitSphericalRepresentation
from astropy import units as u

# Import existing kinematic database and velocity calculator from Step 003
from scripts.steps.step_003_closure_delays_final import (
    PULSAR_PARAMS,
    circular_mean_and_rbar,
    rayleigh_test,
    v_test_circular,
    circular_bootstrap_ci,
)

RESULTS_DIR = PROJECT_ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


# =============================================================================
# Pulsar sky positions for tangent-plane kinematics
# =============================================================================
_PULSAR_RADECS_ICRS = {
    # ATNF PSRCAT / standard timing ephemerides (ICRS J2000)
    # J0437-4715: RA 04:37:15.9, Dec -47:15:09
    "J0437-4715": ("04:37:15.9", "-47:15:09"),
    # J1603-7202: RA 16:03:35.7, Dec -72:02:32.7
    "J1603-7202": ("16:03:35.7", "-72:02:32.7"),
}


def _east_north_basis_icrs(pulsar_name: str) -> tuple[np.ndarray, np.ndarray]:
    """Return ICRS unit vectors (east, north) at the pulsar sky position.

    East points toward increasing right ascension at fixed declination.
    North points toward increasing declination.
    """
    if pulsar_name not in _PULSAR_RADECS_ICRS:
        raise KeyError(
            f"No ICRS RA/Dec registered for {pulsar_name}. "
            "Add it to _PULSAR_RADECS_ICRS for tangent-plane kinematics."
        )
    ra_str, dec_str = _PULSAR_RADECS_ICRS[pulsar_name]
    from astropy.coordinates import SkyCoord

    c = SkyCoord(ra=ra_str, dec=dec_str, unit=(u.hourangle, u.deg), frame="icrs")
    ra = float(c.ra.to_value(u.rad))
    dec = float(c.dec.to_value(u.rad))

    # Equatorial-cartesian basis in ICRS: x toward RA=0,Dec=0; y toward RA=90°,Dec=0; z toward Dec=+90°.
    east = np.array([-math.sin(ra), math.cos(ra), 0.0], dtype=float)
    north = np.array(
        [-math.cos(ra) * math.sin(dec), -math.sin(ra) * math.sin(dec), math.cos(dec)],
        dtype=float,
    )
    east = east / (np.linalg.norm(east) + 1e-15)
    north = north / (np.linalg.norm(north) + 1e-15)
    return east, north


def _project_to_tangent_en(v_icrs_3d: np.ndarray, pulsar_name: str) -> np.ndarray:
    """Project a 3D ICRS vector into the pulsar tangent plane (east, north)."""
    east, north = _east_north_basis_icrs(pulsar_name)
    v = np.asarray(v_icrs_3d, dtype=float)
    return np.array([float(np.dot(v, east)), float(np.dot(v, north))], dtype=float)


def calculate_velocity_vector_tangent_plane(
    mjd: float,
    pulsar_name: str,
    bulk_icrs_3d: Optional[np.ndarray] = None,
    verbose: bool = False,
) -> np.ndarray:
    """Effective transverse velocity in the pulsar sky tangent plane (east, north).

    This is a more literal kinematic construction than the legacy Step 003
    ICRS (vx, vy) reduction: Earth barycentric velocity and any bulk vector
    are projected into the pulsar's local tangent plane before forming v_eff.
    """
    p = PULSAR_PARAMS[pulsar_name]
    weight_pulsar = (1.0 - p["s"]) / p["s"]

    # Proper motion terms are treated as already in tangent-plane components.
    v_pm_e = 4.74 * p["pm_ra"] * p["dist"] / 1000.0
    v_pm_n = 4.74 * p["pm_dec"] * p["dist"] / 1000.0
    v_pm = np.array([v_pm_e, v_pm_n], dtype=float)

    # Binary reflex velocity: same tangent (east, north) basis as proper motion (Step 003)
    v_orbit = np.array([0.0, 0.0], dtype=float)
    if p["pb"] > 0:
        if pulsar_name == "J0437-4715":
            v_orbit = reflex_binary_transverse_velocity_kms(
                float(mjd),
                float(J0437_PB_DAYS),
                float(J0437_T0_MJD),
                float(J0437_A1_LC),
                float(J0437_ECC),
                float(J0437_OM_DEG),
                float(J0437_INC_DEG),
                float(J0437_KOM_DEG),
                float(C_LIGHT_KM_S),
                float(J0437_RA_RAD),
                float(J0437_DEC_RAD),
            ).astype(float)
        else:
            raise NotImplementedError(
                f"Tangent-plane binary velocity only implemented for J0437-4715; got {pulsar_name!r}."
            )

    # Earth barycentric velocity projected into the pulsar tangent plane
    v_earth_3d = get_earth_barycentric_velocity_kms_icrs(mjd)
    v_earth = _project_to_tangent_en(v_earth_3d, pulsar_name)

    v_eff = weight_pulsar * (v_pm + v_orbit) + v_earth

    if bulk_icrs_3d is not None:
        v_bulk = _project_to_tangent_en(np.asarray(bulk_icrs_3d, dtype=float), pulsar_name)
        v_eff = v_eff + v_bulk

    if verbose:
        print_status(
            f"  Tangent-plane v_eff({pulsar_name}) = {v_eff} km/s (east,north), bulk={'yes' if bulk_icrs_3d is not None else 'no'}",
            "PHYSICS",
        )
    return v_eff


def _ols_coeffs_rss_stderr_t(X: np.ndarray, y: np.ndarray) -> tuple:
    """Ordinary least squares: coefficients, RSS, stderr, t-values.

    Coefficients use ``pinv(X) @ y`` (minimum-norm least squares) so rank-
    deficient or ill-conditioned designs do not yield NaNs from ``lstsq``.
    Standard errors use ``pinv(X.T @ X)`` on the residual variance scale.
    """
    n, p = X.shape
    if not np.isfinite(X).all() or not np.isfinite(y).all():
        raise ValueError(
            "OLS inputs contain non-finite values; check per-epoch summaries and predictors."
        )
    # Suppress benign fp warnings from pinv SVD and downstream matmul on stiff designs.
    with np.errstate(divide="ignore", invalid="ignore", over="ignore", under="ignore"):
        coeffs = np.linalg.pinv(X, rcond=1e-12) @ y
        pred = X @ coeffs
        resid = y - pred
        rss = float(np.dot(resid, resid))
        dof = max(1, n - p)
        mse = rss / dof
        gram_inv = np.linalg.pinv(X.T @ X, rcond=1e-12)
        var = mse * np.clip(np.diag(gram_inv), 0.0, None)
        se = np.sqrt(var)
        tvals = coeffs / (se + 1e-15)
    if not np.isfinite(coeffs).all() or not np.isfinite(tvals).all():
        raise ValueError(
            "OLS produced non-finite coefficients or t-statistics; inspect the design matrix."
        )
    return coeffs, rss, se, tvals


# =============================================================================
# CMB DIPOLE PARAMETERS (Planck 2018)
# =============================================================================
CMB_DIPOLE_V_KMS = 369.82  # km/s
CMB_DIPOLE_L_DEG = 264.0211  # Galactic longitude
CMB_DIPOLE_B_DEG = 48.2533  # Galactic latitude


def get_cmb_dipole_icrs_cartesian() -> np.ndarray:
    """
    Compute the CMB dipole velocity vector in ICRS Cartesian coordinates.

    Full vector in ICRS Cartesian coordinates (km/s). Step 003 uses only
    (vx, vy) when forming effective velocities; vz is retained for Earth–CMB
    projection diagnostics in this step.
    """
    cmb_gal = Galactic(
        l=CMB_DIPOLE_L_DEG * u.deg,
        b=CMB_DIPOLE_B_DEG * u.deg,
    )
    cmb_icrs = cmb_gal.transform_to(ICRS())
    # Unit direction in ICRS
    unit_dir = cmb_icrs.represent_as(UnitSphericalRepresentation)
    vx = CMB_DIPOLE_V_KMS * np.cos(unit_dir.lat.rad) * np.cos(unit_dir.lon.rad)
    vy = CMB_DIPOLE_V_KMS * np.cos(unit_dir.lat.rad) * np.sin(unit_dir.lon.rad)
    vz = CMB_DIPOLE_V_KMS * np.sin(unit_dir.lat.rad)
    # Existing code uses only [vx, vy]; return all three for completeness
    return np.array([vx, vy, vz])


# Cache the CMB dipole vector (ICRS Cartesian, km/s) and unit direction
_CMB_DIPOLE_VELOCITY_ICRS = get_cmb_dipole_icrs_cartesian()
_CMB_DIPOLE_UNIT_ICRS = _CMB_DIPOLE_VELOCITY_ICRS / (
    np.linalg.norm(_CMB_DIPOLE_VELOCITY_ICRS) + 1e-15
)


def calculate_velocity_vector_cmb(
    mjd: float,
    pulsar_name: str = "J0437-4715",
    kinematics_model: str = "legacy_xy",
    verbose: bool = False,
) -> np.ndarray:
    """
    Sensitivity-model effective velocity: Step 003 v_eff (SSB kinematics,
    2D ICRS) plus the CMB dipole (vx, vy) used consistently in Step 003.
    """
    if kinematics_model == "tangent_plane":
        v_eff_cmb = calculate_velocity_vector_tangent_plane(
            mjd,
            pulsar_name=pulsar_name,
            bulk_icrs_3d=_CMB_DIPOLE_VELOCITY_ICRS,
            verbose=verbose,
        )
        v_eff_ssb = calculate_velocity_vector_tangent_plane(
            mjd, pulsar_name=pulsar_name, bulk_icrs_3d=None, verbose=verbose
        )
        v_cmb_xy = v_eff_cmb - v_eff_ssb
    elif kinematics_model == "legacy_xy":
        from scripts.steps.step_003_closure_delays_final import calculate_velocity_vector

        v_eff_ssb = calculate_velocity_vector(
            mjd, pulsar_name=pulsar_name, verbose=verbose
        )
        v_cmb_xy = _CMB_DIPOLE_VELOCITY_ICRS[:2]
        v_eff_cmb = v_eff_ssb + v_cmb_xy
    else:
        raise ValueError(f"Unknown kinematics_model: {kinematics_model}")
    if verbose:
        print_status(
            f"  CMB bulk on SSB v_eff: v_ssb={v_eff_ssb}, v_cmb_xy={v_cmb_xy}, v_sum={v_eff_cmb}",
            "PHYSICS",
        )
    return v_eff_cmb


def velocity_projection(v_eff: np.ndarray, psi_deg: float) -> float:
    """Project effective velocity onto the scattering anisotropy axis."""
    rad_psi = np.radians(psi_deg)
    axis = np.array([np.cos(rad_psi), np.sin(rad_psi)])
    return float(np.dot(v_eff, axis))


def _per_epoch_json_candidates(pulsar_name: str) -> List[Path]:
    """Ordered candidate Step 003 per-epoch JSON paths for a pulsar name.

    J0437 and J1603 are matched by **prefix** (``J0437…``, ``J1603…``) only, so
    unrelated names that happen to contain the digit substring ``0437`` (for
    example ``B0437+…``) are not mis-routed to the Parkes J0437 file.
    """
    norm = pulsar_name.strip()
    out: List[Path] = []
    u = norm.upper().replace(" ", "")
    if u.startswith("J0437"):
        out.append(RESULTS_DIR / "step_003_closure_final_per_epoch.json")
    elif u.startswith("J1603"):
        out.append(RESULTS_DIR / "step_003_closure_final_per_epoch_j1603.json")
    else:
        prefix = norm.split("+")[0].split("-")[0]
        out.append(RESULTS_DIR / f"step_003_closure_final_per_epoch_{prefix}.json")
    seen: set[str] = set()
    unique: List[Path] = []
    for p in out:
        key = str(p.resolve())
        if key not in seen:
            seen.add(key)
            unique.append(p)
    return unique


def load_per_epoch_results(pulsar_name: str) -> Optional[List[Dict]]:
    """Load Step 003 per-epoch closure results for the given pulsar."""
    for f in _per_epoch_json_candidates(pulsar_name):
        if f.exists():
            with open(f) as fh:
                return json.load(fh)
    return None


def list_pulsars_with_step003_per_epoch() -> List[str]:
    """Canonical PULSAR_PARAMS names that have a matching Step 003 per-epoch JSON."""
    found: List[str] = []
    for name in PULSAR_PARAMS:
        if load_per_epoch_results(name) is not None:
            found.append(name)
    return sorted(found)


def recompute_epoch_statistics(
    per_epoch_data: List[Dict],
    pulsar_name: str,
    frame: str = "cmb",
    custom_vel_fn=None,
    kinematics_model: str = "legacy_xy",
    include_per_epoch_series: bool = False,
) -> Dict[str, Any]:
    """
    Recompute all closure statistics for a given pulsar in the specified frame.

    Parameters
    ----------
    per_epoch_data : list
        Step 003 per-epoch results (list of dicts with 'mjd', 'triplets').
    pulsar_name : str
    frame : {"hel", "cmb", "ssb"}
        "hel" / "ssb" use Step 003 velocities only (SSB / ICRS convention).
        "cmb" adds the Planck 2018 dipole (vx, vy) to Step 003 v_eff.
    custom_vel_fn : optional callable(mjd, pulsar_name=..., verbose=...)
        Overrides frame-based velocity choice when provided.
    include_per_epoch_series : bool
        If True, attach long per-epoch arrays to the returned dict (large JSON).

    Returns
    -------
    summary : dict
        Statistic summary aligned with Step 003 closure summaries.
    """
    p = PULSAR_PARAMS[pulsar_name]
    psi_deg = p["psi"]

    # Choose velocity calculator
    if custom_vel_fn is not None:
        vel_fn = custom_vel_fn
    else:
        if frame == "cmb":
            def vel_fn(mjd, pulsar_name=pulsar_name, verbose=False):
                return calculate_velocity_vector_cmb(
                    mjd,
                    pulsar_name=pulsar_name,
                    kinematics_model=kinematics_model,
                    verbose=verbose,
                )
        else:
            if kinematics_model == "tangent_plane":
                def vel_fn(mjd, pulsar_name=pulsar_name, verbose=False):
                    return calculate_velocity_vector_tangent_plane(
                        mjd,
                        pulsar_name=pulsar_name,
                        bulk_icrs_3d=None,
                        verbose=verbose,
                    )
            elif kinematics_model == "legacy_xy":
                from scripts.steps.step_003_closure_delays_final import calculate_velocity_vector

                vel_fn = calculate_velocity_vector
            else:
                raise ValueError(f"Unknown kinematics_model: {kinematics_model}")

    # Per-epoch aggregates
    epoch_psi_means = []
    epoch_abs_delta_means = []
    epoch_delta_means = []
    epoch_weights = []
    epoch_geom_signs = []
    epoch_noise_floors_us = []

    n_total_triplets = 0
    n_epochs_used = 0

    per_epoch_v_projs = []
    per_epoch_delta_means_ns = []
    per_epoch_abs_delta_means_ns = []
    per_epoch_psi_means_rad = []

    for epoch in per_epoch_data:
        triplets = epoch.get("triplets", [])
        if len(triplets) < 5:
            continue

        mjd = float(epoch.get("mjd", 0.0))
        v_eff = vel_fn(mjd, pulsar_name=pulsar_name, verbose=False)
        v_proj = velocity_projection(v_eff, psi_deg)
        v_weight = v_proj / 50.0  # Same characteristic scale as Step 003

        # Recompute geometric_delta_us for each triplet
        deltas_cmb = []
        phases = []
        geom_signs = []

        for t in triplets:
            raw_du = t.get("delta_us")
            if raw_du is None:
                raise ValueError(
                    "recompute_epoch_statistics: triplet missing 'delta_us' "
                    f"(epoch={epoch.get('epoch')!r})"
                )
            delta_us_raw = float(raw_du)
            gs = t.get("geom_sign")
            if gs is None:
                raise ValueError(
                    "recompute_epoch_statistics: triplet missing 'geom_sign' "
                    f"(epoch={epoch.get('epoch')!r})"
                )
            geom_sign = float(gs)
            phase = t.get("phase_closure_rad")
            if phase is None:
                raise ValueError(
                    "recompute_epoch_statistics: triplet missing 'phase_closure_rad' "
                    f"(epoch={epoch.get('epoch')!r})"
                )
            geometric_delta = delta_us_raw * geom_sign * v_weight
            deltas_cmb.append(geometric_delta)
            geom_signs.append(geom_sign)
            phases.append(float(phase))

        if not deltas_cmb:
            continue

        n_total_triplets += len(triplets)
        n_epochs_used += 1

        ep_delta = np.array(deltas_cmb)
        ep_abs_delta = np.abs(ep_delta)
        ep_geom_signs = np.array(geom_signs)

        # Phase Closure (frame-independent): match Step 003 SNR^2-weighted epoch mean
        ep_phases = np.array(phases)
        snrs = np.array([float(t["snr"]) for t in triplets], dtype=float)
        w_psi = np.square(np.maximum(snrs, 1e-6))
        ep_psi_mean, _ = circular_mean_and_rbar(ep_phases, w_psi)
        ep_psi_mean = float(ep_psi_mean)

        ep_delta_mean = float(np.mean(ep_delta))
        ep_abs_delta_mean = float(np.mean(ep_abs_delta))

        # Weight by inverse variance
        ep_delta_std = float(np.std(ep_delta, ddof=1)) if len(ep_delta) > 1 else 1e-3
        ep_sem = ep_delta_std / np.sqrt(len(ep_delta))
        weight = 1.0 / (ep_sem**2 + (0.001) ** 2)

        # Majority geometric orientation
        majority_geom = float(np.sign(np.sum(ep_geom_signs)))

        epoch_psi_means.append(ep_psi_mean)
        epoch_abs_delta_means.append(ep_abs_delta_mean)
        epoch_delta_means.append(ep_delta_mean)
        epoch_weights.append(weight)
        epoch_geom_signs.append(majority_geom)
        med_nf = float(np.median(ep_delta))
        mad_nf = float(np.median(np.abs(ep_delta - med_nf)))
        sigma_nf = 1.4826 * mad_nf
        if (not np.isfinite(sigma_nf)) or sigma_nf <= 0.0:
            sigma_nf = float(np.std(ep_delta, ddof=1))
        epoch_noise_floors_us.append(sigma_nf * np.sqrt(2.0 / np.pi))

        per_epoch_v_projs.append(v_proj)
        per_epoch_delta_means_ns.append(ep_delta_mean * 1e3)
        per_epoch_abs_delta_means_ns.append(ep_abs_delta_mean * 1e3)
        per_epoch_psi_means_rad.append(ep_psi_mean)

    if n_epochs_used == 0:
        return {
            "error": (
                f"No qualifying epochs for {pulsar_name} in {frame} frame "
                "(need >=5 triplets with phase_closure_rad on every triplet)"
            )
        }

    epoch_psi_means = np.array(epoch_psi_means)
    epoch_abs_delta_means = np.array(epoch_abs_delta_means)
    epoch_delta_means = np.array(epoch_delta_means)
    epoch_weights = np.array(epoch_weights)
    epoch_geom_signs = np.array(epoch_geom_signs)

    n_independent = n_epochs_used
    sum_w = np.sum(epoch_weights)

    # --- Phase Closure (Primary) ---
    # Weighted mean (same as Step 003)
    psi_mean, psi_rbar = circular_mean_and_rbar(epoch_psi_means, epoch_weights)
    # Unweighted mean (truly frame-invariant, since individual triplet phases
    # do not depend on velocity; only the delay-domain weights do)
    psi_mean_unw, psi_rbar_unw = circular_mean_and_rbar(epoch_psi_means)

    if psi_rbar > 1e-12:
        psi_circ_std = float(np.sqrt(-2.0 * np.log(psi_rbar)))
        psi_circ_se = float(psi_circ_std / np.sqrt(n_independent))
    else:
        psi_circ_std = float(np.pi)
        psi_circ_se = float(np.pi)

    psi_std_err = np.std(epoch_psi_means, ddof=1) / np.sqrt(n_independent)
    psi_t = psi_mean / psi_std_err if psi_std_err > 0 else 0.0

    rayleigh_z, rayleigh_p = rayleigh_test(epoch_psi_means, epoch_weights)
    v_stat, v_p = v_test_circular(epoch_psi_means, mu0=0.0, weights=epoch_weights)
    # Unweighted Rayleigh/V-test (frame-invariant)
    rayleigh_z_unw, rayleigh_p_unw = rayleigh_test(epoch_psi_means)
    v_stat_unw, v_p_unw = v_test_circular(epoch_psi_means, mu0=0.0)
    psi_boot_se, psi_boot_ci = circular_bootstrap_ci(
        epoch_psi_means, epoch_weights, n_boot=10000
    )
    # Unweighted bootstrap
    psi_boot_se_unw, psi_boot_ci_unw = circular_bootstrap_ci(
        epoch_psi_means, n_boot=10000
    )

    # --- Monopole / Bipole ---
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

    psi_monopole = float(((pos_angle + neg_angle) / 2.0 + np.pi) % (2 * np.pi) - np.pi)
    psi_bipole = float(((pos_angle - neg_angle) / 2.0 + np.pi) % (2 * np.pi) - np.pi)

    angular_sep = float(abs(pos_angle - neg_angle))
    if angular_sep > np.pi:
        angular_sep = 2 * np.pi - angular_sep

    bipole_ratio = (
        float(abs(psi_bipole) / abs(psi_monopole))
        if abs(psi_monopole) > 0.1
        else float("inf")
    )

    # --- |H| Magnitude ---
    H_mean_us = (
        np.sum(epoch_weights * epoch_abs_delta_means) / sum_w if sum_w > 0 else 0.0
    )
    H_sem_us = float(np.sqrt(1.0 / sum_w)) if sum_w > 0 else 0.0
    H_sem_between_epoch_us = (
        float(np.std(epoch_abs_delta_means, ddof=1) / np.sqrt(n_independent))
        if n_independent > 1
        else 0.0
    )
    H_t = H_mean_us / H_sem_us if H_sem_us > 0 else 0.0
    H_p = float(2.0 * stats.norm.sf(abs(H_t))) if sum_w > 0 else 1.0

    # --- Signed Mean ---
    H_signed_mean_us = (
        np.sum(epoch_weights * epoch_delta_means) / sum_w if sum_w > 0 else 0.0
    )
    H_signed_sem_us = H_sem_us
    H_signed_t = (
        H_signed_mean_us / H_signed_sem_us if H_signed_sem_us > 0 else 0.0
    )

    epoch_noise_floors_us = np.array(epoch_noise_floors_us, dtype=float)
    if len(epoch_noise_floors_us) == 0 or np.any(np.isnan(epoch_noise_floors_us)):
        H_noise_bias_ns = float("nan")
        H_excess_ns = 0.0
        H_excess_t = 0.0
    else:
        H_noise_bias_us = float(np.sum(epoch_weights * epoch_noise_floors_us) / sum_w)
        H_noise_bias_ns = float(H_noise_bias_us * 1e3)
        H_excess_ns = float(max(0.0, H_mean_us * 1e3 - H_noise_bias_ns))
        H_excess_t = (
            H_excess_ns / (H_sem_us * 1e3)
            if H_excess_ns > 0 and H_sem_us > 0
            else 0.0
        )

    # --- Kinematics ---
    v_pm_ra_kms = 4.74 * p["pm_ra"] * p["dist"] / 1000.0
    v_pm_dec_kms = 4.74 * p["pm_dec"] * p["dist"] / 1000.0
    v_transverse_kms = float(np.sqrt(v_pm_ra_kms**2 + v_pm_dec_kms**2))
    dist_pc = float(p["dist"])
    sigma_circ_sq = psi_circ_std**2
    dv_ratio = dist_pc / v_transverse_kms if v_transverse_kms > 0 else 0.0

    ci_excludes_zero = (
        bool(not (psi_boot_ci[0] <= 0 <= psi_boot_ci[1]))
        if np.isfinite(psi_boot_ci[0]) and np.isfinite(psi_boot_ci[1])
        else False
    )
    detected_3sigma = bool(
        (rayleigh_p < 0.003 or v_p < 0.003) and ci_excludes_zero
    )
    detected_5sigma = bool(
        (rayleigh_p < 1e-6 or v_p < 1e-6) and ci_excludes_zero
    )

    out: Dict[str, Any] = {
        "pulsar": pulsar_name,
        "frame": frame,
        "n_epochs": len(per_epoch_data),
        "n_independent_samples": n_independent,
        "n_total_triplets": n_total_triplets,
        "mean_geometric_closure_us": float(H_mean_us),
        "sem_geometric_closure_us": float(H_sem_us),
        "H_magnitude_ns": float(H_mean_us * 1e3),
        "H_sem_ns": float(H_sem_us * 1e3),
        "H_sem_between_epoch_unweighted_ns": float(H_sem_between_epoch_us * 1e3),
        "H_t_statistic": float(H_t),
        "H_p_value": float(H_p),
        "H_noise_bias_ns": H_noise_bias_ns,
        "H_noise_floor_method": "mad_median_folded_normal_E_abs",
        "phase_closure_epoch_mean_weighting": "triplet_snr_squared_circular_mean",
        "H_excess_ns": H_excess_ns,
        "H_excess_t_statistic": float(H_excess_t),
        "H_signed_mean_ns": float(H_signed_mean_us * 1e3),
        "H_signed_sem_ns": float(H_signed_sem_us * 1e3),
        "H_signed_t_statistic": float(H_signed_t),
        "phase_closure_mean_rad": float(psi_mean),
        "phase_closure_t_statistic": float(psi_t),
        "phase_closure_rbar": float(psi_rbar),
        "phase_closure_circ_std_rad": float(psi_circ_std),
        "phase_closure_circ_se_rad": float(psi_circ_se),
        "phase_closure_rayleigh_z": float(rayleigh_z),
        "phase_closure_rayleigh_p": float(rayleigh_p),
        "phase_closure_v_stat": float(v_stat),
        "phase_closure_v_p": float(v_p),
        "phase_closure_bootstrap_ci_95_lower_rad": (
            float(psi_boot_ci[0]) if np.isfinite(psi_boot_ci[0]) else None
        ),
        "phase_closure_bootstrap_ci_95_upper_rad": (
            float(psi_boot_ci[1]) if np.isfinite(psi_boot_ci[1]) else None
        ),
        # Unweighted (phase-only) Phase Closure
        "phase_closure_mean_unweighted_rad": float(psi_mean_unw),
        "phase_closure_rbar_unweighted": float(psi_rbar_unw),
        "phase_closure_rayleigh_z_unweighted": float(rayleigh_z_unw),
        "phase_closure_rayleigh_p_unweighted": float(rayleigh_p_unw),
        "phase_closure_v_stat_unweighted": float(v_stat_unw),
        "phase_closure_v_p_unweighted": float(v_p_unw),
        "phase_closure_bootstrap_ci_95_lower_unweighted_rad": (
            float(psi_boot_ci_unw[0]) if np.isfinite(psi_boot_ci_unw[0]) else None
        ),
        "phase_closure_bootstrap_ci_95_upper_unweighted_rad": (
            float(psi_boot_ci_unw[1]) if np.isfinite(psi_boot_ci_unw[1]) else None
        ),
        "phase_closure_circ_var_rad2": float(sigma_circ_sq),
        "phase_closure_dv_ratio_pc_per_kms": float(dv_ratio),
        "phase_closure_monopole_rad": float(psi_monopole),
        "phase_closure_bipole_rad": float(psi_bipole),
        "phase_closure_bipole_ratio": (
            float(bipole_ratio) if np.isfinite(bipole_ratio) else None
        ),
        "phase_closure_angular_sep_deg": float(np.degrees(angular_sep)),
        "phase_closure_pos_orientation_n": pos_n,
        "phase_closure_pos_orientation_psi_rad": float(pos_angle),
        "phase_closure_neg_orientation_n": neg_n,
        "phase_closure_neg_orientation_psi_rad": float(neg_angle),
        "pulsar_v_transverse_kms": v_transverse_kms,
        "pulsar_dist_pc": dist_pc,
        "detected_3sigma": detected_3sigma,
        "detected_5sigma": detected_5sigma,
    }
    if include_per_epoch_series:
        out["per_epoch_v_proj_kms"] = [float(x) for x in per_epoch_v_projs]
        out["per_epoch_delta_mean_ns"] = [float(x) for x in per_epoch_delta_means_ns]
        out["per_epoch_abs_delta_mean_ns"] = [
            float(x) for x in per_epoch_abs_delta_means_ns
        ]
        out["per_epoch_psi_mean_rad"] = [float(x) for x in per_epoch_psi_means_rad]
    return out


def compute_velocity_comparison(
    pulsar_name: str, kinematics_model: str = "legacy_xy"
) -> Dict[str, Any]:
    """
    Compare Step-003 (SSB) effective velocity vs the same vector plus CMB bulk (vx, vy).

    Uses a representative epoch (median list index). Projection is onto the
    pulsar scattering anisotropy axis psi from PULSAR_PARAMS.
    """
    per_epoch = load_per_epoch_results(pulsar_name)
    if not per_epoch:
        return {"error": "No per-epoch data found"}

    # Use the middle epoch as representative
    mid_epoch = per_epoch[len(per_epoch) // 2]
    mjd = float(mid_epoch.get("mjd", 0.0))

    p = PULSAR_PARAMS[pulsar_name]
    psi_deg = p["psi"]

    if kinematics_model == "tangent_plane":
        v_ssb = calculate_velocity_vector_tangent_plane(
            mjd, pulsar_name=pulsar_name, bulk_icrs_3d=None, verbose=False
        )
        v_cmb = calculate_velocity_vector_tangent_plane(
            mjd,
            pulsar_name=pulsar_name,
            bulk_icrs_3d=_CMB_DIPOLE_VELOCITY_ICRS,
            verbose=False,
        )
        v_cmb_added = v_cmb - v_ssb
    elif kinematics_model == "legacy_xy":
        from scripts.steps.step_003_closure_delays_final import calculate_velocity_vector

        v_ssb = calculate_velocity_vector(mjd, pulsar_name=pulsar_name, verbose=False)
        v_cmb = calculate_velocity_vector_cmb(
            mjd,
            pulsar_name=pulsar_name,
            kinematics_model="legacy_xy",
            verbose=False,
        )
        v_cmb_added = _CMB_DIPOLE_VELOCITY_ICRS[:2]
    else:
        raise ValueError(f"Unknown kinematics_model: {kinematics_model}")

    v_proj_ssb = velocity_projection(v_ssb, psi_deg)
    v_proj_cmb = velocity_projection(v_cmb, psi_deg)

    # Angle between Step-003 vector and CMB-offset vector (2D ICRS components)
    dot = np.dot(v_ssb, v_cmb)
    norm_ssb = np.linalg.norm(v_ssb)
    norm_cmb = np.linalg.norm(v_cmb)
    angle_between = (
        float(np.degrees(np.arccos(np.clip(dot / (norm_ssb * norm_cmb), -1, 1))))
        if norm_ssb > 0 and norm_cmb > 0
        else 0.0
    )

    return {
        "pulsar": pulsar_name,
        "mjd": mjd,
        "kinematics_model": kinematics_model,
        "v_eff_ssb_kms": [float(v_ssb[0]), float(v_ssb[1])],
        "v_eff_cmb_kms": [float(v_cmb[0]), float(v_cmb[1])],
        "v_eff_hel_kms": [float(v_ssb[0]), float(v_ssb[1])],
        "v_proj_ssb_kms": float(v_proj_ssb),
        "v_proj_cmb_kms": float(v_proj_cmb),
        "v_proj_hel_kms": float(v_proj_ssb),
        "v_proj_sign_ssb": float(np.sign(v_proj_ssb)),
        "v_proj_sign_cmb": float(np.sign(v_proj_cmb)),
        "v_proj_sign_hel": float(np.sign(v_proj_ssb)),
        "v_cmb_added_kms": [float(v_cmb_added[0]), float(v_cmb_added[1])],
        "angle_between_ssb_and_cmb_bulk_deg": angle_between,
        "angle_between_hel_and_cmb_deg": angle_between,
        "norm_ssb_kms": float(norm_ssb),
        "norm_cmb_kms": float(norm_cmb),
        "norm_hel_kms": float(norm_ssb),
    }


# =============================================================================
# ENHANCED CMB FRAME TESTS (matching rigor of Papers 2, 15, 17)
# =============================================================================

def get_earth_barycentric_velocity_kms_icrs(mjd: float) -> np.ndarray:
    """Earth barycentric velocity in ICRS Cartesian (vx, vy, vz) km/s."""
    from astropy.time import Time
    from astropy.coordinates import get_body_barycentric_posvel

    t = Time(mjd, format="mjd")
    _, v_earth_bary = get_body_barycentric_posvel("earth", t)
    return np.asarray(v_earth_bary.xyz.value, dtype=float) * 149597870.7 / 86400.0


def cmb_velocity_projection_3d(v_3d: np.ndarray, cmb_unit_3d: np.ndarray) -> float:
    """Scalar component of v along the unit CMB dipole direction (km/s)."""
    return float(np.dot(np.asarray(v_3d, dtype=float), np.asarray(cmb_unit_3d, dtype=float)))


def get_perpendicular_cmb_vectors() -> Dict[str, np.ndarray]:
    """Return alternative CMB direction vectors for specificity testing."""
    cmb_vec = _CMB_DIPOLE_VELOCITY_ICRS
    cmb_unit = _CMB_DIPOLE_UNIT_ICRS
    z = np.array([0.0, 0.0, 1.0])
    x = np.array([1.0, 0.0, 0.0])
    perp1_3d = np.cross(cmb_unit, z)
    if np.linalg.norm(perp1_3d) < 1e-6:
        perp1_3d = np.cross(cmb_unit, x)
    perp1_3d = perp1_3d / (np.linalg.norm(perp1_3d) + 1e-15)
    perp1 = perp1_3d * CMB_DIPOLE_V_KMS
    perp2_3d = np.cross(cmb_unit, perp1_3d)
    perp2_3d = perp2_3d / (np.linalg.norm(perp2_3d) + 1e-15)
    perp2 = perp2_3d * CMB_DIPOLE_V_KMS
    return {
        "true_cmb": cmb_vec,
        "anti_cmb": -cmb_vec,
        "perpendicular_1": perp1,
        "perpendicular_2": perp2,
    }


def calculate_velocity_vector_custom_cmb(
    mjd: float, pulsar_name: str, cmb_vec: np.ndarray, kinematics_model: str = "legacy_xy"
) -> np.ndarray:
    """Effective velocity plus an arbitrary bulk vector under a chosen kinematics model."""
    if kinematics_model == "tangent_plane":
        return calculate_velocity_vector_tangent_plane(
            mjd, pulsar_name=pulsar_name, bulk_icrs_3d=cmb_vec, verbose=False
        )
    if kinematics_model == "legacy_xy":
        from scripts.steps.step_003_closure_delays_final import calculate_velocity_vector

        v_eff_ssb = calculate_velocity_vector(mjd, pulsar_name=pulsar_name, verbose=False)
        return v_eff_ssb + cmb_vec[:2]
    raise ValueError(f"Unknown kinematics_model: {kinematics_model}")


def make_custom_vel_fn(
    pulsar_name: str, cmb_vec: np.ndarray, kinematics_model: str = "legacy_xy"
):
    """Factory for velocity functions with custom bulk vectors."""
    def _vel_fn(mjd, pulsar_name=pulsar_name, verbose=False):
        return calculate_velocity_vector_custom_cmb(
            mjd,
            pulsar_name=pulsar_name,
            cmb_vec=cmb_vec,
            kinematics_model=kinematics_model,
        )
    return _vel_fn


def recompute_with_custom_cmb(
    per_epoch_data: List[Dict],
    pulsar_name: str,
    cmb_vec: np.ndarray,
    kinematics_model: str = "legacy_xy",
) -> Dict[str, Any]:
    """Recompute statistics with a custom CMB vector."""
    vel_fn = make_custom_vel_fn(pulsar_name, cmb_vec, kinematics_model=kinematics_model)
    return recompute_epoch_statistics(per_epoch_data, pulsar_name, frame="ssb", custom_vel_fn=vel_fn)


def compute_per_epoch_data(
    per_epoch_data: List[Dict],
    pulsar_name: str,
    cmb_vec: np.ndarray,
    kinematics_model: str = "legacy_xy",
) -> Dict[str, np.ndarray]:
    """
    Per-epoch predictors for regression tests.

    v_parallel is the component of Earth's ICRS barycentric velocity along the
    unit vector of ``cmb_vec`` (3D). v_perp is the norm of the perpendicular
    component in 3D.

    ``orbital_phase_rad`` is a simple annual phase proxy (MJD modulo 365.25),
    not Keplerian true anomaly.
    """
    p = PULSAR_PARAMS[pulsar_name]
    psi_deg = p["psi"]
    cmb_unit = np.asarray(cmb_vec, dtype=float)
    cmb_unit = cmb_unit / (np.linalg.norm(cmb_unit) + 1e-15)

    mjds = []
    v_projs = []
    v_projs_sq = []
    v_parallels = []
    v_perps = []
    orbital_phases = []
    abs_deltas = []
    signed_deltas = []
    psi_means = []

    for epoch in per_epoch_data:
        triplets = epoch.get("triplets", [])
        if len(triplets) < 5:
            continue

        mjd = float(epoch.get("mjd", 0.0))

        # Effective velocity for Stokes weight under the chosen kinematics model
        v_eff = calculate_velocity_vector_custom_cmb(
            mjd,
            pulsar_name=pulsar_name,
            cmb_vec=cmb_vec,
            kinematics_model=kinematics_model,
        )
        v_proj = velocity_projection(v_eff, psi_deg)

        # Earth barycentric velocity: parallel / perpendicular to CMB unit (3D ICRS)
        v_earth_3d = get_earth_barycentric_velocity_kms_icrs(mjd)
        v_parallel = cmb_velocity_projection_3d(v_earth_3d, cmb_unit)
        v_perp_vec = v_earth_3d - v_parallel * cmb_unit

        # Dependent variables
        v_weight = v_proj / 50.0
        deltas = []
        phases = []
        for t in triplets:
            raw_du = t.get("delta_us")
            if raw_du is None:
                raise ValueError(
                    "compute_per_epoch_data: triplet missing 'delta_us' "
                    f"(epoch={epoch.get('epoch')!r})"
                )
            delta_us_raw = float(raw_du)
            gs = t.get("geom_sign")
            if gs is None:
                raise ValueError(
                    "compute_per_epoch_data: triplet missing 'geom_sign' "
                    f"(epoch={epoch.get('epoch')!r})"
                )
            geom_sign = float(gs)
            phase_val = t.get("phase_closure_rad")
            if phase_val is None:
                raise ValueError(
                    "compute_per_epoch_data: triplet missing 'phase_closure_rad' "
                    f"(epoch={epoch.get('epoch')!r})"
                )
            deltas.append(delta_us_raw * geom_sign * v_weight)
            phases.append(float(phase_val))

        mjds.append(mjd)
        v_projs.append(v_proj)
        v_projs_sq.append(v_proj**2)
        v_perps.append(float(np.linalg.norm(v_perp_vec)))
        v_parallels.append(v_parallel)

        # Approximate annual phase (calendar proxy; confounder control)
        phase = 2.0 * np.pi * ((mjd - 51544.5) % 365.25) / 365.25
        orbital_phases.append(phase)

        abs_deltas.append(float(np.mean(np.abs(deltas))))
        signed_deltas.append(float(np.mean(deltas)))

        ep_phases = np.array(phases)
        ep_psi_vec = np.mean(np.exp(1j * ep_phases))
        psi_means.append(float(np.angle(ep_psi_vec)))

    orb = np.array(orbital_phases, dtype=float)
    return {
        "mjd": np.array(mjds, dtype=float),
        "v_proj_kms": np.array(v_projs, dtype=float),
        "v_proj_sq_kms2": np.array(v_projs_sq, dtype=float),
        "v_parallel_kms": np.array(v_parallels, dtype=float),
        "v_perp_kms": np.array(v_perps, dtype=float),
        "orbital_phase_rad": orb,
        "sin_orbital": np.sin(orb),
        "cos_orbital": np.cos(orb),
        "abs_delta_ns": np.array(abs_deltas) * 1e3,
        "signed_delta_ns": np.array(signed_deltas) * 1e3,
        "psi_rad": np.array(psi_means, dtype=float),
        "sin_psi": np.sin(np.array(psi_means, dtype=float)),
        "cos_psi": np.cos(np.array(psi_means, dtype=float)),
    }


def test_annual_cmb_projection_modulation(
    per_epoch_data: List[Dict], pulsar_name: str, kinematics_model: str = "legacy_xy"
) -> Dict[str, Any]:
    """
    Test 1: Annual modulation of Earth velocity along the CMB dipole.

    v_parallel(t) is the component of Earth's ICRS barycentric velocity along
    the Planck 2018 dipole unit vector (3D). It varies annually as the Earth
    orbits the SSB. Epoch |H| (after applying the CMB-augmented Stokes model)
    is regressed on v_parallel and on a simple annual harmonic for
    confounder control.
    """
    cmb_vec = _CMB_DIPOLE_VELOCITY_ICRS
    d = compute_per_epoch_data(
        per_epoch_data, pulsar_name, cmb_vec, kinematics_model=kinematics_model
    )
    n = len(d["v_parallel_kms"])
    if n < 10:
        return {"error": f"Insufficient epochs for regression (n={n})"}

    v_par = d["v_parallel_kms"]

    # |H| vs v_parallel
    slope_h, intercept_h, r_h, p_h, se_h = stats.linregress(v_par, d["abs_delta_ns"])

    # Signed mean vs v_parallel
    slope_s, intercept_s, r_s, p_s, se_s = stats.linregress(v_par, d["signed_delta_ns"])

    # Unweighted psi (frame-invariant) — exclude epochs without phase
    m_psi = np.isfinite(d["psi_rad"])
    n_psi = int(np.count_nonzero(m_psi))
    if n_psi >= 10:
        v_psi = v_par[m_psi]
        slope_sin, _, r_sin, p_sin, _ = stats.linregress(v_psi, d["sin_psi"][m_psi])
        slope_cos, _, r_cos, p_cos, _ = stats.linregress(v_psi, d["cos_psi"][m_psi])
    else:
        slope_sin = slope_cos = float("nan")
        r_sin = r_cos = float("nan")
        p_sin = p_cos = float("nan")

    # Joint model: |H| = b0 + b1 * v_parallel + b2 * sin(orbital) + b3 * cos(orbital)
    X_joint = np.column_stack([np.ones(n), v_par, d["sin_orbital"], d["cos_orbital"]])
    coeffs_joint, _, se_joint, t_joint = _ols_coeffs_rss_stderr_t(X_joint, d["abs_delta_ns"])

    return {
        "n_epochs": n,
        "psi_regression_n_epochs": n_psi,
        "v_parallel_range_kms": [float(np.min(v_par)), float(np.max(v_par))],
        "v_parallel_mean_kms": float(np.mean(v_par)),
        "v_parallel_std_kms": float(np.std(v_par)),
        "H_vs_vparallel_slope_ns_per_kms": float(slope_h),
        "H_vs_vparallel_intercept_ns": float(intercept_h),
        "H_vs_vparallel_r": float(r_h),
        "H_vs_vparallel_p": float(p_h),
        "H_vs_vparallel_se": float(se_h),
        "signed_vs_vparallel_slope_ns_per_kms": float(slope_s),
        "signed_vs_vparallel_intercept_ns": float(intercept_s),
        "signed_vs_vparallel_r": float(r_s),
        "signed_vs_vparallel_p": float(p_s),
        "signed_vs_vparallel_se": float(se_s),
        "psi_sin_vs_vparallel_slope": float(slope_sin),
        "psi_sin_vs_vparallel_p": float(p_sin),
        "psi_cos_vs_vparallel_slope": float(slope_cos),
        "psi_cos_vs_vparallel_p": float(p_cos),
        "joint_model_intercept_ns": float(coeffs_joint[0]),
        "joint_model_vparallel_coeff_ns_per_kms": float(coeffs_joint[1]),
        "joint_model_vparallel_t": float(t_joint[1]),
        "joint_model_sin_coeff": float(coeffs_joint[2]),
        "joint_model_sin_t": float(t_joint[2]),
        "joint_model_cos_coeff": float(coeffs_joint[3]),
        "joint_model_cos_t": float(t_joint[3]),
    }


def test_directional_specificity(
    per_epoch_data: List[Dict], pulsar_name: str, kinematics_model: str = "legacy_xy"
) -> Dict[str, Any]:
    """
    Test 2: Directional specificity.

    Recompute statistics with the CMB dipole rotated to alternative
    directions: anti-CMB, two perpendiculars. The unweighted Phase
    Closure should be invariant; velocity-dependent stats should
    show the expected geometric pattern.
    """
    alt_vecs = get_perpendicular_cmb_vectors()
    results = {}

    for label, cmb_vec in alt_vecs.items():
        summary = recompute_with_custom_cmb(
            per_epoch_data,
            pulsar_name,
            cmb_vec,
            kinematics_model=kinematics_model,
        )
        if summary.get("error"):
            return {
                "error": (
                    f"directional_specificity: recomputation failed for '{label}': "
                    f"{summary['error']}"
                ),
            }
        results[label] = {
            "weighted_psi_rad": summary.get("phase_closure_mean_rad"),
            "weighted_rbar": summary.get("phase_closure_rbar"),
            "weighted_rayleigh_p": summary.get("phase_closure_rayleigh_p"),
            "unweighted_psi_rad": summary.get("phase_closure_mean_unweighted_rad"),
            "unweighted_rbar": summary.get("phase_closure_rbar_unweighted"),
            "unweighted_rayleigh_p": summary.get("phase_closure_rayleigh_p_unweighted"),
            "H_magnitude_ns": summary.get("H_magnitude_ns"),
            "H_signed_mean_ns": summary.get("H_signed_mean_ns"),
        }

    true = results["true_cmb"]
    anti = results["anti_cmb"]
    perp1 = results["perpendicular_1"]
    perp2 = results["perpendicular_2"]

    # Unweighted phase-only mean: must match across bulk-vector directions (same phases)
    ref = float(true["unweighted_psi_rad"])
    uw_deltas = [
        abs(ref - float(anti["unweighted_psi_rad"])),
        abs(ref - float(perp1["unweighted_psi_rad"])),
        abs(ref - float(perp2["unweighted_psi_rad"])),
    ]
    uw_max_delta = max(uw_deltas)
    uw_invariant = uw_max_delta < 1e-9

    H_values = {k: v["H_magnitude_ns"] for k, v in results.items()}

    return {
        "alternatives": results,
        "unweighted_psi_invariant_all_directions": bool(uw_invariant),
        "unweighted_psi_max_abs_delta_rad": float(uw_max_delta),
        "H_magnitudes_by_direction_ns": H_values,
    }


def test_joint_model(
    per_epoch_data: List[Dict], pulsar_name: str, kinematics_model: str = "legacy_xy"
) -> Dict[str, Any]:
    """
    Test 3: Joint model with annual proxy confounders.

    Fit: |H| = b0 + b1 * v_parallel + b2 * v_perp + b3 * sin(orbital)
         + b4 * cos(orbital) + epsilon
    Tests whether Earth velocity components relative to the CMB dipole add
    explanatory power beyond a simple annual harmonic.
    """
    cmb_vec = _CMB_DIPOLE_VELOCITY_ICRS
    d = compute_per_epoch_data(
        per_epoch_data, pulsar_name, cmb_vec, kinematics_model=kinematics_model
    )
    n = len(d["v_parallel_kms"])
    if n < 10:
        return {"error": f"Insufficient epochs (n={n})"}

    y = d["abs_delta_ns"]
    v_par = d["v_parallel_kms"]
    v_perp = d["v_perp_kms"]
    sin_orb = d["sin_orbital"]
    cos_orb = d["cos_orbital"]

    # Full model
    X_full = np.column_stack([np.ones(n), v_par, v_perp, sin_orb, cos_orb])
    coeffs_full, ss_full, se_full, t_full = _ols_coeffs_rss_stderr_t(X_full, y)

    # Reduced model (no CMB velocity terms)
    X_red = np.column_stack([np.ones(n), sin_orb, cos_orb])
    coeffs_red, ss_red, _, _ = _ols_coeffs_rss_stderr_t(X_red, y)

    # F-test for adding v_par and v_perp
    df1 = 2
    df2 = n - 5
    f_stat = ((ss_red - ss_full) / df1) / (ss_full / df2) if ss_full > 0 else 0.0
    f_p = float(stats.f.sf(f_stat, df1, df2))

    # AIC comparison
    aic_full = n * np.log(ss_full / n) + 2 * 5
    aic_red = n * np.log(ss_red / n) + 2 * 3

    return {
        "n_epochs": n,
        "full_model": {
            "intercept_ns": float(coeffs_full[0]),
            "v_parallel_coeff_ns_per_kms": float(coeffs_full[1]),
            "v_parallel_t": float(t_full[1]),
            "v_perp_coeff_ns_per_kms": float(coeffs_full[2]),
            "v_perp_t": float(t_full[2]),
            "sin_orbital_coeff": float(coeffs_full[3]),
            "sin_orbital_t": float(t_full[3]),
            "cos_orbital_coeff": float(coeffs_full[4]),
            "cos_orbital_t": float(t_full[4]),
            "RSS": float(ss_full),
            "AIC": float(aic_full),
        },
        "reduced_model": {
            "intercept_ns": float(coeffs_red[0]),
            "sin_orbital_coeff": float(coeffs_red[1]),
            "cos_orbital_coeff": float(coeffs_red[2]),
            "RSS": float(ss_red),
            "AIC": float(aic_red),
        },
        "f_test_adding_cmb_terms": float(f_stat),
        "f_test_p_value": f_p,
        "delta_AIC": float(aic_red - aic_full),
    }


def test_multipole(
    per_epoch_data: List[Dict], pulsar_name: str, kinematics_model: str = "legacy_xy"
) -> Dict[str, Any]:
    """
    Test 4: Check for strong quadratic dependence of epoch |H| on v_proj.

    Dependent variable is mean(|geometric_delta|) per epoch in the CMB-bulk
    Stokes model; v_proj enters the weights linearly, but the mapping to this
    summary can be nonlinear. A significant v_proj^2 term flags residual
    curvature beyond a line in (v_proj, |H|) for this diagnostic.
    """
    cmb_vec = _CMB_DIPOLE_VELOCITY_ICRS
    d = compute_per_epoch_data(
        per_epoch_data, pulsar_name, cmb_vec, kinematics_model=kinematics_model
    )
    n = len(d["v_parallel_kms"])
    if n < 10:
        return {"error": f"Insufficient epochs (n={n})"}

    y = d["abs_delta_ns"]
    v_proj = d["v_proj_kms"]
    v_proj_sq = d["v_proj_sq_kms2"]

    # Linear and quadratic nested models (RSS from OLS residuals)
    X_lin = np.column_stack([np.ones(n), v_proj])
    coeffs_lin, ss_lin, _, _ = _ols_coeffs_rss_stderr_t(X_lin, y)

    X_quad = np.column_stack([np.ones(n), v_proj, v_proj_sq])
    coeffs_quad, ss_quad, _, t_quad = _ols_coeffs_rss_stderr_t(X_quad, y)

    # F-test for quadratic term
    df1 = 1
    df2 = n - 3
    f_stat = ((ss_lin - ss_quad) / df1) / (ss_quad / df2) if ss_quad > 0 else 0.0
    f_p = float(stats.f.sf(f_stat, df1, df2))

    aic_lin = n * np.log(ss_lin / n) + 2 * 2
    aic_quad = n * np.log(ss_quad / n) + 2 * 3

    return {
        "n_epochs": n,
        "linear_model": {
            "intercept_ns": float(coeffs_lin[0]),
            "v_proj_coeff_ns_per_kms": float(coeffs_lin[1]),
            "RSS": float(ss_lin),
            "AIC": float(aic_lin),
        },
        "quadratic_model": {
            "intercept_ns": float(coeffs_quad[0]),
            "v_proj_coeff_ns_per_kms": float(coeffs_quad[1]),
            "v_proj_sq_coeff_ns_per_kms2": float(coeffs_quad[2]),
            "v_proj_sq_t": float(t_quad[2]),
            "RSS": float(ss_quad),
            "AIC": float(aic_quad),
        },
        "f_test_quadratic": float(f_stat),
        "f_test_p_value": f_p,
        "delta_AIC": float(aic_lin - aic_quad),
        "conclusion": "dipole_only" if f_p > 0.05 else "quadrupole_significant",
    }


def test_random_direction_permutation(
    per_epoch_data: List[Dict],
    pulsar_name: str,
    n_random: int = 50,
    seed: int = 42,
    kinematics_model: str = "legacy_xy",
) -> Dict[str, Any]:
    """
    Test 5: Random direction control on the sphere.

    Draws N random bulk vectors with the same magnitude as the kinematic CMB
    dipole. Unweighted phase-only ψ must be identical for every draw (same
    stored triplet phases). Weighted ψ and |H| vary with direction; this test
    summarizes how typical the true dipole direction is for those
    velocity-sensitive statistics.
    """
    rng = np.random.default_rng(seed)

    true_summary = recompute_with_custom_cmb(
        per_epoch_data,
        pulsar_name,
        _CMB_DIPOLE_VELOCITY_ICRS,
        kinematics_model=kinematics_model,
    )
    if true_summary.get("error"):
        return {"error": true_summary["error"]}
    true_uw_psi = float(true_summary.get("phase_closure_mean_unweighted_rad", 0.0))
    true_w_psi = float(true_summary.get("phase_closure_mean_rad", 0.0))
    true_H = float(true_summary.get("H_magnitude_ns", 0.0))
    true_w_rbar = float(true_summary.get("phase_closure_rbar", 0.0))
    true_w_rayleigh_z = float(true_summary.get("phase_closure_rayleigh_z", 0.0))

    random_uw_psis: List[float] = []
    random_w_psis: List[float] = []
    random_Hs: List[float] = []
    random_w_rbars: List[float] = []
    random_w_rayleigh_zs: List[float] = []

    for _ in range(n_random):
        cos_theta = 2.0 * rng.random() - 1.0
        sin_theta = float(np.sqrt(max(0.0, 1.0 - cos_theta**2)))
        phi = 2.0 * np.pi * rng.random()
        rand_vec = np.array(
            [
                CMB_DIPOLE_V_KMS * sin_theta * np.cos(phi),
                CMB_DIPOLE_V_KMS * sin_theta * np.sin(phi),
                CMB_DIPOLE_V_KMS * cos_theta,
            ],
            dtype=float,
        )

        summary = recompute_with_custom_cmb(
            per_epoch_data,
            pulsar_name,
            rand_vec,
            kinematics_model=kinematics_model,
        )
        if summary.get("error"):
            return {"error": f"random_permutation failed mid-run: {summary['error']}"}
        random_uw_psis.append(float(summary.get("phase_closure_mean_unweighted_rad", 0.0)))
        random_w_psis.append(float(summary.get("phase_closure_mean_rad", 0.0)))
        random_Hs.append(float(summary.get("H_magnitude_ns", 0.0)))
        random_w_rbars.append(float(summary.get("phase_closure_rbar", 0.0)))
        random_w_rayleigh_zs.append(float(summary.get("phase_closure_rayleigh_z", 0.0)))

    random_uw_psis_arr = np.asarray(random_uw_psis, dtype=float)
    random_w_psis_arr = np.asarray(random_w_psis, dtype=float)
    random_Hs_arr = np.asarray(random_Hs, dtype=float)
    random_w_rbars_arr = np.asarray(random_w_rbars, dtype=float)
    random_w_rayleigh_zs_arr = np.asarray(random_w_rayleigh_zs, dtype=float)

    uw_max_abs_diff = float(np.max(np.abs(random_uw_psis_arr - true_uw_psi)))
    uw_all_identical = bool(np.all(np.abs(random_uw_psis_arr - true_uw_psi) < 1e-10))

    def _circ_abs_diff(a: float, b: float) -> float:
        d = (a - b + math.pi) % (2.0 * math.pi) - math.pi
        return abs(d)

    w_diffs = np.array([_circ_abs_diff(true_w_psi, x) for x in random_w_psis_arr])
    frac_H_ge_true = float(np.mean(random_Hs_arr >= true_H))
    frac_H_le_true = float(np.mean(random_Hs_arr <= true_H))
    # Two-sided rarity: distance from median rank
    rank_pct = float(stats.percentileofscore(random_Hs_arr, true_H, kind="rank"))
    rbar_rank_pct = float(
        stats.percentileofscore(random_w_rbars_arr, true_w_rbar, kind="rank")
    )
    rayleigh_z_rank_pct = float(
        stats.percentileofscore(random_w_rayleigh_zs_arr, true_w_rayleigh_z, kind="rank")
    )

    return {
        "n_random": n_random,
        "seed": seed,
        "true_unweighted_psi_rad": true_uw_psi,
        "random_unweighted_psi_max_abs_delta_from_true_rad": uw_max_abs_diff,
        "unweighted_psi_identical_all_random_directions": uw_all_identical,
        "true_weighted_psi_rad": true_w_psi,
        "random_weighted_psi_mean_rad": float(np.mean(random_w_psis_arr)),
        "random_weighted_psi_std_rad": float(np.std(random_w_psis_arr)),
        "weighted_psi_mean_abs_circular_diff_from_true_rad": float(np.mean(w_diffs)),
        "weighted_psi_max_abs_circular_diff_from_true_rad": float(np.max(w_diffs)),
        "true_weighted_rbar": true_w_rbar,
        "random_weighted_rbar_mean": float(np.mean(random_w_rbars_arr)),
        "random_weighted_rbar_std": float(np.std(random_w_rbars_arr)),
        "weighted_rbar_percentile_rank_true_among_random": rbar_rank_pct,
        "true_weighted_rayleigh_z": true_w_rayleigh_z,
        "random_weighted_rayleigh_z_mean": float(np.mean(random_w_rayleigh_zs_arr)),
        "random_weighted_rayleigh_z_std": float(np.std(random_w_rayleigh_zs_arr)),
        "weighted_rayleigh_z_percentile_rank_true_among_random": rayleigh_z_rank_pct,
        "true_H_ns": true_H,
        "random_H_mean_ns": float(np.mean(random_Hs_arr)),
        "random_H_std_ns": float(np.std(random_Hs_arr)),
        "random_H_min_ns": float(np.min(random_Hs_arr)),
        "random_H_max_ns": float(np.max(random_Hs_arr)),
        "H_percentile_rank_true_among_random": rank_pct,
        "H_empirical_p_ge_true": float((1.0 + np.sum(random_Hs_arr >= true_H)) / (n_random + 1.0)),
        "H_empirical_p_le_true": float((1.0 + np.sum(random_Hs_arr <= true_H)) / (n_random + 1.0)),
        "fraction_random_H_ge_true": frac_H_ge_true,
        "fraction_random_H_le_true": frac_H_le_true,
    }


def main(verbose: bool = True) -> bool:
    parser = argparse.ArgumentParser(
        description="CMB dipole frame analysis for TEP closure delays"
    )
    parser.add_argument(
        "--pulsars",
        nargs="+",
        default=["J0437-4715", "J1603-7202"],
        help="Pulsars to analyze (must match keys in PULSAR_PARAMS)",
    )
    parser.add_argument(
        "--all-pulsars-with-data",
        action="store_true",
        help="Run every pulsar in PULSAR_PARAMS that has a Step 003 per-epoch JSON",
    )
    parser.add_argument(
        "--include-per-epoch-series",
        action="store_true",
        help="Attach long per-epoch arrays to summaries (large JSON output)",
    )
    parser.add_argument(
        "--kinematics-model",
        choices=["legacy_xy", "tangent_plane"],
        default="legacy_xy",
        help=(
            "Velocity-geometry model for projecting Earth/CMB bulk onto the 2D "
            "scintillation plane. legacy_xy reproduces Step 003 ICRS (vx,vy) "
            "behavior; tangent_plane projects 3D vectors into the pulsar-local "
            "east/north tangent plane before applying psi."
        ),
    )
    parser.add_argument(
        "--n-random",
        type=int,
        default=50,
        metavar="N",
        help="Number of random bulk-vector draws for Test 5 (default: 50)",
    )
    args = parser.parse_args()

    def msg(m, level=None):
        if verbose or level in ("ERROR", "WARNING"):
            print_status(m, level)

    if args.all_pulsars_with_data:
        pulsar_list = list_pulsars_with_step003_per_epoch()
        if not pulsar_list:
            msg(
                "No Step 003 per-epoch JSON files found under results/. Run Step 003 first.",
                "ERROR",
            )
            return False
    else:
        pulsar_list = list(args.pulsars)

    msg("=" * 70, "TITLE")
    msg("STEP 048: CMB DIPOLE FRAME ANALYSIS", "TITLE")
    msg("=" * 70, "TITLE")

    cmb_vec = _CMB_DIPOLE_VELOCITY_ICRS
    msg(
        f"CMB dipole (Planck 2018): v={CMB_DIPOLE_V_KMS:.2f} km/s, "
        f"l={CMB_DIPOLE_L_DEG:.4f} deg, b={CMB_DIPOLE_B_DEG:.4f} deg",
        "INFO",
    )
    msg(
        f"  ICRS Cartesian: vx={cmb_vec[0]:.2f}, vy={cmb_vec[1]:.2f}, vz={cmb_vec[2]:.2f} km/s",
        "INFO",
    )
    msg(f"  Pulsars: {', '.join(pulsar_list)}", "INFO")
    msg("", "INFO")

    all_results: Dict[str, Any] = {
        "methodology_note": (
            "Step 003 uses SSB barycentric Earth velocity in ICRS; this step adds a "
            "constant Planck 2018 kinematic dipole vector to the same (vx, vy) effective "
            "velocity used for Stokes weighting. It is a sensitivity model, not a full "
            "scattering-frame transformation. Unweighted Phase Closure uses only stored "
            "triplet phases and is invariant under bulk-vector substitution."
        ),
        "cmb_dipole_parameters": {
            "v_kms": CMB_DIPOLE_V_KMS,
            "l_deg": CMB_DIPOLE_L_DEG,
            "b_deg": CMB_DIPOLE_B_DEG,
            "icrs_vx_kms": float(cmb_vec[0]),
            "icrs_vy_kms": float(cmb_vec[1]),
            "icrs_vz_kms": float(cmb_vec[2]),
            "icrs_unit_x": float(_CMB_DIPOLE_UNIT_ICRS[0]),
            "icrs_unit_y": float(_CMB_DIPOLE_UNIT_ICRS[1]),
            "icrs_unit_z": float(_CMB_DIPOLE_UNIT_ICRS[2]),
            "reference": "Planck Collaboration 2018 (kinematic dipole parameters)",
        },
        "pulsars": {},
        "random_direction_permutation_n": int(args.n_random),
        "kinematics_model": str(args.kinematics_model),
    }

    n_ok = 0
    for pulsar_name in pulsar_list:
        msg(f"\n--- {pulsar_name} ---", "TITLE")

        if pulsar_name not in PULSAR_PARAMS:
            msg(
                f"Skipping {pulsar_name}: not in PULSAR_PARAMS (Step 003 kinematic DB)",
                "WARNING",
            )
            continue

        per_epoch = load_per_epoch_results(pulsar_name)
        if per_epoch is None:
            msg(
                f"Skipping {pulsar_name}: no Step 003 per-epoch JSON matched "
                f"{[str(p) for p in _per_epoch_json_candidates(pulsar_name)]}",
                "WARNING",
            )
            continue

        inc = args.include_per_epoch_series

        vel_comp = compute_velocity_comparison(
            pulsar_name, kinematics_model=args.kinematics_model
        )
        if vel_comp.get("error"):
            msg(f"Velocity comparison failed: {vel_comp['error']}", "WARNING")
            all_results["pulsars"][pulsar_name] = {"error": vel_comp["error"]}
            continue

        msg(
            f"SSB (Step 003) v_eff: {vel_comp['v_eff_ssb_kms']}, "
            f"v_proj={vel_comp['v_proj_ssb_kms']:+.2f} km/s",
            "INFO",
        )
        msg(
            f"CMB bulk model v_eff: {vel_comp['v_eff_cmb_kms']}, "
            f"v_proj={vel_comp['v_proj_cmb_kms']:+.2f} km/s",
            "INFO",
        )
        _plane = (
            "tangent east–north plane"
            if args.kinematics_model == "tangent_plane"
            else "Step 003 ICRS (vx, vy) plane"
        )
        msg(
            f"Angle SSB vs CMB-bulk vector ({_plane}): "
            f"{vel_comp['angle_between_ssb_and_cmb_bulk_deg']:.2f} deg",
            "INFO",
        )
        msg(
            f"Projection sign: SSB={vel_comp['v_proj_sign_ssb']:+.0f}, "
            f"CMB_bulk={vel_comp['v_proj_sign_cmb']:+.0f}",
            "INFO",
        )

        summary_ssb = recompute_epoch_statistics(
            per_epoch,
            pulsar_name,
            frame="ssb",
            kinematics_model=args.kinematics_model,
            include_per_epoch_series=inc,
        )
        if summary_ssb.get("error"):
            msg(f"SSB recomputation: {summary_ssb['error']}", "WARNING")
            all_results["pulsars"][pulsar_name] = {
                "velocity_comparison": vel_comp,
                "error": summary_ssb["error"],
            }
            continue

        msg(
            f"SSB ψ (weighted) = {summary_ssb['phase_closure_mean_rad']:+.4f} rad, "
            f"R_bar={summary_ssb['phase_closure_rbar']:.4f}, "
            f"Rayleigh p={summary_ssb['phase_closure_rayleigh_p']:.3e}",
            "INFO",
        )
        msg(
            f"SSB ψ (unweighted) = {summary_ssb['phase_closure_mean_unweighted_rad']:+.4f} rad, "
            f"R_bar={summary_ssb['phase_closure_rbar_unweighted']:.4f}, "
            f"Rayleigh p={summary_ssb['phase_closure_rayleigh_p_unweighted']:.3e}",
            "INFO",
        )
        msg(
            f"SSB |H| = {summary_ssb['H_magnitude_ns']:.3f} ± "
            f"{summary_ssb['H_sem_ns']:.3f} ns, signed={summary_ssb['H_signed_mean_ns']:+.3f} ns",
            "INFO",
        )

        summary_cmb = recompute_epoch_statistics(
            per_epoch,
            pulsar_name,
            frame="cmb",
            kinematics_model=args.kinematics_model,
            include_per_epoch_series=inc,
        )
        if summary_cmb.get("error"):
            msg(f"CMB recomputation: {summary_cmb['error']}", "WARNING")
            all_results["pulsars"][pulsar_name] = {
                "velocity_comparison": vel_comp,
                "ssb_frame_summary": summary_ssb,
                "error": summary_cmb["error"],
            }
            continue

        msg(
            f"CMB bulk ψ (weighted) = {summary_cmb['phase_closure_mean_rad']:+.4f} rad, "
            f"R_bar={summary_cmb['phase_closure_rbar']:.4f}, "
            f"Rayleigh p={summary_cmb['phase_closure_rayleigh_p']:.3e}",
            "INFO",
        )
        msg(
            f"CMB bulk ψ (unweighted) = {summary_cmb['phase_closure_mean_unweighted_rad']:+.4f} rad, "
            f"R_bar={summary_cmb['phase_closure_rbar_unweighted']:.4f}, "
            f"Rayleigh p={summary_cmb['phase_closure_rayleigh_p_unweighted']:.3e}",
            "INFO",
        )
        msg(
            f"CMB bulk |H| = {summary_cmb['H_magnitude_ns']:.3f} ± "
            f"{summary_cmb['H_sem_ns']:.3f} ns, signed={summary_cmb['H_signed_mean_ns']:+.3f} ns",
            "INFO",
        )

        sign_unchanged = vel_comp["v_proj_sign_ssb"] == vel_comp["v_proj_sign_cmb"]
        status = "PASS" if sign_unchanged else "WARNING"
        msg(
            f"Velocity projection sign unchanged (SSB vs CMB bulk): {sign_unchanged}",
            status,
        )

        msg("\n  --- Enhanced CMB Frame Tests ---", "TITLE")

        msg("  Test 1: Earth velocity along CMB dipole (annual)", "INFO")
        mod_result = test_annual_cmb_projection_modulation(
            per_epoch, pulsar_name, kinematics_model=args.kinematics_model
        )
        if "error" in mod_result:
            msg(f"    {mod_result['error']}", "WARNING")
        else:
            msg(
                f"    v_parallel range: [{mod_result['v_parallel_range_kms'][0]:+.2f}, "
                f"{mod_result['v_parallel_range_kms'][1]:+.2f}] km/s",
                "INFO",
            )
            msg(
                f"    |H| vs v_parallel: r={mod_result['H_vs_vparallel_r']:+.4f}, "
                f"p={mod_result['H_vs_vparallel_p']:.4f}",
                "INFO",
            )
            msg(
                f"    Joint model v_parallel t={mod_result['joint_model_vparallel_t']:+.2f}",
                "INFO",
            )

        msg("  Test 2: Directional specificity", "INFO")
        dir_result = test_directional_specificity(
            per_epoch, pulsar_name, kinematics_model=args.kinematics_model
        )
        if "error" in dir_result:
            msg(f"    {dir_result['error']}", "WARNING")
        else:
            msg(
                f"    Unweighted psi invariant: {dir_result['unweighted_psi_invariant_all_directions']}, "
                f"max |Δ|={dir_result.get('unweighted_psi_max_abs_delta_rad', 0):.2e} rad",
                "INFO",
            )
            for k, v in dir_result["H_magnitudes_by_direction_ns"].items():
                msg(f"    |H| {k}: {v:.3f} ns", "INFO")

        msg("  Test 3: Joint model (annual + Earth||CMB + Earth⊥CMB)", "INFO")
        joint_result = test_joint_model(
            per_epoch, pulsar_name, kinematics_model=args.kinematics_model
        )
        if "error" in joint_result:
            msg(f"    {joint_result['error']}", "WARNING")
        else:
            msg(
                f"    v_parallel t={joint_result['full_model']['v_parallel_t']:+.2f}",
                "INFO",
            )
            msg(
                f"    F-test for Earth CMB-geometry terms: F={joint_result['f_test_adding_cmb_terms']:.2f}, "
                f"p={joint_result['f_test_p_value']:.4f}",
                "INFO",
            )
            msg(
                f"    delta AIC (reduced - full): {joint_result['delta_AIC']:+.1f}",
                "INFO",
            )

        msg("  Test 4: Multipole (v_proj^2) check on |H|", "INFO")
        multi_result = test_multipole(
            per_epoch, pulsar_name, kinematics_model=args.kinematics_model
        )
        if "error" in multi_result:
            msg(f"    {multi_result['error']}", "WARNING")
        else:
            msg(
                f"    v_proj^2 t={multi_result['quadratic_model']['v_proj_sq_t']:+.2f}, "
                f"p={multi_result['f_test_p_value']:.4f}",
                "INFO",
            )
            msg(f"    Conclusion: {multi_result['conclusion']}", "INFO")

        msg(
            f"  Test 5: Random direction control (n={args.n_random})",
            "INFO",
        )
        perm_result = test_random_direction_permutation(
            per_epoch,
            pulsar_name,
            n_random=args.n_random,
            kinematics_model=args.kinematics_model,
        )
        if "error" in perm_result:
            msg(f"    {perm_result['error']}", "WARNING")
        else:
            msg(
                f"    Unweighted ψ identical all draws: "
                f"{perm_result['unweighted_psi_identical_all_random_directions']}, "
                f"max |Δ|={perm_result['random_unweighted_psi_max_abs_delta_from_true_rad']:.2e} rad",
                "INFO",
            )
            msg(
                f"    |H| true={perm_result['true_H_ns']:.3f} ns, random mean="
                f"{perm_result['random_H_mean_ns']:.3f} ± {perm_result['random_H_std_ns']:.3f} ns, "
                f"percentile={perm_result['H_percentile_rank_true_among_random']:.1f}",
                "INFO",
            )

        all_results["pulsars"][pulsar_name] = {
            "velocity_comparison": vel_comp,
            "ssb_frame_summary": summary_ssb,
            "cmb_frame_summary": summary_cmb,
            "sign_consistency": {
                "sign_unchanged": sign_unchanged,
                "status": status,
            },
            "enhanced_tests": {
                "annual_cmb_modulation": mod_result,
                "directional_specificity": dir_result,
                "joint_model": joint_result,
                "multipole": multi_result,
                "random_permutation": perm_result,
            },
        }
        n_ok += 1

    out_path = RESULTS_DIR / "step_048_cmb_dipole_frame_analysis.json"
    all_results["n_pulsars_analyzed"] = n_ok
    all_results["n_pulsars_requested"] = len(pulsar_list)
    with open(out_path, "w") as fh:
        json.dump(all_results, fh, cls=NpEncoder, indent=2)
    msg(f"\nResults saved to {out_path}", "SUCCESS")

    if n_ok == 0:
        msg("No pulsars completed Step 048 (missing data or errors).", "ERROR")
        return False
    return all_results


def step_main(logger=None, verbose: bool = True):
    """Pipeline entry point (run_pipeline.py)."""
    if logger is not None:
        set_step_logger(logger)
    return main(verbose=verbose)


if __name__ == "__main__":
    ok = main()
    raise SystemExit(0 if ok else 1)
