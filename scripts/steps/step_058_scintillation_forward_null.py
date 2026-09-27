#!/usr/bin/env python3
r"""
================================================================================
STEP 058: STANDARD-SCINTILLATION FORWARD NULL FOR THE CLOSURE-PHASE STATISTIC
================================================================================

Motivation
----------
The flagship phase-closure detection (psi-bar_uw = 1.120 rad, R-bar = 0.304,
Z_uw = 202, p = 1.3e-44; epoch-weighted psi-bar = 0.984 rad) tests the null
hypothesis of *circular uniformity* of arc-triplet closure phases
psi = phi_01 + phi_12 - phi_02 measured at cross-term peaks of the complex
secondary spectrum. A referee-grade objection is that standard scintillation
physics may already violate uniformity mundanely: multi-screen propagation,
non-quadratic arc curvature, discrete arclet families, and persistent
(long-lived) arclet structure can concentrate closure phases without any
non-factorizable transport.

This step builds that forward null directly. Synthetic thin-screen dynamic
spectra I(nu,t) = |sum_k g_k exp[i 2 pi (nu tau_k - fD_k t) + i phi_k]|^2 are
generated under standard scintillation physics (scalar additive image phases),
pushed through the UNCHANGED production chain:

    step_002.compute_secondary_spectrum  ->  step_002.detect_arcs
    -> step_002.detect_arclets / detect_peaks_direct
    -> step_003.find_best_triplets       ->  per-triplet phase_closure_rad

including the Hough arc gate, the two-arc (cross-screen) triplet requirement,
sub-pixel cross-term measurement, 3x3 complex-patch phase extraction, and all
SNR cuts. Any psi clustering the simulators produce is therefore a property of
standard physics plus the identical measurement chain.

Epoch dimensions, sampling, and bands are drawn from the empirical
distribution of the real J0437 processed dynamic spectra so the null inherits
the same (f_D, tau) grids the data used.

Models
------
noise_only             pure Gaussian noise (absolute null)
single_screen          one parabolic arc; tests the cross-screen gate itself
two_screens            two arcs, fresh image populations per epoch (the correct
                       ISS null: arc curvatures are persistent, image content
                       decorrelates on the scintillation timescale)
curved_arc             two screens with non-quadratic curvature
                       tau = eta fD^2 (1 + gamma fD^2)
arclet_families        discrete arclet families (clustered images per screen)
interior_fill          arcs plus sub-arc power filling the parabola interior
three_screens          multi-screen propagation (three arcs)
persistent_arclets     30% of images fixed across epochs (documented
                       long-lived J0437 arclets), remainder fresh
fully_persistent       image geometry and phases fixed; only noise varies
                       (maximally persistent screen -- tests whether a fixed
                       screen realization yields a fixed non-zero psi-bar)
two_screens_scrambled  two_screens with arg F randomized at fixed |F| --
                       isolates amplitude-structure-only contributions

Outputs
-------
results/step_058_scintillation_forward_null.json
    per-model pooled and epoch-level circular statistics, sign-pattern
    decomposition of the measured-delay sign combinations, psi histograms,
    and a verdict block comparing each model against the observed headline
    values loaded from results/step_003_closure_final_summary_j0437.json and
    results/step_048_cmb_dipole_frame_analysis.json.

Usage
-----
    PYTHONPATH=. python scripts/steps/step_058_scintillation_forward_null.py
        [--epochs 1500] [--geometries 8] [--geom-epochs 60] [--workers N]
"""

from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "steps"))

from scripts.utils.json_numpy import NpEncoder
from scripts.utils.logger import TEPLogger, print_status, set_step_logger
from scripts.utils.parallel_workers import configure_blas_thread_env, worker_count

RESULTS_DIR = PROJECT_ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

RNG_SEED = 20260927

# ---------------------------------------------------------------------------
# Empirical epoch-configuration distribution (measured from
# data/processed/j0437/*.npz, every 20th file): (n_time, n_freq, f_min, f_max,
# dt_s, weight)
# ---------------------------------------------------------------------------
EPOCH_CONFIGS: List[Tuple[int, int, float, float, float, int]] = [
    (64, 1024, 1241.0, 1497.0, 60.0, 54),
    (64, 1024, 2589.0, 3612.0, 60.0, 32),
    (121, 512, 1182.0, 1582.0, 32.0, 25),
    (58, 128, 653.0, 717.0, 67.0, 21),
    (58, 128, 1373.0, 1437.0, 67.0, 20),
    (58, 128, 1309.0, 1373.0, 67.0, 12),
    (120, 512, 1182.0, 1582.0, 32.0, 12),
    (58, 128, 3224.0, 3288.0, 67.0, 11),
    (29, 128, 653.0, 717.0, 67.0, 9),
    (58, 128, 3096.0, 3160.0, 67.0, 8),
    (16, 2048, 1241.0, 1497.0, 60.0, 8),
    (121, 256, 628.0, 828.0, 32.0, 7),
    (121, 256, 620.0, 820.0, 32.0, 7),
    (10, 1024, 1241.0, 1497.0, 60.0, 7),
    (64, 512, 700.0, 764.0, 60.0, 7),
]

# Arc-curvature support for the J0437 screens (us/mHz^2). Hough detections on
# real epochs return eta ~ 0.004-0.010 with the two detected arcs often closely
# spaced; the physical screen pair (89.8 pc / 124 pc) gives eta ratio ~1.5.
ETA_MIN, ETA_MAX = 0.0030, 0.0100


def _draw_epoch_config(rng: np.random.RandomState) -> Tuple[int, int, np.ndarray, float]:
    """Draw (n_time, n_freq, freq_MHz, dt_s) from the empirical distribution."""
    weights = np.array([c[5] for c in EPOCH_CONFIGS], dtype=float)
    weights /= weights.sum()
    idx = rng.choice(len(EPOCH_CONFIGS), p=weights)
    n_time, n_freq, fmin, fmax, dt_s, _ = EPOCH_CONFIGS[idx]
    freq_MHz = np.linspace(fmin, fmax, n_freq)
    return int(n_time), int(n_freq), freq_MHz, float(dt_s)


def _draw_eta_pair(rng: np.random.RandomState) -> Tuple[float, float]:
    """Two-screen curvatures: eta1 uniform in range, eta2 = eta1 * U(1.02,1.6)."""
    eta1 = rng.uniform(ETA_MIN, ETA_MAX)
    eta2 = eta1 * rng.uniform(1.02, 1.6)
    return eta1, min(eta2, ETA_MAX * 1.6)


def _draw_images_on_arc(
    rng: np.random.RandomState,
    eta: float,
    fmax_mHz: float,
    n_img: int,
    curve_gamma: float = 0.0,
    interior_frac: float = 0.0,
    n_families: int = 0,
) -> List[Tuple[float, float, float, float]]:
    """Draw (tau_us, fD_mHz, amplitude, phase) image list for one arc.

    Image fD values are drawn over the observable range (Gaussian sigma ~4-10
    mHz truncated at 0.85*f_Nyquist), matching the concentration of detected
    arclets at low |fD| in the real catalogue. tau = eta fD^2 (optionally with
    a quartic correction); a fraction may fill the parabola interior.
    """
    terms: List[Tuple[float, float, float, float]] = []
    if n_families > 0:
        n_fam = int(n_families)
        cents = rng.uniform(-0.55 * fmax_mHz, 0.55 * fmax_mHz, n_fam)
        per = max(2, n_img // n_fam)
        fd_all = np.concatenate(
            [c + rng.normal(0.0, 1.5, per) for c in cents]
        )
    else:
        sigma_f = rng.uniform(4.0, 10.0)
        fd_all = rng.normal(0.0, sigma_f, n_img)
    fd_all = fd_all[np.abs(fd_all) < 0.85 * fmax_mHz]
    for fk in fd_all:
        tk = eta * fk * fk
        if curve_gamma != 0.0:
            tk *= 1.0 + curve_gamma * fk * fk
        if interior_frac > 0.0 and rng.uniform() < interior_frac:
            tk *= rng.uniform(0.15, 1.0)
        terms.append(
            (float(tk), float(fk), float(rng.uniform(0.1, 0.7)), float(rng.uniform(0, 2 * np.pi)))
        )
    return terms


def _dynspec_from_terms(
    terms: List[Tuple[float, float, float, float]],
    n_time: int,
    freq_MHz: np.ndarray,
    dt_s: float,
    rng: np.random.RandomState,
    noise_frac: float,
) -> np.ndarray:
    """Standard additive-image dynamic spectrum I = |sum_k g_k e^{i theta_k}|^2."""
    t = np.arange(n_time) * dt_s
    nu_Hz = freq_MHz * 1e6
    E = np.zeros((n_time, len(freq_MHz)), dtype=np.complex128)
    for tk, fk, ak, pk in terms:
        E += ak * np.exp(
            1j
            * (
                2 * np.pi * (nu_Hz[None, :] * tk * 1e-6 - fk * 1e-3 * t[:, None])
                + pk
            )
        )
    I = np.abs(E) ** 2
    if noise_frac > 0.0:
        I = I + rng.normal(0.0, noise_frac * float(I.mean()), I.shape)
    return I


def _screen_terms(
    rng: np.random.RandomState,
    model: str,
    fmax_mHz: float,
    geom: Optional[List[Tuple[float, float, float, float]]] = None,
    persist_frac: float = 0.0,
    persist_terms: Optional[List[Tuple[float, float, float, float]]] = None,
) -> Tuple[List[Tuple[float, float, float, float]], List[Tuple[float, float, float, float]]]:
    """Return (terms, persistent_terms) for one epoch under the given model."""
    if geom is not None:
        # fully persistent geometry: identical image list every epoch
        return list(geom), []

    n_img = int(rng.randint(25, 61))
    carrier_amp = float(rng.uniform(1.5, 2.5))
    noise_imgs: List[Tuple[float, float, float, float]] = []

    def _screen_block(eta: float) -> List[Tuple[float, float, float, float]]:
        block = [(0.0, 0.0, carrier_amp, float(rng.uniform(0, 2 * np.pi)))]
        if model == "curved_arc":
            gamma = float(rng.uniform(-0.02, 0.02))
            block += _draw_images_on_arc(rng, eta, fmax_mHz, n_img, curve_gamma=gamma)
        elif model == "arclet_families":
            block += _draw_images_on_arc(
                rng, eta, fmax_mHz, n_img, n_families=int(rng.randint(2, 6))
            )
        elif model == "interior_fill":
            block += _draw_images_on_arc(
                rng, eta, fmax_mHz, n_img, interior_frac=0.6
            )
        else:
            block += _draw_images_on_arc(rng, eta, fmax_mHz, n_img)
        return block

    if model == "single_screen":
        eta1 = rng.uniform(ETA_MIN, ETA_MAX)
        terms = _screen_block(eta1)
    elif model == "three_screens":
        eta1, eta2 = _draw_eta_pair(rng)
        eta3 = min(eta2 * rng.uniform(1.05, 1.4), ETA_MAX * 2.0)
        terms = _screen_block(eta1) + _screen_block(eta2) + _screen_block(eta3)
    else:
        eta1, eta2 = _draw_eta_pair(rng)
        terms = _screen_block(eta1) + _screen_block(eta2)

    if persist_frac > 0.0:
        # Split off a persistent subset: identical (tau,fD,amp,phase) each epoch.
        if persist_terms is None or len(persist_terms) == 0:
            n_p = max(2, int(persist_frac * len(terms)))
            idx = rng.choice(len(terms), size=n_p, replace=False)
            persist_terms = [terms[i] for i in idx]
        fresh_terms = [t for t in terms if t not in persist_terms]
        terms = fresh_terms + list(persist_terms)
    return terms, (persist_terms or [])


def _simulate_epoch(
    model: str,
    seed: int,
    geom_terms: Optional[List[Tuple[float, float, float, float]]] = None,
    persist_terms: Optional[List[Tuple[float, float, float, float]]] = None,
) -> Dict[str, Any]:
    """Simulate one epoch and run the identical production chain.

    Returns dict with per-triplet psi, sign pattern, snr, and diagnostics; or
    {'viable': False} when no cross-screen triplets survive (as in production).
    """
    from scripts.steps.step_002_secondary_spectra import (
        compute_secondary_spectrum,
        detect_arcs,
        detect_arclets,
        detect_peaks_direct,
    )
    from scripts.steps.step_003_closure_delays_final import find_best_triplets

    rng = np.random.RandomState(seed)
    n_time, n_freq, freq_MHz, dt_s = _draw_epoch_config(rng)
    fmax_mHz = 1.0 / (2.0 * dt_s) * 1e3

    if model == "noise_only":
        I = rng.normal(0.0, 1.0, (n_time, n_freq)) ** 2
        persist_out: List = []
    else:
        noise_frac = float(rng.uniform(0.03, 0.10))
        persist_frac = 0.30 if model == "persistent_arclets" else 0.0
        terms, persist_out = _screen_terms(
            rng, model, fmax_mHz,
            geom=geom_terms if model == "fully_persistent" else None,
            persist_frac=persist_frac,
            persist_terms=persist_terms,
        )
        I = _dynspec_from_terms(terms, n_time, freq_MHz, dt_s, rng, noise_frac)

    ss = compute_secondary_spectrum(I, dt_s, freq_MHz)
    S = ss["secondary"]
    F = ss["secondary_complex"]
    tau_us = ss["tau_us"]
    fD_mHz = ss["fD_mHz"]

    if model == "two_screens_scrambled":
        rng_s = np.random.RandomState(seed + 5_000_000)
        F = (np.abs(F) * np.exp(1j * rng_s.uniform(0, 2 * np.pi, F.shape))).astype(
            np.complex64
        )
        S = np.abs(F) ** 2

    arcs = detect_arcs(S, tau_us, fD_mHz)
    arclets = detect_arclets(S, tau_us, fD_mHz, arcs, min_snr=1.5)
    if len(arclets) < 3:
        dp = detect_peaks_direct(S, tau_us, fD_mHz, min_snr=1.5)
        if len(dp) > len(arclets):
            arclets = dp

    eta1 = arcs[0]["eta"] if len(arcs) >= 1 else 0.0
    eta2 = arcs[1]["eta"] if len(arcs) >= 2 else 0.0

    res = find_best_triplets(
        S,
        F,
        tau_us,
        fD_mHz,
        arclets,
        np.array([-80.0, -40.0]),  # v_eff; psi is independent of it
        eta1=eta1,
        eta2=eta2,
        mjd=55000.0,
        nu_ref_mhz=float(np.mean(freq_MHz)),
        pulsar_name="J0437-4715",
    )

    if not res:
        return {"viable": False, "n_arcs": len(arcs), "persist_terms": persist_out if model == "persistent_arclets" else []}

    psis = np.array([r.phase_closure_rad for r in res], dtype=float)
    snrs = np.array([r.snr for r in res], dtype=float)
    signs = np.array(
        [
            [
                int(np.sign(r.tau_01.tau_meas)),
                int(np.sign(r.tau_12.tau_meas)),
                int(np.sign(r.tau_02.tau_meas)),
            ]
            for r in res
        ],
        dtype=np.int8,
    )
    return {
        "viable": True,
        "n_arcs": len(arcs),
        "eta1": float(eta1),
        "eta2": float(eta2),
        "n_arclets": int(len(arclets)),
        "psis": psis,
        "snrs": snrs,
        "signs": signs,
        "persist_terms": persist_out if model == "persistent_arclets" else [],
    }


def _circular_mean_rbar(angles: np.ndarray, weights: Optional[np.ndarray] = None) -> Tuple[float, float]:
    if weights is None:
        weights = np.ones(len(angles))
    z = np.sum(weights * np.exp(1j * angles))
    w = np.sum(weights)
    if w <= 0:
        return 0.0, 0.0
    return float(np.angle(z / w)), float(np.abs(z) / w)


def _rayleigh(angles: np.ndarray) -> Tuple[float, float]:
    from scipy import stats as _st

    n = len(angles)
    if n < 3:
        return 0.0, 1.0
    _, r = _circular_mean_rbar(angles)
    z = 2.0 * n * r * r
    return float(z), float(_st.chi2.sf(z, 2))


def _worker(task: Tuple[str, int, Optional[list], Optional[list]]) -> Dict[str, Any]:
    model, seed, geom, persist = task
    try:
        return _simulate_epoch(model, seed, geom_terms=geom, persist_terms=persist)
    except Exception as e:  # noqa: BLE001 - keep ensemble robust
        return {"viable": False, "error": f"{type(e).__name__}: {e}"}


def _run_model(
    model: str,
    n_epochs: int,
    seed0: int,
    workers: int,
    n_geometries: int = 1,
    geom_epochs: int = 60,
) -> Dict[str, Any]:
    """Run the forward null for one model and return full statistics."""
    all_psis: List[np.ndarray] = []
    all_snrs: List[np.ndarray] = []
    all_signs: List[np.ndarray] = []
    epoch_dirs: List[float] = []
    n_arc_hist: List[int] = []
    n_viable = 0
    n_err = 0

    tasks: List[Tuple[str, int, Optional[list], Optional[list]]] = []

    if model == "fully_persistent":
        # G independent fixed geometries x E epochs each: the per-geometry
        # pooled psi-bar gives the direction scatter across screen realizations.
        geom_records = []
        for g in range(n_geometries):
            grng = np.random.RandomState(seed0 + 1000 * g)
            gterms: List[Tuple[float, float, float, float]] = []
            eta1, eta2 = _draw_eta_pair(grng)
            n_img = int(grng.randint(25, 61))
            for eta in (eta1, eta2):
                gterms.append((0.0, 0.0, float(grng.uniform(1.5, 2.5)), float(grng.uniform(0, 2 * np.pi))))
                gterms += _draw_images_on_arc(grng, eta, 12.0, n_img)
            geom_records.append({"seed": seed0 + 1000 * g, "terms": gterms})
            for e in range(geom_epochs):
                tasks.append((model, seed0 + 1000 * g + e + 1, gterms, None))
        results_by_geom: Dict[int, List[Dict[str, Any]]] = {g: [] for g in range(n_geometries)}
    elif model == "persistent_arclets":
        # One persistent arclet set shared by all epochs; fresh remainder.
        prng = np.random.RandomState(seed0 - 1)
        pterms: List[Tuple[float, float, float, float]] = []
        eta1, eta2 = _draw_eta_pair(prng)
        for eta in (eta1, eta2):
            pterms.append((0.0, 0.0, float(prng.uniform(1.5, 2.5)), float(prng.uniform(0, 2 * np.pi))))
            pterms += _draw_images_on_arc(prng, eta, 12.0, 12)
        for e in range(n_epochs):
            tasks.append((model, seed0 + e + 1, None, pterms))
    else:
        for e in range(n_epochs):
            tasks.append((model, seed0 + e + 1, None, None))

    print_status(f"  {model}: submitting {len(tasks)} simulated epochs", "INFO")

    if model == "fully_persistent":
        # collect per-geometry then aggregate
        geom_out: Dict[int, Dict[str, Any]] = {}
        geom_epoch_dirs: Dict[int, List[float]] = {g: [] for g in range(n_geometries)}
        geom_psis: Dict[int, List[np.ndarray]] = {g: [] for g in range(n_geometries)}
        geom_signs: Dict[int, List[np.ndarray]] = {g: [] for g in range(n_geometries)}
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futs = {}
            for t_ in tasks:
                g_idx = (t_[1] - seed0 - 1) // 1000 if False else None
                # reconstruct geometry index from stored seed base
                g_idx = int((t_[1] - 1 - seed0) // 1000)
                futs[pool.submit(_worker, t_)] = g_idx
            for fut in as_completed(futs):
                g_idx = futs[fut]
                r = fut.result()
                if r.get("error"):
                    n_err += 1
                    continue
                if not r.get("viable"):
                    continue
                geom_psis[g_idx].append(r["psis"])
                geom_signs[g_idx].append(r["signs"])
                m_e, _ = _circular_mean_rbar(r["psis"])
                geom_epoch_dirs[g_idx].append(m_e)
                all_psis.append(r["psis"])
                all_signs.append(r["signs"])
                epoch_dirs.append(m_e)
                n_viable += 1
        per_geom = []
        for g in range(n_geometries):
            if geom_psis[g]:
                arr = np.concatenate(geom_psis[g])
                m_, r_ = _circular_mean_rbar(arr)
                per_geom.append({"geometry": int(g), "n_triplets": int(len(arr)), "psi_uw": m_, "rbar": r_})
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futs = [pool.submit(_worker, t_) for t_ in tasks]
            for fut in as_completed(futs):
                r = fut.result()
                if r.get("error"):
                    n_err += 1
                    continue
                if not r.get("viable"):
                    continue
                all_psis.append(r["psis"])
                all_snrs.append(r["snrs"])
                all_signs.append(r["signs"])
                m_e, _ = _circular_mean_rbar(r["psis"])
                epoch_dirs.append(m_e)
                n_arc_hist.append(int(r.get("n_arclets", 0)))
                n_viable += 1

    out: Dict[str, Any] = {
        "model": model,
        "n_epochs_submitted": len(tasks),
        "n_viable_epochs": n_viable,
        "yield_fraction": float(n_viable / max(len(tasks), 1)),
        "n_errors": n_err,
    }

    if not all_psis:
        out["note"] = "no viable epochs - model produced no cross-screen triplets"
        return out

    pooled = np.concatenate(all_psis)
    signs_arr = np.concatenate(all_signs)
    psi_uw, rbar_uw = _circular_mean_rbar(pooled)
    z_uw, p_uw = _rayleigh(pooled)
    epoch_dirs_arr = np.array(epoch_dirs)
    dir_mean, dir_rbar = _circular_mean_rbar(epoch_dirs_arr)
    z_dir, p_dir = _rayleigh(epoch_dirs_arr)

    hist, edges = np.histogram(pooled, bins=12, range=(-np.pi, np.pi))
    sign_table = {}
    for pat in [(1, 1, 1), (-1, 1, 1), (1, -1, 1), (-1, -1, 1), (1, 1, -1), (1, -1, -1), (-1, 1, -1), (-1, -1, -1)]:
        mask = np.all(signs_arr == np.array(pat, dtype=np.int8), axis=1)
        if mask.sum() >= 10:
            mm, rr = _circular_mean_rbar(pooled[mask])
            sign_table["".join("+" if s > 0 else "-" for s in pat)] = {
                "n": int(mask.sum()),
                "psi_uw": mm,
                "rbar": rr,
            }

    out.update(
        {
            "n_triplets": int(len(pooled)),
            "psi_uw_rad": psi_uw,
            "rbar_uw": rbar_uw,
            "rayleigh_z_uw": z_uw,
            "rayleigh_p_uw": p_uw,
            "epoch_direction_mean_rad": dir_mean,
            "epoch_direction_rbar": dir_rbar,
            "epoch_direction_rayleigh_z": z_dir,
            "epoch_direction_rayleigh_p": p_dir,
            "n_epoch_means": int(len(epoch_dirs_arr)),
            "psi_hist_12bin": [int(h) for h in hist],
            "sign_pattern_breakdown": sign_table,
            "median_n_arclets": float(np.median(n_arc_hist)) if n_arc_hist else None,
        }
    )
    if model == "fully_persistent":
        out["per_geometry_pooled"] = per_geom
        out["geometry_direction_rbar"] = (
            _circular_mean_rbar(np.array([g["psi_uw"] for g in per_geom]))[1]
            if per_geom
            else None
        )
    return out


def _load_observed() -> Dict[str, Any]:
    obs: Dict[str, Any] = {}
    f3 = RESULTS_DIR / "step_003_closure_final_summary_j0437.json"
    f48 = RESULTS_DIR / "step_048_cmb_dipole_frame_analysis.json"
    if f3.exists():
        d3 = json.load(open(f3))
        obs["epoch_weighted"] = {
            "psi_mean_rad": d3.get("phase_closure_mean_rad"),
            "rbar": d3.get("phase_closure_rbar"),
            "rayleigh_z": d3.get("phase_closure_rayleigh_z"),
            "rayleigh_p": d3.get("phase_closure_rayleigh_p"),
            "n_independent_samples": d3.get("n_independent_samples"),
        }
    if f48.exists():
        d48 = json.load(open(f48))
        s = d48.get("pulsars", {}).get("J0437-4715", {}).get("ssb_frame_summary", {})
        obs["pooled_unweighted"] = {
            "psi_mean_rad": s.get("phase_closure_mean_unweighted_rad"),
            "rbar": s.get("phase_closure_rbar_unweighted"),
            "rayleigh_z": s.get("phase_closure_rayleigh_z_unweighted"),
            "rayleigh_p": s.get("phase_closure_rayleigh_p_unweighted"),
            "n_total_triplets": s.get("n_total_triplets"),
        }
    return obs


def run(
    n_epochs: int = 1500,
    n_geometries: int = 8,
    geom_epochs: int = 60,
    workers: Optional[int] = None,
    models: Optional[List[str]] = None,
) -> Dict[str, Any]:
    configure_blas_thread_env()
    if workers is None:
        workers = worker_count(role="cpu_bound", reserve=2)

    all_models = [
        "noise_only",
        "single_screen",
        "two_screens",
        "curved_arc",
        "arclet_families",
        "interior_fill",
        "three_screens",
        "persistent_arclets",
        "fully_persistent",
        "two_screens_scrambled",
    ]
    if models:
        all_models = [m for m in all_models if m in models]

    observed = _load_observed()

    print_status("=" * 70, "TITLE")
    print_status("STEP 058: STANDARD-SCINTILLATION FORWARD NULL", "TITLE")
    print_status("=" * 70, "TITLE")
    print_status(f"Models: {all_models}", "INFO")
    print_status(f"Epochs/model: {n_epochs} (fully_persistent: {n_geometries} geometries x {geom_epochs})", "INFO")
    print_status(f"Workers: {workers}", "INFO")

    per_model: Dict[str, Any] = {}
    for mi, model in enumerate(all_models):
        print_status(f"--- model {mi+1}/{len(all_models)}: {model} ---", "TITLE")
        seed0 = RNG_SEED + 100_000 * mi
        per_model[model] = _run_model(
            model,
            n_epochs=n_epochs,
            seed0=seed0,
            workers=workers,
            n_geometries=n_geometries,
            geom_epochs=geom_epochs,
        )
        pm = per_model[model]
        if pm.get("n_triplets"):
            print_status(
                f"    viable {pm['n_viable_epochs']}/{pm['n_epochs_submitted']}, "
                f"triplets {pm['n_triplets']}, psi_uw={pm['psi_uw_rad']:+.3f}, "
                f"Rbar={pm['rbar_uw']:.3f}, Z={pm['rayleigh_z_uw']:.1f}, "
                f"dir_Rbar={pm['epoch_direction_rbar']:.3f}",
                "RESULT",
            )

    # Verdicts: a model 'explains' the observation only if it reproduces BOTH
    # the pooled concentration magnitude AND the epoch-direction persistence.
    obs_pooled = observed.get("pooled_unweighted", {})
    obs_epoch = observed.get("epoch_weighted", {})
    obs_rbar = obs_pooled.get("rbar") or 0.304
    obs_dir_rbar = None
    # epoch-direction persistence of the real sample: rbar of per-epoch means
    # (read from per-epoch file when present)
    f_ep = RESULTS_DIR / "step_003_closure_final_per_epoch_j0437.json"
    if f_ep.exists():
        ep = json.load(open(f_ep))
        means = []
        for e in ep:
            ps = [t["phase_closure_rad"] for t in e.get("triplets", []) if t.get("phase_closure_rad") is not None]
            if len(ps) >= 5:
                means.append(_circular_mean_rbar(np.array(ps))[0])
        if means:
            obs_dir_rbar = _circular_mean_rbar(np.array(means))[1]

    for model, pm in per_model.items():
        if not pm.get("n_triplets"):
            pm["verdict"] = "no_triplets"
            continue
        mag_ok = pm["rbar_uw"] >= 0.5 * obs_rbar  # within factor ~2 of observed
        dir_ok = obs_dir_rbar is not None and pm["epoch_direction_rbar"] >= 0.5 * obs_dir_rbar
        pm["verdict"] = {
            "reproduces_concentration_magnitude": bool(mag_ok),
            "reproduces_direction_persistence": bool(dir_ok),
            "rbar_ratio_model_to_observed": float(pm["rbar_uw"] / obs_rbar),
            "dir_rbar_ratio_model_to_observed": (
                float(pm["epoch_direction_rbar"] / obs_dir_rbar) if obs_dir_rbar else None
            ),
        }

    result = {
        "step": "step_058_scintillation_forward_null",
        "description": (
            "Forward null for the phase-closure statistic under standard "
            "thin-screen scintillation: additive-image dynamic spectra pushed "
            "through the unchanged step_002 -> step_003 production chain "
            "(Hough arcs, arclet detection, cross-screen triplet requirement, "
            "sub-pixel cross-term measurement, complex-patch phase extraction)."
        ),
        "rng_seed": RNG_SEED,
        "n_epochs_per_model": n_epochs,
        "n_geometries_persistent": n_geometries,
        "geom_epochs_persistent": geom_epochs,
        "workers": workers,
        "epoch_config_table": "empirical (n_time,n_freq,band,dt) from data/processed/j0437",
        "observed": observed,
        "observed_epoch_direction_rbar": obs_dir_rbar,
        "models": per_model,
    }

    out_path = RESULTS_DIR / "step_058_scintillation_forward_null.json"
    with open(out_path, "w") as fh:
        json.dump(result, fh, indent=2, cls=NpEncoder)
    print_status(f"Wrote {out_path}", "SUCCESS")
    return result


def step_main(logger=None, verbose=True, n_epochs: int = 1500, **kw):
    if logger:
        set_step_logger(logger)
    return run(n_epochs=n_epochs, **kw)


def main():
    ap = argparse.ArgumentParser(description="Standard-scintillation forward null")
    ap.add_argument("--epochs", type=int, default=1500)
    ap.add_argument("--geometries", type=int, default=8)
    ap.add_argument("--geom-epochs", type=int, default=60)
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--models", type=str, default=None, help="comma-separated subset")
    args = ap.parse_args()
    models = args.models.split(",") if args.models else None
    run(
        n_epochs=args.epochs,
        n_geometries=args.geometries,
        geom_epochs=args.geom_epochs,
        workers=args.workers,
        models=models,
    )


if __name__ == "__main__":
    main()
