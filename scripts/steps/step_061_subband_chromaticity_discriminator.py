#!/usr/bin/env python3
"""
Step 061: Sub-band Chromaticity and Arc-Geometry Discriminator
=============================================================

PURPOSE: The J0437 phase-closure detection leaves the transport origin
ambiguous between a non-exact (TEP holonomy) channel and residual
scintillation-optics systematics. This step registers the plasma-side
discriminators on the two within-source sub-bands of the 20-cm parent
band (sb0: ~1433 MHz upper half; sb1: ~1305 MHz lower half; nu0/nu1 =
1.098):

  T1. psi chromaticity — are the sub-band circular means equal?
      Epoch-level permutation test on the angular separation of the
      sb0 vs sb1 circular means, plus a pooled triplet-level
      Watson-Williams F-test.

  T2. Delay-domain frequency scaling — the epoch-mean |phase_delta_ns|
      ratio r = sb1/sb0 compared against the three template predictions:
      achromatic delay (r=1), phase-fixed holonomy (r = nu0/nu1 = 1.098),
      plasma dispersive delay nu^{-2} (r = (nu0/nu1)^2 = 1.206).
      Uncertainty by epoch-block bootstrap.

  T3. Arc-geometry correlations — per-epoch circular-mean psi regressed
      against epoch-mean |f_D| and the arc-model delay mismatch
      |tau - tau_pred| / sigma_tau (Spearman, pooled and per band).
      A scintillation-optics artifact should co-vary with arc
      morphology; a transport holonomy should not.

Inputs are the existing Step 003 per-epoch triplet products; no new
arclet detection is performed.

Author: TEP-J0437 Analysis Pipeline
"""

import json
import os
import numpy as np
from datetime import datetime
from pathlib import Path
from scipy.stats import spearmanr

SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"
DATA_DIR = PROJECT_ROOT / "data" / "processed" / "j0437"
OUTPUT_FILE = RESULTS_DIR / "step_061_subband_chromaticity_discriminator.json"

N_PERM = 20000
N_BOOT = 10000
RNG_SEED = 20261001


def print_status(msg, status="INFO"):
    symbols = {"INFO": "i", "SUCCESS": "+", "WARNING": "!", "ERROR": "x", "WORKING": "~"}
    print(f"[{symbols.get(status, '-')}] {msg}")


def circ_mean(angles, weights=None):
    """Weighted circular mean; returns (mean_rad, Rbar)."""
    a = np.asarray(angles, dtype=float)
    w = np.ones_like(a) if weights is None else np.asarray(weights, dtype=float)
    s = np.sum(w * np.sin(a))
    c = np.sum(w * np.cos(a))
    return float(np.arctan2(s, c)), float(np.hypot(s, c) / np.sum(w))


def angular_sep(a, b):
    return float(np.abs(np.angle(np.exp(1j * (a - b)))))


def epoch_psi(epoch):
    """SNR^2-weighted circular mean of triplet psi within one epoch."""
    tr = epoch["triplets"]
    psi = np.array([t["phase_closure_rad"] for t in tr])
    snr = np.array([max(t["closure_snr"], 1e-3) for t in tr])
    m, rbar = circ_mean(psi, snr ** 2)
    return m, rbar, len(tr)


def epoch_nu(epoch, tag):
    """Sub-band centre frequency from the parent npz frequency axis."""
    base = epoch["epoch"].replace(f"_{tag}_secondary", "") + ".npz"
    p = DATA_DIR / base
    if not p.exists():
        return None
    f = np.load(p, allow_pickle=True)["freq_MHz"]
    n = len(f)
    # axis is descending; sb0 = first (upper) half, sb1 = second (lower) half
    half = f[: n // 2] if tag == "sb0" else f[n // 2:]
    return float(np.mean(half))


def load_band(tag):
    eps = json.loads((RESULTS_DIR /
                      f"step_003_closure_final_per_epoch_j0437_{tag}.json").read_text())
    recs = []
    for e in eps:
        nu = epoch_nu(e, tag)
        m, rbar, ntr = epoch_psi(e)
        tr = e["triplets"]
        recs.append({
            "epoch": e["epoch"], "mjd": e["mjd"], "nu_mhz": nu,
            "psi_mean": m, "rbar": rbar, "n_triplets": ntr,
            "mean_abs_fD": float(np.mean([np.mean([abs(t["fD_01"]),
                                                   abs(t["fD_12"]),
                                                   abs(t["fD_02"])])
                                          for t in tr])),
            "mean_arc_mismatch": float(np.mean([
                np.mean([abs(t["tau_01"] - t["tau_01_pred"]),
                         abs(t["tau_12"] - t["tau_12_pred"]),
                         abs(t["tau_02"] - t["tau_02_pred"])])
                / max(t["sigma_us"], 1e-9) for t in tr])),
            "mean_abs_phase_delta_ns": float(np.mean([abs(t["phase_delta_ns"])
                                                      for t in tr])),
            "triplets": tr,
        })
    return recs


def main():
    print_status("Step 061: Sub-band chromaticity / arc-geometry discriminator", "WORKING")
    rng = np.random.default_rng(RNG_SEED)

    bands = {t: load_band(t) for t in ("sb0", "sb1")}
    for t, recs in bands.items():
        print_status(f"{t}: {len(recs)} epochs, "
                     f"nu_ref ~ {np.nanmean([r['nu_mhz'] for r in recs]):.1f} MHz",
                     "SUCCESS")

    nu0 = float(np.nanmean([r["nu_mhz"] for r in bands["sb0"]]))
    nu1 = float(np.nanmean([r["nu_mhz"] for r in bands["sb1"]]))

    # ---------- T1: psi chromaticity ----------
    psi0 = np.array([r["psi_mean"] for r in bands["sb0"]])
    psi1 = np.array([r["psi_mean"] for r in bands["sb1"]])
    m0, r0 = circ_mean(psi0)
    m1, r1 = circ_mean(psi1)
    sep_obs = angular_sep(m0, m1)

    pooled = np.concatenate([psi0, psi1])
    n0, n1 = len(psi0), len(psi1)
    perm_sep = np.empty(N_PERM)
    for i in range(N_PERM):
        idx = rng.permutation(n0 + n1)
        a, b = pooled[idx[:n0]], pooled[idx[n0:]]
        perm_sep[i] = angular_sep(circ_mean(a)[0], circ_mean(b)[0])
    p_chrom = float((np.sum(perm_sep >= sep_obs) + 1) / (N_PERM + 1))

    # Triplet-level Watson-Williams (pooled triplets, unweighted)
    t0 = np.array([t["phase_closure_rad"] for r in bands["sb0"] for t in r["triplets"]])
    t1 = np.array([t["phase_closure_rad"] for r in bands["sb1"] for t in r["triplets"]])
    N = len(t0) + len(t1)
    R0 = np.hypot(np.sum(np.cos(t0)), np.sum(np.sin(t0)))
    R1 = np.hypot(np.sum(np.cos(t1)), np.sum(np.sin(t1)))
    Rp = np.hypot(np.sum(np.cos(np.concatenate([t0, t1]))),
                  np.sum(np.sin(np.concatenate([t0, t1]))))
    # Watson-Williams F (common mean vs different means)
    if R0 + R1 < N and (N - 2) > 0:
        F_ww = ((N - 2) * (R0 + R1 - Rp)) / (N - R0 - R1)
        from scipy.stats import f as fdist
        p_ww = float(1 - fdist.cdf(F_ww, 1, N - 2))
    else:
        F_ww, p_ww = np.nan, np.nan

    # ---------- T2: delay-domain frequency scaling ----------
    d0 = np.array([r["mean_abs_phase_delta_ns"] for r in bands["sb0"]])
    d1 = np.array([r["mean_abs_phase_delta_ns"] for r in bands["sb1"]])
    ratio_obs = float(np.mean(d1) / np.mean(d0))
    boot_ratio = np.empty(N_BOOT)
    for i in range(N_BOOT):
        b0 = rng.choice(d0, size=len(d0), replace=True)
        b1 = rng.choice(d1, size=len(d1), replace=True)
        boot_ratio[i] = np.mean(b1) / np.mean(b0)
    ratio_ci = np.quantile(boot_ratio, [0.025, 0.975]).tolist()
    ratio_se = float(np.std(boot_ratio))

    predictions = {
        "achromatic_delay": 1.0,
        "phase_fixed_holonomy_nu^-1": nu0 / nu1,
        "plasma_dispersive_nu^-2": (nu0 / nu1) ** 2,
    }
    tensions = {k: (ratio_obs - v) / ratio_se if ratio_se > 0 else np.nan
                for k, v in predictions.items()}

    # ---------- T3: arc-geometry correlations ----------
    arc = {}
    for label, recs in [("sb0", bands["sb0"]), ("sb1", bands["sb1"]),
                        ("pooled", bands["sb0"] + bands["sb1"])]:
        psi_e = np.array([r["psi_mean"] for r in recs])
        fd = np.array([r["mean_abs_fD"] for r in recs])
        mm = np.array([r["mean_arc_mismatch"] for r in recs])
        mask = np.isfinite(fd) & np.isfinite(mm) & np.isfinite(psi_e)
        if mask.sum() >= 4:
            r_fd = spearmanr(np.abs(psi_e[mask]), fd[mask])
            r_mm = spearmanr(np.abs(psi_e[mask]), mm[mask])
            arc[label] = {
                "n_epochs": int(mask.sum()),
                "spearman_abspsi_vs_absfD": {"rho": float(r_fd.statistic),
                                             "p": float(r_fd.pvalue)},
                "spearman_abspsi_vs_arc_mismatch": {"rho": float(r_mm.statistic),
                                                      "p": float(r_mm.pvalue)},
            }

    out = {
        "step": "061",
        "title": "Sub-band chromaticity and arc-geometry discriminator",
        "timestamp": datetime.now().isoformat() + "Z",
        "bands": {
            "sb0": {"n_epochs": n0, "nu_ref_mhz_mean": nu0,
                    "psi_circ_mean_rad": m0, "rbar": r0,
                    "mean_abs_phase_delta_ns": float(np.mean(d0))},
            "sb1": {"n_epochs": n1, "nu_ref_mhz_mean": nu1,
                    "psi_circ_mean_rad": m1, "rbar": r1,
                    "mean_abs_phase_delta_ns": float(np.mean(d1))},
            "nu_ratio_nu0_over_nu1": nu0 / nu1,
        },
        "T1_psi_chromaticity": {
            "angular_separation_deg": float(np.degrees(sep_obs)),
            "epoch_permutation_p": p_chrom,
            "n_permutations": N_PERM,
            "triplet_watson_williams_F": float(F_ww),
            "triplet_watson_williams_p": p_ww,
        },
        "T2_delay_frequency_scaling": {
            "ratio_sb1_over_sb0": ratio_obs,
            "bootstrap_se": ratio_se,
            "bootstrap_ci95": ratio_ci,
            "template_predictions": predictions,
            "signed_tensions_sigma": tensions,
            "note": ("ratio near 1.0 = achromatic delay; near nu0/nu1 = "
                     "1.098 = phase-fixed (psi constant, delay ~ 1/nu); "
                     "near 1.206 = plasma dispersive nu^-2"),
        },
        "T3_arc_geometry_correlations": arc,
        "scope": ("Sub-band chromaticity on the 20-cm parent band only; "
                  "per-epoch DM series and secondary-spectrum arc-curvature "
                  "catalogues are not carried in-repo, so DM scaling and "
                  "arc-curvature regression remain registered open items."),
    }

    # Headline
    hl = []
    hl.append(f"Sub-band psi means separated by {np.degrees(sep_obs):.0f} deg "
              f"(epoch permutation p={p_chrom:.3f}; Watson-Williams p={p_ww:.3g})")
    hl.append(f"Delay-domain band ratio {ratio_obs:.3f} +/- {ratio_se:.3f}: "
              f"achromatic tension {abs(tensions['achromatic_delay']):.2f} sigma, "
              f"nu^-2 tension {abs(tensions['plasma_dispersive_nu^-2']):.2f} sigma, "
              f"nu^-1 tension {abs(tensions['phase_fixed_holonomy_nu^-1']):.2f} sigma")
    out["headline"] = hl

    OUTPUT_FILE.write_text(json.dumps(out, indent=2))
    print_status(f"psi separation {np.degrees(sep_obs):.1f} deg, perm p={p_chrom:.4f}, WW p={p_ww:.3g}", "SUCCESS")
    print_status(f"delay ratio {ratio_obs:.3f} ± {ratio_se:.3f} CI {ratio_ci}", "SUCCESS")
    print_status(f"Results written to {OUTPUT_FILE.name}", "SUCCESS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
