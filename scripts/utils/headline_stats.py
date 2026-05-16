#!/usr/bin/env python3
"""
Canonical headline statistics from frozen pipeline JSON outputs.

Manuscript HTML and replication checks should match these values after each
full pipeline run.
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"

SUMMARY_J0437 = RESULTS_DIR / "step_003_closure_final_summary_j0437.json"
SUMMARY_FALLBACK = RESULTS_DIR / "step_003_closure_final_summary.json"
SUMMARY_J1603 = RESULTS_DIR / "step_003_closure_final_summary_j1603.json"
FRAME_ANALYSIS = RESULTS_DIR / "step_048_cmb_dipole_frame_analysis.json"
STEP_007 = RESULTS_DIR / "step_007_independent_statistical_validation_results.json"


def load_results_json(name: str) -> Dict[str, Any]:
    path = RESULTS_DIR / name
    if not path.exists():
        raise FileNotFoundError(f"Missing results file: {path}")
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


@lru_cache(maxsize=1)
def load_j0437_step003() -> Dict[str, Any]:
    path = SUMMARY_J0437 if SUMMARY_J0437.exists() else SUMMARY_FALLBACK
    return load_results_json(path.name)


@lru_cache(maxsize=1)
def load_j1603_step003() -> Dict[str, Any]:
    return load_results_json(SUMMARY_J1603.name)


@lru_cache(maxsize=1)
def load_step007() -> Dict[str, Any]:
    return load_results_json(STEP_007.name)


@lru_cache(maxsize=1)
def load_frame_analysis() -> Dict[str, Any]:
    return load_results_json(FRAME_ANALYSIS.name)


def frame_summary(pulsar: str, frame: str) -> Dict[str, Any]:
    data = load_frame_analysis()
    key = f"{frame}_frame_summary"
    return data["pulsars"][pulsar][key]


def j0437_headlines() -> Dict[str, Any]:
    """Headline bundle for J0437-4715 (weighted Step 003 + unweighted Step 048 SSB)."""
    s3 = load_j0437_step003()
    s48 = frame_summary("J0437-4715", "ssb")

    psi = float(s3["phase_closure_mean_rad"])
    psi_se = float(s3["phase_closure_circ_se_rad"])
    ci_lo = float(s3["phase_closure_bootstrap_ci_95_lower_rad"])
    ci_hi = float(s3["phase_closure_bootstrap_ci_95_upper_rad"])
    psi_uw = float(s48["phase_closure_mean_unweighted_rad"])

    return {
        "n_epochs": int(s3["n_epochs"]),
        "n_triplets": int(s3["n_total_triplets"]),
        "n_independent": int(s3["n_independent_samples"]),
        "psi_weighted_rad": psi,
        "psi_weighted_se_rad": psi_se,
        "psi_weighted_str": f"{psi:.3f}",
        "psi_weighted_se_str": f"{psi_se:.3f}",
        "rbar": float(s3["phase_closure_rbar"]),
        "rayleigh_z": float(s3["phase_closure_rayleigh_z"]),
        "rayleigh_p": float(s3["phase_closure_rayleigh_p"]),
        "v_p": float(s3["phase_closure_v_p"]),
        "ci_lo": ci_lo,
        "ci_hi": ci_hi,
        "psi_unweighted_rad": psi_uw,
        "rbar_unweighted": float(s48["phase_closure_rbar_unweighted"]),
        "rayleigh_p_unweighted": float(s48["phase_closure_rayleigh_p_unweighted"]),
        "H_mag_ns": float(s3["H_magnitude_ns"]),
        "H_sem_ns": float(s3["H_sem_ns"]),
        "H_noise_ns": float(s3["H_noise_bias_ns"]),
        "H_excess_ns": float(s3["H_excess_ns"]),
        "H_excess_t": float(s3["H_excess_t_statistic"]),
        "H_trim_ns": float(s3["H_trim_magnitude_ns"]),
        "H_trim_sem_ns": float(s3["H_trim_sem_ns"]),
        "H_trim_t": float(s3["H_trim_t_statistic"]),
        "H_signed_ns": float(s3["H_signed_mean_ns"]),
        "H_signed_t": abs(float(s3["H_signed_t_statistic"])),
        "circ_var_rad2": float(s3["phase_closure_circ_var_rad2"]),
        "noise_floor_method": s3.get("H_noise_floor_method", "mad_median_folded_normal_E_abs"),
    }


def j1603_headlines() -> Dict[str, Any]:
    s3 = load_j1603_step003()
    s48 = frame_summary("J1603-7202", "ssb")
    return {
        "n_epochs": int(s3["n_epochs"]),
        "n_triplets": int(s3["n_total_triplets"]),
        "n_independent": int(s3["n_independent_samples"]),
        "H_mag_ns": float(s3["H_magnitude_ns"]),
        "H_noise_ns": float(s3["H_noise_bias_ns"]),
        "H_excess_ns": float(s3["H_excess_ns"]),
        "H_excess_t": float(s3["H_excess_t_statistic"]),
        "rayleigh_p": float(s3["phase_closure_rayleigh_p"]),
        "rbar": float(s3["phase_closure_rbar"]),
        "circ_var_rad2": float(s3["phase_closure_circ_var_rad2"]),
        "dv_ratio": float(s3["phase_closure_dv_ratio_pc_per_kms"]),
        "psi_unweighted_rad": float(s48["phase_closure_mean_unweighted_rad"]),
        "rbar_unweighted": float(s48["phase_closure_rbar_unweighted"]),
        "rayleigh_p_unweighted": float(s48["phase_closure_rayleigh_p_unweighted"]),
        "psi_weighted_rad": float(s3["phase_closure_mean_rad"]),
    }
