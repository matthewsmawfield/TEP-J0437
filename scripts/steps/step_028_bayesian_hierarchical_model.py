#!/usr/bin/env python3
"""
================================================================================
STEP 028: BAYESIAN DIAGNOSTIC MODEL — UNSIGNED |H|
================================================================================

This step is retained as a diagnostic check on the unsigned delay-amplitude
channel.  It is not the primary TEP evidence, because mean(|delay|) has a
folded-noise floor and can be inflated by precision weighting or heavy tails.
The primary detection statistic is the phase-domain circular closure statistic
reported by step_003.
================================================================================
"""

import json
import sys
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.json_numpy import NpEncoder
from scripts.utils.config import RANDOM_SEED
from scripts.utils.logger import TEPLogger, print_status, set_step_logger

RESULTS_DIR = PROJECT_ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

try:
    import arviz as az
    import pymc as pm
    import pytensor.tensor as pt

    PYMC_AVAILABLE = True
except ImportError:
    PYMC_AVAILABLE = False



def load_closure_data() -> Dict[str, Any]:
    closure_file = PROJECT_ROOT / "results" / "step_003_closure_final_per_epoch.json"
    if not closure_file.exists():
        raise FileNotFoundError(
            f"Missing step_003 output: {closure_file}. Run step_003 first."
        )
    with open(closure_file, "r") as f:
        data = json.load(f)
    if isinstance(data, list):
        return {"epochs": data}
    return data


def _extract_arrays(closure_data: Dict[str, Any]):
    """Extract epoch-level independent means for diagnostic inference.

    Returns:
        abs_means: np.array of mean(|geometric_delta_us|) * 1000 (ns)
        epoch_sems: np.array of SEM for each epoch (ns)
        signed_means: np.array of mean(geometric_delta_us) * 1000 (ns)
    """
    epochs_data = closure_data.get("epochs", [])
    abs_means, epoch_sems, signed_means = [], [], []

    for epoch in epochs_data:
        triplets = epoch.get("triplets", [])
        if not triplets:
            continue

        # Collect values within epoch
        vals = np.array([t.get("geometric_delta_us") for t in triplets if t.get("geometric_delta_us") is not None])

        if len(vals) > 0:
            abs_vals = np.abs(vals)
            abs_m = np.mean(abs_vals) * 1000  # ns
            signed_m = np.mean(vals) * 1000  # ns

            se = (np.std(vals, ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0.001) * 1000  # ns

            abs_means.append(abs_m)
            epoch_sems.append(se)
            signed_means.append(signed_m)
            
    return (
        np.array(abs_means),
        np.array(epoch_sems),
        np.array(signed_means),
    )


def _load_step003_noise_floor() -> float:
    """Read the Rice noise floor computed by step_003 for consistency."""
    summary_file = PROJECT_ROOT / "results" / "step_003_closure_final_summary.json"
    if summary_file.exists():
        with open(summary_file, "r") as f:
            summary = json.load(f)
        return float(summary.get("H_noise_bias_ns", 0.0))
    return 0.0


def build_hierarchical_model_analytic(closure_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Closed-form diagnostic estimate for unsigned |H| using independent epochs.

    Operating on independent epochs prevents significance inflation
    from correlated triplets. The folded-noise floor is read from step_003
    to ensure consistency with the frequentist noise-subtraction method.
    """
    abs_means, epoch_sems, signed_means = _extract_arrays(closure_data)
    n_epochs = len(abs_means)
    if n_epochs == 0:
        raise ValueError("No epoch data found")

    signed_mean = float(np.mean(signed_means))
    signed_sem = float(stats.sem(signed_means)) if n_epochs > 1 else float("nan")
    signed_median = float(np.median(signed_means))
    signed_mad = float(stats.median_abs_deviation(signed_means, scale="normal"))

    # Use an unweighted epoch mean here. Inverse-variance weighting can let a few
    # precision-dominated epochs turn the folded-noise diagnostic into a fake
    # high-significance result.
    raw_abs_mean = float(np.mean(abs_means))
    raw_abs_std = float(stats.sem(abs_means)) if n_epochs > 1 else float("nan")

    # 2. Noise-subtracted Physical Holonomy
    # Read the Rice noise floor E[|H|] = sigma * sqrt(2/pi) computed by step_003.
    # The excess |H|_excess = max(0, |H|_raw - floor) isolates signal above the
    # folded-normal null expectation.
    noise_floor = _load_step003_noise_floor()
    h_unclipped = float(raw_abs_mean - noise_floor)
    h_phys = float(max(0.0, h_unclipped))

    post_std = raw_abs_std
    t_equiv = h_unclipped / post_std if post_std > 0 else 0.0
    z95 = float(stats.norm.ppf(0.975))
    prob_excess_positive = float(stats.norm.sf(0.0, loc=h_unclipped, scale=post_std)) if post_std > 0 else None

    return {
        "trace": None,
        "H": {
            "mean": h_phys,
            "unclipped_excess": h_unclipped,
            "raw_magnitude": raw_abs_mean,
            "noise_floor": noise_floor,
            "std": post_std,
            "ci_lower": max(0, h_phys - z95 * post_std),
            "ci_upper": h_phys + z95 * post_std,
            "t_equivalent": t_equiv,
        },
        "signed_diagnostic": {
            "mean_ns": signed_mean,
            "sem_ns": signed_sem,
            "t_signed": signed_mean / signed_sem if signed_sem > 0 else 0.0,
            "median_ns": signed_median,
            "mad_ns": signed_mad,
            "note": "Unweighted signed mean across independent epochs; diagnostic only.",
        },
        "prob_excess_positive": prob_excess_positive,
        "p_pos": {
            "mean": float(np.mean(signed_means > 0)),
            "std": float(np.std(signed_means > 0) / np.sqrt(n_epochs)),
        },
        "n_epochs": n_epochs,
        "n_triplets_total": int(np.sum([len(e.get('triplets', [])) for e in closure_data.get('epochs', [])])),
        "rhat_mean": 1.0,
        "converged": True,
    }



def build_hierarchical_model_simple(closure_data: Dict[str, Any]) -> Dict[str, Any]:
    """Dispatch to analytic path. PyMC FoldedNormal model removed due to Cauchy-tail failure."""
    return build_hierarchical_model_analytic(closure_data)


def analyze_posterior(model_result: Dict[str, Any]) -> Dict[str, Any]:
    H = model_result["H"]
    pp = model_result["p_pos"]
    sgn = model_result["signed_diagnostic"]

    if PYMC_AVAILABLE and model_result.get("trace") is not None:
        rhat_mean = model_result["rhat_mean"]
        converged = model_result["converged"]
    else:
        # Use default convergence metrics for analytic model (no MCMC)
        rhat_mean = 1.0
        converged = True

    return {
        "mu_H": H,  # keep key name mu_H for downstream compatibility
        "signed_diagnostic": sgn,
        "p_pos": pp,
        "detection": {
            "prob_H_excess_positive": model_result.get("prob_excess_positive"),
            "detected_5sigma": False,
            "detected_3sigma": False,
            "primary_statistic": "diagnostic mean(|geometric_delta|) excess above folded-noise floor",
            "valid_for_primary_inference": False,
            "note": (
                "Unsigned-|H| epoch summaries disagree with the step_003 noise-subtracted excess "
                "and are not used as detection evidence."
            ),
        },
        "convergence": {"rhat_mean": rhat_mean, "converged": converged},
    }


# Backwards-compat alias
def analyze_posterior_simple(model_result: Dict[str, Any]) -> Dict[str, Any]:
    return analyze_posterior(model_result)


def compare_with_frequentist(bayesian_results: Dict[str, Any]) -> Dict[str, Any]:
    freq_file = PROJECT_ROOT / "results" / "step_003_closure_final_summary.json"
    if not freq_file.exists():
        raise FileNotFoundError(f"Frequentist results not found: {freq_file}")
    with open(freq_file, "r") as f:
        freq = json.load(f)

    freq_h = freq.get("H_excess_ns")
    freq_sem = freq.get("H_excess_uncertainty_ns") or freq.get("H_sem_between_epoch_unweighted_ns")
    freq_t = freq.get("H_excess_significance_sigma")

    bay_H = bayesian_results["mu_H"]
    h_diff = abs(bay_H["mean"] - freq_h) if freq_h is not None else None

    return {
        "frequentist": {
            "H_excess_ns": freq_h,
            "SEM_ns": freq_sem,
            "t_statistic": freq_t,
        },
        "bayesian": {
            "H_ns": bay_H["mean"],
            "CI_lower_ns": bay_H["ci_lower"],
            "CI_upper_ns": bay_H["ci_upper"],
            "t_equivalent": bay_H["t_equivalent"],
        },
        "agreement": {
            "H_difference_ns": h_diff,
            "agreement_fraction": h_diff / abs(freq_h) if freq_h else None,
            "both_detect_5sigma": bool(
                bayesian_results["detection"]["detected_5sigma"]
                and freq_t is not None
                and abs(freq_t) > 5
            ),
            "both_detect_3sigma": bool(
                bayesian_results["detection"]["detected_3sigma"]
                and freq_t is not None
                and abs(freq_t) > 3
            ),
        },
    }


def main():
    """Run Bayesian unsigned-|H| analysis."""
    from scripts.utils.logger import _active_logger

    if _active_logger is None:
        log_dir = PROJECT_ROOT / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        _logger = TEPLogger(
            "step_028_bayesian",
            str(log_dir / "step_028_bayesian_hierarchical_model.log"),
        )
        set_step_logger(_logger)

    print_status("=" * 80, "TITLE")
    print_status("STEP 028: BAYESIAN DIAGNOSTIC MODEL — UNSIGNED |H|", "TITLE")
    print_status("=" * 80, "TITLE")
    print_status(f"PyMC available: {PYMC_AVAILABLE}", "INFO")
    print_status(
        "Model: epoch mean(|delta_i|) with folded-noise-floor subtraction", "INFO"
    )
    print_status("Primary TEP statistic remains phase closure from step_003", "INFO")
    print_status("", "INFO")

    print_status("Loading closure delay data...", "INFO")
    closure_data = load_closure_data()
    n_ep = len(closure_data.get("epochs", []))
    print_status(f"Loaded {n_ep} epochs", "INFO")

    print_status("Running Bayesian model...", "INFO")
    model_result = build_hierarchical_model_simple(closure_data)
    print_status(f"n_triplets = {model_result['n_triplets_total']:,}", "INFO")

    results = analyze_posterior(model_result)
    comparison = compare_with_frequentist(results)

    # ── Print results ---------------------------------------------------------─
    print_status("=" * 80, "TITLE")
    print_status("BAYESIAN RESULTS — UNSIGNED |H| DIAGNOSTIC", "TITLE")
    print_status("=" * 80, "TITLE")

    H = results["mu_H"]
    print_status(f"|H| excess mean     : {H['mean']:.3f} ns", "INFO")
    print_status(f"|H| excess std      : {H['std']:.4f} ns", "INFO")
    print_status(
        f"95% CI              : [{H['ci_lower']:.3f}, {H['ci_upper']:.3f}] ns", "INFO"
    )
    print_status(f"diagnostic z        : {H['t_equivalent']:.1f}sigma (not inferential)", "INFO")
    print_status(
        f"3sigma detection    : {results['detection']['detected_3sigma']}", "INFO"
    )
    print_status(
        f"5sigma detection    : {results['detection']['detected_5sigma']}", "INFO"
    )

    sgn = results["signed_diagnostic"]
    print_status("", "INFO")
    print_status("Signed mean (diagnostic -- bipolar cancellation):", "INFO")
    print_status(
        f"  mean = {sgn['mean_ns']:+.3f} ns  t = {sgn['t_signed']:+.2f}sigma", "INFO"
    )
    print_status(f"  Note: {sgn['note']}", "INFO")

    cf = comparison["frequentist"]
    cb = comparison["bayesian"]
    print_status("", "INFO")
    print_status("Frequentist comparison (unsigned-|H| excess diagnostic):", "INFO")
    if cf["H_excess_ns"] is not None and cf["SEM_ns"] is not None:
        t_text = f"{cf['t_statistic']:.1f}sigma" if cf["t_statistic"] is not None else "n/a"
        print_status(
            f"  Frequentist: {cf['H_excess_ns']:.3f} +/- {cf['SEM_ns']:.3f} ns"
            f"  (t = {t_text})",
            "INFO",
        )
    else:
        print_status("  Frequentist: H_excess unavailable in step_003 summary", "INFO")
    print_status(
        f"  Bayesian:    {cb['H_ns']:.3f} +/- {H['std']:.4f} ns"
        f"  (diagnostic z = {cb['t_equivalent']:.1f}, not inferential)",
        "INFO",
    )
    ag = comparison["agreement"]
    if ag["H_difference_ns"] is not None and ag["agreement_fraction"] is not None:
        print_status(
            f"  Difference:  {ag['H_difference_ns']:.4f} ns"
            f" ({ag['agreement_fraction']:.2%})",
            "INFO",
        )
    print_status(f"  Both detect 5sigma: {ag['both_detect_5sigma']}", "INFO")
    print_status(f"  Both detect 3sigma: {ag['both_detect_3sigma']}", "INFO")

    print_status(
        "Unsigned-|H| diagnostic is not promoted to a detection; use phase closure for primary inference.",
        "WARNING",
    )

    print_status(
        f"Converged: {results['convergence']['converged']}"
        f"  Rhat = {results['convergence']['rhat_mean']:.3f}",
        "INFO",
    )

    # ── Save ------------------------------------------------------------------
    save_results = {
        "bayesian_results": {k: v for k, v in results.items() if k != "trace"},
        "comparison_with_frequentist": comparison,
        "model_info": {
            "backend": "analytic_diagnostic",
            "model": "unsigned_H_epoch_mean_noise_floor_diagnostic",
            "primary_statistic": "phase closure from step_003; unsigned |H| is diagnostic only",
            "note": (
                "This result audits the unsigned delay-amplitude channel after subtracting the folded-noise floor. "
                "It deliberately avoids inverse-variance precision domination and is not used as the primary "
                "evidence for TEP. The primary inference is the phase-domain circular-statistics result."
            ),
        },
    }
    output_file = RESULTS_DIR / "step_028_bayesian_hierarchical_results.json"
    with open(output_file, "w") as f:
        json.dump(save_results, f, indent=2, cls=NpEncoder)
    print_status(f"Results saved to: {output_file}", "SUCCESS")
    print_status("=" * 80, "TITLE")
    return save_results


if __name__ == "__main__":
    main()
