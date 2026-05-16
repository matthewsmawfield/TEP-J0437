#!/usr/bin/env python3
"""
================================================================================
STEP 039: SYSTEMATIC ERROR ANALYSIS
================================================================================

Purpose: Quantify systematic uncertainties in the TEP measurement.

Systematic Effects to Analyze:
-----------------------------
- Instrumental drift (temporal stability already tested in step_034)
- ISM variations (correlation with scintillation parameters)
- Selection effects (already tested in step_038)
- Seasonal/orbital effects (already tested in step_034)
- Data quality dependence (SNR, dynamic range)
- Arclet number dependence

Expected Outcomes:
----------------
Quantify each systematic effect and add to total uncertainty budget.
If all systematic effects are small compared to statistical uncertainty,
the detection is robust.

================================================================================
"""

import json
import sys
from pathlib import Path
from typing import Optional, Any, Dict, List

import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
RESULTS_DIR = PROJECT_ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

from scripts.utils.json_numpy import NpEncoder
from scripts.utils.config import RANDOM_SEED


def load_closure_data() -> List[Dict[str, Any]]:
    """Load closure delay data."""
    closure_file = PROJECT_ROOT / "results" / "step_003_closure_final_per_epoch.json"

    if not closure_file.exists():
        raise FileNotFoundError(f"Closure delay results not found: {closure_file}")

    with open(closure_file, "r") as f:
        data = json.load(f)

    if isinstance(data, list):
        epochs = data
    elif isinstance(data, dict) and "epochs" in data:
        epochs = data["epochs"]
    else:
        raise ValueError(f"Unexpected format in {closure_file}")

    print(f"Loaded {len(epochs)} epochs with closure delay data")
    return epochs


def analyze_arclet_number_dependence(epochs: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Analyze if H depends on the number of arclets."""

    results_by_arclets = {}
    for epoch in epochs:
        n_arclets = epoch.get("n_arclets", 0)
        triplets = epoch.get("triplets", [])

        if len(triplets) < 1:
            continue

        closures = [
            triplet.get("geometric_delta_us") * 1000
            for triplet in triplets
            if triplet.get("geometric_delta_us") is not None
        ]
        if not closures:
            continue
        epoch_H = np.mean(closures)

        if n_arclets not in results_by_arclets:
            results_by_arclets[n_arclets] = []
        results_by_arclets[n_arclets].append(epoch_H)

    # Compute statistics for each arclet number
    arclet_stats = {}
    for n_arclets, H_values in sorted(results_by_arclets.items()):
        if len(H_values) < 3:
            continue
        arclet_stats[n_arclets] = {
            "n_epochs": len(H_values),
            "mean_H_ns": float(np.mean(H_values)),
            "std_H_ns": float(np.std(H_values, ddof=1)),
            "sem_ns": float(np.std(H_values, ddof=1) / np.sqrt(len(H_values))),
        }

    # Test for correlation
    if len(arclet_stats) >= 3:
        arclet_numbers = np.array(list(arclet_stats.keys()))
        H_means = np.array([s["mean_H_ns"] for s in arclet_stats.values()])
        if len(arclet_numbers) > 1:
            corr, p_value = stats.pearsonr(arclet_numbers, H_means)
        else:
            raise ValueError(
                "Cannot compute arclet number correlation: insufficient data (need > 1 epoch). "
                "This test requires multiple epochs with varying arclet counts."
            )
    else:
        corr, p_value = 0.0, 1.0

    return {
        "arclet_number_stats": arclet_stats,
        "correlation_with_arclets": float(corr),
        "p_value": float(p_value),
        "significant": bool(p_value < 0.05),
        "interpretation": "H depends on arclet number"
        if p_value < 0.05
        else "H independent of arclet number",
    }


def analyze_snr_dependence(epochs: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Analyze if H depends on SNR.

    NOTE: The triplet `snr` field in step_003 is an independent cross-term SNR,
    not the circular closure-SNR statistic. However, dedicated SNR analyses are
    already performed in steps 025 and 028 using arclet-level SNR observables,
    so this systematic-error summary defers to those authoritative outputs.

    Instead, the analysis should correlate H with truly independent variables like:
    - Number of arclets
    - Number of triplets
    - Dedicated arclet-level SNR diagnostics from steps 025/028
    """

    print(f"   NOTE: Triplet `snr` is independent cross-term SNR, not closure SNR")
    print(f"   NOTE: Dedicated SNR diagnostics already exist in steps 025 and 028")
    print(f"   NOTE: Skipping duplicate SNR dependence analysis in this summary step")

    return {
        "note": "Triplet snr in step_003 is independent cross-term SNR; dedicated SNR diagnostics are handled in steps 025 and 028",
        "reason": "Avoid duplicating a more rigorous dedicated SNR analysis in the systematic-error summary",
        "correlation_with_snr": None,
        "interpretation": "Analysis deferred to steps 025 and 028",
    }


def analyze_triplet_number_dependence(epochs: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Analyze if H depends on the number of triplets."""

    epoch_H_values = []
    epoch_n_triplets = []

    for epoch in epochs:
        triplets = epoch.get("triplets", [])
        if len(triplets) < 1:
            continue

        closures = [
            triplet.get("geometric_delta_us") * 1000
            for triplet in triplets
            if triplet.get("geometric_delta_us") is not None
        ]
        if not closures:
            continue
        epoch_H = np.mean(closures)

        epoch_H_values.append(epoch_H)
        epoch_n_triplets.append(len(triplets))

    epoch_H_values = np.array(epoch_H_values)
    epoch_n_triplets = np.array(epoch_n_triplets)

    # Test for correlation
    corr, p_value = stats.pearsonr(epoch_n_triplets, epoch_H_values)

    return {
        "correlation_with_n_triplets": float(corr),
        "p_value": float(p_value),
        "significant": bool(p_value < 0.05),
        "interpretation": "H depends on number of triplets"
        if p_value < 0.05
        else "H independent of number of triplets",
    }


def analyze_pixel_discretization_error(epochs):
    """
    Calculate systematic from tau-axis pixel discretisation.

    Each cross-term peak is located to sub-pixel precision via parabolic
    interpolation. The tau grid spacing d_tau is derived from sigma_us measurements.
    The discretization error is bounded by the theoretical limit of parabolic interpolation.
    """
    sigma_vals = []
    for epoch in epochs:
        for t in epoch.get("triplets", []):
            sv = t.get("sigma_us", None)
            if sv is not None and sv > 0:
                sigma_vals.append(sv * 1000)  # ns

    if not sigma_vals:
        return {"pixel_discretization_ns": 0.0, "note": "No sigma_us data"}

    # sigma_us = d_tau * sqrt(3)  (from step_003 error model)
    # d_tau = sigma_us / sqrt(3)
    median_sigma = float(np.median(sigma_vals))
    d_tau_ns = median_sigma / np.sqrt(3)

    # Discretisation error: theoretical upper bound for parabolic interpolation
    # Error ≤ 0.25 * f''(ξ) * h² where h is grid spacing (numerical analysis theory)
    # For parabolic fit to symmetric peaks, this bounds to ≤ 0.1 * d_tau
    discretization_ns = 0.1 * d_tau_ns

    return {
        "median_sigma_us_ns": median_sigma,
        "d_tau_ns": d_tau_ns,
        "pixel_discretization_ns": discretization_ns,
        "note": "Conservative theoretical upper bound for parabolic sub-pixel interpolation: ≤ 0.1 × d_tau. Actual contribution to phase-angle variance is subsumed in the empirical circular SE.",
    }


def analyze_thermal_noise_contribution(epochs):
    """
    Estimate systematic contribution from thermal noise floor.

    The per-triplet measurement uncertainty sigma_us captures thermal noise.
    The systematic contribution to the MEAN closure delay is sigma_us / sqrt(N),
    where N is the total number of triplets. This is already captured by the
    statistical SEM; here we estimate the ADDITIONAL systematic from any
    epoch-dependent noise-floor variation.
    """
    epoch_sigmas = []
    for epoch in epochs:
        triplets = epoch.get("triplets", [])
        if not triplets:
            continue
        sigs = [t.get("sigma_us", 0.003) * 1000 for t in triplets]
        epoch_sigmas.append(np.mean(sigs))

    if len(epoch_sigmas) < 2:
        return {"thermal_noise_systematic_ns": 0.0, "note": "Insufficient epochs"}

    epoch_sigmas = np.array(epoch_sigmas)
    # Epoch-to-epoch variation in the noise floor
    sigma_noise_variation = float(np.std(epoch_sigmas, ddof=1))
    # Systematic contribution: noise-floor variation propagated through mean
    thermal_systematic = sigma_noise_variation / np.sqrt(len(epoch_sigmas))

    return {
        "mean_sigma_per_triplet_ns": float(np.mean(epoch_sigmas)),
        "std_sigma_across_epochs_ns": sigma_noise_variation,
        "thermal_noise_systematic_ns": thermal_systematic,
        "n_epochs": len(epoch_sigmas),
        "note": "Epoch-to-epoch noise-floor variation propagated to mean",
    }


def analyze_bandpass_calibration_error(epochs):
    """
    Calculate systematic from bandpass calibration uncertainty.

    The dynamic spectra are bandpass-normalised by median per-channel flux.
    Calibration systematic is bounded by instrumental specifications.
    Parkes/UWL calibration accuracy: < 1% residual after bandpass correction.
    """
    all_H = []
    for epoch in epochs:
        for t in epoch.get("triplets", []):
            gd = t.get("geometric_delta_us")
            if gd is not None:
                all_H.append(abs(gd * 1000))  # ns

    if not all_H:
        return {"bandpass_calibration_ns": 0.0, "note": "No data"}

    median_H = float(np.median(all_H))
    # Calibration systematic: 1% of median |H| (Parkes/UWL specification)
    calibration_systematic = 0.01 * median_H

    return {
        "median_abs_H_ns": median_H,
        "calibration_fraction": 0.01,
        "bandpass_calibration_ns": calibration_systematic,
        "note": "1% of |median H|, consistent with Parkes bandpass accuracy. Closure observables largely cancel common-mode bandpass residuals; this is a conservative upper bound.",
    }


def compute_circular_psi_stats(epochs: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Compute circular statistics on epoch-level Phase Closure ψ angles.

    Each epoch contributes its mean phase_closure_rad across triplets.
    Returns circular mean, R_bar, Rayleigh Z, V-test, and bootstrap CI.
    """
    epoch_psi = []
    for epoch in epochs:
        triplets = epoch.get("triplets", [])
        psi_vals = [
            t.get("phase_closure_rad", None)
            for t in triplets
            if t.get("phase_closure_rad") is not None
        ]
        if len(psi_vals) >= 3:
            # Epoch-level circular mean
            z = np.sum(np.exp(1j * np.array(psi_vals)))
            epoch_psi.append(float(np.angle(z / len(psi_vals))))

    if len(epoch_psi) < 3:
        return {
            "n_epochs": len(epoch_psi),
            "psi_mean_rad": 0.0,
            "r_bar": 0.0,
            "rayleigh_z": 0.0,
            "rayleigh_p": 1.0,
            "v_stat": 0.0,
            "v_p": 1.0,
            "bootstrap_se": float("inf"),
            "bootstrap_ci_95": (float("-inf"), float("inf")),
        }

    angles = np.array(epoch_psi)
    n = len(angles)
    z_sum = np.sum(np.exp(1j * angles))
    psi_mean = float(np.angle(z_sum / n))
    r_bar = float(np.abs(z_sum) / n)
    rayleigh_z = float(2.0 * n * r_bar**2)
    rayleigh_p = float(stats.chi2.sf(rayleigh_z, 2))
    v_stat = float(np.sqrt(2.0 * n) * r_bar * np.cos(psi_mean))
    v_p = float(2.0 * stats.norm.sf(abs(v_stat)))

    # Bootstrap CI
    rng = np.random.RandomState(RANDOM_SEED)
    boot_means = []
    for _ in range(10000):
        idx = rng.choice(n, size=n, replace=True)
        z_boot = np.sum(np.exp(1j * angles[idx]))
        boot_means.append(float(np.angle(z_boot / n)))
    boot_means = np.array(boot_means)
    centered = (boot_means - psi_mean + np.pi) % (2.0 * np.pi) - np.pi
    boot_se = float(np.std(centered, ddof=1))
    ci_low = float((psi_mean + np.percentile(centered, 2.5) + np.pi) % (2.0 * np.pi) - np.pi)
    ci_high = float((psi_mean + np.percentile(centered, 97.5) + np.pi) % (2.0 * np.pi) - np.pi)

    return {
        "n_epochs": n,
        "psi_mean_rad": psi_mean,
        "r_bar": r_bar,
        "rayleigh_z": rayleigh_z,
        "rayleigh_p": rayleigh_p,
        "v_stat": v_stat,
        "v_p": v_p,
        "bootstrap_se": boot_se,
        "bootstrap_ci_95": (ci_low, ci_high),
    }


def compile_systematic_error_budget(
    arclet_dep: Dict,
    snr_dep: Dict,
    triplet_dep: Dict,
    selection_bias: Optional[Dict],
    pixel_disc=None,
    thermal_noise=None,
    bandpass_cal=None,
) -> Dict[str, Any]:
    """Compile systematic error budget from all analyses.

    The budget is split into:
    (a) Core metrology systematics — physical uncertainties in delay measurement
        (pixel discretization, thermal noise floor, bandpass calibration).
    (b) Threshold-robustness allowance — variation in signed-mean across SNR/arclet
        cuts from step_041. This is a diagnostic quantity expected near zero due
        to bipolar cancellation, so its spread reflects estimator noise rather
        than a bias in the primary Phase Closure ψ.
    """

    # Load standard statistical error
    summary_file = PROJECT_ROOT / "results" / "step_003_closure_final_summary.json"
    with open(summary_file, "r") as f:
        summary = json.load(f)

    statistical_error = summary["H_sem_ns"]

    # --- Core metrology systematics ---
    core_systematic_errors = []
    core_systematic_sources = []

    if arclet_dep["significant"]:
        arclet_stats = arclet_dep["arclet_number_stats"]
        H_values = [s["mean_H_ns"] for s in arclet_stats.values()]
        core_systematic_errors.append(np.std(H_values))
        core_systematic_sources.append("arclet_number")

    if triplet_dep["significant"]:
        triplet_epoch_H = []
        if selection_bias is not None:
            for epoch in selection_bias.get("arclet_threshold_analysis", []):
                if "mean_H_ns" in epoch:
                    triplet_epoch_H.append(epoch["mean_H_ns"])
        if len(triplet_epoch_H) >= 2:
            core_systematic_errors.append(float(np.std(triplet_epoch_H, ddof=1)))
        else:
            core_systematic_errors.append(0.0)
        core_systematic_sources.append("triplet_number")

    # Physics-motivated systematics — always included regardless of significance guards
    if pixel_disc is not None and pixel_disc.get("pixel_discretization_ns", 0) > 0:
        core_systematic_errors.append(pixel_disc["pixel_discretization_ns"])
        core_systematic_sources.append("pixel_discretization")

    if (
        thermal_noise is not None
        and thermal_noise.get("thermal_noise_systematic_ns", 0) > 0
    ):
        core_systematic_errors.append(thermal_noise["thermal_noise_systematic_ns"])
        core_systematic_sources.append("thermal_noise_floor")

    if bandpass_cal is not None and bandpass_cal.get("bandpass_calibration_ns", 0) > 0:
        core_systematic_errors.append(bandpass_cal["bandpass_calibration_ns"])
        core_systematic_sources.append("bandpass_calibration")

    core_systematic = np.sqrt(sum(e**2 for e in core_systematic_errors)) if core_systematic_errors else 0.0

    # --- Threshold-robustness allowance (signed-mean diagnostic) ---
    threshold_means = []
    if selection_bias is not None:
        for block_name in ["snr_threshold_analysis", "arclet_threshold_analysis"]:
            for result in selection_bias.get(block_name, []):
                if "mean_H_ns" in result:
                    threshold_means.append(float(result["mean_H_ns"]))

    threshold_conclusions = []
    if selection_bias is not None:
        threshold_conclusions = selection_bias.get("conclusions", [])
    threshold_evidence_mixed = any(
        "limited or mixed" in conclusion.lower() or "false" in conclusion.lower()
        for conclusion in threshold_conclusions
    )
    threshold_allowance = 0.0
    if threshold_evidence_mixed and len(threshold_means) >= 2:
        threshold_allowance = float(np.std(threshold_means, ddof=1))

    # Total systematic including threshold allowance
    total_systematic = np.sqrt(core_systematic**2 + threshold_allowance**2)
    total_error = np.sqrt(statistical_error**2 + total_systematic**2)

    return {
        "statistical_error_ns": float(statistical_error),
        "core_systematic_error_ns": float(core_systematic),
        "core_systematic_sources": core_systematic_sources,
        "core_systematic_errors_ns": [float(e) for e in core_systematic_errors],
        "threshold_robustness_ns": float(threshold_allowance),
        "total_systematic_error_ns": float(total_systematic),
        "total_error_ns": float(total_error),
        "systematic_fraction": float(total_systematic / statistical_error)
        if statistical_error > 0
        else 0.0,
        "core_fraction": float(core_systematic / statistical_error)
        if statistical_error > 0
        else 0.0,
        "interpretation": "Core metrology systematic errors are small"
        if core_systematic < statistical_error
        else "Core metrology systematic errors exceed statistical error",
    }


def main():
    """Run systematic error analysis."""
    import sys

    sys.path.insert(0, str(PROJECT_ROOT))
    from scripts.utils.logger import (
        TEPLogger,
        _active_logger,
        print_status,
        set_step_logger,
    )

    if _active_logger is None:
        log_dir = PROJECT_ROOT / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        _logger = TEPLogger("step_042", str(log_dir / "step_042_systematic_errors.log"))
        set_step_logger(_logger)

    print_status("=" * 80, "INFO")
    print_status("STEP 042: SYSTEMATIC ERROR ANALYSIS", "TITLE")
    print_status("=" * 80, "INFO")
    print_status("", "INFO")
    print_status(
        "Purpose: Quantify systematic uncertainties in TEP measurement", "INFO"
    )
    print_status("", "INFO")

    # Load data
    print_status("1. LOADING CLOSURE DATA:", "PROCESS")
    epochs = load_closure_data()
    print_status("", "INFO")

    # Analyze arclet number dependence
    print_status("2. ARCLET NUMBER DEPENDENCE:", "PROCESS")
    arclet_dep = analyze_arclet_number_dependence(epochs)
    print_status(
        f"   Correlation: {arclet_dep['correlation_with_arclets']:.3f}", "INFO"
    )
    print_status(f"   Significant: {arclet_dep['significant']}", "INFO")
    print_status(f"   {arclet_dep['interpretation']}", "INFO")
    print_status("", "INFO")

    # Analyze SNR dependence
    print_status("3. SNR DEPENDENCE:", "PROCESS")
    snr_dep = analyze_snr_dependence(epochs)
    if snr_dep["correlation_with_snr"] is not None:
        print_status(f"   Correlation: {snr_dep['correlation_with_snr']:.3f}", "INFO")
        print_status(f"   Significant: {snr_dep['significant']}", "INFO")
    print_status(f"   {snr_dep['interpretation']}", "INFO")
    print_status("", "INFO")

    # Analyze triplet number dependence
    print_status("4. TRIPLET NUMBER DEPENDENCE:", "PROCESS")
    triplet_dep = analyze_triplet_number_dependence(epochs)
    print_status(
        f"   Correlation: {triplet_dep['correlation_with_n_triplets']:.3f}", "INFO"
    )
    print_status(f"   Significant: {triplet_dep['significant']}", "INFO")
    print_status(f"   {triplet_dep['interpretation']}", "INFO")
    print_status("", "INFO")

    # Load selection bias results
    print_status("5. SELECTION BIAS (from step_041):", "PROCESS")
    # step_041 writes to results/step_041_selection_bias_results.json
    selection_bias_file = (
        PROJECT_ROOT / "results" / "step_041_selection_bias_results.json"
    )
    # Allow graceful skip if selection bias analysis didn't run
    if not selection_bias_file.exists():
        print_status("   Selection bias results not found - skipping this section.", "WARNING")
        print_status("   Run step_041_selection_bias_analysis.py first to generate selection bias analysis.", "INFO")
        selection_bias = None
    else:
        with open(selection_bias_file, "r") as f:
            selection_bias = json.load(f)
        for conclusion in selection_bias.get("conclusions", []):
            print_status(f"   {conclusion}", "INFO")
        print_status("", "INFO")

    # Physics-motivated systematics
    print_status("7. PIXEL DISCRETIZATION ERROR:", "PROCESS")
    pixel_disc = analyze_pixel_discretization_error(epochs)
    print_status(
        f"   d_tau = {pixel_disc.get('d_tau_ns', 0):.3f} ns, discretization systematic = {pixel_disc.get('pixel_discretization_ns', 0):.4f} ns",
        "INFO",
    )
    print_status("", "INFO")

    print_status("8. THERMAL NOISE CONTRIBUTION:", "PROCESS")
    thermal_noise = analyze_thermal_noise_contribution(epochs)
    print_status(
        f"   systematic = {thermal_noise.get('thermal_noise_systematic_ns', 0):.4f} ns",
        "INFO",
    )
    print_status("", "INFO")

    print_status("9. BANDPASS CALIBRATION ERROR:", "PROCESS")
    bandpass_cal = analyze_bandpass_calibration_error(epochs)
    print_status(
        f"   systematic = {bandpass_cal.get('bandpass_calibration_ns', 0):.4f} ns",
        "INFO",
    )
    print_status("", "INFO")

    # Compile error budget
    print_status("6. SYSTEMATIC ERROR BUDGET:", "PROCESS")
    error_budget = compile_systematic_error_budget(
        arclet_dep,
        snr_dep,
        triplet_dep,
        selection_bias,
        pixel_disc=pixel_disc,
        thermal_noise=thermal_noise,
        bandpass_cal=bandpass_cal,
    )
    print_status(
        f"   Statistical error: {error_budget['statistical_error_ns']:.3f} ns", "INFO"
    )
    print_status(
        f"   Systematic error: {error_budget['total_systematic_error_ns']:.3f} ns",
        "INFO",
    )
    print_status(f"   Total error: {error_budget['total_error_ns']:.3f} ns", "INFO")
    print_status(
        f"   Systematic fraction: {error_budget['systematic_fraction']:.3f}", "INFO"
    )
    print_status(f"   {error_budget['interpretation']}", "INFO")
    print_status("", "INFO")

    # Circular-statistics check on ψ
    print_status("10. PHASE CLOSURE ψ CIRCULAR STATISTICS:", "PROCESS")
    psi_stats = compute_circular_psi_stats(epochs)
    print_status(
        f"   Epochs with ψ: {psi_stats['n_epochs']}, ψ_mean = {psi_stats['psi_mean_rad']:.3f} rad, "
        f"R_bar = {psi_stats['r_bar']:.3f}", "INFO"
    )
    print_status(
        f"   Rayleigh Z = {psi_stats['rayleigh_z']:.2f} (p = {psi_stats['rayleigh_p']:.2e}), "
        f"V = {psi_stats['v_stat']:.2f} (p = {psi_stats['v_p']:.2e})", "INFO"
    )
    print_status(
        f"   Bootstrap SE = {psi_stats['bootstrap_se']:.3f} rad, "
        f"95% CI = [{psi_stats['bootstrap_ci_95'][0]:.3f}, {psi_stats['bootstrap_ci_95'][1]:.3f}] rad",
        "INFO",
    )
    print_status("   NOTE: ψ is extracted from complex cross-term phases prior to SNR weighting; threshold cuts do not directly affect phase measurement", "INFO")
    print_status("", "INFO")

    # Total systematic error is sqrt(sum(e^2))
    # Total error is sqrt(stat^2 + sys^2)
    sys_fraction = error_budget["systematic_fraction"]
    core_fraction = error_budget["core_fraction"]

    # Compile results
    full_results = {
        "arclet_number_dependence": arclet_dep,
        "snr_dependence": snr_dep,
        "triplet_number_dependence": triplet_dep,
        "selection_bias": selection_bias,
        "pixel_discretization": pixel_disc,
        "thermal_noise": thermal_noise,
        "bandpass_calibration": bandpass_cal,
        "psi_circular_statistics": psi_stats,
        "error_budget": error_budget,
        "conclusions": [
            f"Arclet number dependence: {arclet_dep['interpretation']}",
            f"SNR dependence: {snr_dep['interpretation']}",
            f"Triplet number dependence: {triplet_dep['interpretation']}",
            f"Core metrology systematics: {error_budget['interpretation']}",
            f"Phase Closure ψ: Rayleigh Z = {psi_stats['rayleigh_z']:.2f} (p = {psi_stats['rayleigh_p']:.2e}), circular SE = {psi_stats['bootstrap_se']:.3f} rad",
            "Threshold-robustness concerns apply to signed-mean diagnostic (expected near zero by bipolar cancellation), not to primary Phase Closure ψ",
        ],
        "implications": {
            "core_systematic_ratio": f"{core_fraction:.2f}",
            "total_systematic_ratio": f"{sys_fraction:.2f}",
            "core_systematic_ns": f"{error_budget['core_systematic_error_ns']:.3f}",
            "threshold_robustness_ns": f"{error_budget['threshold_robustness_ns']:.3f}",
            "psi_circular_se_rad": f"{psi_stats['bootstrap_se']:.3f}",
            "measurement_status": (
                "Core metrology systematic exceeds the delay-domain statistical error; "
                "this weakens delay-amplitude claims but does not directly control the "
                "phase-closure ψ circular-statistics inference."
                if core_fraction >= 1.0
                else "Core metrology systematic is sub-statistical; total systematic including threshold allowance is larger but does not affect circular-statistics inference on ψ"
            ),
        },
    }


    # Save results
    output_file = RESULTS_DIR / "step_042_systematic_error_results.json"
    with open(output_file, "w") as f:
        json.dump(full_results, f, indent=2, cls=NpEncoder)

    print_status("=" * 80, "INFO")
    print_status("CONCLUSIONS:", "TITLE")
    print_status("=" * 80, "INFO")
    for conclusion in full_results["conclusions"]:
        print_status(f"  * {conclusion}", "INFO")
    print_status("", "INFO")
    print_status("=" * 80, "INFO")
    print_status("IMPLICATIONS:", "TITLE")
    print_status("=" * 80, "INFO")
    for key, value in full_results["implications"].items():
        print_status(f"  {key}: {value}", "INFO")
    print_status("", "INFO")
    print_status("=" * 80, "INFO")
    print_status(f"Results saved to: {output_file}", "SUCCESS")
    print_status("=" * 80, "INFO")

    return full_results


if __name__ == "__main__":
    main()
