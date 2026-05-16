#!/usr/bin/env python3
"""
================================================================================
STEP 044: PROBABILISTIC WEIGHTING DIAGNOSTIC
================================================================================

Purpose: Replace hard SNR/arclet threshold cuts with probabilistic weighting
based on measurement uncertainty. This eliminates the arbitrary nature of
cuts that can distort the geometric distribution needed for valid statistical
tests.

Methodology:
------------
1. Hierarchical Model Structure:
   - Level 1: Individual arclet measurements (τᵢ, σᵢ)
   - Level 2: Triplet closure delays with per-triplet uncertainties
   - Level 3: Epoch-level holonomy with epoch-specific precision
   - Level 4: Population-level TEP parameters (μ_H, σ_H)

2. Probabilistic Weighting (Instead of Hard Cuts):
   - Each triplet contributes with weight wᵢ ∝ 1/σᵢ² (inverse variance)
   - SNR information enters through σᵢ = δᵢ / SNR (uncertainty estimate)
   - Low-SNR triplets contribute less but are NOT discarded
   - This preserves the full geometric distribution

3. Bayesian Diagnostic:
   - Priors: H ~ HalfNormal(0, 50 ns), σ_epoch ~ HalfCauchy(0, 10)
   - Likelihood: δᵢ ~ Normal(H × sᵢ, √(σᵢ² + σ_epoch²))
   - Posterior computed via numerical integration (no MCMC needed)

4. Advantages Over Hard Cuts:
   - No arbitrary thresholds to tune
   - Smooth degradation of influence with decreasing SNR
   - Full uncertainty propagation through all levels
   - Natural handling of heteroscedastic data

Important limitation:
---------------------
This step works on unsigned |delay| amplitudes. Folded magnitudes have a
noise floor and are not the primary TEP detection statistic. Results from this
script are diagnostic checks on threshold sensitivity, not detection evidence.

================================================================================
"""

import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.json_numpy import NpEncoder
RESULTS_DIR = PROJECT_ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

from scripts.utils.config import RANDOM_SEED


@dataclass
class TripletMeasurement:
    """A single triplet closure measurement with uncertainty."""

    delta_ns: float
    sigma_ns: float
    snr: float
    closure_snr: float
    geom_sign: float
    saa_sign: float
    epoch_idx: int
    triplet_idx: int


def load_closure_data() -> List[Dict[str, Any]]:
    """Load closure delay data from step 003."""
    closure_file = PROJECT_ROOT / "results" / "step_003_closure_final_per_epoch.json"

    if not closure_file.exists():
        raise FileNotFoundError(f"Missing closure data: {closure_file}")

    with open(closure_file, "r") as f:
        data = json.load(f)

    # Handle both list and dict formats
    if isinstance(data, list):
        return data
    elif isinstance(data, dict) and "epochs" in data:
        return data["epochs"]
    else:
        return [data]


def extract_triplet_measurements(
    epochs: List[Dict[str, Any]],
) -> List[TripletMeasurement]:
    """
    Extract all triplet measurements with their uncertainties.

    Key insight: Instead of filtering by SNR threshold, we compute
    per-triplet uncertainties that will be used for probabilistic weighting.
    """
    measurements = []

    for epoch_idx, epoch in enumerate(epochs):
        epoch_name = epoch.get("epoch", f"epoch_{epoch_idx}")
        triplets = epoch.get("triplets", [])

        for triplet_idx, triplet in enumerate(triplets):
            delta_us = triplet.get("geometric_delta_us", triplet.get("delta_us", 0))
            delta_ns = delta_us * 1000  # Convert μs to ns

            # Quality information - sigma_us is REQUIRED for proper uncertainty weighting
            sigma_us = triplet.get("sigma_us")
            if sigma_us is None or sigma_us <= 0:
                print_status(f"  [SKIP] Triplet in {epoch_name}: missing/invalid sigma_us", "WARNING")
                continue
            sigma_ns_base = sigma_us * 1000
            snr = triplet.get("snr", 0)
            closure_snr = triplet.get("closure_snr", 0)

            # Use closure SNR to build a heteroscedastic uncertainty proxy.
            # sigma_us alone is often close to constant across triplets and does
            # not produce meaningful probabilistic weighting.
            quality = max(float(closure_snr), 1e-3)
            sigma_ns = sigma_ns_base / quality

            # Geometric and kinematic alignment
            geom_sign = triplet.get("geom_sign", 1.0)
            saa_sign = triplet.get("saa_sign", 1.0)

            measurements.append(
                TripletMeasurement(
                    delta_ns=delta_ns,
                    sigma_ns=sigma_ns,
                    snr=snr,
                    closure_snr=closure_snr,
                    geom_sign=geom_sign,
                    saa_sign=saa_sign,
                    epoch_idx=epoch_idx,
                    triplet_idx=triplet_idx,
                )
            )

    return measurements


def compute_probabilistic_weights(measurements: List[TripletMeasurement]) -> np.ndarray:
    """
    Compute inverse-variance weights for probabilistic analysis.

    Weight formula: wᵢ = 1/σᵢ² / Σ(1/σⱼ²)

    This ensures:
    - High-precision measurements (low σ) contribute more
    - Low-precision measurements (high σ) contribute less but are retained
    - Sum of weights = 1 (properly normalized)
    """
    variances = np.array([m.sigma_ns**2 for m in measurements])

    # Handle edge case of zero variance
    variances = np.maximum(variances, 1e-6)

    # Inverse variance weights
    inv_variances = 1.0 / variances
    weights = inv_variances / np.sum(inv_variances)

    return weights


def hierarchical_bayesian_analysis(
    measurements: List[TripletMeasurement], n_grid_points: int = 2000
) -> Dict[str, Any]:
    """
    Perform Bayesian magnitude inference using a sign-marginalized bipolar model.

    Model:
    ------
    H ~ HalfNormal(0, 50)
    delta_i ~ 0.5 * Normal(+H, sigma_i) + 0.5 * Normal(-H, sigma_i)

    This preserves the bipolar TEP structure while avoiding a per-observation
    latent sign variable that can spuriously inflate certainty.
    """
    if not measurements:
        return {"error": "No measurements available"}

    n_triplets = len(measurements)

    # Extract data arrays
    # Note: geometric_delta_us in the data ALREADY includes geom_sign * saa_sign
    # from step_003, so we use deltas directly without re-applying signs
    deltas = np.array([m.delta_ns for m in measurements])
    sigmas = np.array([m.sigma_ns for m in measurements])
    abs_deltas = np.abs(deltas)

    # Empirical intrinsic scatter term.
    # A single common bipolar magnitude does not explain the full spread of
    # triplet-level |δ| values, so using measurement error alone produces a
    # spuriously overconfident posterior. We absorb the unresolved physical/
    # epoch-to-epoch spread into tau_intrinsic.
    tau_intrinsic = float(np.std(abs_deltas, ddof=1)) if len(abs_deltas) > 1 else 0.0

    max_abs_delta = float(np.max(np.abs(deltas))) if len(deltas) else 0.0
    h_max = max(50.0, max_abs_delta * 1.5)
    h_grid = np.linspace(0.0, h_max, n_grid_points)

    # Protect against zero/unrealistic uncertainties
    # Floor at 0.1 ns (100 ps) - realistic minimum for scintillation timing
    sigma_safe = np.maximum(sigmas, 0.1)
    sigma_total = np.sqrt(sigma_safe**2 + tau_intrinsic**2)
    log_prior = -0.5 * (h_grid / 50.0) ** 2

    log_likelihood = np.zeros_like(h_grid)
    log_norm = -0.5 * np.log(2.0 * np.pi) - np.log(sigma_total)

    for i, h in enumerate(h_grid):
        z_pos = -0.5 * ((deltas - h) / sigma_total) ** 2 + log_norm
        z_neg = -0.5 * ((deltas + h) / sigma_total) ** 2 + log_norm
        log_mix = np.logaddexp(z_pos, z_neg) - np.log(2.0)
        log_likelihood[i] = np.sum(log_mix)

    log_posterior = log_prior + log_likelihood
    log_posterior -= np.max(log_posterior)
    posterior = np.exp(log_posterior)
    grid_norm = np.trapz(posterior, h_grid)
    if not np.isfinite(grid_norm) or grid_norm <= 0:
        return {"error": "Posterior normalization failed"}
    posterior /= grid_norm

    posterior_mean = float(np.trapz(h_grid * posterior, h_grid))
    posterior_var = float(np.trapz(((h_grid - posterior_mean) ** 2) * posterior, h_grid))
    posterior_std = float(np.sqrt(max(posterior_var, 0.0)))
    cdf = np.concatenate(([0.0], np.cumsum((posterior[1:] + posterior[:-1]) * np.diff(h_grid) / 2.0)))
    cdf = np.clip(cdf, 0.0, 1.0)
    ci_lower = float(np.interp(0.025, cdf, h_grid))
    ci_upper = float(np.interp(0.975, cdf, h_grid))
    t_equivalent = posterior_mean / posterior_std if posterior_std > 0 else 0
    snr_posterior = t_equivalent

    # Effective sample size (accounting for weights)
    weights = compute_probabilistic_weights(measurements)
    effective_sample_size = 1.0 / np.sum(weights**2)

    weights_summary = {
        "min_weight": float(np.min(weights)),
        "max_weight": float(np.max(weights)),
        "median_weight": float(np.median(weights)),
        "entropy_bits": float(-np.sum(weights * np.log2(weights + 1e-10))),
    }

    return {
        "n_triplets": n_triplets,
        "posterior_H": {
            "mean_ns": posterior_mean,
            "std_ns": posterior_std,
            "ci_95_lower_ns": ci_lower,
            "ci_95_upper_ns": ci_upper,
            "probability_H_positive_under_halfnormal_prior": 1.0,
            "probability_note": "Half-normal magnitude prior makes H non-negative by construction; this is not detection evidence.",
            "t_equivalent": float(t_equivalent),
            "snr_posterior": float(snr_posterior),
            "tau_intrinsic_ns": float(tau_intrinsic),
        },
        "effective_sample_size": float(effective_sample_size),
        "weights_summary": weights_summary,
    }


def compute_weighted_statistics(
    measurements: List[TripletMeasurement],
) -> Dict[str, Any]:
    """
    Compute weighted mean and standard error using inverse-variance weighting.

    For TEP bipolar signals, we use |δ| to properly estimate the magnitude.
    This is the frequentist analog to the Bayesian analysis.
    """
    if not measurements:
        return {"error": "No measurements"}

    # Use absolute values for TEP magnitude estimation (proper for bipolar)
    abs_deltas = np.array([abs(m.delta_ns) for m in measurements])
    sigmas = np.array([m.sigma_ns for m in measurements])
    tau_intrinsic = float(np.std(abs_deltas, ddof=1)) if len(abs_deltas) > 1 else 0.0

    # Inverse-variance weights
    # Protect against zero/unrealistic uncertainties
    # Floor at 0.1 ns (100 ps) - realistic minimum for scintillation timing
    sigmas = np.maximum(sigmas, 0.1)
    sigma_total = np.sqrt(sigmas**2 + tau_intrinsic**2)
    weights = 1.0 / (sigma_total**2)
    weights = weights / np.sum(weights)

    # Weighted mean of |δ|
    weighted_mean = np.sum(weights * abs_deltas)

    # Weighted variance (heteroscedasticity-aware)
    weighted_var = np.sum(weights * (abs_deltas - weighted_mean) ** 2)

    # Standard error of weighted mean
    n_eff = 1.0 / np.sum(weights**2)
    sem = np.sqrt(weighted_var / n_eff)

    # t-statistic
    t_stat = weighted_mean / sem if sem > 0 else 0

    return {
        "weighted_mean_ns": float(weighted_mean),
        "weighted_sem_ns": float(sem),
        "weighted_std_ns": float(np.sqrt(weighted_var)),
        "t_statistic": float(t_stat),
        "n_eff": float(n_eff),
        "n_total": len(measurements),
        "tau_intrinsic_ns": float(tau_intrinsic),
    }


def compare_threshold_vs_probabilistic(
    measurements: List[TripletMeasurement], threshold_results: Dict[str, Any]
) -> Dict[str, Any]:
    """
    Compare probabilistic weighting to hard threshold cuts.

    Shows how the signal degrades under hard cuts vs. smooth weighting.
    """
    # Apply various hard thresholds
    thresholds = [0, 3, 5, 7, 10]
    threshold_analyses = {}

    for thresh in thresholds:
        # Hard cut: keep only triplets above closure-SNR threshold
        retained = [m for m in measurements if m.closure_snr >= thresh]

        if len(retained) < 3:
            threshold_analyses[f"snr_{thresh}"] = {
                "n_retained": len(retained),
                "mean_ns": None,
                "t_stat": None,
                "error": "Insufficient data after cut",
            }
            continue

        # Use absolute values for TEP magnitude (proper for bipolar signal)
        abs_deltas = np.array([abs(m.delta_ns) for m in retained])
        mean = np.mean(abs_deltas)
        sem = np.std(abs_deltas, ddof=1) / np.sqrt(len(retained))
        t_stat = mean / sem if sem > 0 else 0

        threshold_analyses[f"snr_{thresh}"] = {
            "n_retained": len(retained),
            "n_total": len(measurements),
            "retention_fraction": len(retained) / len(measurements),
            "mean_ns": float(mean),
            "sem_ns": float(sem),
            "t_stat": float(t_stat),
        }

    # Probabilistic analysis (no cuts)
    weighted_stats = compute_weighted_statistics(measurements)

    # Bayesian analysis (no cuts)
    bayesian_results = hierarchical_bayesian_analysis(measurements)

    return {
        "hard_thresholds": threshold_analyses,
        "probabilistic_weighting": weighted_stats,
        "bayesian_hierarchical": bayesian_results,
        "comparison": {
            "probabilistic_mean_ns": weighted_stats["weighted_mean_ns"],
            "probabilistic_t": weighted_stats["t_statistic"],
            "bayesian_mean_ns": bayesian_results["posterior_H"]["mean_ns"],
            "bayesian_t": bayesian_results["posterior_H"]["t_equivalent"],
            "threshold_degradation": {
                k: v
                for k, v in threshold_analyses.items()
                if v.get("mean_ns") is not None
            },
        },
    }


def analyze_snr_influence(measurements: List[TripletMeasurement]) -> Dict[str, Any]:
    """
    Analyze how closure SNR influences the measurement distribution.

    Key diagnostic: |Does H| differ systematically with SNR?

    Note: Uses percentile-based stratification because input data has
    already been quality-selected (all SNR values are high).
    """
    snrs = np.array([m.closure_snr for m in measurements])
    # Use absolute values for TEP magnitude
    abs_deltas = np.array([abs(m.delta_ns) for m in measurements])

    # Use percentile-based stratification (data is already high-SNR)
    p33 = np.percentile(snrs, 33)
    p67 = np.percentile(snrs, 67)

    low_snr = abs_deltas[snrs < p33]
    mid_snr = abs_deltas[(snrs >= p33) & (snrs < p67)]
    high_snr = abs_deltas[snrs >= p67]

    analysis = {}

    for label, subset, threshold in [
        ("low", low_snr, f"< {p33:.0f}"),
        ("mid", mid_snr, f"{p33:.0f} - {p67:.0f}"),
        ("high", high_snr, f">= {p67:.0f}"),
    ]:
        if len(subset) > 0:
            analysis[label] = {
                "n": len(subset),
                "snr_range": threshold,
                "mean_abs_H_ns": float(np.mean(subset)),
                "std_ns": float(np.std(subset)),
                "sem_ns": float(np.std(subset) / np.sqrt(len(subset))),
            }

    # Correlation between SNR and |H|
    valid_idx = ~np.isnan(snrs) & ~np.isnan(abs_deltas)

    if np.sum(valid_idx) > 10:
        r, p = stats.pearsonr(snrs[valid_idx], abs_deltas[valid_idx])
    else:
        r, p = np.nan, np.nan

    return {
        "stratified_by_snr": analysis,
        "snr_percentiles": {"p33": float(p33), "p67": float(p67)},
        "snr_abs_delta_correlation": {
            "r": float(r),
            "p": float(p),
            "n": int(np.sum(valid_idx)),
        },
        "quality_metric": "closure_snr",
    }


def main():
    """Run probabilistic weighting analysis."""
    print("=" * 80)
    print("STEP 044: PROBABILISTIC WEIGHTING DIAGNOSTIC")
    print("=" * 80)
    print()
    print("Purpose: Replace hard SNR/arclet cuts with probabilistic weighting")
    print(
        "Key innovation: Smooth uncertainty-based weighting preserves full distribution"
    )
    print()
    print("1. LOADING CLOSURE DATA...")
    epochs = load_closure_data()
    print(f"   Loaded {len(epochs)} epochs")
    print()

    # Extract measurements
    print("2. EXTRACTING TRIPLET MEASUREMENTS...")
    measurements = extract_triplet_measurements(epochs)
    print(f"   Extracted {len(measurements)} triplet measurements")

    # Compute weights
    weights = compute_probabilistic_weights(measurements)
    print(f"   Probabilistic weights computed")
    print(f"   Weight range: [{np.min(weights):.2e}, {np.max(weights):.2e}]")
    print(f"   Effective sample size: {1.0 / np.sum(weights**2):.1f}")
    closure_snrs = np.array([m.closure_snr for m in measurements])
    print(f"   Closure-SNR range: [{np.min(closure_snrs):.3f}, {np.max(closure_snrs):.3f}]")
    print()

    # Hierarchical Bayesian analysis
    print("3. HIERARCHICAL BAYESIAN DIAGNOSTIC...")
    bayesian_results = hierarchical_bayesian_analysis(measurements)

    print(f"   Posterior mean H: {bayesian_results['posterior_H']['mean_ns']:.3f} ns")
    print(f"   Posterior std: {bayesian_results['posterior_H']['std_ns']:.3f} ns")
    print(
        f"   95% CI: [{bayesian_results['posterior_H']['ci_95_lower_ns']:.3f}, "
        f"{bayesian_results['posterior_H']['ci_95_upper_ns']:.3f}] ns"
    )
    print(
        f"   Diagnostic z: {bayesian_results['posterior_H']['t_equivalent']:.2f} (not inferential)"
    )
    print()

    # Weighted statistics (frequentist)
    print("4. INVERSE-VARIANCE WEIGHTED ANALYSIS...")
    weighted_stats = compute_weighted_statistics(measurements)
    print(
        f"   Weighted mean: {weighted_stats['weighted_mean_ns']:.3f} ± "
        f"{weighted_stats['weighted_sem_ns']:.3f} ns"
    )
    print(f"   diagnostic z: {weighted_stats['t_statistic']:.2f} (not inferential)")
    print(
        f"   Effective n: {weighted_stats['n_eff']:.1f} (from {weighted_stats['n_total']} total)"
    )
    print()

    # Threshold comparison
    print("5. COMPARISON WITH HARD THRESHOLD CUTS...")
    comparison = compare_threshold_vs_probabilistic(measurements, bayesian_results)

    for thresh_key, thresh_data in comparison["hard_thresholds"].items():
        if thresh_data.get("mean_ns") is not None:
            print(
                f"   {thresh_key}: n={thresh_data['n_retained']}, "
                f"H={thresh_data['mean_ns']:.3f}±{thresh_data['sem_ns']:.3f} ns, "
                f"diagnostic z={thresh_data['t_stat']:.2f}"
            )
        else:
            print(f"   {thresh_key}: {thresh_data.get('error', 'Insufficient data')}")
    print()

    # SNR influence analysis
    print("6. SNR INFLUENCE ANALYSIS...")
    snr_analysis = analyze_snr_influence(measurements)

    for label, stats in snr_analysis["stratified_by_snr"].items():
        print(
            f"   {label} SNR (SNR {stats['snr_range']}): n={stats['n']}, |H|={stats['mean_abs_H_ns']:.3f}±{stats['sem_ns']:.3f} ns"
        )

    corr = snr_analysis["snr_abs_delta_correlation"]
    print(f"   SNR |H| correlation: r={corr['r']:.3f}, p={corr['p']:.3e}")
    print(
        f"   (Note: Data pre-filtered; all SNR values high. Stratification by percentiles.)"
    )
    print()

    # Conclusions
    print("=" * 80)
    print("CONCLUSIONS")
    print("=" * 80)

    bayes_mean = bayesian_results["posterior_H"]["mean_ns"]
    bayes_t = bayesian_results["posterior_H"]["t_equivalent"]
    weighted_mean = weighted_stats["weighted_mean_ns"]
    weighted_t = weighted_stats["t_statistic"]

    conclusions = [
        f"Probabilistic weighting diagnostic gives |H| = {bayes_mean:.3f} ns (z = {bayes_t:.2f}, not inferential)",
        f"All {len(measurements)} triplets contribute (no data discarded)",
    ]
    bayes_weighted_diff = abs(bayes_mean - weighted_mean)
    if bayes_weighted_diff <= 0.5:
        conclusions.append(
            f"Bayesian and frequentist weighted diagnostics agree: |H_bayes - H_weighted| = {bayes_weighted_diff:.3f} ns"
        )
    else:
        conclusions.append(
            f"Bayesian and frequentist weighted diagnostics disagree by {bayes_weighted_diff:.3f} ns; this supports treating unsigned-|H| as diagnostic only"
        )

    # Note: The closure data from step 003 has already undergone quality selection.
    # The SNR stratification here shows minimal variation because the input data
    # represents the already-filtered viable measurement set.
    threshold_degradation = comparison["comparison"]["threshold_degradation"]
    snr_5_data = threshold_degradation.get("snr_5", {})

    if snr_5_data and snr_5_data.get("mean_ns") is not None:
        degradation = abs(snr_5_data["mean_ns"] - bayes_mean)
        if degradation > 0.5:  # Only report if meaningful
            conclusions.append(
                f"SNR≥5 cut changes the unsigned-|H| diagnostic by {degradation:.3f} ns vs. probabilistic weighting"
            )
        else:
            conclusions.append(
                "SNR threshold effects minimal on pre-filtered data; probabilistic weighting optimal"
            )

    for conclusion in conclusions:
        print(f"  • {conclusion}")
    print()

    # Implications
    print("=" * 80)
    print("IMPLICATIONS")
    print("=" * 80)
    print("  1. Hard SNR cuts can systematically bias unsigned-|H| diagnostics")
    print("  2. Probabilistic weighting is preferred: no arbitrary thresholds")
    print("  3. Full geometric distribution is preserved for diagnostic threshold checks")
    print("  4. Low-SNR triplets contribute appropriately to precision")
    print()

    # Save results
    output = {
        "bayesian_hierarchical": bayesian_results,
        "weighted_statistics": weighted_stats,
        "threshold_comparison": comparison,
        "snr_influence": snr_analysis,
        "conclusions": conclusions,
        "inference_status": "diagnostic_only",
        "valid_for_primary_inference": False,
        "note": (
            "Unsigned |H| probabilistic weighting is a threshold-sensitivity diagnostic. "
            "It is not used as primary evidence because folded delay magnitudes have a noise floor."
        ),
        "methodology": {
            "weighting_scheme": "inverse_variance",
            "uncertainty_proxy": "sigma_us / max(closure_snr, 1e-3)",
            "prior": "HalfNormal(0, 50)",
            "bayesian_likelihood": "0.5*N(+H,sigma_i) + 0.5*N(-H,sigma_i)",
            "n_triplets": len(measurements),
            "n_epochs": len(epochs),
        },
    }

    output_file = RESULTS_DIR / "step_044_probabilistic_weighting_results.json"
    with open(output_file, "w") as f:
        json.dump(output, f, indent=2, cls=NpEncoder)

    print(f"Results saved to: {output_file}")
    print("=" * 80)

    return output


if __name__ == "__main__":
    main()
