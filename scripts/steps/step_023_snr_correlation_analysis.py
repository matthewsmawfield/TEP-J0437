#!/usr/bin/env python3
"""
================================================================================
STEP 025: SNR CORRELATION ANALYSIS
================================================================================

Purpose: Investigate the extreme SNR-weighted outlier (|H| = 26.0 ns vs standard
8.9 ns, 2.9x larger). This addresses the concern that the signal correlates
with measurement precision, which would be characteristic of systematic errors.

This step performs:
1. Plot |H| vs SNR for all triplets
2. Fit linear regression: |H| = a + b x SNR
3. Test if b ≠ 0 (SNR-dependent systematic)
4. Compare SNR-weighted vs SNR-independent estimators
5. Identify which triplets drive the SNR-weighted outlier

Expected Outcomes:
- If b ≈ 0: No SNR correlation -> use standard |H|
- If b > 0: Positive correlation -> SNR-dependent systematic
- If b < 0: Negative correlation -> high-SNR underestimates

If significant correlation found:
- Use SNR-corrected or SNR-independent estimator
- Report both SNR-weighted and SNR-independent results
- Investigate instrumental sources of SNR dependence

================================================================================
"""

import json
import numpy as np
import sys
from pathlib import Path
from scipy import stats
from datetime import datetime
from typing import Optional, Dict, Any, List, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.json_numpy import NpEncoder

from scripts.utils.config import RANDOM_SEED
RESULTS_DIR = PROJECT_ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def load_triplet_data() -> List[Dict]:
    """Load all triplet-level closure delay data with arclet SNRs."""
    per_epoch_file = PROJECT_ROOT / "results" / "step_003_closure_final_per_epoch.json"
    
    with open(per_epoch_file, 'r') as f:
        epochs = json.load(f)
    
    # Extract all triplets with arclet SNRs
    all_triplets = []
    for epoch in epochs:
        for triplet in epoch.get("triplets", []):
            triplet["mjd"] = epoch.get("mjd", 0)
            # Use mean of arclet SNRs as the independent SNR measurement
            arclet_snr_list = triplet.get("arclet_snrs", [])
            if arclet_snr_list:
                triplet["arclet_snr_mean"] = np.mean(arclet_snr_list)
            else:
                triplet["arclet_snr_mean"] = np.nan
            all_triplets.append(triplet)
    
    return all_triplets


def compute_snr_correlation(triplets: List[Dict]) -> Dict[str, Any]:
    """
    Compute correlation between |H| and SNR at triplet level.
    
    NOTE: Uses arclet_snr_mean (mean of 3 arclet SNRs) as the independent SNR
    measurement. The "snr" field in the triplet data is derived from |H| itself
    (|H|/1.708), so it cannot be used for correlation analysis.
    
    Parameters
    ----------
    triplets : list of dict
        Triplet data with delta_us and arclet_snr_mean
    
    Returns
    -------
    dict
        Correlation analysis results
    """
    # Extract |H| and SNR (using arclet_snr_mean)
    abs_delays = np.array([abs(t.get("geometric_delta_us", t.get("delta_us", 0))) for t in triplets])
    snrs = np.array([t.get("arclet_snr_mean", np.nan) for t in triplets])
    
    # Convert |H| to ns
    abs_delay_ns = abs_delays * 1e3
    
    # Remove NaN values
    valid_idx = ~np.isnan(snrs)
    snrs = snrs[valid_idx]
    abs_delay_ns = abs_delay_ns[valid_idx]
    
    # Pearson correlation
    r, p = stats.pearsonr(snrs, abs_delay_ns)
    
    # Spearman correlation (rank-based, robust to outliers)
    r_spearman, p_spearman = stats.spearmanr(snrs, abs_delay_ns)
    
    # Linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(snrs, abs_delay_ns)
    
    # Predicted values and residuals
    predicted = intercept + slope * snrs
    residuals = abs_delay_ns - predicted
    
    return {
        "n_triplets": len(triplets),
        "pearson_r": float(r),
        "pearson_p": float(p),
        "spearman_r": float(r_spearman),
        "spearman_p": float(p_spearman),
        "regression": {
            "slope": float(slope),
            "intercept": float(intercept),
            "r_value": float(r_value),
            "p_value": float(p_value),
            "std_err": float(std_err)
        },
        "residuals": {
            "mean": float(np.mean(residuals)),
            "std": float(np.std(residuals)),
            "max_residual": float(np.max(np.abs(residuals)))
        },
        "snr_stats": {
            "mean": float(np.mean(snrs)),
            "std": float(np.std(snrs)),
            "min": float(np.min(snrs)),
            "max": float(np.max(snrs))
        },
        "abs_delay_stats_ns": {
            "mean": float(np.mean(abs_delay_ns)),
            "std": float(np.std(abs_delay_ns)),
            "min": float(np.min(abs_delay_ns)),
            "max": float(np.max(abs_delay_ns))
        }
    }


def snr_stratified_analysis(triplets: List[Dict]) -> Dict[str, Any]:
    """
    Stratify triplets by SNR and compare |H|.
    
    Parameters
    ----------
    triplets : list of dict
        Triplet data
    
    Returns
    -------
    dict
        SNR-stratified analysis
    """
    abs_delays = np.array([abs(t.get("geometric_delta_us", t.get("delta_us", 0))) for t in triplets])
    snrs = np.array([t.get("arclet_snr_mean", np.nan) for t in triplets])
    abs_delay_ns = abs_delays * 1e3
    
    # Remove NaN
    valid_idx = ~np.isnan(snrs)
    snrs = snrs[valid_idx]
    abs_delay_ns = abs_delay_ns[valid_idx]
    
    # Divide into SNR quartiles
    q1 = np.percentile(snrs, 25)
    q2 = np.percentile(snrs, 50)
    q3 = np.percentile(snrs, 75)
    
    groups = {
        "Q1_lowest_snr": [i for i, s in enumerate(snrs) if s <= q1],
        "Q2": [i for i, s in enumerate(snrs) if q1 < s <= q2],
        "Q3": [i for i, s in enumerate(snrs) if q2 < s <= q3],
        "Q4_highest_snr": [i for i, s in enumerate(snrs) if s > q3]
    }
    
    results = {}
    for group_name, indices in groups.items():
        group_abs_delay = abs_delay_ns[indices]
        group_snr = snrs[indices]
        
        results[group_name] = {
            "n_triplets": len(indices),
            "mean_abs_delay_ns": float(np.mean(group_abs_delay)),
            "std_abs_delay_ns": float(np.std(group_abs_delay)),
            "median_abs_delay_ns": float(np.median(group_abs_delay)),
            "mean_snr": float(np.mean(group_snr)),
            "t_statistic": float(np.mean(group_abs_delay) / (np.std(group_abs_delay) / np.sqrt(len(indices)))) if len(indices) > 0 else None
        }
    
    return results


def compute_snr_weighted_vs_unweighted(triplets: List[Dict]) -> Dict[str, Any]:
    """
    Compare SNR-weighted and unweighted |H| estimators.
    
    Parameters
    ----------
    triplets : list of dict
        Triplet data
    
    Returns
    -------
    dict
        Comparison of estimators
    """
    abs_delays = np.array([abs(t.get("geometric_delta_us", t.get("delta_us", 0))) for t in triplets])
    snrs = np.array([t.get("arclet_snr_mean", 1.0) for t in triplets])
    abs_delay_ns = abs_delays * 1e3
    
    # Unweighted (standard)
    H_unweighted = np.mean(abs_delay_ns)
    sem_unweighted = np.std(abs_delay_ns, ddof=1) / np.sqrt(len(abs_delay_ns))
    
    # SNR-weighted
    H_weighted = np.average(abs_delay_ns, weights=snrs)
    
    # Compute weighted SEM
    # Use bootstrap for weighted SEM
    n_boot = 1000
    bootstrap_weighted = []
    for _ in range(n_boot):
        boot_idx = np.random.choice(len(abs_delay_ns), size=len(abs_delay_ns), replace=True)
        boot_H = abs_delay_ns[boot_idx]
        boot_snr = snrs[boot_idx]
        bootstrap_weighted.append(np.average(boot_H, weights=boot_snr))
    
    sem_weighted = np.std(bootstrap_weighted, ddof=1)
    
    # SNR-independent (equal weight per SNR bin)
    # Normalize SNR to unit weight
    H_snr_independent = np.mean(abs_delay_ns / snrs * np.mean(snrs))
    sem_snr_independent = np.std(abs_delay_ns / snrs * np.mean(snrs), ddof=1) / np.sqrt(len(abs_delay_ns))
    
    return {
        "unweighted": {
            "abs_delay_ns": float(H_unweighted),
            "sem_ns": float(sem_unweighted),
            "t_statistic": float(H_unweighted / sem_unweighted) if sem_unweighted > 0 else 0.0
        },
        "snr_weighted": {
            "abs_delay_ns": float(H_weighted),
            "sem_ns": float(sem_weighted),
            "t_statistic": float(H_weighted / sem_weighted) if sem_weighted > 0 else 0.0
        },
        "snr_independent": {
            "abs_delay_ns": float(H_snr_independent),
            "sem_ns": float(sem_snr_independent),
            "t_statistic": float(H_snr_independent / sem_snr_independent) if sem_snr_independent > 0 else 0.0
        },
        "ratio_weighted_unweighted": float(H_weighted / H_unweighted),
        "ratio_independent_unweighted": float(H_snr_independent / H_unweighted)
    }


def identify_outlier_triplets(triplets: List[Dict], correlation: Dict) -> Dict[str, Any]:
    """
    Identify which triplets drive the SNR-weighted outlier.
    
    Parameters
    ----------
    triplets : list of dict
        Triplet data
    correlation : dict
        Correlation analysis results
    
    Returns
    -------
    dict
        Outlier triplet analysis
    """
    abs_delays = np.array([abs(t.get("geometric_delta_us", t.get("delta_us", 0))) for t in triplets])
    snrs = np.array([t.get("arclet_snr_mean", 1.0) for t in triplets])
    abs_delay_ns = abs_delays * 1e3
    
    # Compute residuals from regression
    slope = correlation["regression"]["slope"]
    intercept = correlation["regression"]["intercept"]
    predicted = intercept + slope * snrs
    residuals = abs_delay_ns - predicted
    
    # Identify outliers (|residual| > 2sigma)
    residual_std = np.std(residuals)
    outlier_threshold = 2 * residual_std
    outlier_indices = np.where(np.abs(residuals) > outlier_threshold)[0]
    
    # Get outlier triplets
    outliers = [triplets[i] for i in outlier_indices]
    
    # Compute statistics for outliers
    outlier_abs_delay = abs_delay_ns[outlier_indices]
    outlier_snr = snrs[outlier_indices]
    
    return {
        "n_outliers": len(outlier_indices),
        "outlier_fraction": float(len(outlier_indices) / len(triplets)),
        "outlier_threshold_ns": float(outlier_threshold),
        "outlier_mean_abs_delay_ns": float(np.mean(outlier_abs_delay)) if len(outlier_abs_delay) > 0 else None,
        "outlier_mean_snr": float(np.mean(outlier_snr)) if len(outlier_snr) > 0 else None,
        "outlier_indices": outlier_indices.tolist()
    }


def main():
    print("=" * 80)
    print("STEP 025: SNR CORRELATION ANALYSIS")
    print("=" * 80)
    print("\nPurpose: Investigate SNR-weighted outlier (|H| = 26.0 ns vs 8.9 ns)")
    print("Concern: Signal may correlate with measurement precision")
    print("Expected for systematics: Strong correlation with SNR")
    print("Expected for genuine signal: Weak or no correlation with SNR")
    
    # Load data
    print("\nLoading triplet data...")
    triplets = load_triplet_data()
    print(f"Loaded {len(triplets)} triplets")
    
    # Correlation analysis
    print("\nComputing SNR correlation...")
    correlation = compute_snr_correlation(triplets)
    
    print(f"\nPearson correlation:")
    print(f"  r = {correlation['pearson_r']:.3f}")
    print(f"  p = {correlation['pearson_p']:.2e}")
    print(f"  Significant: {correlation['pearson_p'] < 0.05}")
    
    print(f"\nSpearman correlation (rank-based):")
    print(f"  r = {correlation['spearman_r']:.3f}")
    print(f"  p = {correlation['spearman_p']:.2e}")
    print(f"  Significant: {correlation['spearman_p'] < 0.05}")
    
    print(f"\nLinear regression: |H| = a + b x SNR")
    print(f"  Slope (b): {correlation['regression']['slope']:.3f} ns per SNR unit")
    print(f"  Intercept (a): {correlation['regression']['intercept']:.3f} ns")
    print(f"  R²: {correlation['regression']['r_value']**2:.3f}")
    print(f"  p-value: {correlation['regression']['p_value']:.2e}")
    
    # SNR-stratified analysis
    print("\nSNR-stratified analysis:")
    strat_snr = snr_stratified_analysis(triplets)
    for group, data in strat_snr.items():
        print(f"  {group}: |Delta| = {data['mean_abs_delay_ns']:.1f} ns (SNR = {data['mean_snr']:.1f})")
    
    # Compare estimators
    print("\nEstimator comparison:")
    estimators = compute_snr_weighted_vs_unweighted(triplets)
    print(f"  Unweighted: |Delta| = {estimators['unweighted']['abs_delay_ns']:.1f} ns ({estimators['unweighted']['t_statistic']:.1f}sigma)")
    print(f"  SNR-weighted: |Delta| = {estimators['snr_weighted']['abs_delay_ns']:.1f} ns ({estimators['snr_weighted']['t_statistic']:.1f}sigma)")
    print(f"  SNR-independent: |Delta| = {estimators['snr_independent']['abs_delay_ns']:.1f} ns ({estimators['snr_independent']['t_statistic']:.1f}sigma)")
    print(f"\n  Ratio (weighted/unweighted): {estimators['ratio_weighted_unweighted']:.2f}x")
    print(f"  Ratio (independent/unweighted): {estimators['ratio_independent_unweighted']:.2f}x")
    
    # Outlier analysis
    print("\nOutlier triplet analysis:")
    outliers = identify_outlier_triplets(triplets, correlation)
    print(f"  Outliers (|residual| > 2sigma): {outliers['n_outliers']}/{correlation['n_triplets']} ({outliers['outlier_fraction']*100:.1f}%)")
    if outliers['outlier_mean_abs_delay_ns'] is not None:
        print(f"  Outlier |Delta|_mean: {outliers['outlier_mean_abs_delay_ns']:.1f} ns")
        print(f"  Outlier mean SNR: {outliers['outlier_mean_snr']:.1f}")
    
    # Interpretation
    print("\n" + "=" * 80)
    print("INTERPRETATION")
    print("=" * 80)
    
    if correlation['pearson_p'] < 0.05 and abs(correlation['pearson_r']) > 0.1:
        print("\n[WARN] SIGNIFICANT SNR CORRELATION DETECTED")
        print(f"   Pearson r = {correlation['pearson_r']:.3f}, p = {correlation['pearson_p']:.2e}")
        
        if correlation['regression']['slope'] > 0:
            print("   POSITIVE correlation: higher SNR -> higher |H|")
            print("   Suggests SNR-dependent systematic error")
            print("   Recommendation: Use SNR-independent estimator")
        else:
            print("   NEGATIVE correlation: higher SNR -> |H|")
            print("   Suggests high-SNR measurements underestimate |H|")
            print("   Recommendation: Investigate high-SNR processing")
        
        print(f"\n   SNR-weighted outlier explained by correlation:")
        print(f"   - Slope = {correlation['regression']['slope']:.3f} ns/SNR")
        print(f"   - At SNR = {correlation['snr_stats']['max']:.1f}: predicted excess = {correlation['regression']['slope'] * (correlation['snr_stats']['max'] - correlation['snr_stats']['mean']):.1f} ns")
    else:
        print("\n[OK] NO SIGNIFICANT SNR CORRELATION")
        print(f"   Pearson r = {correlation['pearson_r']:.3f}, p = {correlation['pearson_p']:.2e}")
        print("   SNR-weighted outlier (26.0 ns) may be statistical fluke")
        print("   Recommendation: Use standard unweighted estimator")
    
    if abs(estimators['ratio_weighted_unweighted'] - 1.0) > 0.5:
        print(f"\n[WARN] LARGE DISCREPANCY BETWEEN ESTIMATORS")
        print(f"   Weighted/unweighted ratio = {estimators['ratio_weighted_unweighted']:.2f}x")
        print("   This confirms SNR-dependent bias")
    else:
        print(f"\n[OK] ESTIMATORS AGREE")
        print(f"   Weighted/unweighted ratio = {estimators['ratio_weighted_unweighted']:.2f}x")
    
    # Save results
    results = {
        "validation_type": "SNR Correlation Analysis",
        "validation_date": datetime.now().isoformat(),
        "correlation_analysis": correlation,
        "snr_stratified": strat_snr,
        "estimator_comparison": estimators,
        "outlier_analysis": outliers
    }
    
    output_file = RESULTS_DIR / "step_023_snr_correlation_analysis.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2, cls=NpEncoder)
    
    print(f"\nResults saved to: {output_file}")
    
    print("\n" + "=" * 80)
    print("STEP 025 COMPLETED")
    print("=" * 80)


if __name__ == "__main__":
    main()
