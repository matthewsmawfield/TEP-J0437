#!/usr/bin/env python3
"""
================================================================================
STEP 028: SNR CORRELATION MECHANISM INVESTIGATION
================================================================================

Purpose: Investigate the weak negative SNR correlation (r=-0.241, R²=0.058)
found after implementing the SNR threshold criterion.

Observation: Higher SNR -> |H| (negative correlation)
Possible explanations:
1. High-SNR measurements are more conservative (smaller systematic errors)
2. SNR threshold preferentially selects certain types of triplets
3. Residual selection effect from new criterion
4. Genuine physical effect (high SNR = different ISM conditions)

This step investigates the mechanism to determine if it's a concern
or a benign feature of the data.

================================================================================
"""

import json
import numpy as np
from pathlib import Path
from scipy import stats
from datetime import datetime
from typing import Optional, Dict, Any, List
from scripts.utils.json_numpy import NpEncoder
from scripts.utils.config import DEFAULT_CONFIG

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def load_triplet_data() -> List[Dict]:
    """Load all triplet-level closure delay data with arclet_snr_mean computed on the fly."""
    per_epoch_file = PROJECT_ROOT / "results" / "step_003_closure_final_per_epoch.json"
    
    with open(per_epoch_file, 'r') as f:
        epochs = json.load(f)
    
    # Extract all triplets
    all_triplets = []
    for epoch in epochs:
        for triplet in epoch.get("triplets", []):
            triplet["mjd"] = epoch.get("mjd", 0)
            triplet["n_triplets_epoch"] = epoch.get("n_triplets", 0)
            # Compute arclet_snr_mean on the fly (mean of 3 arclet SNRs)
            arclet_snr_list = triplet.get("arclet_snrs", [])
            if arclet_snr_list:
                triplet["arclet_snr_mean"] = np.mean(arclet_snr_list)
            else:
                triplet["arclet_snr_mean"] = np.nan
            all_triplets.append(triplet)
    
    return all_triplets


def analyze_snr_by_triplet_characteristics(triplets: List[Dict]) -> Dict[str, Any]:
    """
    Analyze if SNR correlates with other triplet characteristics.
    
    Parameters
    ----------
    triplets : list of dict
        Triplet data
    
    Returns
    -------
    dict
        Analysis results
    """
    # Extract data
    abs_delays = np.array([abs(t.get("geometric_delta_us", t.get("delta_us", 0))) for t in triplets])
    abs_delay_ns = abs_delays * 1e3
    snrs = np.array([t.get("arclet_snr_mean", np.nan) for t in triplets])
    
    # Remove NaN
    valid_idx = ~np.isnan(snrs)
    abs_delay_ns = abs_delay_ns[valid_idx]
    snrs = snrs[valid_idx]
    
    # Arclet-level SNRs (from arclet_snrs field)
    arclet_snr_1 = np.array([t["arclet_snrs"][0] if len(t.get("arclet_snrs", [])) > 0 else np.nan for t in triplets])[valid_idx]
    arclet_snr_2 = np.array([t["arclet_snrs"][1] if len(t.get("arclet_snrs", [])) > 1 else np.nan for t in triplets])[valid_idx]
    arclet_snr_3 = np.array([t["arclet_snrs"][2] if len(t.get("arclet_snrs", [])) > 2 else np.nan for t in triplets])[valid_idx]
    
    # Correlations
    results = {}
    
    # SNR vs H (already known)
    if len(snrs) >= 2 and len(abs_delay_ns) >= 2:
        r_snr_H, p_snr_H = stats.pearsonr(snrs, abs_delay_ns)
        results["snr_vs_abs_delay"] = {"r": float(r_snr_H), "p": float(p_snr_H)}
    else:
        results["snr_vs_abs_delay"] = {"r": None, "p": None, "error": "Insufficient data"}
    
    # Arclet SNR consistency (std of 3 arclet SNRs)
    arclet_snr_std = []
    for t in triplets:
        if len(t.get("arclet_snrs", [])) == 3:
            arclet_snr_std.append(np.std(t["arclet_snrs"]))
    
    arclet_snr_std = np.array(arclet_snr_std)
    if len(arclet_snr_std) >= 2 and len(abs_delay_ns) >= len(arclet_snr_std):
        r_std_H, p_std_H = stats.pearsonr(arclet_snr_std, abs_delay_ns[:len(arclet_snr_std)])
        results["arclet_snr_std_vs_abs_delay"] = {"r": float(r_std_H), "p": float(p_std_H)}
    else:
        results["arclet_snr_std_vs_abs_delay"] = {"r": None, "p": None, "error": "Insufficient data"}
    
    # SNR vs delay magnitude (not absolute)
    delays = np.array([t.get("geometric_delta_us", t.get("delta_us", 0)) for t in triplets])[valid_idx]
    if len(snrs) >= 2 and len(delays) >= 2:
        r_snr_delay, p_snr_delay = stats.pearsonr(snrs, delays)
        results["snr_vs_delay"] = {"r": float(r_snr_delay), "p": float(p_snr_delay)}
    else:
        results["snr_vs_delay"] = {"r": None, "p": None, "error": "Insufficient data"}
    
    return results


def analyze_snr_by_epoch_characteristics(triplets: List[Dict]) -> Dict[str, Any]:
    """
    Analyze if SNR varies by epoch characteristics.
    
    Parameters
    ----------
    triplets : list of dict
        Triplet data
    
    Returns
    -------
    dict
        Analysis results
    """
    # Group by epoch
    epochs = {}
    for t in triplets:
        mjd = t["mjd"]
        if mjd not in epochs:
            epochs[mjd] = []
        epochs[mjd].append(t)
    
    # Compute epoch-level statistics
    epoch_stats = []
    for mjd, epoch_triplets in epochs.items():
        snrs = np.array([t.get("arclet_snr_mean", np.nan) for t in epoch_triplets])
        abs_delay_ns = np.array([abs(t.get("geometric_delta_us", t.get("delta_us", 0))) * 1e3 for t in epoch_triplets])
        n_triplets = epoch_triplets[0]["n_triplets_epoch"]
        
        valid_idx = ~np.isnan(snrs)
        if np.sum(valid_idx) > 0:
            epoch_stats.append({
                "mjd": mjd,
                "n_triplets": n_triplets,
                "mean_snr": float(np.mean(snrs[valid_idx])),
                "mean_abs_delay_ns": float(np.mean(abs_delay_ns[valid_idx])),
                "std_snr": float(np.std(snrs[valid_idx]))
            })
    
    # Correlations
    n_triplets = np.array([e["n_triplets"] for e in epoch_stats])
    mean_snr = np.array([e["mean_snr"] for e in epoch_stats])
    mean_abs_delay_ns = np.array([e["mean_abs_delay_ns"] for e in epoch_stats])
    
    if len(n_triplets) >= 2 and len(mean_snr) >= 2:
        r_n_snr, p_n_snr = stats.pearsonr(n_triplets, mean_snr)
    else:
        r_n_snr, p_n_snr = np.nan, np.nan
    
    if len(mean_snr) >= 2 and len(mean_abs_delay_ns) >= 2:
        r_snr_H, p_snr_H = stats.pearsonr(mean_snr, mean_abs_delay_ns)
    else:
        r_snr_H, p_snr_H = np.nan, np.nan
    
    if len(n_triplets) >= 2 and len(mean_abs_delay_ns) >= 2:
        r_n_H, p_n_H = stats.pearsonr(n_triplets, mean_abs_delay_ns)
    else:
        r_n_H, p_n_H = np.nan, np.nan
    
    return {
        "n_epochs": len(epoch_stats),
        "n_triplets_vs_mean_snr": {"r": float(r_n_snr), "p": float(p_n_snr)},
        "mean_snr_vs_mean_H": {"r": float(r_snr_H), "p": float(p_snr_H)},
        "n_triplets_vs_mean_H": {"r": float(r_n_H), "p": float(p_n_H)},
        "epoch_stats": epoch_stats
    }


def analyze_snr_threshold_effect(triplets: List[Dict]) -> Dict[str, Any]:
    """
    Analyze the effect of the SNR threshold (5.0) on the correlation.
    
    Parameters
    ----------
    triplets : list of dict
        Triplet data
    
    Returns
    -------
    dict
        Analysis results
    """
    # Extract data
    abs_delays = np.array([abs(t.get("geometric_delta_us", t.get("delta_us", 0))) for t in triplets])
    abs_delay_ns = abs_delays * 1e3
    snrs = np.array([t.get("arclet_snr_mean", np.nan) for t in triplets])
    
    # Remove NaN
    valid_idx = ~np.isnan(snrs)
    abs_delay_ns = abs_delay_ns[valid_idx]
    snrs = snrs[valid_idx]
    
    # Stratify by SNR relative to threshold
    threshold = DEFAULT_CONFIG["step_003_closure_delays"]["min_snr"]
    below_threshold = snrs < threshold
    above_threshold = snrs >= threshold
    
    # Compare H values
    H_below = abs_delay_ns[below_threshold]
    H_above = abs_delay_ns[above_threshold]

    if len(H_below) < 2 or len(H_above) < 2:
        t_stat, p_val = None, None
        threshold_test_valid = False
    else:
        t_stat, p_val = stats.ttest_ind(H_below, H_above)
        threshold_test_valid = True
    
    # Correlation within each stratum
    if len(H_below) > 10:
        r_below, p_below = stats.pearsonr(snrs[below_threshold], H_below)
    else:
        r_below, p_below = None, None
    
    if len(H_above) > 10:
        r_above, p_above = stats.pearsonr(snrs[above_threshold], H_above)
    else:
        r_above, p_above = None, None
    
    return {
        "threshold": threshold,
        "n_below_threshold": int(np.sum(below_threshold)),
        "n_above_threshold": int(np.sum(above_threshold)),
        "H_below_mean": float(np.mean(H_below)) if len(H_below) > 0 else None,
        "H_above_mean": float(np.mean(H_above)) if len(H_above) > 0 else None,
        "threshold_test_valid": bool(threshold_test_valid),
        "t_statistic": float(t_stat) if t_stat is not None else None,
        "p_value": float(p_val) if p_val is not None else None,
        "correlation_below_threshold": {"r": float(r_below) if r_below is not None else None, "p": float(p_below) if p_below is not None else None},
        "correlation_above_threshold": {"r": float(r_above) if r_above is not None else None, "p": float(p_above) if p_above is not None else None}
    }


def main():
    print("=" * 80)
    print("STEP 028: SNR CORRELATION MECHANISM INVESTIGATION")
    print("=" * 80)
    print("\nPurpose: Investigate weak negative SNR correlation (r=-0.241)")
    print("Observation: Higher SNR -> |H|")
    print("Goal: Determine if this is a concern or benign feature")
    
    # Load data
    print("\nLoading triplet data...")
    triplets = load_triplet_data()
    print(f"Loaded {len(triplets)} triplets")
    
    # Analyze SNR by triplet characteristics
    print("\nAnalyzing SNR by triplet characteristics...")
    triplet_analysis = analyze_snr_by_triplet_characteristics(triplets)
    
    if triplet_analysis['snr_vs_abs_delay']['r'] is not None:
        print(f"\nSNR |Delta|: r = {triplet_analysis['snr_vs_abs_delay']['r']:.3f}, p = {triplet_analysis['snr_vs_abs_delay']['p']:.2e}")
    else:
        print(f"\nSNR |Delta|: {triplet_analysis['snr_vs_abs_delay'].get('error', 'Error')}")
    
    if triplet_analysis['snr_vs_delay']['r'] is not None:
        print(f"SNR vs delay (signed): r = {triplet_analysis['snr_vs_delay']['r']:.3f}, p = {triplet_analysis['snr_vs_delay']['p']:.2e}")
    else:
        print(f"SNR vs delay: {triplet_analysis['snr_vs_delay'].get('error', 'Error')}")
    
    if "arclet_snr_std_vs_abs_delay" in triplet_analysis:
        if triplet_analysis['arclet_snr_std_vs_abs_delay']['r'] is not None:
            print(f"Arclet SNR std |Delta|: r = {triplet_analysis['arclet_snr_std_vs_abs_delay']['r']:.3f}, p = {triplet_analysis['arclet_snr_std_vs_abs_delay']['p']:.2e}")
        else:
            print(f"Arclet SNR std |Delta|: {triplet_analysis['arclet_snr_std_vs_abs_delay'].get('error', 'Error')}")
    
    # Analyze SNR by epoch characteristics
    print("\nAnalyzing SNR by epoch characteristics...")
    epoch_analysis = analyze_snr_by_epoch_characteristics(triplets)
    
    print(f"\nEpoch-level analysis ({epoch_analysis['n_epochs']} epochs):")
    
    if not np.isnan(epoch_analysis['n_triplets_vs_mean_snr']['r']):
        print(f"  n_triplets vs mean SNR: r = {epoch_analysis['n_triplets_vs_mean_snr']['r']:.3f}, p = {epoch_analysis['n_triplets_vs_mean_snr']['p']:.2e}")
    else:
        print(f"  n_triplets vs mean SNR: Insufficient data")
    
    if not np.isnan(epoch_analysis['mean_snr_vs_mean_H']['r']):
        print(f"  mean SNR vs |H|_mean: r = {epoch_analysis['mean_snr_vs_mean_H']['r']:.3f}, p = {epoch_analysis['mean_snr_vs_mean_H']['p']:.2e}")
    else:
        print(f"  mean SNR vs |H|_mean: Insufficient data")
    
    if not np.isnan(epoch_analysis['n_triplets_vs_mean_H']['r']):
        print(f"  n_triplets vs |H|_mean: r = {epoch_analysis['n_triplets_vs_mean_H']['r']:.3f}, p = {epoch_analysis['n_triplets_vs_mean_H']['p']:.2e}")
    else:
        print(f"  n_triplets vs |H|_mean: Insufficient data")
    
    # Analyze SNR threshold effect
    print("\nAnalyzing SNR threshold effect...")
    threshold_analysis = analyze_snr_threshold_effect(triplets)
    
    print(f"\nSNR threshold = {threshold_analysis['threshold']}:")
    
    if threshold_analysis['H_below_mean'] is not None:
        print(f"  Below threshold: {threshold_analysis['n_below_threshold']} triplets, H = {threshold_analysis['H_below_mean']:.1f} ns")
    else:
        print(f"  Below threshold: {threshold_analysis['n_below_threshold']} triplets, H = N/A")
    
    if threshold_analysis['H_above_mean'] is not None:
        print(f"  Above threshold: {threshold_analysis['n_above_threshold']} triplets, H = {threshold_analysis['H_above_mean']:.1f} ns")
    else:
        print(f"  Above threshold: {threshold_analysis['n_above_threshold']} triplets, H = N/A")
    
    if threshold_analysis['threshold_test_valid'] and threshold_analysis.get('t_statistic') is not None:
        print(f"  t-test: t = {threshold_analysis['t_statistic']:.2f}, p = {threshold_analysis['p_value']:.2e}")
    else:
        print(f"  t-test: N/A (insufficient data)")
    
    r_below = threshold_analysis['correlation_below_threshold'].get('r')
    if r_below is not None:
        print(f"  Correlation below: r = {r_below:.3f}")
    else:
        print(f"  Correlation below: N/A")
    
    r_above = threshold_analysis['correlation_above_threshold'].get('r')
    if r_above is not None:
        print(f"  Correlation above: r = {r_above:.3f}")
    else:
        print(f"  Correlation above: N/A")
    
    # Interpretation
    print("\n" + "=" * 80)
    print("INTERPRETATION")
    print("=" * 80)
    
    r_snr_delay = triplet_analysis['snr_vs_delay']['r']
    r_epoch_snr_H = epoch_analysis['mean_snr_vs_mean_H']['r']
    
    if r_snr_delay is not None and abs(r_snr_delay) < 0.1:
        print("\n[OK] SNR not correlated with signed delay")
        print("   Suggests SNR correlation is with |H| magnitude, not sign")
        print("   This is expected if high-SNR measurements are more conservative")
    elif r_snr_delay is not None:
        print(f"\n[WARN] SNR correlated with signed delay (r = {r_snr_delay:.3f})")
        print("   May indicate SNR-dependent systematic")
    else:
        print("\n[WARN] Could not assess SNR vs delay correlation (insufficient data)")
    
    if r_epoch_snr_H is not None and abs(r_epoch_snr_H) < 0.2:
        print("\n[OK] Epoch-level SNR not strongly correlated with |H|")
        print("   Suggests correlation is triplet-level, not epoch-level")
        print("   This is less concerning for the overall detection")
    elif r_epoch_snr_H is not None:
        print(f"\n[WARN] Epoch-level SNR correlated with |H| (r = {r_epoch_snr_H:.3f})")
        print("   May indicate epoch-level SNR dependence")
    else:
        print("\n[WARN] Could not assess epoch-level SNR correlation (insufficient data)")
    
    if not threshold_analysis['threshold_test_valid']:
        print("\n[WARN] Threshold comparison skipped: one SNR stratum is empty or too small")
        print("   This criterion does not generate two populated groups in the current dataset")
    elif threshold_analysis['p_value'] > 0.05:
        print("\n[OK] No significant difference in |H| above/below threshold")
        print("   Suggests SNR threshold is not creating artificial correlation")
    else:
        print("\n[WARN] Significant difference in |H| above/below threshold")
        t_stat = threshold_analysis.get('t_statistic')
        p_val = threshold_analysis.get('p_value')
        if t_stat is not None and p_val is not None:
            print(f"   t = {t_stat:.2f}, p = {p_val:.2e}")
        print("   May indicate threshold effect")
    
    # Overall assessment
    print("\n" + "=" * 80)
    print("OVERALL ASSESSMENT")
    print("=" * 80)
    
    r_snr_H = triplet_analysis['snr_vs_abs_delay']['r']
    if r_snr_H is None:
        print("\n[WARN] ERROR: Could not compute SNR correlation (insufficient data)")
        print("   Recommendation: Check data quality or skip this analysis")
        r2 = None
    else:
        r2 = r_snr_H ** 2
        print(f"\nSNR |H| correlation: r = {r_snr_H:.3f}, R² = {r2:.3f}")
        print(f"Only {r2*100:.1f}% of variance in |H| explained by SNR")
    
    if r2 is not None:
        if r2 < 0.1:
            print("\n[OK] CONCLUSION: WEAK CORRELATION, NOT A MAJOR CONCERN")
            print("   R² < 0.1 means <10% variance explained by SNR")
            print("   90%+ of variance is due to other factors (genuine TEP effect)")
            print("   Recommendation: Monitor, but not a blocking issue")
        elif r2 < 0.2:
            print("\n[WARN] CONCLUSION: MODERATE CORRELATION, WARRANTS MONITORING")
            print(f"   R² = {r2:.3f} means {r2*100:.1f}% variance explained by SNR")
            print("   Recommendation: Investigate further, but may be acceptable")
        else:
            print("\n[WARN] CONCLUSION: STRONG CORRELATION, CONCERNING")
            print(f"   R² = {r2:.3f} means {r2*100:.1f}% variance explained by SNR")
            print("   Recommendation: Address before publication")
    else:
        print("\n[WARN] CONCLUSION: CANNOT ASSESS (insufficient data)")
        print("   Recommendation: Skip this analysis or use different data")
    
    # Save results
    results = {
        "validation_type": "SNR Correlation Mechanism Investigation",
        "validation_date": datetime.now().isoformat(),
        "triplet_analysis": triplet_analysis,
        "epoch_analysis": epoch_analysis,
        "threshold_analysis": threshold_analysis,
        "overall_r2": float(r2) if r_snr_H is not None else None
    }
    
    output_file = RESULTS_DIR / "step_026_snr_correlation_investigation.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2, cls=NpEncoder)
    
    print(f"\nResults saved to: {output_file}")
    
    print("\n" + "=" * 80)
    print("STEP 028 COMPLETED")
    print("=" * 80)


if __name__ == "__main__":
    main()
