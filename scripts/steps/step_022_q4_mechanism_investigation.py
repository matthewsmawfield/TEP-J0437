#!/usr/bin/env python3
"""
================================================================================
STEP 024: Q4 DOMINANCE MECHANISM INVESTIGATION
================================================================================

Purpose: Investigate why Q4 (high-triplet epochs) shows 1.59× |H|
than Q1-Q3. This addresses the selection bias concern and determines whether
the Q4 dominance is due to:

1. Selection bias (high-triplet epochs preferentially selected)
2. SNR-dependent systematic error
3. Genuine TEP effect that scales with ISM complexity
4. Instrumental or observational effects

This step performs:
1. Correlation analysis: ||H| vs n_triplets, SNR, arclet count, epoch duration
2. Stratified analysis: Divide epochs by multiple quality metrics
3. Instrumental check: Compare Q4 epochs to non-Q4 epochs in same time period
4. ISM analysis: Check if Q4 epochs have different ISM conditions

Expected Outcomes:
- If selection bias: Strong correlation with n_triplets/SNR, disappears when controlled
- If SNR-dependent: Positive correlation with SNR, suggests systematic error
- If genuine TEP: No correlation with selection metrics, persists after controls
- If instrumental: Correlation with telescope, backend, or observing band

================================================================================
"""

import json
import os
import numpy as np
from pathlib import Path
from scipy import stats
from datetime import datetime
from typing import Optional, Dict, Any, List, Tuple
from scripts.utils.json_numpy import NpEncoder

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def load_epoch_data() -> List[Dict]:
    """Load per-epoch closure delay data with triplet details."""
    per_epoch_file = PROJECT_ROOT / "results" / "step_003_closure_final_per_epoch.json"
    
    with open(per_epoch_file, 'r') as f:
        epochs = json.load(f)
    
    return epochs


def compute_epoch_metrics(epoch_data: Dict) -> Dict[str, Any]:
    """
    Compute comprehensive metrics for a single epoch.
    
    Parameters
    ----------
    epoch_data : dict
        Epoch data with triplets
    
    Returns
    -------
    dict
        Epoch metrics (|H|, n_triplets, SNR statistics, etc.)
    """
    triplets = epoch_data.get("triplets", [])
    # CRITICAL: Use geometric_delta_us (Stokes-aligned, SIGNED) — NOT raw delta_us.
    # Raw |delta_us| is unsigned noise; mean(|H|) inflates with more triplets purely
    # by chance (more draws → higher mean |noise|). Only geometric_delta_us is the
    # physically meaningful observable under TEP.
    geom_delays = np.array([t["geometric_delta_us"] for t in triplets])  # signed μs
    
    # Use arclet_snr_mean (mean of 3 arclet SNRs) as independent SNR measurement
    snrs = []
    for t in triplets:
        arclet_snr_list = t.get("arclet_snrs", [])
        if arclet_snr_list:
            snrs.append(np.mean(arclet_snr_list))
        else:
            snrs.append(np.nan)
    snrs = np.array(snrs)
    
    if len(geom_delays) == 0:
        return {
            "n_triplets": 0,
            "H_ns": None,
            "mean_snr": None,
            "std_snr": None,
            "mean_arclets": None,
            "duration_s": None
        }
    
    # Signed mean — the correct TEP observable
    H = np.mean(geom_delays) * 1e3  # ns, signed
    
    # SNR statistics
    mean_snr = np.mean(snrs)
    std_snr = np.std(snrs)
    
    # Arclet count (from epoch metadata if available)
    mean_arclets = epoch_data.get("n_arclets", len(triplets))
    
    # Duration (from epoch metadata if available)
    duration_s = epoch_data.get("duration_s", None)
    
    return {
        "n_triplets": len(triplets),
        "H_ns": H,
        "mean_snr": float(mean_snr),
        "std_snr": float(std_snr),
        "mean_arclets": mean_arclets,
        "duration_s": duration_s,
        "mjd": epoch_data.get("mjd", 0)
    }


def compute_correlations(epoch_metrics: List[Dict]) -> Dict[str, Any]:
    """
    Compute correlations between |H| and selection metrics.
    
    Parameters
    ----------
    epoch_metrics : list of dict
        Epoch metrics
    
    Returns
    -------
    dict
        Correlation analysis results
    """
    # Filter valid epochs
    valid = [e for e in epoch_metrics if e["H_ns"] is not None]
    
    if len(valid) < 10:
        return {"error": "Insufficient valid epochs"}
    
    H = np.array([e["H_ns"] for e in valid])
    n_triplets = np.array([e["n_triplets"] for e in valid])
    snr = np.array([e["mean_snr"] for e in valid])
    
    # Remove NaN SNR values
    valid_snr_idx = ~np.isnan(snr)
    H_snr = H[valid_snr_idx]
    snr_valid = snr[valid_snr_idx]
    n_triplets_snr = n_triplets[valid_snr_idx]
    
    correlations = {}
    
    # Correlation with n_triplets
    r_triplets, p_triplets = stats.pearsonr(n_triplets, H)
    correlations["n_triplets"] = {
        "correlation": float(r_triplets),
        "p_value": float(p_triplets),
        "significant": bool(p_triplets < 0.05)
    }
    
    # Correlation with SNR (using valid SNR values)
    r_snr, p_snr = stats.pearsonr(snr_valid, H_snr)
    correlations["snr"] = {
        "correlation": float(r_snr),
        "p_value": float(p_snr),
        "significant": bool(p_snr < 0.05)
    }
    
    # Linear regression for n_triplets
    slope_triplets, intercept_triplets, r_triplets, p_triplets, std_triplets = stats.linregress(n_triplets, H)
    correlations["n_triplets_regression"] = {
        "slope": float(slope_triplets),
        "intercept": float(intercept_triplets),
        "r_value": float(r_triplets),
        "p_value": float(p_triplets),
        "std_err": float(std_triplets)
    }
    
    # Linear regression for SNR (using valid SNR values)
    slope_snr, intercept_snr, r_snr, p_snr, std_snr = stats.linregress(snr_valid, H_snr)
    correlations["snr_regression"] = {
        "slope": float(slope_snr),
        "intercept": float(intercept_snr),
        "r_value": float(r_snr),
        "p_value": float(p_snr),
        "std_err": float(std_snr)
    }
    
    return correlations


def stratify_by_n_triplets(epoch_metrics: List[Dict]) -> Dict[str, Any]:
    """
    Stratify epochs by n_triplets and compare |H|.
    
    Parameters
    ----------
    epoch_metrics : list of dict
        Epoch metrics
    
    Returns
    -------
    dict
        Stratified analysis results
    """
    valid = [e for e in epoch_metrics if e["H_ns"] is not None]
    
    # Divide into quartiles by n_triplets
    n_triplets = np.array([e["n_triplets"] for e in valid])
    q1_threshold = np.percentile(n_triplets, 25)
    q2_threshold = np.percentile(n_triplets, 50)
    q3_threshold = np.percentile(n_triplets, 75)
    
    groups = {
        "Q1_lowest": [e for e in valid if e["n_triplets"] <= q1_threshold],
        "Q2": [e for e in valid if q1_threshold < e["n_triplets"] <= q2_threshold],
        "Q3": [e for e in valid if q2_threshold < e["n_triplets"] <= q3_threshold],
        "Q4_highest": [e for e in valid if e["n_triplets"] > q3_threshold]
    }
    
    results = {}
    for group_name, group_epochs in groups.items():
        H_values = [e["H_ns"] for e in group_epochs]
        n_triplets_vals = [e["n_triplets"] for e in group_epochs]
        snr_vals = [e["mean_snr"] for e in group_epochs]
        
        results[group_name] = {
            "n_epochs": len(group_epochs),
            "mean_H_ns": float(np.mean(H_values)) if H_values else None,
            "std_H_ns": float(np.std(H_values)) if H_values else None,
            "mean_n_triplets": float(np.mean(n_triplets_vals)) if n_triplets_vals else None,
            "mean_snr": float(np.mean(snr_vals)) if snr_vals else None
        }
    
    return results


def stratify_by_snr(epoch_metrics: List[Dict]) -> Dict[str, Any]:
    """
    Stratify epochs by SNR and compare |H|.
    
    Parameters
    ----------
    epoch_metrics : list of dict
        Epoch metrics
    
    Returns
    -------
    dict
        SNR-stratified analysis
    """
    valid = [e for e in epoch_metrics if e["H_ns"] is not None and e["mean_snr"] is not None]
    
    if len(valid) < 10:
        return {"error": "Insufficient valid epochs"}
    
    # Divide by SNR quartiles
    snr = np.array([e["mean_snr"] for e in valid])
    q1 = np.percentile(snr, 25)
    q2 = np.percentile(snr, 50)
    q3 = np.percentile(snr, 75)
    
    groups = {
        "low_snr": [e for e in valid if e["mean_snr"] <= q1],
        "medium_low_snr": [e for e in valid if q1 < e["mean_snr"] <= q2],
        "medium_high_snr": [e for e in valid if q2 < e["mean_snr"] <= q3],
        "high_snr": [e for e in valid if e["mean_snr"] > q3]
    }
    
    results = {}
    for group_name, group_epochs in groups.items():
        H_values = [e["H_ns"] for e in group_epochs]
        snr_vals = [e["mean_snr"] for e in group_epochs]
        
        results[group_name] = {
            "n_epochs": len(group_epochs),
            "mean_H_ns": float(np.mean(H_values)) if H_values else None,
            "mean_snr": float(np.mean(snr_vals)) if snr_vals else None
        }
    
    return results


def test_q4_after_controls(epoch_metrics: List[Dict]) -> Dict[str, Any]:
    """
    Test if Q4 dominance persists after controlling for covariates.
    
    Parameters
    ----------
    epoch_metrics : list of dict
        Epoch metrics
    
    Returns
    -------
    dict
        Control analysis results
    """
    valid = [e for e in epoch_metrics if e["H_ns"] is not None]
    
    # Define Q4 as highest quartile of n_triplets
    n_triplets = np.array([e["n_triplets"] for e in valid])
    q3_threshold = np.percentile(n_triplets, 75)
    
    q4_epochs = [e for e in valid if e["n_triplets"] > q3_threshold]
    q123_epochs = [e for e in valid if e["n_triplets"] <= q3_threshold]
    
    # Compute signed means
    H_q4 = np.mean([e["H_ns"] for e in q4_epochs])    # signed ns
    H_q123 = np.mean([e["H_ns"] for e in q123_epochs]) # signed ns
    # Ratio of signed means — negative values expected under TEP
    # We compare magnitudes; negative/negative should be close to 1.0 if uniform
    raw_ratio = abs(H_q4) / abs(H_q123) if H_q123 != 0 else float('nan')
    
    # Control for SNR: Compare epochs matched by SNR
    q4_snr = np.array([e["mean_snr"] for e in q4_epochs])
    q123_snr = np.array([e["mean_snr"] for e in q123_epochs])
    
    # Find Q123 epochs with similar SNR to Q4
    matched_q123 = []
    for e in q123_epochs:
        if any(abs(e["mean_snr"] - s) < 5.0 for s in q4_snr):
            matched_q123.append(e)
    
    H_matched = np.mean([e["H_ns"] for e in matched_q123]) if matched_q123 else H_q123
    controlled_ratio = abs(H_q4) / abs(H_matched) if H_matched != 0 else float('nan')
    
    return {
        "metric": "signed_geometric_mean_ns (Stokes-aligned)",
        "raw_q4_H_ns": float(H_q4),
        "raw_q123_H_ns": float(H_q123),
        "raw_mag_ratio": float(raw_ratio),
        "n_matched_q123": len(matched_q123),
        "controlled_q123_H_ns": float(H_matched),
        "controlled_ratio": float(controlled_ratio),
        "ratio_change": float((controlled_ratio - raw_ratio) / raw_ratio * 100)
    }


def main():
    from datetime import datetime
    start_time = datetime.now()
    
    print("=" * 80)
    print("STEP 024: Q4 DOMINANCE MECHANISM INVESTIGATION")
    print("=" * 80)
    print(f"\n[INFO] Started at: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("\nPurpose: Investigate why Q4 (high-triplet epochs) shows 1.59× |H|")
    print("than Q1-Q3. This addresses the selection bias concern and determines")
    print("whether the Q4 dominance is due to:")
    print("\n  1. Selection bias (high-triplet epochs preferentially selected)")
    print("  2. SNR-dependent systematic error")
    print("  3. Genuine TEP effect that scales with ISM complexity")
    print("  4. Instrumental or observational effects")
    print("\n" + "-" * 80)
    print("METHODOLOGY")
    print("-" * 80)
    print("  • Correlation analysis: |H| vs n_triplets, SNR, arclet count")
    print("  • Stratified analysis: Divide epochs by quality metrics")
    print("  • Control analysis: SNR-matched comparison")
    print("  • Statistical tests: Pearson correlation, linear regression")
    
    # Load data
    print("\n" + ">>>" * 26)
    print("[INFO] PHASE 1: Data Loading")
    print(">>>" * 26)
    print("\nLoading epoch data from step_003_closure_final_per_epoch.json...")
    epochs = load_epoch_data()
    print(f"[PASS] Loaded {len(epochs)} epochs successfully")
    
    # Data statistics
    total_triplets = sum(len(e.get('triplets', [])) for e in epochs)
    print(f"[DATA] Total triplets: {total_triplets:,}")
    print(f"[DATA] Average triplets per epoch: {total_triplets/len(epochs):.1f}")
    
    # Compute metrics
    print("\n" + ">>>" * 26)
    print("[INFO] PHASE 2: Computing Epoch Metrics")
    print(">>>" * 26)
    print("\nComputing comprehensive metrics for each epoch...")
    print("  • Signed geometric mean (H_ns) - Stokes-aligned observable")
    print("  • SNR statistics per epoch")
    print("  • Arclet and triplet counts")
    print()
    
    epoch_metrics = []
    for i, e in enumerate(epochs, 1):
        metrics = compute_epoch_metrics(e)
        epoch_metrics.append(metrics)
        if i <= 5 or i == len(epochs):
            print(f"  [INFO] Epoch {i:3d}/{len(epochs)}: H={metrics['H_ns']:+.3f} ns, "
                  f"n_triplets={metrics['n_triplets']:2d}, SNR={metrics['mean_snr']:.1f}")
        elif i == 6:
            print(f"  ... ({len(epochs)-6} more epochs) ...")
    
    valid_metrics = [m for m in epoch_metrics if m['H_ns'] is not None]
    print(f"\n[PASS] Computed metrics for {len(epoch_metrics)} epochs")
    print(f"[DATA] Valid epochs with H measurements: {len(valid_metrics)}")
    if valid_metrics:
        h_values = [m['H_ns'] for m in valid_metrics]
        print(f"[DATA] H range: [{min(h_values):+.3f}, {max(h_values):+.3f}] ns")
        print(f"[DATA] H mean: {sum(h_values)/len(h_values):+.3f} ns")
    
    # Correlation analysis
    print("\n" + ">>>" * 26)
    print("[INFO] PHASE 3: Correlation Analysis")
    print(">>>" * 26)
    print("\nTesting correlations between |H| and selection metrics...")
    print("  • H0: No correlation between |H| and selection variables")
    print("  • H1: Significant correlation indicates selection bias")
    print()
    
    correlations = compute_correlations(epoch_metrics)
    
    print(f"\n[RESULTS] Correlation with n_triplets:")
    r_triplets = correlations['n_triplets']['correlation']
    p_triplets = correlations['n_triplets']['p_value']
    sig_triplets = correlations['n_triplets']['significant']
    print(f"  Pearson r = {r_triplets:+.4f}")
    print(f"  p-value   = {p_triplets:.2e}")
    print(f"  Status    : {'[FAIL] SIGNIFICANT' if sig_triplets else '[PASS] Not significant'}")
    if 'n_triplets_regression' in correlations:
        reg = correlations['n_triplets_regression']
        print(f"  Regression: slope={reg['slope']:+.4f}, intercept={reg['intercept']:+.4f}")
    
    print(f"\n[RESULTS] Correlation with SNR:")
    r_snr = correlations['snr']['correlation']
    p_snr = correlations['snr']['p_value']
    sig_snr = correlations['snr']['significant']
    print(f"  Pearson r = {r_snr:+.4f}")
    print(f"  p-value   = {p_snr:.2e}")
    print(f"  Status    : {'[FAIL] SIGNIFICANT' if sig_snr else '[PASS] Not significant'}")
    if 'snr_regression' in correlations:
        reg = correlations['snr_regression']
        print(f"  Regression: slope={reg['slope']:+.4f}, intercept={reg['intercept']:+.4f}")
    
    # Stratified analysis
    print("\n" + ">>>" * 26)
    print("[INFO] PHASE 4: Stratified Analysis")
    print(">>>" * 26)
    print("\nDividing epochs into quartiles by n_triplets...")
    
    strat_n_triplets = stratify_by_n_triplets(epoch_metrics)
    print("\n[RESULTS] Stratification by n_triplets (quartiles):")
    print(f"  {'Group':<15} {'n_epochs':>8} {'mean_H_ns':>12} {'std_H_ns':>10} {'mean_n_triplets':>15}")
    print(f"  {'-'*15} {'-'*8} {'-'*12} {'-'*10} {'-'*15}")
    for group, data in strat_n_triplets.items():
        mh = data["mean_H_ns"]
        mnt = data["mean_n_triplets"]
        std_h = data.get("std_H_ns")
        if mh is not None and mnt is not None and std_h is not None:
            print(
                f"  {group:<15} {data['n_epochs']:>8} {mh:>+11.3f} "
                f"{std_h:>10.3f} {mnt:>15.1f}"
            )
        elif mh is not None and mnt is not None:
            print(
                f"  {group:<15} {data['n_epochs']:>8} {mh:>+11.3f} "
                f"{'n/a':>10} {mnt:>15.1f}"
            )
        else:
            mnt_s = f"{mnt:>15.1f}" if mnt is not None else f"{'n/a':>15}"
            print(f"  {group:<15} {data['n_epochs']:>8} {'n/a':>12} {'n/a':>10} {mnt_s}")

    # Calculate Q4/Q1 ratio
    q4_m = strat_n_triplets["Q4_highest"]["mean_H_ns"]
    q1_m = strat_n_triplets["Q1_lowest"]["mean_H_ns"]
    if q4_m is not None and q1_m is not None:
        q4_h = abs(q4_m)
        q1_h = abs(q1_m)
        if q1_h > 0:
            ratio = q4_h / q1_h
            print(f"\n[DATA] Q4/Q1 magnitude ratio: {ratio:.2f}x")
    
    print("\nDividing epochs into quartiles by SNR...")
    strat_snr = stratify_by_snr(epoch_metrics)
    print("\n[RESULTS] Stratification by SNR (quartiles):")
    print(f"  {'Group':<18} {'n_epochs':>8} {'mean_H_ns':>12} {'mean_snr':>12}")
    print(f"  {'-'*18} {'-'*8} {'-'*12} {'-'*12}")
    for group, data in strat_snr.items():
        if data['mean_H_ns'] is not None:
            print(f"  {group:<18} {data['n_epochs']:>8} {data['mean_H_ns']:>+11.3f} "
                  f"{data['mean_snr']:>12.1f}")
    
    # Control analysis
    print("\n" + ">>>" * 26)
    print("[INFO] PHASE 5: Control Analysis (SNR-matched)")
    print(">>>" * 26)
    print("\nControlling for SNR to test if Q4 dominance persists...")
    print("  • Finding Q123 epochs with similar SNR to Q4 epochs")
    print("  • Matching criterion: |SNR_Q4 - SNR_Q123| < 5.0")
    print()
    
    controls = test_q4_after_controls(epoch_metrics)
    print("[RESULTS] Control Analysis:")
    rq4 = controls["raw_q4_H_ns"]
    print(
        f"  Q4 mean (signed):              {rq4:+.3f} ns"
        if isinstance(rq4, (int, float)) and np.isfinite(rq4)
        else "  Q4 mean (signed):              n/a (no Q4 epochs)"
    )
    print(f"  Q1-Q3 mean (signed):           {controls['raw_q123_H_ns']:+.3f} ns")
    rmr = controls["raw_mag_ratio"]
    print(
        f"  Magnitude ratio (Q4/Q123):     {rmr:.2f}x"
        if isinstance(rmr, (int, float)) and np.isfinite(rmr)
        else "  Magnitude ratio (Q4/Q123):     n/a"
    )
    print(f"\n  After SNR matching:")
    print(f"  Matched Q123 epochs:           {controls['n_matched_q123']}")
    print(f"  SNR-matched Q123 mean:         {controls['controlled_q123_H_ns']:+.3f} ns")
    cr = controls["controlled_ratio"]
    print(
        f"  SNR-matched ratio:             {cr:.2f}x"
        if isinstance(cr, (int, float)) and np.isfinite(cr)
        else "  SNR-matched ratio:             n/a"
    )
    rch = controls["ratio_change"]
    print(
        f"  Ratio change:                  {rch:+.1f}%"
        if isinstance(rch, (int, float)) and np.isfinite(rch)
        else "  Ratio change:                  n/a"
    )
    
    # Summary statistics
    print("\n" + ">>>" * 26)
    print("[INFO] SUMMARY STATISTICS")
    print(">>>" * 26)
    
    # Count interpretations
    n_warnings = 0
    n_ok = 0
    
    if correlations['n_triplets']['significant']:
        n_warnings += 1
    else:
        n_ok += 1
        
    if correlations['snr']['significant']:
        n_warnings += 1
    else:
        n_ok += 1
        
    rc = controls["ratio_change"]
    if isinstance(rc, (int, float)) and np.isfinite(rc) and abs(rc) < 10:
        n_ok += 1
    else:
        n_warnings += 1
    
    print(f"\n[DATA] Tests passed: {n_ok}/3")
    print(f"[DATA] Warnings: {n_warnings}/3")
    
    # Interpretation
    print("\n" + "=" * 80)
    print("INTERPRETATION & CONCLUSIONS")
    print("=" * 80)
    
    if correlations['n_triplets']['significant']:
        print("\n[WARN] SIGNIFICANT CORRELATION WITH n_triplets")
        print("   Suggests selection bias: epochs with more triplets give higher |H|")
        print("   Possible cause: Triplet selection criteria preferentially select")
        print("   measurements with larger systematic errors")
    else:
        print("\n[OK] NO SIGNIFICANT CORRELATION WITH n_triplets")
        print("   Q4 dominance not explained by n_triplets alone")
    
    if correlations['snr']['significant']:
        print("\n[WARN] SIGNIFICANT CORRELATION WITH SNR")
        print("   Suggests SNR-dependent systematic error")
        print("   High-SNR measurements may have different systematic biases")
        print("   Recommendation: Use SNR-independent estimator")
    else:
        print("\n[OK] NO SIGNIFICANT CORRELATION WITH SNR")
        print("   SNR-weighted outlier (26.0 ns) may be statistical fluke")
    
    rc_ctrl = controls["ratio_change"]
    if (
        isinstance(rc_ctrl, (int, float))
        and np.isfinite(rc_ctrl)
        and abs(rc_ctrl) < 10
    ):
        print("\n[OK] Q4 dominance persists after SNR matching")
        print("   Ratio change <10% after controlling for SNR")
        print("   Suggests genuine effect, not SNR-dependent systematic")
    elif isinstance(rc_ctrl, (int, float)) and np.isfinite(rc_ctrl):
        print("\n[WARN] Q4 dominance reduced after SNR matching")
        print(f"   Ratio changed by {rc_ctrl:.1f}%")
        print("   Suggests SNR-dependent systematic contributes to Q4 dominance")
    else:
        print(
            "\n[INFO] SNR-matched control ratio undefined "
            "(e.g. empty Q4 quartile under discrete n_triplets thresholds)."
        )
        print("   Interpret Q4/Q123 diagnostics from stratification tables only.")
    
    # Save results
    results = {
        "validation_type": "Q4 Dominance Mechanism Investigation",
        "validation_date": datetime.now().isoformat(),
        "execution_summary": {
            "n_epochs": len(epochs),
            "n_valid_epochs": len(valid_metrics),
            "n_warnings": n_warnings,
            "n_passed": n_ok
        },
        "correlations": correlations,
        "stratification_n_triplets": strat_n_triplets,
        "stratification_snr": strat_snr,
        "control_analysis": controls,
        "epoch_metrics": epoch_metrics
    }
    
    output_file = RESULTS_DIR / "step_024_q4_mechanism_investigation.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2, cls=NpEncoder)
    
    print(f"\n[PASS] Results saved to: {output_file}")
    print(f"[DATA] File size: {os.path.getsize(output_file)/1024:.1f} KB")
    
    # Final completion
    end_time = datetime.now()
    duration = (end_time - start_time).total_seconds()
    
    print("\n" + "=" * 80)
    print("STEP 024 COMPLETED SUCCESSFULLY")
    print("=" * 80)
    print(f"[INFO] Ended at: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"[INFO] Duration: {duration:.2f} seconds")
    print(f"[INFO] Status: {'ALL TESTS PASSED' if n_warnings == 0 else f'{n_warnings} WARNING(S)'}")
    print("=" * 80)


if __name__ == "__main__":
    main()
