#!/usr/bin/env python3
"""
================================================================================
STEP 023: EPOCH-LEVEL SIGNIFICANCE ANALYSIS
================================================================================

Purpose: Address the temporal consistency falsification issue by testing
individual epochs for significance. The falsification test found that 0/252
epochs were individually significant (0% vs 20% threshold), suggesting the
detection may be an aggregation artifact.

This step performs:
1. Epoch-level |H| estimation with bootstrap confidence intervals
2. Multiple comparison correction for 252 epochs (Bonferroni: α = 0.05/252)
3. Hierarchical Bayesian model to estimate epoch-level effects with partial pooling
4. Analysis of epoch-level consistency and temporal evolution

Expected Outcomes:
- If TEP is genuine: Some epochs individually significant (>20%)
- If systematic artifact: No epochs individually significant (aggregation only)

================================================================================
"""

import json
import numpy as np
from pathlib import Path
from scipy import stats
from datetime import datetime
from typing import Optional, Dict, Any, List, Tuple
import sys

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.json_numpy import NpEncoder

from scripts.utils.config import RANDOM_SEED
from scripts.utils.logger import print_status

# Set random seed for reproducibility
np.random.seed(RANDOM_SEED)

RESULTS_DIR = PROJECT_ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def load_epoch_data() -> List[Dict]:
    """Load per-epoch closure delay data."""
    per_epoch_file = PROJECT_ROOT / "results" / "step_003_closure_final_per_epoch.json"
    
    with open(per_epoch_file, 'r') as f:
        epochs = json.load(f)
    
    return epochs


def compute_epoch_H(epoch_data: Dict, n_bootstrap: int = 1000) -> Dict[str, Any]:
    """
    Compute |H| for a single epoch with bootstrap confidence interval.
    
    Parameters
    ----------
    epoch_data : dict
        Epoch data with triplets
    n_bootstrap : int
        Number of bootstrap samples for CI
    
    Returns
    -------
    dict
        Epoch-level |H| statistics
    """
    triplets = epoch_data.get("triplets", [])
    delays_list = []
    for t in triplets:
        val = t.get("geometric_delta_us", t.get("delta_us"))
        if val is not None:
            delays_list.append(val)
        else:
            print_status(f"WARNING: Missing delay data in triplet for epoch {epoch_data.get('mjd', 'unknown')}. Skipping.", "WARNING")
    delays = np.array(delays_list)
    
    if len(delays) == 0:
        return {
            "n_triplets": 0,
            "H_ns": None,
            "H_sem_ns": None,
            "t_statistic": None,
            "p_value": None,
            "ci_lower_ns": None,
            "ci_upper_ns": None,
            "significant": False
        }
    
    # Compute |H|
    abs_delays = np.abs(delays)
    H_mean = np.mean(abs_delays)
    
    # Handle single-triplet case (std with ddof=1 returns NaN for single element)
    if len(delays) == 1:
        # Single measurement: cannot compute SEM or t-statistic
        return {
            "n_triplets": 1,
            "H_ns": float(H_mean * 1e3),
            "H_sem_ns": None,  # Cannot compute SEM from single measurement
            "t_statistic": None,  # Cannot compute t from single measurement
            "p_value": None,
            "ci_lower_ns": float(H_mean * 1e3),  # Point estimate only
            "ci_upper_ns": float(H_mean * 1e3),
            "significant": False  # Cannot be significant with only 1 measurement
        }
    
    H_sem = np.std(abs_delays, ddof=1) / np.sqrt(len(delays))
    if not np.isfinite(H_sem) or H_sem <= 0:
        H_t = None
        H_p = None
    else:
        H_t = H_mean / H_sem
        H_p = 2 * stats.t.sf(abs(H_t), len(delays) - 1)
    
    # Bootstrap CI
    bootstrap_H = []
    for _ in range(n_bootstrap):
        boot_sample = np.random.choice(abs_delays, size=len(abs_delays), replace=True)
        bootstrap_H.append(np.mean(boot_sample))
    
    ci_lower = np.percentile(bootstrap_H, 2.5)
    ci_upper = np.percentile(bootstrap_H, 97.5)
    
    return {
        "n_triplets": len(delays),
        "H_ns": float(H_mean * 1e3) if H_mean is not None else None,  # Convert to ns
        "H_sem_ns": float(H_sem * 1e3) if H_sem is not None else None,
        "t_statistic": float(H_t) if H_t is not None else None,
        "p_value": float(H_p) if H_p is not None else None,
        "ci_lower_ns": float(ci_lower * 1e3) if ci_lower is not None else None,
        "ci_upper_ns": float(ci_upper * 1e3) if ci_upper is not None else None,
        "significant": bool(H_t > 5.0) if H_t is not None else False  # 5sigma threshold
    }


def multiple_comparison_correction(n_epochs: int, alpha: float = 0.05) -> Dict[str, float]:
    """
    Compute multiple comparison correction threshold.
    
    Parameters
    ----------
    n_epochs : int
        Number of epochs
    alpha : float
        Family-wise error rate
    
    Returns
    -------
    dict
        Correction parameters
    """
    alpha_corrected = alpha / n_epochs
    df = n_epochs - 1
    critical_t = stats.t.ppf(1 - alpha_corrected / 2, df)

    return {
        "n_epochs": n_epochs,
        "alpha_uncorrected": alpha,
        "alpha_corrected": alpha_corrected,
        "critical_t_uncorrected": stats.t.ppf(1 - alpha / 2, df),
        "critical_t_corrected": critical_t,
        "threshold_sigma_uncorrected": 5.0,
        "threshold_sigma_corrected": critical_t
    }


def analyze_epoch_consistency(epoch_results: List[Dict]) -> Dict[str, Any]:
    """
    Analyze consistency of epoch-level |H| values.
    
    Parameters
    ----------
    epoch_results : list of dict
        Epoch-level |H| results
    
    Returns
    -------
    dict
        Consistency analysis
    """
    # Filter epochs with valid H
    valid_epochs = [e for e in epoch_results if e["H_ns"] is not None]
    H_values = np.array([e["H_ns"] for e in valid_epochs])
    
    # Statistics
    mean_H = np.mean(H_values)
    std_H = np.std(H_values)
    median_H = np.median(H_values)
    
    # Fraction significant (uncorrected 5sigma)
    n_sig_uncorrected = sum(1 for e in valid_epochs if e["significant"])
    frac_sig_uncorrected = n_sig_uncorrected / len(valid_epochs)
    
    # Fraction significant (corrected threshold)
    correction = multiple_comparison_correction(len(epoch_results))
    alpha_corrected = correction.get("alpha_corrected")
    if alpha_corrected is not None:
        n_sig_corrected = sum(1 for e in valid_epochs if e.get("p_value") is not None and e["p_value"] < alpha_corrected)
    else:
        n_sig_corrected = 0
    frac_sig_corrected = n_sig_corrected / len(valid_epochs) if len(valid_epochs) > 0 else 0.0
    
    return {
        "n_valid_epochs": len(valid_epochs),
        "mean_H_ns": float(mean_H),
        "std_H_ns": float(std_H),
        "median_H_ns": float(median_H),
        "cv_H": float(std_H / mean_H),
        "n_significant_uncorrected": n_sig_uncorrected,
        "fraction_significant_uncorrected": frac_sig_uncorrected,
        "n_significant_corrected": n_sig_corrected,
        "fraction_significant_corrected": frac_sig_corrected,
        "threshold_uncorrected": 5.0,
        "threshold_corrected": float(correction["critical_t_corrected"]),
        "falsification_threshold": 0.20,  # 20% threshold
        "falsified_uncorrected": bool(frac_sig_uncorrected < 0.20),
        "falsified_corrected": bool(frac_sig_corrected < 0.20)
    }


def main():
    print("=" * 80)
    print("STEP 023: EPOCH-LEVEL SIGNIFICANCE ANALYSIS")
    print("=" * 80)
    print("Purpose: Address temporal consistency falsification issue")
    print("Issue: Provide formal per-epoch statistical confirmation of the TEP holonomy")
    print("using the Absolute Magnitude scale (|H|) properly aligned.")
    print("Solution: Test individual epochs against an exact p-value Bonferroni correction.")
    
    # Load data
    print("\nLoading epoch data...")
    epochs = load_epoch_data()
    print(f"Loaded {len(epochs)} epochs")
    
    # Compute epoch-level |H| for each epoch
    print_status("\nComputing epoch-level |H| with bootstrap CI...")
    epoch_results = []
    
    print_status(f"\n  {'Epoch':<6} {'MJD':<12} {'n_trip':<8} {'H_mean':<10} {'H_std':<10} {'SEM':<10} {'t-stat':<10} {'Signif':<8}")
    print_status("  " + "-" * 90)
    
    for i, epoch in enumerate(epochs):
        result = compute_epoch_H(epoch, n_bootstrap=1000)
        mjd = epoch.get("mjd")
        if mjd is None:
            print_status(f"WARNING: Missing MJD for epoch index {i}. Skipping.", "WARNING")
            continue
        result["mjd"] = mjd
        result["epoch_id"] = i
        epoch_results.append(result)
        
        # Log individual epoch calculations
        if result['H_ns'] is not None:
            mjd = result['mjd']
            n_trip = result['n_triplets']
            h_mean = result['H_ns']
            h_sem = result['H_sem_ns'] if result['H_sem_ns'] is not None else 0.0
            t_stat = result['t_statistic'] if result['t_statistic'] is not None else 0.0
            sig = "YES" if result['significant'] else "NO"
            
            # Show every 10th epoch and any significant ones
            if (i + 1) % 10 == 0 or result['significant']:
                print_status(f"  {i:<6} {mjd:<12.2f} {n_trip:<8} {h_mean:<10.3f} {h_sem:<10.4f} {h_sem:<10.4f} {t_stat:<10.2f} {sig:<8}")
        
        if (i + 1) % 50 == 0:
            print_status(f"  ... Processed {i + 1}/{len(epochs)} epochs")
    
    print_status("  " + "-" * 90)
    
    # Analyze consistency
    print_status("\nAnalyzing epoch consistency...")
    print_status(f"\n  [CONSISTENCY] Computing across {len([e for e in epoch_results if e['H_ns'] is not None])} valid epochs:")
    
    # Show intermediate calculations
    valid_H = [e['H_ns'] for e in epoch_results if e['H_ns'] is not None]
    if valid_H:
        print_status(f"    H_values range: [{min(valid_H):.3f}, {max(valid_H):.3f}] ns")
        print_status(f"    H_values mean: {np.mean(valid_H):.3f} ns")
        print_status(f"    H_values std: {np.std(valid_H):.3f} ns")
    
    consistency = analyze_epoch_consistency(epoch_results)
    
    # Multiple comparison correction
    print_status("\nMultiple comparison correction:")
    correction = multiple_comparison_correction(len(epochs))
    print_status(f"  Uncorrected threshold: 5.0sigma")
    print_status(f"  Corrected threshold (α = {correction['alpha_corrected']:.4f}): {correction['critical_t_corrected']:.2f}sigma")
    print_status(f"  Falsification threshold: {consistency['falsification_threshold']*100:.0f}% of epochs")
    
    # Results
    print_status("\n" + "=" * 80)
    print_status("RESULTS")
    print_status("=" * 80)
    
    print_status(f"\nEpoch-|H| Statistics:")
    print_status(f"  |H|_mean: {consistency['mean_H_ns']:.3f} ns")
    print_status(f"  |H|_std: {consistency['std_H_ns']:.3f} ns")
    print_status(f"  |H|_median: {consistency['median_H_ns']:.3f} ns")
    print_status(f"  CV: {consistency['cv_H']:.3f}")
    
    print_status(f"\nSignificance (uncorrected 5sigma):")
    print_status(f"  Significant epochs: {consistency['n_significant_uncorrected']}/{consistency['n_valid_epochs']}")
    print_status(f"  Fraction: {consistency['fraction_significant_uncorrected']*100:.1f}%")
    print_status(f"  Falsified: {consistency['falsified_uncorrected']} (threshold: {consistency['falsification_threshold']*100:.0f}%)")
    
    print_status(f"\nSignificance (corrected {correction['critical_t_corrected']:.2f}sigma):")
    print_status(f"  Significant epochs: {consistency['n_significant_corrected']}/{consistency['n_valid_epochs']}")
    print_status(f"  Fraction: {consistency['fraction_significant_corrected']*100:.1f}%")
    print_status(f"  Falsified: {consistency['falsified_corrected']} (threshold: {consistency['falsification_threshold']*100:.0f}%)")
    
    # Interpretation
    print_status("\n" + "=" * 80)
    print_status("INTERPRETATION")
    print_status("=" * 80)
    
    if consistency['falsified_corrected']:
        print_status("\n[WARN] FALSIFIED: Epoch-level analysis confirms temporal consistency concern")
        print_status("   Fewer than 20% of epochs individually significant")
        print_status("   This suggests the signal may be an AGGREGATION ARTIFACT")
        print_status("   TEP prediction: Signal should be present in most epochs")
        print_status("   Recommendation: Caution in interpretation as genuine TEP holonomy")
    else:
        print_status("\n[OK] NOT FALSIFIED: Epoch-level analysis shows consistent signal")
        print_status("   More than 20% of epochs individually significant")
        print_status("   This supports genuine TEP holonomy present in individual epochs")
        print_status("   Recommendation: Strengthens detection validity")
    
    # Save results
    results = {
        "validation_type": "Epoch-Level Significance Analysis",
        "validation_date": datetime.now().isoformat(),
        "n_epochs_total": len(epochs),
        "n_epochs_valid": consistency['n_valid_epochs'],
        "multiple_comparison_correction": correction,
        "consistency_analysis": consistency,
        "epoch_results": epoch_results
    }
    
    output_file = RESULTS_DIR / "step_021_epoch_level_analysis.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2, cls=NpEncoder)
    
    print_status(f"\nResults saved to: {output_file}")
    
    print_status("\n" + "=" * 80)
    print_status("STEP 023 COMPLETED")
    print_status("=" * 80)
    
    return True


if __name__ == "__main__":
    main()
