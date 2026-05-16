#!/usr/bin/env python3
"""
================================================================================
STEP 027: RANDOM TRIPLET SUBSET ANALYSIS
================================================================================

Purpose: Test if n_triplets correlation (r=0.409) is due to selection bias
or genuine physical sampling effect.

Method:
1. For each epoch, randomly select N triplets (where N = min_triplets)
2. Compute |H| using only these N triplets
3. Repeat 1000 times to get bootstrap distribution
4. Compare with original |H| using all available triplets

Expected outcomes:
- If selection bias: H_ratio < 1 for high-triplet epochs (|H decreases)
- If physical effect: H_ratio ≈ 1 for all epochs (|H consistent)
- If statistical artifact: H_ratio < 1 due to larger SEM, but not mean

This directly tests whether the n_triplets correlation is a selection artifact
or a genuine physical effect.

================================================================================
"""

import json
import numpy as np
from pathlib import Path
from scipy import stats
from datetime import datetime
from typing import Optional, Dict, Any, List
from scripts.utils.json_numpy import NpEncoder

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def load_epoch_data() -> List[Dict]:
    """Load per-epoch closure delay data."""
    per_epoch_file = PROJECT_ROOT / "results" / "step_003_closure_final_per_epoch.json"
    
    with open(per_epoch_file, 'r') as f:
        epochs = json.load(f)
    
    return epochs


def compute_H_from_triplets(triplets: List[Dict]) -> float:
    """Compute |H| from a list of triplets."""
    delays = np.array([t.get("geometric_delta_us", t.get("delta_us", 0)) for t in triplets])
    H = abs(np.mean(delays)) * 1e3  # ns
    return H


def random_subset_analysis(epochs: List[Dict], n_bootstrap: int = 1000) -> List[Dict]:
    """
    Perform random triplet subset analysis for each epoch.
    
    Parameters
    ----------
    epochs : list of dict
        Epoch data with triplets
    n_bootstrap : int
        Number of bootstrap samples
    
    Returns
    -------
    list of dict
        Epoch results with subset analysis
    """
    # Find minimum number of triplets across all epochs
    min_triplets = min(e["n_triplets"] for e in epochs)
    print(f"Minimum triplets per epoch: {min_triplets}")
    print(f"Maximum triplets per epoch: {max(e['n_triplets'] for e in epochs)}")
    
    results = []
    
    for i, epoch in enumerate(epochs):
        triplets = epoch.get("triplets", [])
        n_triplets = len(triplets)
        
        if n_triplets == 0:
            continue
        
        # Original |H| using all triplets
        H_original = compute_H_from_triplets(triplets)
        
        # Bootstrap random subsets
        subset_H = []
        subset_size = min(min_triplets, n_triplets)
        
        for _ in range(n_bootstrap):
            if subset_size >= n_triplets:
                # If epoch has fewer triplets than min, use all
                subset = triplets
            else:
                # Random subset
                subset = np.random.choice(triplets, size=subset_size, replace=False)
            
            H_subset = compute_H_from_triplets(subset)
            subset_H.append(H_subset)
        
        # Statistics
        H_subset_mean = np.mean(subset_H)
        H_subset_std = np.std(subset_H)
        H_subset_ci_lower = np.percentile(subset_H, 2.5)
        H_subset_ci_upper = np.percentile(subset_H, 97.5)
        
        # Ratio
        H_ratio = H_subset_mean / H_original if H_original > 0 else 1.0
        
        results.append({
            "epoch_id": i,
            "mjd": epoch.get("mjd", 0),
            "n_triplets": n_triplets,
            "H_original_ns": H_original,
            "H_subset_mean_ns": H_subset_mean,
            "H_subset_std_ns": H_subset_std,
            "H_subset_ci_lower_ns": H_subset_ci_lower,
            "H_subset_ci_upper_ns": H_subset_ci_upper,
            "H_ratio": H_ratio,
            "subset_size": subset_size
        })
        
        if (i + 1) % 50 == 0:
            print(f"  Processed {i + 1}/{len(epochs)} epochs")
    
    return results


def analyze_subset_results(results: List[Dict]) -> Dict[str, Any]:
    """
    Analyze the random subset results.
    
    Parameters
    ----------
    results : list of dict
        Epoch subset analysis results
    
    Returns
    -------
    dict
        Analysis results
    """
    # Extract data
    n_triplets = np.array([r["n_triplets"] for r in results])
    H_original = np.array([r["H_original_ns"] for r in results])
    H_subset_mean = np.array([r["H_subset_mean_ns"] for r in results])
    H_ratio = np.array([r["H_ratio"] for r in results])
    
    # Correlation analysis
    r_n_triplets_original, p_n_triplets_original = stats.pearsonr(n_triplets, H_original)
    r_n_triplets_subset, p_n_triplets_subset = stats.pearsonr(n_triplets, H_subset_mean)
    r_n_triplets_ratio, p_n_triplets_ratio = stats.pearsonr(n_triplets, H_ratio)
    
    # Stratify by n_triplets quartiles
    q1_threshold = np.percentile(n_triplets, 25)
    q2_threshold = np.percentile(n_triplets, 50)
    q3_threshold = np.percentile(n_triplets, 75)
    
    groups = {
        "Q1_lowest": [r for r in results if r["n_triplets"] <= q1_threshold],
        "Q2": [r for r in results if q1_threshold < r["n_triplets"] <= q2_threshold],
        "Q3": [r for r in results if q2_threshold < r["n_triplets"] <= q3_threshold],
        "Q4_highest": [r for r in results if r["n_triplets"] > q3_threshold]
    }
    
    group_stats = {}
    for group_name, group_results in groups.items():
        group_H_ratio = np.array([r["H_ratio"] for r in group_results])
        group_n_triplets = np.array([r["n_triplets"] for r in group_results])
        
        group_stats[group_name] = {
            "n_epochs": len(group_results),
            "mean_n_triplets": float(np.mean(group_n_triplets)),
            "mean_H_ratio": float(np.mean(group_H_ratio)),
            "std_H_ratio": float(np.std(group_H_ratio)),
            "median_H_ratio": float(np.median(group_H_ratio)),
            "fraction_H_ratio_lt_1": float(np.mean(group_H_ratio < 1.0))
        }
    
    return {
        "correlations": {
            "n_triplets_vs_H_original": {
                "r": float(r_n_triplets_original),
                "p": float(p_n_triplets_original)
            },
            "n_triplets_vs_H_subset": {
                "r": float(r_n_triplets_subset),
                "p": float(p_n_triplets_subset)
            },
            "n_triplets_vs_H_ratio": {
                "r": float(r_n_triplets_ratio),
                "p": float(p_n_triplets_ratio)
            }
        },
        "group_statistics": group_stats,
        "overall_statistics": {
            "mean_H_ratio": float(np.mean(H_ratio)),
            "std_H_ratio": float(np.std(H_ratio)),
            "median_H_ratio": float(np.median(H_ratio)),
            "fraction_H_ratio_lt_1": float(np.mean(H_ratio < 1.0))
        }
    }


def main():
    print("=" * 80)
    print("STEP 027: RANDOM TRIPLET SUBSET ANALYSIS")
    print("=" * 80)
    print("\nPurpose: Test if n_triplets correlation is selection bias or physical effect")
    print("Method: Randomly select subsets of triplets and |H|")
    
    # Load data
    print("\nLoading epoch data...")
    epochs = load_epoch_data()
    print(f"Loaded {len(epochs)} epochs")
    
    # Perform random subset analysis
    print("\nPerforming random subset analysis (1000 bootstrap samples)...")
    results = random_subset_analysis(epochs, n_bootstrap=1000)
    print(f"Analyzed {len(results)} epochs")
    
    # Analyze results
    print("\nAnalyzing subset results...")
    analysis = analyze_subset_results(results)
    
    # Print results
    print("\n" + "=" * 80)
    print("RESULTS")
    print("=" * 80)
    
    print(f"\nCorrelation analysis:")
    print(f"  n_triplets vs H_original: r = {analysis['correlations']['n_triplets_vs_H_original']['r']:.3f}, p = {analysis['correlations']['n_triplets_vs_H_original']['p']:.2e}")
    print(f"  n_triplets vs H_subset: r = {analysis['correlations']['n_triplets_vs_H_subset']['r']:.3f}, p = {analysis['correlations']['n_triplets_vs_H_subset']['p']:.2e}")
    print(f"  n_triplets vs H_ratio: r = {analysis['correlations']['n_triplets_vs_H_ratio']['r']:.3f}, p = {analysis['correlations']['n_triplets_vs_H_ratio']['p']:.2e}")
    
    print(f"\nOverall statistics:")
    print(f"  Mean H_ratio: {analysis['overall_statistics']['mean_H_ratio']:.3f}")
    print(f"  Std H_ratio: {analysis['overall_statistics']['std_H_ratio']:.3f}")
    print(f"  Median H_ratio: {analysis['overall_statistics']['median_H_ratio']:.3f}")
    print(f"  Fraction H_ratio < 1: {analysis['overall_statistics']['fraction_H_ratio_lt_1']*100:.1f}%")
    
    print(f"\nStratified by n_triplets:")
    for group_name, stats in analysis["group_statistics"].items():
        print(f"  {group_name}:")
        print(f"    n_epochs: {stats['n_epochs']}")
        print(f"    mean_n_triplets: {stats['mean_n_triplets']:.1f}")
        print(f"    mean_H_ratio: {stats['mean_H_ratio']:.3f}")
        print(f"    fraction_H_ratio_lt_1: {stats['fraction_H_ratio_lt_1']*100:.1f}%")
    
    # Interpretation
    print("\n" + "=" * 80)
    print("INTERPRETATION")
    print("=" * 80)
    
    r_ratio = analysis['correlations']['n_triplets_vs_H_ratio']['r']
    p_ratio = analysis['correlations']['n_triplets_vs_H_ratio']['p']
    
    q4_fraction = analysis["group_statistics"]["Q4_highest"]["fraction_H_ratio_lt_1"]
    
    if r_ratio < -0.2 and p_ratio < 0.05:
        print("\n[WARN] SELECTION BIAS DETECTED")
        print(f"   Negative correlation between n_triplets and H_ratio (r = {r_ratio:.3f})")
        print("   High-triplet epochs show H_ratio < 1 (|H| decreases with subsets)")
        print("   Suggests current selection criteria preferentially select larger delays")
        print("   Recommendation: Change selection criteria to be SNR-independent or random")
    elif r_ratio > 0.2 and p_ratio < 0.05:
        print("\n[OK] PHYSICAL EFFECT DETECTED")
        print(f"   Positive correlation between n_triplets and H_ratio (r = {r_ratio:.3f})")
        print("   High-triplet epochs show H_ratio > 1 (|H| increases with subsets)")
        print("   Suggests more triplets = better ISM sampling = more accurate |H|")
        print("   Recommendation: Report as genuine physical effect")
    else:
        print("\n[OK] NO SIGNIFICANT CORRELATION")
        print(f"   n_triplets vs H_ratio: r = {r_ratio:.3f}, p = {p_ratio:.2e}")
        print("   H_ratio consistent across all epochs regardless of n_triplets")
        print("   Suggests n_triplets correlation is not due to selection bias")
        print("   Recommendation: Original n_triplets correlation may be physical or statistical")
    
    if q4_fraction > 0.7:
        print(f"\n[WARN] Q4 EPOCHS SHOW REDUCTION")
        print(f"   {q4_fraction*100:.1f}% of Q4 epochs have H_ratio < 1")
        print("   Suggests Q4 dominance is due to selection bias")
    elif q4_fraction < 0.3:
        print(f"\n[OK] Q4 EPOCHS CONSISTENT")
        print(f"   Only {q4_fraction*100:.1f}% of Q4 epochs have H_ratio < 1")
        print("   Suggests Q4 dominance is not due to selection bias")
    
    # Save results
    output = {
        "validation_type": "Random Triplet Subset Analysis",
        "validation_date": datetime.now().isoformat(),
        "n_bootstrap": 1000,
        "analysis": analysis,
        "epoch_results": results
    }
    
    output_file = RESULTS_DIR / "step_027_random_triplet_subset_analysis.json"
    with open(output_file, 'w') as f:
        json.dump(output, f, indent=2, cls=NpEncoder)
    
    print(f"\nResults saved to: {output_file}")
    
    print("\n" + "=" * 80)
    print("STEP 027 COMPLETED")
    print("=" * 80)


if __name__ == "__main__":
    main()
