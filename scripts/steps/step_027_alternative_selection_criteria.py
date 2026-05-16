#!/usr/bin/env python3
"""
================================================================================
STEP 029: ALTERNATIVE SELECTION CRITERIA TEST
================================================================================

Purpose: Test sensitivity to selection criteria after Q4 dominance was de-artifacted.
(Note: The strongest Q4 effect was identified as a bug in unsigned-delay analysis.)

Method:
1. Define alternative triplet selection criteria:
   - Max 10 triplets per epoch (instead of 20)
   - Random selection (instead of best SNR)
   - Equal SNR bins (instead of top SNR)
   - Minimum SNR threshold (instead of max n_triplets)
2. Re-compute signed geometric H with each criterion
3. Compare Q4 dominance across criteria

Expected outcomes:
- Post de-artifacting: Signed geometric H should remain stable across criteria
- Simple triplet-count bias was not supported (r = −0.006, p = 0.92)
- Threshold sensitivity remains mixed; this test characterizes selection robustness

This evaluates selection criteria sensitivity after Q4 dominance was resolved as
an unsigned-delay artifact (not a genuine physical effect).

================================================================================
"""

import json
import sys
import numpy as np
from scripts.utils.config import RANDOM_SEED
from scripts.utils.logger import print_status
from pathlib import Path
from scipy import stats
from datetime import datetime
from typing import Dict, Any, List, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.json_numpy import NpEncoder

RESULTS_DIR = PROJECT_ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def load_epoch_data() -> List[Dict]:
    """Load per-epoch closure delay data."""
    per_epoch_file = PROJECT_ROOT / "results" / "step_003_closure_final_per_epoch.json"
    
    with open(per_epoch_file, 'r') as f:
        epochs = json.load(f)
    
    return epochs


def compute_H_from_triplets(triplets: List[Dict]) -> float:
    """Compute the authoritative signed geometric H statistic from a list of triplets."""
    delays = np.array([t.get("geometric_delta_us", t.get("delta_us", 0)) for t in triplets])
    H = np.mean(delays) * 1e3  # ns
    return H


def apply_criterion_max_n_triplets(epoch: Dict, max_triplets: int) -> List[Dict]:
    """Select top N triplets by SNR."""
    triplets = epoch.get("triplets", [])
    
    # Sort by SNR
    sorted_triplets = sorted(triplets, key=lambda t: t.get("snr", 0), reverse=True)
    
    # Select top N
    selected = sorted_triplets[:max_triplets]
    
    return selected


def apply_criterion_random_selection(epoch: Dict, n_triplets: int, rng: np.random.Generator) -> List[Dict]:
    """Randomly select N triplets."""
    triplets = epoch.get("triplets", [])
    
    if len(triplets) <= n_triplets:
        return triplets
    
    # Random selection with fixed seed for reproducibility
    selected_indices = rng.choice(len(triplets), size=n_triplets, replace=False)
    selected = [triplets[int(idx)] for idx in selected_indices]
    
    return selected


def apply_criterion_snr_threshold(epoch: Dict, min_snr: float) -> List[Dict]:
    """Select all triplets above SNR threshold."""
    triplets = epoch.get("triplets", [])
    
    # Filter by SNR threshold
    selected = [t for t in triplets if t.get("snr", 0) >= min_snr]
    
    return selected


def apply_criterion_equal_snr_bins(epoch: Dict, n_bins: int, n_per_bin: int, rng: np.random.Generator) -> List[Dict]:
    """Select equal numbers from each SNR bin."""
    triplets = epoch.get("triplets", [])
    
    if len(triplets) == 0:
        return []
    
    # Sort by SNR
    sorted_triplets = sorted(triplets, key=lambda t: t.get("snr", 0))
    
    # Divide into bins
    bin_size = len(sorted_triplets) // n_bins
    selected = []
    
    for i in range(n_bins):
        start = i * bin_size
        end = start + bin_size
        bin_triplets = sorted_triplets[start:end]
        
        # Select n_per_bin from this bin (or all if fewer)
        if len(bin_triplets) <= n_per_bin:
            selected.extend(bin_triplets)
        else:
            # Random selection from bin with fixed seed for reproducibility
            selected_indices = rng.choice(len(bin_triplets), size=n_per_bin, replace=False)
            bin_selected = [bin_triplets[int(idx)] for idx in selected_indices]
            selected.extend(bin_selected)
    
    return selected


def evaluate_criterion(epochs: List[Dict], criterion_func, criterion_name: str, **kwargs) -> Dict[str, Any]:
    """
    Evaluate a selection criterion across all epochs.
    
    Parameters
    ----------
    epochs : list of dict
        Epoch data
    criterion_func : function
        Function to apply selection criterion
    criterion_name : str
        Name of the criterion
    **kwargs : dict
        Arguments for criterion function
    
    Returns
    -------
    dict
        Criterion evaluation results
    """
    epoch_results = []
    pooled_triplets = []
    rng = np.random.default_rng(kwargs.pop("random_seed", RANDOM_SEED))
    
    for epoch in epochs:
        # Apply criterion
        criterion_kwargs = dict(kwargs)
        if criterion_func in (apply_criterion_random_selection, apply_criterion_equal_snr_bins):
            criterion_kwargs["rng"] = rng
        selected_triplets = criterion_func(epoch, **criterion_kwargs)
        
        if len(selected_triplets) == 0:
            continue
        pooled_triplets.extend(selected_triplets)
        
        # Compute epoch-level signed H
        H = compute_H_from_triplets(selected_triplets)
        
        epoch_results.append({
            "mjd": epoch.get("mjd", 0),
            "n_triplets_original": epoch.get("n_triplets", 0),
            "n_triplets_selected": len(selected_triplets),
            "H_ns": H
        })
    
    if len(epoch_results) == 0:
        return {"error": "No epochs with selected triplets"}

    pooled_H_ns = compute_H_from_triplets(pooled_triplets)
    pooled_delays_ns = np.array([t.get("geometric_delta_us", t.get("delta_us", 0)) * 1e3 for t in pooled_triplets])
    pooled_sem_ns = float(np.std(pooled_delays_ns, ddof=1) / np.sqrt(len(pooled_delays_ns))) if len(pooled_delays_ns) > 1 else 0.0
    pooled_t = float(abs(pooled_H_ns) / pooled_sem_ns) if pooled_sem_ns > 0 else 0.0
    
    # Compute Q4 dominance against original epoch richness rather than the
    # selected-triplet count, which can become degenerate for fixed-cap criteria.
    n_triplets = np.array([e["n_triplets_original"] for e in epoch_results])
    H_values = np.array([e["H_ns"] for e in epoch_results])
    
    # Quartiles by n_triplets
    q3_threshold = np.percentile(n_triplets, 75)
    
    q4_indices = n_triplets > q3_threshold
    q123_indices = n_triplets <= q3_threshold
    
    H_q4 = np.mean(H_values[q4_indices]) if np.any(q4_indices) else 0
    H_q123 = np.mean(H_values[q123_indices]) if np.any(q123_indices) else 0
    q4_ratio = abs(H_q4) / abs(H_q123) if H_q123 != 0 else float('nan')
    
    # Correlation with n_triplets
    r, p = stats.pearsonr(n_triplets, H_values)
    
    return {
        "criterion": criterion_name,
        "n_epochs": len(epoch_results),
        "n_triplets_pooled": len(pooled_triplets),
        "pooled_H_ns": float(pooled_H_ns),
        "pooled_sem_ns": pooled_sem_ns,
        "pooled_t_statistic": pooled_t,
        "mean_H_ns": float(np.mean(H_values)),
        "std_H_ns": float(np.std(H_values)),
        "q4_H_ns": float(H_q4),
        "q123_H_ns": float(H_q123),
        "q4_q123_ratio": float(q4_ratio),
        "n_triplets_correlation": float(r),
        "n_triplets_p_value": float(p),
        "epoch_results": epoch_results
    }


def main():
    print_status("===" * 80)
    print("STEP 029: ALTERNATIVE SELECTION CRITERIA TEST")
    print_status("===" * 80)
    print_status("Purpose: Test if Q4 dominance is due to current selection criteria")
    print("Method: Compare signed geometric H and Q4 richness diagnostics across different selection criteria")
    print(f"Randomized criteria use fixed seed = {RANDOM_SEED}")
    
    # Load data
    print_status("Loading epoch data...")
    epochs = load_epoch_data()
    print_status(f"Loaded {len(epochs)} epochs")
    
    # Define criteria to test
    criteria = [
        {
            "name": "Original (max 20, top SNR)",
            "func": apply_criterion_max_n_triplets,
            "kwargs": {"max_triplets": 20}
        },
        {
            "name": "Max 10 triplets (top SNR)",
            "func": apply_criterion_max_n_triplets,
            "kwargs": {"max_triplets": 10}
        },
        {
            "name": "Max 5 triplets (top SNR)",
            "func": apply_criterion_max_n_triplets,
            "kwargs": {"max_triplets": 5}
        },
        {
            "name": "Random 10 triplets",
            "func": apply_criterion_random_selection,
            "kwargs": {"n_triplets": 10}
        },
        {
            "name": "Random 5 triplets",
            "func": apply_criterion_random_selection,
            "kwargs": {"n_triplets": 5}
        },
        {
            "name": "SNR threshold 5.0",
            "func": apply_criterion_snr_threshold,
            "kwargs": {"min_snr": 5.0}
        },
        {
            "name": "SNR threshold 3.0",
            "func": apply_criterion_snr_threshold,
            "kwargs": {"min_snr": 3.0}
        },
        {
            "name": "Equal SNR bins (3 bins, 3 per bin)",
            "func": apply_criterion_equal_snr_bins,
            "kwargs": {"n_bins": 3, "n_per_bin": 3}
        }
    ]
    
    # Evaluate each criterion
    print_status("Evaluating criteria...")
    results = []
    
    for criterion in criteria:
        print(f"\n  Testing: {criterion['name']}")
        result = evaluate_criterion(epochs, criterion["func"], criterion["name"], **criterion["kwargs"])
        
        if "error" in result:
            print(f"    Error: {result['error']}")
            continue
        
        print(f"    Pooled H: {result['pooled_H_ns']:+.2f} +/- {result['pooled_sem_ns']:.2f} ns (t = {result['pooled_t_statistic']:.2f})")
        print(f"    Mean epoch H: {result['mean_H_ns']:+.1f} ns")
        print(f"    Q4/Q123 ratio: {result['q4_q123_ratio']:.2f}x")
        print(f"    n_triplets correlation: r = {result['n_triplets_correlation']:.3f}")
        
        results.append(result)
    
    # Compare results
    print_status("" + "=" * 80)
    print("CRITERION COMPARISON")
    print_status("===" * 80)
    
    print(f"\n{'Criterion':<40} {'Pooled H':<18} {'Q4/Q123':<10} {'n_triplets r':<12}")
    print("-" * 80)
    
    for result in results:
        pooled_summary = f"{result['pooled_H_ns']:+.2f}+/-{result['pooled_sem_ns']:.2f}"
        print(f"{result['criterion']:<40} {pooled_summary:<18} {result['q4_q123_ratio']:<10.2f} {result['n_triplets_correlation']:<12.3f}")
    
    # Analysis
    print_status("" + "=" * 80)
    print("ANALYSIS")
    print_status("===" * 80)
    
    q4_ratios = [r["q4_q123_ratio"] for r in results]
    correlations = [r["n_triplets_correlation"] for r in results]
    pooled_H_values = [r["pooled_H_ns"] for r in results]
    
    q4_ratio_range = max(q4_ratios) - min(q4_ratios)
    correlation_range = max(correlations) - min(correlations)
    pooled_H_range = max(pooled_H_values) - min(pooled_H_values)
    
    print(f"\nQ4/Q123 ratio range: {q4_ratio_range:.2f}x")
    print(f"n_triplets correlation range: {correlation_range:.3f}")
    print(f"Pooled signed H range: {pooled_H_range:.2f} ns")
    
    if q4_ratio_range < 0.5:
        print("\n[OK] Q4 dominance CONSISTENT across criteria")
        print("   Suggests epoch-richness contrast is not strongly criterion-dependent")
    elif q4_ratio_range > 1.0:
        print("\n[WARN] Q4 dominance VARIES widely across criteria")
        print("   Suggests criterion choice materially affects the richness contrast")
    else:
        print("\n[WARN] Q4 dominance MODERATE variation across criteria")
        print("   Suggests partial criterion sensitivity")
    
    if correlation_range < 0.1:
        print("\n[OK] n_triplets correlation CONSISTENT across criteria")
        print("   Suggests the triplet-richness diagnostic is not highly criterion-specific")
    elif correlation_range > 0.3:
        print("\n[WARN] n_triplets correlation VARIES across criteria")
        print("   Suggests the triplet-richness diagnostic is criterion-sensitive")
    else:
        print("\n[WARN] n_triplets correlation MODERATE variation")

    if pooled_H_range < 0.5:
        print("\n[OK] Pooled signed H remains broadly stable across tested criteria")
    else:
        print("\n[WARN] Pooled signed H changes materially across tested criteria")
        print("   This indicates non-negligible selection sensitivity in the pooled estimator")
    
    # Identify best criterion (minimizes Q4 dominance and n_triplets correlation)
    print_status("" + "=" * 80)
    print("BEST CRITERION (minimizes Q4 dominance and n_triplets correlation)")
    print_status("===" * 80)
    
    # Score each criterion (lower is better)
    scores = []
    for result in results:
        score = abs(result["q4_q123_ratio"] - 1.0) + abs(result["n_triplets_correlation"])
        scores.append((score, result))
    
    scores.sort(key=lambda item: (item[0], item[1]["criterion"]))
    best_score, best_result = scores[0]
    
    print(f"\nBest criterion: {best_result['criterion']}")
    print_status(f"Score: {best_score:.3f}")
    print_status(f"Q4/Q123 ratio: {best_result['q4_q123_ratio']:.2f}x")
    print_status(f"n_triplets correlation: r = {best_result['n_triplets_correlation']:.3f}")
    
    # Save results
    output = {
        "validation_type": "Alternative Selection Criteria Test",
        "validation_date": datetime.now().isoformat(),
        "criteria_tested": len(criteria),
        "results": results,
        "analysis": {
            "q4_ratio_range": float(q4_ratio_range),
            "correlation_range": float(correlation_range),
            "pooled_H_range_ns": float(pooled_H_range),
            "best_criterion": best_result["criterion"],
            "best_criterion_score": float(best_score)
        },
        "random_seed": RANDOM_SEED
    }
    
    output_file = RESULTS_DIR / "step_029_alternative_selection_criteria.json"
    with open(output_file, 'w') as f:
        json.dump(output, f, indent=2, cls=NpEncoder)
    
    print_status(f"Results saved to: {output_file}")
    
    print_status("" + "=" * 80)
    print("STEP 029 COMPLETED")
    print_status("===" * 80)


if __name__ == "__main__":
    main()
