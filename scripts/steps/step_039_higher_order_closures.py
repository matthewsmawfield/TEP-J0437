#!/usr/bin/env python3
"""
================================================================================
STEP 035: HIGHER-ORDER CLOSURE QUANTITIES
================================================================================

Purpose: Analyze higher-order closure quantities (quadruplets, etc.) to test
if TEP signatures persist beyond simple triplet closures.

Higher-Order Closures:
--------------------
- Triplet closure (3rd order): Δ = τ₀₁ + τ₁₂ - τ₀₂
- Quadruplet closure (4th order): Q = Δ₀₁₂ + Δ₁₂₃ - Δ₀₂₃
- General n-th order closure: Recursive construction

Theoretical Prediction:
----------------------
If TEP is a genuine phase transport effect:
- Higher-order closures should show consistent TEP signatures
- Magnitude should scale predictably with order
- Bipolar structure should persist across orders

If TEP is an artifact:
- Higher-order closures may show different behavior
- Signal may not persist beyond triplets
- Inconsistent patterns across orders

================================================================================
"""

import json
import numpy as np
import sys
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple
from itertools import combinations
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.json_numpy import NpEncoder

from scripts.utils.config import RANDOM_SEED
RESULTS_DIR = PROJECT_ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

def get_base_h():
    freq_file = PROJECT_ROOT / "results" / "step_003_closure_final_summary.json"
    if freq_file.exists():
        with open(freq_file, 'r') as f:
            data = json.load(f)
            return abs(data.get("H_trim_magnitude_ns", 1.07))
    return 1.07

BASE_H = get_base_h()


def generate_arclet_triplets(n_arclets: int) -> List[Tuple[int, int, int]]:
    """Generate all possible triplets from n arclets."""
    return list(combinations(range(n_arclets), 3))


def compute_triplet_closure(
    tau_values: np.ndarray,
    triplet: Tuple[int, int, int]
) -> float:
    """
    Compute triplet closure delay.
    
    Δ = τ_ij + τ_jk - τ_ik
    """
    i, j, k = triplet
    return tau_values[i] + tau_values[j] - tau_values[k]


def compute_quadruplet_closure(
    tau_values: np.ndarray,
    quadruplet: Tuple[int, int, int, int]
) -> float:
    """
    Compute quadruplet closure delay.
    
    Q = Δ_012 + Δ_123 - Δ_023
    
    Where each Δ is a triplet closure.
    """
    i, j, k, l = quadruplet
    
    # Compute component triplet closures
    delta_ijk = compute_triplet_closure(tau_values, (i, j, k))
    delta_jkl = compute_triplet_closure(tau_values, (j, k, l))
    delta_ikl = compute_triplet_closure(tau_values, (i, k, l))
    
    # Quadruplet closure
    Q = delta_ijk + delta_jkl - delta_ikl
    return Q


def analyze_triplet_closures(tau_values: np.ndarray) -> Dict[str, Any]:
    """Analyze triplet closures."""
    n_arclets = len(tau_values)
    triplets = generate_arclet_triplets(n_arclets)
    
    closures = []
    for triplet in triplets:
        closure = compute_triplet_closure(tau_values, triplet)
        closures.append(closure)
    
    closures = np.array(closures)
    abs_closures = np.abs(closures)
    
    # Statistics
    mean_H = np.mean(abs_closures)
    # Handle single-closure case (std with ddof=1 returns NaN)
    if len(closures) > 1:
        sem_H = np.std(abs_closures, ddof=1) / np.sqrt(len(closures))
        t_stat = mean_H / sem_H if sem_H > 0 else 0.0
    else:
        sem_H = None  # Cannot compute SEM from single measurement
        t_stat = 0.0
    
    # Bipolar structure
    neg_closures = closures[closures < 0]
    pos_closures = closures[closures > 0]
    
    if len(neg_closures) > 0 and len(pos_closures) > 0:
        neg_mean = np.mean(np.abs(neg_closures))
        pos_mean = np.mean(np.abs(pos_closures))
        magnitude_ratio = neg_mean / pos_mean
    else:
        neg_mean = pos_mean = 0
        magnitude_ratio = 0
    
    return {
        "order": 3,
        "n_closures": len(closures),
        "mean_H_ns": float(mean_H),
        "sem_ns": float(sem_H) if sem_H is not None else None,
        "t_statistic": float(t_stat),
        "detected": bool(t_stat > 5.0) if t_stat > 0 else False,
        "magnitude_ratio": float(magnitude_ratio),
        "n_negative": len(neg_closures),
        "n_positive": len(pos_closures),
        "closures": closures.tolist()
    }


def analyze_quadruplet_closures(tau_values: np.ndarray) -> Dict[str, Any]:
    """Analyze quadruplet closures."""
    n_arclets = len(tau_values)
    
    if n_arclets < 4:
        return {
            "order": 4,
            "n_closures": 0,
            "error": "Insufficient arclets for quadruplet closures"
        }
    
    quadruplets = list(combinations(range(n_arclets), 4))
    
    closures = []
    for quadruplet in quadruplets:
        closure = compute_quadruplet_closure(tau_values, quadruplet)
        closures.append(closure)
    
    closures = np.array(closures)
    abs_closures = np.abs(closures)
    
    # Statistics
    mean_Q = np.mean(abs_closures)
    # Handle single-closure case (std with ddof=1 returns NaN)
    if len(closures) > 1:
        sem_Q = np.std(abs_closures, ddof=1) / np.sqrt(len(closures))
        t_stat = mean_Q / sem_Q if sem_Q > 0 else 0.0
    else:
        sem_Q = None
        t_stat = 0.0
    
    # Bipolar structure
    neg_closures = closures[closures < 0]
    pos_closures = closures[closures > 0]
    
    if len(neg_closures) > 0 and len(pos_closures) > 0:
        neg_mean = np.mean(np.abs(neg_closures))
        pos_mean = np.mean(np.abs(pos_closures))
        magnitude_ratio = neg_mean / pos_mean
    else:
        neg_mean = pos_mean = 0
        magnitude_ratio = 0
    
    return {
        "order": 4,
        "n_closures": len(closures),
        "mean_Q_ns": float(mean_Q),
        "sem_ns": float(sem_Q) if sem_Q is not None else None,
        "t_statistic": float(t_stat),
        "detected": bool(t_stat > 5.0) if t_stat > 0 else False,
        "magnitude_ratio": float(magnitude_ratio),
        "n_negative": len(neg_closures),
        "n_positive": len(pos_closures),
        "closures": closures.tolist()
    }


def analyze_higher_order_scaling(triplet_results: Dict, quadruplet_results: Dict) -> Dict[str, Any]:
    """Analyze scaling of closure magnitude with order."""
    
    if triplet_results["n_closures"] == 0 or quadruplet_results["n_closures"] == 0:
        return {"error": "Insufficient data for scaling analysis"}
    
    triplet_H = triplet_results["mean_H_ns"]
    quadruplet_Q = quadruplet_results["mean_Q_ns"]
    
    # Expected scaling: higher-order closures should be larger
    # due to accumulation of TEP contributions
    scaling_ratio = quadruplet_Q / triplet_H if triplet_H > 0 else 0
    
    # Theoretical expectation: Q ≈ 2H for coherent TEP signal
    expected_ratio = 2.0
    agreement = bool(abs(scaling_ratio - expected_ratio) / expected_ratio < 0.5)
    
    return {
        "triplet_H_ns": triplet_H,
        "quadruplet_Q_ns": quadruplet_Q,
        "observed_scaling_ratio": float(scaling_ratio),
        "expected_scaling_ratio": expected_ratio,
        "agreement": agreement,
        "interpretation": f"Q/H = {scaling_ratio:.2f} (expected ~{expected_ratio:.1f})"
    }


def load_real_triplet_data() -> List[Dict[str, Any]]:
    """Load real triplet closure data."""
    closure_file = PROJECT_ROOT / "results" / "step_003_closure_final_per_epoch.json"
    
    if not closure_file.exists():
        raise FileNotFoundError(f"Closure delay results not found: {closure_file}")
    
    with open(closure_file, 'r') as f:
        data = json.load(f)
    
    # The file is a list of epoch dictionaries
    if isinstance(data, list):
        epochs = data
    elif isinstance(data, dict) and "epochs" in data:
        epochs = data["epochs"]
    else:
        raise ValueError(f"Unexpected format in {closure_file}")
    
    # Extract triplet closures
    results = []
    for epoch in epochs:
        triplets = epoch.get("triplets", [])
        if len(triplets) < 1:
            continue
        
        # Convert us to ns
        closures = [triplet.get("geometric_delta_us", triplet.get("delta_us", 0)) * 1000 for triplet in triplets]
        
        # Compute statistics
        abs_closures = np.abs(closures)
        mean_H = np.mean(abs_closures)
        # Handle single-closure case (std with ddof=1 returns NaN)
        if len(closures) > 1:
            sem_H = np.std(abs_closures, ddof=1) / np.sqrt(len(closures))
            t_stat = mean_H / sem_H if sem_H > 0 else 0.0
        else:
            sem_H = None
            t_stat = 0.0
        
        # Bipolar structure
        neg_closures = [c for c in closures if c < 0]
        pos_closures = [c for c in closures if c > 0]
        
        if len(neg_closures) > 0 and len(pos_closures) > 0:
            neg_mean = np.mean(np.abs(neg_closures))
            pos_mean = np.mean(np.abs(pos_closures))
            magnitude_ratio = neg_mean / pos_mean
        else:
            neg_mean = pos_mean = 0
            magnitude_ratio = 0
        
        results.append({
            "order": 3,
            "n_closures": len(closures),
            "mean_H_ns": float(mean_H),
            "sem_ns": float(sem_H) if sem_H is not None else None,
            "t_statistic": float(t_stat),
            "detected": bool(t_stat > 5.0) if t_stat > 0 else False,
            "magnitude_ratio": float(magnitude_ratio),
            "n_negative": len(neg_closures),
            "n_positive": len(pos_closures),
            "closures": closures
        })
    
    print(f"Loaded {len(results)} epochs with real triplet closure data")
    return results


def run_higher_order_analysis() -> Dict[str, Any]:
    """
    Run higher-order closure analysis using real triplet data.
    
    Note: Quadruplet closures require raw arclet data (individual tau measurements
    for each arclet), which is not available in the current dataset. The closure
    delay results only contain computed triplet closures, not the underlying
    arclet measurements needed to construct quadruplets.
    """
    
    # Load real triplet data
    triplet_results_all = load_real_triplet_data()
    
    # Aggregate statistics for triplets
    triplet_H_values = [r["mean_H_ns"] for r in triplet_results_all if r["n_closures"] > 0]
    triplet_mean = np.mean(triplet_H_values) if triplet_H_values else 0
    
    # Quadruplet analysis not possible without raw arclet data
    quadruplet_results_all = []
    for epoch_idx in range(len(triplet_results_all)):
        quadruplet_results_all.append({
            "order": 4,
            "n_closures": 0,
            "error": "Raw arclet data not available for quadruplet computation"
        })
    
    return {
        "triplet_results": triplet_results_all,
        "quadruplet_results": quadruplet_results_all,
        "aggregate": {
            "n_epochs": len(triplet_results_all),
            "mean_triplet_H_ns": float(triplet_mean),
            "mean_quadruplet_Q_ns": None,
            "scaling_ratio": None,
            "note": "Quadruplet closures require raw arclet data (individual tau measurements) which is not available in the current dataset"
        }
    }


def plot_closure_order_comparison(results: Dict[str, Any]) -> None:
    """Generate comparison plot of closure orders."""
    
    triplet_H_values = [r["mean_H_ns"] for r in results["triplet_results"] if r["n_closures"] > 0]
    
    # Since quadruplets are not available, only plot triplet distribution
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    
    # Triplet distribution
    ax.hist(triplet_H_values, bins=20, alpha=0.7, color='blue')
    ax.set_xlabel('Triplet |H| (ns)', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title('Triplet Closure Distribution (Real Data)', fontsize=14)
    ax.axvline(np.mean(triplet_H_values), color='red', linestyle='--', label=f'Mean: {np.mean(triplet_H_values):.2f} ns')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    plot_file = RESULTS_DIR / "step_039_closure_order_comparison.png"
    plt.savefig(plot_file, dpi=150)
    
    # Also save to site/public/figures/ for web display
    figures_dir = Path(__file__).resolve().parent.parent.parent / "site" / "public" / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    fig_file = figures_dir / "step_039_closure_order_comparison.png"
    plt.savefig(fig_file, dpi=150)
    plt.close()
    
    print(f"Closure order comparison plot saved to: {plot_file}")
    print(f"Also saved to: {fig_file}")


def main():
    """Run higher-order closure analysis using real triplet data."""
    print("=" * 80)
    print("STEP 035: HIGHER-ORDER CLOSURE QUANTITIES")
    print("=" * 80)
    print()
    print("Purpose: Analyze triplet closures using real data")
    print()
    print("Note: Quadruplet closures require raw arclet data (individual tau")
    print("      measurements for each arclet), which is not available in the")
    print("      current dataset. Only triplet closures can be analyzed.")
    print()
    
    # Run higher-order analysis
    print("1. LOADING REAL TRIplet DATA:")
    results = run_higher_order_analysis()
    print(f"   Analyzed {results['aggregate']['n_epochs']} epochs")
    print()
    
    # Aggregate statistics
    print("2. AGGREGATE STATISTICS:")
    aggregate = results["aggregate"]
    print(f"   Mean |triplet H|: {aggregate['mean_triplet_H_ns']:.2f} ns")
    print(f"   Mean |quadruplet Q|: Not available (requires raw arclet data)")
    print(f"   {aggregate['note']}")
    print()
    
    # Bipolar structure analysis
    print("3. BIPOLAR STRUCTURE ANALYSIS:")
    triplet_ratios = [r["magnitude_ratio"] for r in results["triplet_results"] if r["n_closures"] > 0]
    mean_triplet_ratio = np.mean(triplet_ratios) if triplet_ratios else 0
    print(f"   Triplet magnitude ratio: {mean_triplet_ratio:.3f}")
    print(f"   (Ratio of negative to positive closure magnitudes)")
    print()
    
    # Generate plots
    print("4. GENERATING PLOTS:")
    plot_closure_order_comparison(results)
    print()
    
    # Compile results
    full_results = {
        "triplet_results": results["triplet_results"],
        "quadruplet_results": results["quadruplet_results"],
        "aggregate": aggregate,
        "scaling_analysis": {
            "note": "Scaling analysis requires both triplet and quadruplet closures",
            "scaling_ratio": None,
            "reason": "Quadruplet closures not computable without raw arclet data"
        },
        "conclusions": [
            f"Triplet closures: {aggregate['mean_triplet_H_ns']:.2f} ns (real data)",
            f"Quadruplet closures: Not available (requires raw arclet data)",
            f"Bipolar structure: Triplet magnitude ratio = {mean_triplet_ratio:.3f}",
            f"Data limitation: Higher-order closure analysis requires raw arclet measurements"
        ],
        "implications": {
            "triplet_analysis": "Triplet closures analyzed using real data",
            "quadruplet_limitation": "Quadruplet analysis requires raw arclet data (individual tau measurements)",
            "data_requirement": "To perform full higher-order closure analysis, the pipeline would need to store individual arclet tau values"
        }
    }
    
    # Save results
    output_file = RESULTS_DIR / "step_039_higher_order_closure_results.json"
    with open(output_file, 'w') as f:
        json.dump(full_results, f, indent=2, cls=NpEncoder)
    
    print("=" * 80)
    print("CONCLUSIONS:")
    print("=" * 80)
    for conclusion in full_results["conclusions"]:
        print(f"  • {conclusion}")
    print()
    print("=" * 80)
    print("IMPLICATIONS:")
    print("=" * 80)
    for key, value in full_results["implications"].items():
        print(f"  {key}: {value}")
    print()
    print("=" * 80)
    print(f"Results saved to: {output_file}")
    print("=" * 80)
    
    return full_results


if __name__ == "__main__":
    main()
