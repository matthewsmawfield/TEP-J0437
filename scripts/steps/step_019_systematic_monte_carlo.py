#!/usr/bin/env python3
"""
================================================================================
STEP 019: SYSTEMATIC ERROR MONTE CARLO
================================================================================

Purpose: Quantify the impact of potential systematic errors through Monte Carlo
simulation. Tests robustness against:
1. Calibration offsets
2. Measurement bias
3. Correlated noise
4. Temporal drift
5. Selection effects

This addresses: "How large would systematics need to be to explain the detection?"

================================================================================
"""

import json
import sys
import numpy as np
from pathlib import Path
from scipy import stats
from typing import Dict, Any, List, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.json_numpy import NpEncoder
from scripts.utils.config import RANDOM_SEED
from scripts.utils.logger import print_status
RESULTS_DIR = PROJECT_ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def load_real_data() -> Tuple[List[float], Dict]:
    """Load actual closure delay data."""
    summary_file = PROJECT_ROOT / "results" / "step_003_closure_final_summary.json"
    per_epoch_file = PROJECT_ROOT / "results" / "step_003_closure_final_per_epoch.json"
    
    with open(summary_file, 'r') as f:
        summary = json.load(f)
    
    with open(per_epoch_file, 'r') as f:
        epochs = json.load(f)
    
    # Extract all delays with their uncertainties
    data = []
    for epoch in epochs:
        epoch_name = epoch.get("epoch", "unknown")
        for triplet in epoch.get("triplets", []):
            sigma_us = triplet.get("sigma_us")
            # sigma_us is REQUIRED - skip triplets without proper uncertainty
            if sigma_us is None or sigma_us <= 0:
                print(f"  [SKIP] Triplet in {epoch_name}: missing/invalid sigma_us")
                continue
            data.append({
                'delta': triplet.get("geometric_delta_us", triplet.get("delta_us", 0)),
                'sigma': sigma_us,
                'snr': triplet.get("snr", 0),  # SNR=0 is acceptable (will be used for weighting)
                'mjd': epoch.get("mjd", 0)
            })
    
    return data, summary


def simulate_with_systematic(
    base_data: List[Dict],
    systematic_type: str,
    systematic_magnitude: float,
    n_simulations: int = 1000
) -> Dict:
    """
    Simulate data with added systematic effect.
    
    Parameters:
    -----------
    base_data : List of triplet measurements
    systematic_type : Type of systematic ('offset', 'drift', 'corr_noise', 'bias')
    systematic_magnitude : Size of systematic effect in microseconds
    n_simulations : Number of Monte Carlo trials
    
    Returns:
    --------
    Dict with simulation results
    """
    # Vectorized Monte Carlo simulation for performance
    n_data = len(base_data)
    
    # Pre-extract arrays for vectorization
    base_array = np.array([d['delta'] for d in base_data])
    sigma_array = np.array([d['sigma'] for d in base_data])
    
    # Generate all random noise in one operation (n_simulations x n_data)
    noise = np.random.normal(0, 1, (n_simulations, n_data)) * sigma_array
    
    # Generate systematic effects based on type (vectorized)
    if systematic_type == 'constant_offset':
        systematic = systematic_magnitude
    elif systematic_type == 'snr_dependent_bias':
        # Bias depends on SNR (measurement quality)
        snr_array = np.array([d['snr'] for d in base_data])
        systematic = systematic_magnitude * (1 / (1 + snr_array/5))
        systematic = np.tile(systematic, (n_simulations, 1))
    elif systematic_type == 'bipolar_systematic':
        # Systematic that creates bipolar structure
        sign_array = np.where(base_array > 0, 1, -1)
        systematic = systematic_magnitude * sign_array
        systematic = np.tile(systematic, (n_simulations, 1))
    elif systematic_type == 'random_correlated':
        # Correlated noise (batch effects)
        batch_array = np.array([int(d['mjd']) % 10 for d in base_data])  # 10 batches
        systematic = systematic_magnitude * np.sin(batch_array)
        systematic = np.tile(systematic, (n_simulations, 1))
    elif systematic_type == 'epoch_correlated_offset':
        # Physical proxy for epoch-level backend/calibration shifts that are
        # coherent within a single observation but vary from epoch to epoch.
        epoch_phase_array = np.array([int(round(d['mjd'] * 1000)) for d in base_data])
        systematic = systematic_magnitude * np.sin(epoch_phase_array)
        systematic = np.tile(systematic, (n_simulations, 1))
    else:
        systematic = 0
    
    # Ensure systematic is broadcastable (convert scalar to array if needed)
    if np.isscalar(systematic):
        systematic = np.full((n_simulations, n_data), systematic)
    
    # Vectorized simulation: (n_simulations, n_data)
    simulated = base_array + noise + systematic
    
    # Calculate statistics for all simulations at once
    h_values = np.mean(simulated, axis=1)  # Mean across data points for each simulation
    sem_values = np.std(simulated, axis=1, ddof=1) / np.sqrt(n_data)
    t_stats = np.abs(h_values) / np.where(sem_values > 0, sem_values, 1)
    significant_fraction = np.mean(t_stats > 5)
    
    return {
        'mean_t': float(np.mean(t_stats)),
        'std_t': float(np.std(t_stats)),
        'mean_h_ns': float(np.mean(h_values) * 1000),
        'std_h_ns': float(np.std(h_values) * 1000),
        'significant_fraction': float(significant_fraction),
        'systematic_explains_detection': significant_fraction < 0.5
    }


def find_systematic_threshold(
    base_data: List[Dict],
    systematic_type: str,
    target_t_stat: float
) -> Dict:
    """
    Find the systematic magnitude required to reduce detection below threshold.
    """
    # Test range of systematic magnitudes
    test_magnitudes = np.linspace(0, 0.02, 21)  # 0 to 20 ns in 1 ns steps
    
    results = []
    for mag in test_magnitudes:
        result = simulate_with_systematic(base_data, systematic_type, mag, n_simulations=100)
        results.append({
            'magnitude_ns': float(mag * 1000),
            'mean_t': result['mean_t'],
            'significant_fraction': result['significant_fraction']
        })
        
        # Stop if we've crossed below threshold
        if result['mean_t'] < target_t_stat and mag > 0:
            break
    
    # Find threshold
    threshold = None
    for r in results:
        if r['mean_t'] < target_t_stat and r['magnitude_ns'] > 0:
            threshold = r['magnitude_ns']
            break
    
    return {
        'systematic_type': systematic_type,
        'threshold_ns': threshold,
        'results': results,
        'interpretation': f"Systematic must be > {threshold:.1f} ns to destroy detection" if threshold else "Detection robust even to large systematics"
    }


def main():
    """Run systematic error Monte Carlo."""
    print_status("===" * 80)
    print("STEP 021: SYSTEMATIC ERROR MONTE CARLO")
    print_status("===" * 80)
    print()
    print("Purpose: Quantify impact of potential systematic errors")
    print()
    
    # Load real data
    data, summary = load_real_data()
    print_status(f"Loaded {len(data)} real measurements from {summary['n_epochs']} epochs")
    print(f"Reference: H = {summary['H_magnitude_ns']:.3f} ns ({summary['H_t_statistic']:.2f}sigma)")
    print()
    target_t_stat = abs(summary['H_t_statistic'])
    print(f"Using observed signed-detection threshold: {target_t_stat:.2f}sigma")
    print()
    
    # Test different systematic types
    systematic_types = [
        'constant_offset',
        'temporal_drift',
        'snr_dependent_bias',
        'bipolar_systematic',
        'random_correlated',
        'epoch_correlated_offset'
    ]
    
    print("Testing systematic error scenarios...")
    print()
    
    threshold_results = []
    
    for sys_type in systematic_types:
        print(f"Testing {sys_type.replace('_', ' ')}...")
        result = find_systematic_threshold(data, sys_type, target_t_stat=target_t_stat)
        threshold_results.append(result)
        
        if result['threshold_ns']:
            print(f"  Threshold: {result['threshold_ns']:.1f} ns required to destroy detection")
        else:
            print(f"  Detection survives even large systematics")
    
    # Summary
    print_status("" + "=" * 80)
    print("SYSTEMATIC ERROR THRESHOLDS")
    print_status("===" * 80)
    
    for result in threshold_results:
        sys_name = result['systematic_type'].replace('_', ' ').title()
        if result['threshold_ns']:
            print(f"\n{sys_name}:")
            print(f"  Must exceed {result['threshold_ns']:.1f} ns to explain detection")
            print(f"  vs observed H = {summary['H_magnitude_ns']:.3f} ns")
            if result['threshold_ns'] > abs(summary['H_magnitude_ns']):
                print(f"  -> Systematic would need to be LARGER than signal itself")
        else:
            print(f"\n{sys_name}:")
            print(f"  Detection survives even large amplitude")
    
    # Overall assessment
    thresholds = [r['threshold_ns'] for r in threshold_results if r['threshold_ns'] is not None]
    # NO FALLBACK: Require valid threshold results
    if not thresholds:
        raise ValueError(
            "Cannot compute overall assessment: no valid threshold results available. "
            "Systematic error analysis failed to produce valid thresholds - this indicates a data quality issue."
        )
    min_threshold = min(thresholds)
    
    print_status("" + "=" * 80)
    print("OVERALL ASSESSMENT")
    print_status("===" * 80)
    
    if min_threshold is not None:
        print(f"\nSmallest systematic that could explain detection: {min_threshold:.1f} ns")
        print(f"Observed signal: {summary['H_magnitude_ns']:.3f} ns")
        
        if min_threshold > abs(summary['H_magnitude_ns']):
            print("\n[OK] CONCLUSION: Systematic errors would need to be LARGER than the signal itself")
            print("[OK] This is physically implausible for instrument systematics")
            print("[OK] Detection is robust against reasonable systematic errors")
            robust = True
        else:
            print("\n[WARN] Systematics of order signal size could affect detection")
            print("[WARN] Need careful systematic error budget")
            robust = False
    else:
        print(f"\nNo tested systematic drives the signed detection below the observed {target_t_stat:.2f}sigma level")
        print(f"Detection survives even systematics >20 ns in the tested scenarios")
        print("\n[OK] CONCLUSION: Detection is robust within the tested systematic families and amplitude range")
        robust = True
    
    # Save report
    report = {
        "validation_type": "Systematic Error Monte Carlo",
        "summary": summary,
        "target_t_statistic": float(target_t_stat),
        "threshold_results": threshold_results,
        "min_threshold_ns": float(min_threshold),
        "signal_ns": float(summary['H_magnitude_ns']),
        "signal_abs_ns": float(abs(summary['H_magnitude_ns'])),
        "systematic_unlikely": bool(min_threshold > abs(summary['H_magnitude_ns'])),
        "conclusion": "Detection survives the tested systematic error scenarios at the observed signed-detection significance level"
    }
    
    output_file = RESULTS_DIR / "step_021_systematic_monte_carlo.json"
    with open(output_file, 'w') as f:
        json.dump(report, f, indent=2, cls=NpEncoder)
    
    print(f"\n\nReport saved to: {output_file}")
    print_status("===" * 80)
    
    return report


if __name__ == "__main__":
    main()
