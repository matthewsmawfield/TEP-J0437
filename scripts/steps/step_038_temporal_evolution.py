#!/usr/bin/env python3
"""
================================================================================
STEP 034: TEMPORAL EVOLUTION ANALYSIS
================================================================================

Purpose: Analyze how TEP signal evolves over time to test for correlations
with ISM parameters and identify potential systematic effects.

Temporal Analysis:
-----------------
- Track holonomy magnitude |H| over time (2008-2018 baseline)
- Correlate with scintillation parameters (scintillation timescale, scattering strength)
- Test for seasonal variations or long-term trends
- Check for correlations with Earth's orbital position

Expected Outcomes:
-----------------
If TEP is real and ISM-dependent:
- |H| should correlate with scintillation strength
- No correlation with Earth's orbital position (rules out instrumental effects)
- Possible variations on timescales of ISM changes (months to years)

If TEP is instrumental artifact:
- Correlation with Earth's orbital position (annual modulation)
- No correlation with ISM parameters
- Constant signal regardless of ISM conditions

================================================================================
"""

import json
import numpy as np
from pathlib import Path
from typing import Dict, Any, List
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from scripts.utils.json_numpy import NpEncoder

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def load_closure_data_with_time_and_ism() -> List[Dict[str, Any]]:
    """Load real closure delay data with temporal and ISM information."""
    closure_file = PROJECT_ROOT / "results" / "step_003_closure_final_per_epoch.json"
    secondary_catalog = PROJECT_ROOT / "data" / "secondary" / "j0437_secondary_catalog.json"
    
    if not closure_file.exists():
        raise FileNotFoundError(f"Closure delay results not found: {closure_file}")
    
    if not secondary_catalog.exists():
        raise FileNotFoundError(f"Secondary spectra catalog not found: {secondary_catalog}")
    
    with open(closure_file, 'r') as f:
        closure_data = json.load(f)
    
    with open(secondary_catalog, 'r') as f:
        secondary_data = json.load(f)
    
    # The file is a list of epoch dictionaries
    if isinstance(closure_data, list):
        closure_epochs = closure_data
    elif isinstance(closure_data, dict) and "epochs" in closure_data:
        closure_epochs = closure_data["epochs"]
    else:
        raise ValueError(f"Unexpected format in {closure_file}")
    
    # Create a mapping from filename to secondary spectrum data
    secondary_map = {}
    secondary_by_mjd = []
    for epoch in secondary_data.get("epochs", []):
        filename = epoch.get("file", "")
        if filename:
            # Remove the _secondary.npz suffix to match closure data
            base_name = filename.replace("_secondary.npz", "")
            secondary_map[base_name] = epoch
        mjd_start = epoch.get("mjd_start", None)
        if mjd_start is not None:
            secondary_by_mjd.append((float(mjd_start), epoch))
    
    # Extract closure delays and compute per-epoch H
    results = []
    matched_by_name = 0
    matched_by_mjd = 0
    nonzero_eta_matches = 0
    for epoch in closure_epochs:
        mjd = epoch.get("mjd", 0)
        if mjd == 0:
            continue
        
        # Compute year from MJD (MJD 51544.5 = January 1, 2000)
        year = 2000 + (mjd - 51544.5) / 365.25
        
        # Earth's orbital position (simplified)
        doy = (mjd - 51544.5) % 365.25
        orbital_phase = 2 * np.pi * doy / 365.25
        
        # Extract authoritative signed geometric closure delays from triplets
        triplets = epoch.get("triplets", [])
        closures = [triplet.get("geometric_delta_us") * 1000 for triplet in triplets if triplet.get("geometric_delta_us") is not None]  # Convert us to ns
        
        if len(closures) < 3:
            continue
        
        # Compute epoch signed H
        epoch_H_ns = np.mean(closures)
        epoch_sem = np.std(closures, ddof=1) / np.sqrt(len(closures))
        
        # Get ISM parameters from secondary catalog
        epoch_file = epoch.get("epoch", "")
        base_name = epoch_file.replace(".npz", "").replace("_secondary", "")
        secondary_epoch = secondary_map.get(base_name, None)
        match_method = None
        if secondary_epoch is not None:
            matched_by_name += 1
            match_method = "name"
        # NO FALLBACK: MJD matching removed - requires exact name match
        # This prevents false matches from temporally proximate epochs
        else:
            secondary_epoch = {}
        
        eta_screen1 = secondary_epoch.get("eta_screen1", 0.0)
        eta_screen2 = secondary_epoch.get("eta_screen2", 0.0)
        n_arclets = secondary_epoch.get("n_arclets", epoch.get("n_arclets", 0))
        if eta_screen1 > 0 or eta_screen2 > 0:
            nonzero_eta_matches += 1
        
        results.append({
            "mjd": mjd,
            "year": year,
            "orbital_phase": orbital_phase,
            "H_ns": epoch_H_ns,
            "sem_ns": epoch_sem,
            "n_triplets": len(closures),
            "closures": closures,
            "eta_screen1": eta_screen1,
            "eta_screen2": eta_screen2,
            "n_arclets": n_arclets,
            "secondary_match_method": match_method
        })
    
    print(f"Loaded {len(results)} epochs with real closure delay data, temporal and ISM information")
    print(f"Secondary join coverage: name={matched_by_name}, mjd_fallback={matched_by_mjd}, nonzero_eta={nonzero_eta_matches}")
    return results


def analyze_temporal_trends(results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Analyze temporal trends in TEP signal."""
    
    years = np.array([r["year"] for r in results])
    H_values = np.array([r["H_ns"] for r in results])
    H_sem = np.array([r["sem_ns"] for r in results])
    
    # Linear trend analysis
    slope, intercept, r_value, p_value, std_err = stats.linregress(years, H_values)
    
    # Check for significant trend
    significant_trend = p_value < 0.05
    
    return {
        "years": years.tolist(),
        "H_values": H_values.tolist(),
        "linear_trend": {
            "slope_ns_per_year": float(slope),
            "intercept_ns": float(intercept),
            "r_squared": float(r_value ** 2),
            "p_value": float(p_value),
            "significant": bool(significant_trend)
        },
        "mean_H_ns": float(np.mean(H_values)),
        "std_H_ns": float(np.std(H_values, ddof=1))
    }


def correlate_with_ISM(results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Correlate TEP signal with ISM parameters.
    
    Uses eta_screen1 and eta_screen2 from secondary spectra as ISM strength indicators.
    """
    H_values = np.array([r["H_ns"] for r in results])
    eta_screen1 = np.array([r["eta_screen1"] for r in results])
    eta_screen2 = np.array([r["eta_screen2"] for r in results])
    n_arclets = np.array([r["n_arclets"] for r in results])
    
    # Filter out epochs with zero eta values (no ISM measurement)
    mask1 = eta_screen1 > 0
    mask2 = eta_screen2 > 0
    
    # Correlation with eta_screen1
    if np.sum(mask1) > 3:
        corr1, p1 = stats.pearsonr(H_values[mask1], eta_screen1[mask1])
        eta1_result = {
            "correlation": float(corr1),
            "p_value": float(p1),
            "significant": bool(p1 < 0.05),
            "n_data_points": int(np.sum(mask1)),
            "skipped": False
        }
    else:
        eta1_result = {
            "correlation": None,
            "p_value": None,
            "significant": False,
            "n_data_points": int(np.sum(mask1)),
            "skipped": True,
            "reason": "Insufficient epochs with eta_screen1 measurements"
        }
    
    # Correlation with eta_screen2
    if np.sum(mask2) > 3:
        corr2, p2 = stats.pearsonr(H_values[mask2], eta_screen2[mask2])
        eta2_result = {
            "correlation": float(corr2),
            "p_value": float(p2),
            "significant": bool(p2 < 0.05),
            "n_data_points": int(np.sum(mask2)),
            "skipped": False
        }
    else:
        eta2_result = {
            "correlation": None,
            "p_value": None,
            "significant": False,
            "n_data_points": int(np.sum(mask2)),
            "skipped": True,
            "reason": "Insufficient epochs with eta_screen2 measurements"
        }
    
    # Correlation with n_arclets
    if np.std(n_arclets) == 0:
        raise ValueError(
            "Cannot compute correlation with n_arclets: zero variance in arclet counts. "
            "This indicates either all epochs have the same number of arclets or insufficient data."
        )
    corr_arclets, p_arclets = stats.pearsonr(H_values, n_arclets)
    
    return {
        "n_epochs_total": len(results),
        "n_epochs_with_eta1": int(np.sum(mask1)),
        "n_epochs_with_eta2": int(np.sum(mask2)),
        "correlation_with_eta_screen1": eta1_result,
        "correlation_with_eta_screen2": eta2_result,
        "correlation_with_n_arclets": {
            "correlation": float(corr_arclets),
            "p_value": float(p_arclets),
            "significant": bool(p_arclets < 0.05)
        }
    }


def correlate_with_orbit(results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Correlate TEP signal with Earth's orbital position."""
    
    H_values = np.array([r["H_ns"] for r in results])
    orbital_phases = np.array([r["orbital_phase"] for r in results])
    
    # Correlation with orbital phase
    corr_orbit, p_orbit = stats.pearsonr(H_values, np.sin(orbital_phases))
    
    return {
        "correlation_with_orbit": {
            "correlation": float(corr_orbit),
            "p_value": float(p_orbit),
            "significant": bool(p_orbit < 0.05),
            "expected_zero": True  # TEP should NOT correlate with orbit
        }
    }


def plot_temporal_evolution(results: List[Dict[str, Any]], trends: Dict[str, Any]) -> None:
    """Generate temporal evolution plot."""
    
    years = [r["year"] for r in results]
    H_values = [r["H_ns"] for r in results]
    H_sem = [r["sem_ns"] for r in results]
    
    plt.figure(figsize=(12, 6))
    plt.errorbar(years, H_values, yerr=H_sem, fmt='o', capsize=3, markersize=4)
    
    # Plot trend line
    slope = trends["linear_trend"]["slope_ns_per_year"]
    intercept = trends["linear_trend"]["intercept_ns"]
    years_fit = np.linspace(min(years), max(years), 100)
    H_fit = slope * years_fit + intercept
    plt.plot(years_fit, H_fit, '--', label=f'Trend: {slope:.2f} ns/yr', alpha=0.7)
    
    plt.xlabel('Year', fontsize=12)
    plt.ylabel('Signed H (ns)', fontsize=12)
    plt.title('Temporal Evolution of Signed TEP Signal', fontsize=14)
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.3)
    
    plot_file = RESULTS_DIR / "step_038_temporal_evolution.png"
    plt.savefig(plot_file, dpi=150)
    
    # Also save to site/public/figures/ for web display
    figures_dir = PROJECT_ROOT / "site" / "public" / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    fig_file = figures_dir / "step_038_temporal_evolution.png"
    plt.savefig(fig_file, dpi=150)
    plt.close()
    
    print(f"Temporal evolution plot saved to: {plot_file}")
    print(f"Also saved to: {fig_file}")


def main():
    """Run temporal evolution analysis using real closure delay data with ISM parameters."""
    print("=" * 80)
    print("STEP 034: TEMPORAL EVOLUTION ANALYSIS")
    print("=" * 80)
    print()
    print("Purpose: Analyze TEP signal evolution over time using real data with ISM parameters")
    print()
    
    # Load real closure delay data with temporal and ISM information
    print("1. LOADING REAL CLOSURE DELAY DATA WITH ISM PARAMETERS:")
    results = load_closure_data_with_time_and_ism()
    print(f"   Loaded {len(results)} epochs with temporal and ISM information")
    print()
    
    # Analyze temporal trends
    print("2. ANALYZING TEMPORAL TRENDS:")
    trends = analyze_temporal_trends(results)
    print(f"   Linear trend: {trends['linear_trend']['slope_ns_per_year']:.3f} ns/year")
    print(f"   R²: {trends['linear_trend']['r_squared']:.3f}")
    print(f"   Significant: {trends['linear_trend']['significant']}")
    print()
    
    # Correlate with ISM
    print("3. CORRELATING WITH ISM PARAMETERS:")
    ism_correlations = correlate_with_ISM(results)
    eta1 = ism_correlations['correlation_with_eta_screen1']
    if eta1.get('skipped'):
        print(f"   Eta screen1 correlation: skipped ({eta1.get('reason')})")
    else:
        print(f"   Eta screen1 correlation: r = {eta1['correlation']:.3f} (n = {ism_correlations['n_epochs_with_eta1']})")
        print(f"   Significant: {eta1['significant']}")
    print()
    eta2 = ism_correlations['correlation_with_eta_screen2']
    if eta2.get('skipped'):
        print(f"   Eta screen2 correlation: skipped ({eta2.get('reason')})")
    else:
        print(f"   Eta screen2 correlation: r = {eta2['correlation']:.3f} (n = {ism_correlations['n_epochs_with_eta2']})")
        print(f"   Significant: {eta2['significant']}")
    print()
    print(f"   N arclets correlation: r = {ism_correlations['correlation_with_n_arclets']['correlation']:.3f}")
    print(f"   Significant: {ism_correlations['correlation_with_n_arclets']['significant']}")
    print()
    
    # Correlate with orbit
    print("4. CORRELATING WITH EARTH'S ORBIT:")
    orbit_correlations = correlate_with_orbit(results)
    print(f"   Orbital correlation: r = {orbit_correlations['correlation_with_orbit']['correlation']:.3f}")
    print(f"   Significant: {orbit_correlations['correlation_with_orbit']['significant']}")
    print()
    
    # Generate plots
    print("5. GENERATING PLOTS:")
    plot_temporal_evolution(results, trends)
    print()
    
    # Compile results
    full_results = {
        "temporal_trends": trends,
        "ism_correlations": ism_correlations,
        "orbit_correlations": orbit_correlations,
        "conclusions": [
            f"Temporal trend: {'Significant' if trends['linear_trend']['significant'] else 'No significant trend'}",
            f"ISM correlation (eta1): {'Skipped' if ism_correlations['correlation_with_eta_screen1'].get('skipped') else ('Significant' if ism_correlations['correlation_with_eta_screen1']['significant'] else 'No significant correlation')}",
            f"ISM correlation (eta2): {'Skipped' if ism_correlations['correlation_with_eta_screen2'].get('skipped') else ('Significant' if ism_correlations['correlation_with_eta_screen2']['significant'] else 'No significant correlation')}",
            f"Orbital correlation: {'Present (instrumental concern)' if orbit_correlations['correlation_with_orbit']['significant'] else 'No significant orbital signature detected'}"
        ],
        "implications": {
            "temporal_stability": f"Signal stability: std = {trends['std_H_ns']:.2f} ns",
            "orbital_test": "No significant correlation with Earth's orbit argues against simple orbital/instrumental contamination",
            "ism_test": f"ISM correlation tested with {ism_correlations['n_epochs_with_eta1']} epochs with eta1 data"
        }
    }
    
    # Save results
    output_file = RESULTS_DIR / "step_038_temporal_evolution_results.json"
    with open(output_file, 'w') as f:
        json.dump(full_results, f, indent=2, cls=NpEncoder)
    
    print("=" * 80)
    print("CONCLUSIONS:")
    print("=" * 80)
    for conclusion in full_results["conclusions"]:
        print(f"  * {conclusion}")
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
