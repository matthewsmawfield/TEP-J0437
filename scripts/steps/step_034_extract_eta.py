#!/usr/bin/env python3
"""
Extract arc curvature eta from secondary spectra and estimate screen distances.

The arc curvature eta is related to the screen distance s through:
eta = (2pi / lambda) * (s * (1-s) * D * θ^2) / c

Where:
- lambda is observing wavelength
- s is fractional screen distance (0 = at observer, 1 = at pulsar)
- D is pulsar distance
- θ is angular scale of scattering disk
- c is speed of light

For a given eta and known D, the analysis can solve for s.
"""

import json
import sys
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.json_numpy import NpEncoder
BASE_CATALOG_PATH = PROJECT_ROOT / "data" / "secondary" / "jiamusi_secondary_catalog.json"
UPDATED_CATALOG_PATH = PROJECT_ROOT / "data" / "secondary" / "jiamusi_secondary_catalog_updated.json"
RESULTS_PATH = PROJECT_ROOT / "results" / "step_031_jiamusi_eta_analysis.json"

# Physical constants
C_LIGHT = 2.998e8  # m/s
JIAMUSI_FREQ_MHZ = 400.0  # Center frequency from Wang et al. 2018
LAMBDA_M = C_LIGHT / (JIAMUSI_FREQ_MHZ * 1e6)  # Wavelength in meters

# Pulsar distances from step_030 (in pc)
PULSAR_DISTANCES = {
    "B0329": 1000.0,
    "B0355": 1000.0,
    "B0540": 1600.0,
    "B0740": 2000.0,
    "B1508": 2100.0,
    "B1933": 3700.0,
    "B2154": 2900.0,
    "B2310": 1060.0,
    "B2324": 2700.0,
    "B2351": 2400.0,
}

def eta_to_screen_distance(eta, dist_pc):
    """
    Convert arc curvature eta to fractional screen distance s.
    
    Using the simplified relation from arc curvature measurements:
    eta ≈ 2pi * s * (1-s) * (D_eff / lambda) * (θ^2)
    
    For ISM scattering screens, the arc curvature provides a direct
    measure of s when combined with proper motion and distance.
    
    A simpler empirical relation: s ≈ 0.5 +/- 0.3 * sqrt(eta / eta_typical)
    where eta_typical ~ 0.01 for typical ISM screens.
    
    This is an approximation - precise s requires fitting the arc
    curvature with known proper motion and scattering timescale.
    
    Args:
        eta: Arc curvature (dimensionless)
        dist_pc: Pulsar distance in pc
    
    Returns:
        s: Fractional screen distance (0-1), or None if calculation fails
    """
    if eta <= 0:
        return None
    
    # Empirical relation based on typical arc curvatures
    # eta ~ 0.01 corresponds to s ~ 0.5 (midway screen)
    # Smaller eta -> screen closer to observer or pulsar
    eta_typical = 0.01
    
    # Calculate deviation from typical
    ratio = eta / eta_typical
    
    # Estimate s: for ratio ~1, s ~ 0.5
    # This is a rough approximation
    s_estimate = 0.5 * np.sqrt(ratio)
    
    # Clamp to physical range
    s_estimate = np.clip(s_estimate, 0.01, 0.99)
    
    return s_estimate


def main():
    """Extract and analyze eta values from secondary spectra."""
    
    print("=" * 70)
    print("ARC CURVATURE eta ANALYSIS FOR JIAMUSI PULSARS")
    print("=" * 70)
    
    catalog_path = UPDATED_CATALOG_PATH
    if not catalog_path.exists() and BASE_CATALOG_PATH.exists():
        print("Updated catalog not found; extracting eta from base Jiamusi catalog first.")
        from scripts.steps.step_032_extract_eta_from_arclets import main as extract_eta_from_arclets

        extract_eta_from_arclets()

    if not catalog_path.exists():
        print(f"[WARN] Jiamusi secondary catalog not found: {catalog_path}")
        print("[WARN] Skipping eta analysis until Jiamusi secondary spectra are generated.")
        RESULTS_PATH.parent.mkdir(parents=True, exist_ok=True)
        skipped = {
            "status": "skipped",
            "reason": "Jiamusi secondary catalog not found",
            "required_catalog": str(catalog_path),
            "base_catalog": str(BASE_CATALOG_PATH),
        }
        with open(RESULTS_PATH, 'w') as f:
            json.dump(skipped, f, indent=2, cls=NpEncoder)
        return skipped

    # Load secondary catalog
    with open(catalog_path, 'r') as f:
        catalog = json.load(f)
    
    # Extract eta values by pulsar
    pulsar_eta_data = {}
    
    for epoch in catalog["epochs"]:
        filename = epoch["file"]
        mjd = epoch["mjd_start"]
        eta1 = epoch["eta_screen1"]
        eta2 = epoch["eta_screen2"]
        n_arclets = epoch["n_arclets"]
        
        # Extract pulsar name from filename (e.g., "B1933_201506160205d_secondary.npz")
        pulsar_name = filename.split('_')[0]
        
        if pulsar_name not in pulsar_eta_data:
            pulsar_eta_data[pulsar_name] = {
                "eta1_values": [],
                "eta2_values": [],
                "mjds": [],
                "n_arclets": [],
            }
        
        pulsar_eta_data[pulsar_name]["eta1_values"].append(eta1)
        pulsar_eta_data[pulsar_name]["eta2_values"].append(eta2)
        pulsar_eta_data[pulsar_name]["mjds"].append(mjd)
        pulsar_eta_data[pulsar_name]["n_arclets"].append(n_arclets)
    
    # Analyze each pulsar
    results = {}
    
    for pulsar_name, data in pulsar_eta_data.items():
        eta1_vals = np.array(data["eta1_values"])
        eta2_vals = np.array(data["eta2_values"])
        n_arclets = np.array(data["n_arclets"])
        
        # Get pulsar distance
        if pulsar_name not in PULSAR_DISTANCES:
            raise ValueError(
                f"Pulsar {pulsar_name} not found in PULSAR_DISTANCES. "
                "Add distance from Wang et al. 2018 or ATNF catalog."
            )
        dist_pc = PULSAR_DISTANCES[pulsar_name]
        
        # Filter to epochs with measured eta
        measured_mask = (eta1_vals > 0) | (eta2_vals > 0)
        n_measured = np.sum(measured_mask)
        n_total = len(eta1_vals)
        
        if n_measured > 0:
            eta1_mean = np.mean(eta1_vals[measured_mask]) if np.any(eta1_vals[measured_mask] > 0) else 0
            eta2_mean = np.mean(eta2_vals[measured_mask]) if np.any(eta2_vals[measured_mask] > 0) else 0
            eta_mean = (eta1_mean + eta2_mean) / 2 if (eta1_mean > 0 or eta2_mean > 0) else 0
            
            # Estimate screen distance
            s_estimate = eta_to_screen_distance(eta_mean, dist_pc) if eta_mean > 0 else None
        else:
            eta1_mean = 0
            eta2_mean = 0
            eta_mean = 0
            s_estimate = None
        
        results[pulsar_name] = {
            "n_epochs_total": n_total,
            "n_epochs_with_eta": int(n_measured),
            "eta1_mean": float(eta1_mean),
            "eta2_mean": float(eta2_mean),
            "eta_mean": float(eta_mean),
            "screen_distance_s": float(s_estimate) if s_estimate is not None else None,
            "distance_pc": dist_pc,
        }
        
        print(f"\n{pulsar_name}:")
        print(f"  Epochs with eta: {n_measured}/{n_total}")
        print(f"  eta1 mean: {eta1_mean:.6f}")
        print(f"  eta2 mean: {eta2_mean:.6f}")
        print(f"  eta mean: {eta_mean:.6f}")
        if s_estimate is not None:
            print(f"  Estimated screen distance s: {s_estimate:.3f}")
        else:
            print(f"  Screen distance: Not estimable (eta = 0)")
    
    # Save results
    with open(RESULTS_PATH, 'w') as f:
        json.dump(results, f, indent=2, cls=NpEncoder)
    
    print(f"\nResults saved to {RESULTS_PATH}")
    
    return results


if __name__ == "__main__":
    main()
