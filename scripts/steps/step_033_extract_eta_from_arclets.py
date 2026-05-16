#!/usr/bin/env python3
"""
Extract arc curvature eta from arclets in secondary spectra.

The arc curvature eta is related to the screen distance s through the parabolic
relationship in the secondary spectrum: tau = eta * fD^2

This script extracts eta by fitting this parabolic model to the detected arclets
in each secondary spectrum file, then updates the secondary catalog with the
extracted values.
"""

import json
import numpy as np
from pathlib import Path
from scripts.utils.json_numpy import NpEncoder

PROJECT_ROOT = Path(__file__).parent.parent.parent
SEC_DIR = PROJECT_ROOT / "data" / "secondary" / "jiamusi"
CATALOG_PATH = PROJECT_ROOT / "data" / "secondary" / "jiamusi_secondary_catalog.json"
UPDATED_CATALOG_PATH = PROJECT_ROOT / "data" / "secondary" / "jiamusi_secondary_catalog_updated.json"


def extract_eta_from_arclets(sec_file_path):
    """
    Extract arc curvature eta by fitting parabolic model to arclets.
    
    Args:
        sec_file_path: Path to secondary spectrum .npz file
    
    Returns:
        eta1: Estimated eta from arclets (or 0.0 if extraction fails)
        eta2: Secondary eta estimate (or 0.0)
        n_arclets: Number of arclets detected
    """
    try:
        data = np.load(sec_file_path)
        
        # Check if arclets are available
        if 'arclets' not in data:
            return 0.0, 0.0, 0
        
        arclets = data['arclets']
        n_arclets = len(arclets)
        
        if n_arclets < 2:
            return 0.0, 0.0, n_arclets
        
        # Extract tau and fD from arclets
        # Arclet format: [tau_us, fD_mHz, amplitude]
        tau = arclets[:, 0]  # tau in microseconds
        fD = arclets[:, 1]   # fD in mHz
        
        # Fit parabolic model: tau = eta * fD^2
        # Linear regression: tau = eta * fD_sq
        fD_sq = fD**2
        
        # Weighted fit using amplitude as weight (brighter arclets more reliable)
        if arclets.shape[1] >= 3:
            weights = arclets[:, 2]  # amplitude
            weights = weights / np.mean(weights)  # normalize
        else:
            weights = np.ones_like(tau)
        
        # Weighted least squares fit
        eta_fit = np.sum(weights * tau * fD_sq) / np.sum(weights * fD_sq**2)
        
        # Calculate residuals
        tau_pred = eta_fit * fD_sq
        residuals = tau - tau_pred
        rms_residual = np.sqrt(np.mean(weights * residuals**2))
        
        # Only accept fit if residuals are reasonable
        # RMS residual should be less than the mean tau
        mean_tau = np.mean(tau)
        if rms_residual > mean_tau:
            return 0.0, 0.0, n_arclets
        
        # Secondary estimate: use median of eta calculated from each arclet
        eta_individual = tau / fD_sq
        eta_median = np.median(eta_individual)
        
        return float(eta_fit), float(eta_median), n_arclets
        
    except Exception as e:
        print(f"Error processing {sec_file_path}: {e}")
        return 0.0, 0.0, 0


def main():
    """Extract eta from all Jiamusi secondary spectra and update catalog."""
    
    print("=" * 70)
    print("ARC CURVATURE eta EXTRACTION FROM ARCLETS")
    print("=" * 70)
    
    # Load existing catalog
    with open(CATALOG_PATH, 'r') as f:
        catalog = json.load(f)
    
    # Extract eta from each epoch
    updated_epochs = []
    
    for epoch in catalog["epochs"]:
        filename = epoch["file"]
        sec_file_path = SEC_DIR / filename
        
        if not sec_file_path.exists():
            print(f"  [FAIL] {filename}: File not found")
            updated_epochs.append(epoch)
            continue
        
        eta1, eta2, n_arclets = extract_eta_from_arclets(sec_file_path)
        
        # Update epoch with extracted values
        updated_epoch = epoch.copy()
        if eta1 > 0:
            updated_epoch["eta_screen1"] = eta1
            updated_epoch["eta_screen2"] = eta2
            print(f"  [OK] {filename}: eta1={eta1:.6f}, eta2={eta2:.6f}, n_arclets={n_arclets}")
        else:
            print(f"  [SKIP] {filename}: no arc curvature extracted (eta=0)")
        
        updated_epochs.append(updated_epoch)
    
    # Create updated catalog
    updated_catalog = catalog.copy()
    updated_catalog["epochs"] = updated_epochs
    
    # Save updated catalog
    with open(UPDATED_CATALOG_PATH, 'w') as f:
        json.dump(updated_catalog, f, indent=2, cls=NpEncoder)
    
    print(f"\nUpdated catalog saved to {UPDATED_CATALOG_PATH}")
    
    # Summary statistics
    eta1_values = [e["eta_screen1"] for e in updated_epochs if e["eta_screen1"] > 0]
    if eta1_values:
        print(f"\nSummary:")
        print(f"  Epochs with eta: {len(eta1_values)}/{len(updated_epochs)}")
        print(f"  Mean eta1: {np.mean(eta1_values):.6f}")
        print(f"  Std eta1: {np.std(eta1_values):.6f}")
        print(f"  Range: [{np.min(eta1_values):.6f}, {np.max(eta1_values):.6f}]")
    
    return updated_catalog


if __name__ == "__main__":
    main()
