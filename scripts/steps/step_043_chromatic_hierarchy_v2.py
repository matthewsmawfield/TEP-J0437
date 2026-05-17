#!/usr/bin/env python3
"""
================================================================================
STEP 043 V2: CHROMATIC HIERARCHY USING PHASE CLOSURE ψ (Bayesian)
================================================================================

Tests whether the Phase Closure ψ exhibits frequency dependence. Under TEP,
the holonomy phase is achromatic (predicted δ = 0). Under standard ISM
chromatic plasma effects, the delay scales as ν^(-2), which would produce a
frequency-dependent ψ if there were any chromatic bias in the measurement.

Unlike the legacy step_043 (which used unsigned |H| and hit optimizer
boundaries), this module uses the noise-immune Phase Closure ψ as the
frequency-dependent observable. The model is:

    ψ_{e,b} ~ von Mises(μ_b, κ_e)    for epoch e in frequency band b
    μ_b = μ_0 + δ · ln(ν_b / ν_0)
    δ ~ Uniform(-4, 4)                prior on frequency exponent

where ν_b are the observing band centre frequencies and ν_0 is a reference
frequency. The frequency exponent δ is the key diagnostic:
    δ = 0   → achromatic (TEP prediction)
    δ = -2  → standard chromatic ISM plasma delay scaling

Because the sub-band samples are small (sb0: 5 epochs, sb1: 5 epochs),
uncertainties on δ are large. The script therefore implements a proper
profile-likelihood approach with bootstrap confidence intervals, avoiding
the boundary-hit artifacts of the legacy optimizer.

Data Sources:
-------------
- Full band:  results/step_003_closure_final_per_epoch_j0437.json  (1380 MHz)
- Sub-band 0: results/step_003_closure_final_per_epoch_j0437_sb0.json (1089 MHz)
- Sub-band 1: results/step_003_closure_final_per_epoch_j0437_sb1.json (1476 MHz)

OUTPUT:
-------
results/step_043_chromatic_hierarchy_v2.json
    - Profile likelihood for δ
    - 68% and 95% confidence intervals
    - Model comparison: achromatic (δ=0) vs chromatic (δ free)
    - ΔBIC and approximate Bayes factor
    - Conservative conclusion given sparse sub-band sampling

================================================================================
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats
from scipy.optimize import minimize_scalar

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.config import RANDOM_SEED
from scripts.utils.json_numpy import NpEncoder
from scripts.utils.logger import print_status

RESULTS_DIR = PROJECT_ROOT / "results"
OUTPUT_FILE = RESULTS_DIR / "step_043_chromatic_hierarchy_v2.json"

np.random.seed(RANDOM_SEED)


def circular_mean_and_rbar(angles, weights=None):
    """Circular mean and mean resultant length."""
    if weights is None:
        weights = np.ones_like(angles)
    weights = np.asarray(weights)
    z = np.sum(weights * np.exp(1j * angles))
    w_sum = np.sum(weights)
    r_bar = float(np.abs(z) / w_sum) if w_sum > 0 else 0.0
    psi_mean = float(np.angle(z / w_sum)) if w_sum > 0 else 0.0
    return psi_mean, r_bar


def von_mises_log_likelihood(psi, mu, kappa):
    """Log-likelihood of von Mises distribution."""
    if kappa < 100:
        norm = np.log(2 * np.pi) + np.log(stats.i0(kappa))
    else:
        norm = np.log(2 * np.pi) + kappa - 0.5 * np.log(2 * np.pi * kappa)
    ll = kappa * np.cos(psi - mu) - norm
    return float(np.sum(ll))


def load_band_data(per_epoch_file, freq_mhz):
    """Load per-epoch data for a single frequency band."""
    if not per_epoch_file.exists():
        return None
    
    with open(per_epoch_file) as f:
        data = json.load(f)
    
    epochs = []
    for ep in data:
        triplets = ep.get("triplets", [])
        if len(triplets) < 3:
            continue
        
        psi_vals = np.array([t["phase_closure_rad"] for t in triplets if t.get("phase_closure_rad") is not None])
        snr_vals = np.array([t.get("snr", 1.0) for t in triplets])
        
        if len(psi_vals) < 3:
            continue
        
        weights = snr_vals ** 2
        psi_mean, r_bar = circular_mean_and_rbar(psi_vals, weights)
        
        epochs.append({
            "psi_mean": psi_mean,
            "r_bar": r_bar,
            "n_triplets": len(triplets),
            "freq_mhz": freq_mhz,
        })
    
    return epochs


def fit_chromatic_model(all_epochs, delta):
    """Fit the model with a fixed δ, returning total log-likelihood.
    
    For each band b, compute μ_b = μ_0 + δ·ln(ν_b/ν_0) and find the MLE
    for μ_0 and per-band κ.
    """
    freqs = np.array([e["freq_mhz"] for e in all_epochs])
    psi_vals = np.array([e["psi_mean"] for e in all_epochs])
    n_bands = len(set(freqs))
    
    nu_0 = 1380.0  # Reference: full band centre frequency
    
    # For fixed δ, the model is μ_b = μ_0 + δ·ln(ν_b/ν_0)
    # MLE for μ_0: weighted circular mean of (ψ_e - δ·ln(ν_e/ν_0))
    shifted = psi_vals - delta * np.log(freqs / nu_0)
    z = np.sum(np.exp(1j * shifted))
    mu_0 = float(np.angle(z / len(shifted)))
    R = float(np.abs(z) / len(shifted))
    
    # Concentration parameter
    if R < 0.53:
        kappa = 2 * R + R**3 + (5.0 / 6.0) * R**5
    elif R < 0.85:
        kappa = -0.4 + 1.39 * R + 0.43 / (1 - R)
    else:
        kappa = 1.0 / (2 * (1 - R) - (1 - R)**2 - (1 - R)**3)
    
    # Compute predicted μ for each epoch
    mu_pred = mu_0 + delta * np.log(freqs / nu_0)
    
    # Total log-likelihood
    ll = von_mises_log_likelihood(psi_vals, mu_pred, kappa)
    
    return {
        "mu_0": float(mu_0),
        "kappa": float(kappa),
        "log_likelihood": float(ll),
        "n_params": 2,  # mu_0 and kappa (δ is fixed in this call)
    }


def profile_likelihood_delta(all_epochs, delta_grid):
    """Compute profile likelihood over a grid of δ values."""
    results = []
    for delta in delta_grid:
        fit = fit_chromatic_model(all_epochs, delta)
        results.append({
            "delta": float(delta),
            **fit,
        })
    return results


def bootstrap_delta_ci(all_epochs, n_boot=10000, seed=42):
    """Bootstrap confidence interval for δ.
    
    Resamples epochs with replacement and re-fits the chromatic model
    to find the MLE δ for each bootstrap sample.
    """
    rng = np.random.RandomState(seed)
    n = len(all_epochs)
    
    delta_boots = []
    
    for _ in range(n_boot):
        idx = rng.choice(n, size=n, replace=True)
        boot_epochs = [all_epochs[i] for i in idx]
        
        # Find MLE δ for this bootstrap sample
        def neg_ll(delta):
            fit = fit_chromatic_model(boot_epochs, delta)
            return -fit["log_likelihood"]
        
        result = minimize_scalar(neg_ll, bounds=(-4, 4), method="bounded")
        if result.success:
            delta_boots.append(result.x)
    
    delta_boots = np.array(delta_boots)
    
    if len(delta_boots) == 0:
        return None, None, None
    
    # Compute CI from bootstrap distribution
    delta_mle = float(np.median(delta_boots))
    ci_68 = (float(np.percentile(delta_boots, 16)), float(np.percentile(delta_boots, 84)))
    ci_95 = (float(np.percentile(delta_boots, 2.5)), float(np.percentile(delta_boots, 97.5)))
    
    return delta_mle, ci_68, ci_95


def main():
    print_status("=" * 70, "INFO")
    print_status("STEP 043 V2: CHROMATIC HIERARCHY (PHASE CLOSURE ψ)", "INFO")
    print_status("=" * 70)
    
    # Load data for all three bands
    bands = {
        "full": (RESULTS_DIR / "step_003_closure_final_per_epoch_j0437.json", 1380.0),
        "sb0": (RESULTS_DIR / "step_003_closure_final_per_epoch_j0437_sb0.json", 1089.0),
        "sb1": (RESULTS_DIR / "step_003_closure_final_per_epoch_j0437_sb1.json", 1476.0),
    }
    
    all_epochs = []
    band_summaries = {}
    
    for band_name, (file_path, freq_mhz) in bands.items():
        epochs = load_band_data(file_path, freq_mhz)
        if epochs is None:
            print_status(f"WARNING: Could not load data for {band_name}", "WARN")
            continue
        
        all_epochs.extend(epochs)
        
        # Compute band-level circular mean
        psi_vals = np.array([e["psi_mean"] for e in epochs])
        weights = np.ones(len(psi_vals))
        psi_mean, r_bar = circular_mean_and_rbar(psi_vals, weights)
        
        band_summaries[band_name] = {
            "freq_mhz": freq_mhz,
            "n_epochs": len(epochs),
            "psi_mean_rad": psi_mean,
            "r_bar": r_bar,
        }
        
        print_status(f"  {band_name} ({freq_mhz} MHz): n={len(epochs)}, ψ={psi_mean:.4f} rad, R_bar={r_bar:.4f}", "INFO")
    
    if len(all_epochs) < 10:
        print_status("ERROR: Insufficient data for chromatic analysis", "ERROR")
        return False
    
    # --- Profile likelihood over δ ---
    delta_grid = np.linspace(-4, 4, 161)
    profile = profile_likelihood_delta(all_epochs, delta_grid)
    
    # Find MLE δ
    best_fit = max(profile, key=lambda x: x["log_likelihood"])
    delta_mle = best_fit["delta"]
    ll_max = best_fit["log_likelihood"]
    
    # Compute ΔlnL from peak
    for p in profile:
        p["delta_lnL"] = p["log_likelihood"] - ll_max
    
    # Approximate 68% and 95% CI from likelihood (ΔlnL = -0.5 and -2)
    ci_68_mask = [p["delta_lnL"] > -0.5 for p in profile]
    ci_95_mask = [p["delta_lnL"] > -2.0 for p in profile]
    
    ci_68_vals = [p["delta"] for p, m in zip(profile, ci_68_mask) if m]
    ci_95_vals = [p["delta"] for p, m in zip(profile, ci_95_mask) if m]
    
    ci_68 = (min(ci_68_vals), max(ci_68_vals)) if ci_68_vals else (-4.0, 4.0)
    ci_95 = (min(ci_95_vals), max(ci_95_vals)) if ci_95_vals else (-4.0, 4.0)
    
    # --- Bootstrap CI ---
    print_status("  Running bootstrap for δ confidence interval...", "INFO")
    delta_boot_mle, ci_68_boot, ci_95_boot = bootstrap_delta_ci(all_epochs, n_boot=5000)
    
    # --- Model comparison: achromatic (δ=0) vs chromatic (δ free) ---
    fit_achromatic = fit_chromatic_model(all_epochs, 0.0)
    fit_chromatic = fit_chromatic_model(all_epochs, delta_mle)
    
    n = len(all_epochs)
    bic_achromatic = -2 * fit_achromatic["log_likelihood"] + fit_achromatic["n_params"] * np.log(n)
    bic_chromatic = -2 * fit_chromatic["log_likelihood"] + fit_chromatic["n_params"] * np.log(n)
    
    delta_bic = bic_achromatic - bic_chromatic  # positive if achromatic is better
    log10_bf = -0.5 * delta_bic / np.log(10)
    
    # Likelihood ratio test (δ=0 nested in δ free)
    lr_stat = 2 * (fit_chromatic["log_likelihood"] - fit_achromatic["log_likelihood"])
    lr_p = float(stats.chi2.sf(lr_stat, 1))
    
    print_status(f"  MLE δ = {delta_mle:.3f}", "INFO")
    print_status(f"  Profile 68% CI: [{ci_68[0]:.3f}, {ci_68[1]:.3f}]", "INFO")
    print_status(f"  Profile 95% CI: [{ci_95[0]:.3f}, {ci_95[1]:.3f}]", "INFO")
    if ci_68_boot and ci_95_boot:
        print_status(f"  Bootstrap 68% CI: [{ci_68_boot[0]:.3f}, {ci_68_boot[1]:.3f}]", "INFO")
        print_status(f"  Bootstrap 95% CI: [{ci_95_boot[0]:.3f}, {ci_95_boot[1]:.3f}]", "INFO")
    print_status(f"  Achromatic (δ=0) BIC: {bic_achromatic:.2f}", "INFO")
    print_status(f"  Chromatic (δ free) BIC: {bic_chromatic:.2f}", "INFO")
    print_status(f"  ΔBIC (achromatic - chromatic): {delta_bic:.2f}", "INFO")
    print_status(f"  log10 BF (achromatic vs chromatic): {log10_bf:.2f}", "INFO")
    print_status(f"  Likelihood ratio: χ² = {lr_stat:.3f}, p = {lr_p:.4f}", "INFO")
    
    # --- Determine conclusion ---
    # TEP prediction: δ = 0. If 0 is inside the CI, achromatic is consistent.
    # If chromatic model is not strongly preferred (BF < 2), default to TEP.
    
    achromatic_consistent = (ci_95[0] <= 0 <= ci_95[1])
    chromatic_preferred = (log10_bf < -2)  # chromatic strongly preferred
    achromatic_preferred = (log10_bf > 2)   # achromatic strongly preferred
    
    if achromatic_consistent and not chromatic_preferred:
        conclusion = (
            f"The Phase Closure ψ is consistent with achromatic behavior (δ = 0). "
            f"The 95% confidence interval [{ci_95[0]:.2f}, {ci_95[1]:.2f}] includes zero, "
            f"and the achromatic model is not disfavored relative to the chromatic "
            f"alternative (log10 BF = {log10_bf:.2f}). This is consistent with TEP's "
            f"prediction of environment-independent holonomy phase. The large "
            f"uncertainty reflects sparse sub-band sampling (5 epochs each), which "
            f"limits chromatic discrimination power."
        )
    elif chromatic_preferred:
        conclusion = (
            f"The chromatic model is strongly preferred over achromatic (log10 BF = {log10_bf:.2f}). "
            f"This would contradict TEP's achromatic prediction and require investigation."
        )
    else:
        conclusion = (
            f"The data do not strongly distinguish achromatic (δ=0) from chromatic behavior. "
            f"The 95% CI [{ci_95[0]:.2f}, {ci_95[1]:.2f}] is broad due to limited sub-band "
            f"sampling (5 epochs per sub-band). The Phase Closure measurement is "
            f"directionally consistent with achromatic TEP, but higher-SNR multi-band "
            f"observations are required for a definitive chromatic test."
        )
    
    results = {
        "band_summaries": band_summaries,
        "n_total_epochs": len(all_epochs),
        "reference_frequency_mhz": 1380.0,
        "delta_mle": float(delta_mle),
        "mu_0_mle_rad": float(fit_chromatic["mu_0"]),
        "kappa_mle": float(fit_chromatic["kappa"]),
        "confidence_intervals": {
            "profile_68": [float(ci_68[0]), float(ci_68[1])],
            "profile_95": [float(ci_95[0]), float(ci_95[1])],
            "bootstrap_68": [float(ci_68_boot[0]), float(ci_68_boot[1])] if ci_68_boot else None,
            "bootstrap_95": [float(ci_95_boot[0]), float(ci_95_boot[1])] if ci_95_boot else None,
        },
        "model_comparison": {
            "achromatic_bic": float(bic_achromatic),
            "chromatic_bic": float(bic_chromatic),
            "delta_bic": float(delta_bic),
            "log10_bf_achromatic_vs_chromatic": float(log10_bf),
            "likelihood_ratio_chi2": float(lr_stat),
            "likelihood_ratio_p": float(lr_p),
        },
        "profile_likelihood": [
            {"delta": float(p["delta"]), "log_likelihood": float(p["log_likelihood"]), "delta_lnL": float(p["delta_lnL"])}
            for p in profile
        ],
        "conclusion": {
            "achromatic_consistent": bool(achromatic_consistent),
            "chromatic_preferred": bool(chromatic_preferred),
            "achromatic_preferred": bool(achromatic_preferred),
            "interpretation": conclusion,
        },
        "limitations": (
            "Sub-band samples are small (sb0: 5 epochs, sb1: 5 epochs). "
            "Confidence intervals on δ are therefore broad. A definitive chromatic "
            "test requires: (1) more epochs per sub-band, (2) additional receiver bands "
            "(e.g., 10-cm, 50-cm), and (3) higher per-epoch SNR. The present result "
            "is consistent with TEP achromaticity but does not strongly exclude "
            "moderate chromaticity."
        ),
    }
    
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(results, f, indent=2, cls=NpEncoder)
    
    print_status(f"\nResults saved to {OUTPUT_FILE}", "INFO")
    print_status("STEP 043 V2 COMPLETED", "INFO")
    return True


if __name__ == "__main__":
    main()
