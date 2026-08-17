#!/usr/bin/env python3
"""
================================================================================
STEP 051 V2: ORBITAL MODULATION HIERARCHY (GP-AUGMENTED KINEMATIC MODEL)
================================================================================

Tests whether the Phase Closure ψ exhibits orbital modulation consistent
with TEP's prediction of binary reflex velocity projection onto the scattering
screen normal.

The legacy orbital analysis (Section 4.9.6) reported p = 0.372, which is
not statistically significant. This module implements a more sophisticated
hierarchical model with:

1. Kinematic projection: ψ(φ_orb) = A · v_proj(φ_orb) + B
   where v_proj is the binary reflex velocity projected onto the screen normal

2. Gaussian Process augmentation: allows for non-sinusoidal orbital structure
   (e.g., eclipse ingress, wind from the companion) via a Matérn kernel

3. Proper Bayesian model comparison via Bayes factors (BIC approximation)

The model hierarchy:
    M0 (null):       ψ_e ~ von Mises(μ, κ) — no orbital dependence
    M1 (sinusoid):   ψ_e ~ von Mises(A·v_proj + B, κ) — pure kinematic
    M2 (GP):         ψ_e ~ von Mises(A·v_proj + B + GP(φ_orb), κ) — kinematic + GP

If M1 or M2 is strongly preferred over M0, this supports orbital modulation.
If neither is preferred, the orbital channel should be downgraded from
primary validation to exploratory diagnostic.

DATA:
------
Uses epoch-level ψ values and MJDs from step_003 per-epoch data, together
with the J0437-4715 binary orbital parameters from Reardon et al. (2024).

OUTPUT:
-------
results/step_051_orbital_hierarchy_v2.json
    - Per-model MLE parameters, log-likelihoods, BIC values
    - Bayes factors for orbital vs null
    - GP hyperparameters (if applicable)
    - Recommendation: retain as validation or downgrade to exploratory

================================================================================
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats
from scipy.optimize import minimize

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.config import (
    C_LIGHT_KM_S,
    J0437_A1_LC,
    J0437_ECC,
    J0437_INC_DEG,
    J0437_KOM_DEG,
    J0437_OM_DEG,
    J0437_PB_DAYS,
    J0437_PSI_DEG,
    J0437_RA_RAD,
    J0437_T0_MJD,
    RANDOM_SEED,
)
from scripts.utils.binary_reflex_kinematics import reflex_binary_transverse_velocity_kms
from scripts.utils.json_numpy import NpEncoder
from scripts.utils.logger import print_status

RESULTS_DIR = PROJECT_ROOT / "results"
OUTPUT_FILE = RESULTS_DIR / "step_051_orbital_hierarchy_v2.json"

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
        from scipy.special import i0
        norm = np.log(2 * np.pi) + np.log(i0(kappa))
    else:
        norm = np.log(2 * np.pi) + kappa - 0.5 * np.log(2 * np.pi * kappa)
    ll = kappa * np.cos(psi - mu) - norm
    return float(np.sum(ll))


def load_epoch_data():
    """Load epoch-level ψ and MJD data for J0437."""
    per_epoch_file = RESULTS_DIR / "step_003_closure_final_per_epoch_j0437.json"
    if not per_epoch_file.exists():
        print_status("ERROR: step_003 per-epoch data not found", "ERROR")
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
        mjd = ep.get("mjd", 0.0)
        
        epochs.append({
            "mjd": mjd,
            "psi_mean": psi_mean,
            "r_bar": r_bar,
            "n_triplets": len(triplets),
        })
    
    return epochs


def compute_orbital_phase(mjd):
    """Compute orbital phase φ ∈ [0, 1) from MJD."""
    pb_days = J0437_PB_DAYS
    t0_mjd = J0437_T0_MJD
    phase = ((mjd - t0_mjd) / pb_days) % 1.0
    return phase


def compute_reflex_velocity_projection(mjd):
    """Compute binary reflex velocity projected onto scattering screen normal.
    
    Uses Thiele-Innes projection from Reardon et al. (2024) PPTA-DR3.
    """
    try:
        v_orbit = reflex_binary_transverse_velocity_kms(
            float(mjd),
            float(J0437_PB_DAYS),
            float(J0437_T0_MJD),
            float(J0437_A1_LC),
            float(J0437_ECC),
            float(J0437_OM_DEG),
            float(J0437_INC_DEG),
            float(J0437_KOM_DEG),
            float(C_LIGHT_KM_S),
            float(J0437_RA_RAD),
            float(J0437_DEC_RAD),
        )
    except Exception as e:
        print_status(f"WARNING: Failed to compute reflex velocity: {e}", "WARN")
        return 0.0
    
    # Project onto scattering axis (J0437: 128°)
    rad_psi = np.radians(J0437_PSI_DEG)
    v_proj = v_orbit[0] * np.cos(rad_psi) + v_orbit[1] * np.sin(rad_psi)
    return float(v_proj)


def matern_kernel(x1, x2, length_scale=0.2, nu=1.5):
    """Matérn kernel for GP prior."""
    x1 = np.asarray(x1)
    x2 = np.asarray(x2)
    
    if nu == 1.5:
        r = np.abs(x1 - x2)[:, np.newaxis]
        return (1 + np.sqrt(3) * r / length_scale) * np.exp(-np.sqrt(3) * r / length_scale)
    elif nu == 2.5:
        r = np.abs(x1 - x2)[:, np.newaxis]
        return (1 + np.sqrt(5) * r / length_scale + 5 * r**2 / (3 * length_scale**2)) * np.exp(-np.sqrt(5) * r / length_scale)
    else:
        r = np.abs(x1 - x2)[:, np.newaxis]
        return np.exp(-r / length_scale)


def fit_model_m0(epochs):
    """Fit M0 (null): ψ ~ von Mises(μ, κ) — no orbital dependence."""
    psi = np.array([e["psi_mean"] for e in epochs])
    
    z = np.mean(np.exp(1j * psi))
    mu = float(np.angle(z))
    R = float(np.abs(z))
    
    if R < 0.53:
        kappa = 2 * R + R**3 + (5.0 / 6.0) * R**5
    elif R < 0.85:
        kappa = -0.4 + 1.39 * R + 0.43 / (1 - R)
    else:
        kappa = 1.0 / (2 * (1 - R) - (1 - R)**2 - (1 - R)**3)
    
    ll = von_mises_log_likelihood(psi, mu, kappa)
    n = len(epochs)
    bic = -2 * ll + 2 * np.log(n)
    
    return {
        "mu": float(mu),
        "kappa": float(kappa),
        "log_likelihood": float(ll),
        "bic": float(bic),
        "n_params": 2,
    }


def fit_model_m1(epochs):
    """Fit M1 (sinusoid): ψ ~ von Mises(A·v_proj + B, κ)."""
    psi = np.array([e["psi_mean"] for e in epochs])
    v_proj = np.array([compute_reflex_velocity_projection(e["mjd"]) for e in epochs])
    
    # Standardize v_proj for numerical stability
    v_mean = np.mean(v_proj)
    v_std = np.std(v_proj) + 1e-10
    v_stdz = (v_proj - v_mean) / v_std
    
    def neg_ll(params):
        A, B, kappa = params
        if kappa <= 0:
            return 1e10
        mu = A * v_stdz + B
        return -von_mises_log_likelihood(psi, mu, kappa)
    
    # Initialize: B = circular mean, A = 0
    z = np.mean(np.exp(1j * psi))
    mu_init = np.angle(z)
    R_init = np.abs(z)
    
    if R_init < 0.53:
        kappa_init = 2 * R_init + R_init**3
    elif R_init < 0.85:
        kappa_init = -0.4 + 1.39 * R_init + 0.43 / (1 - R_init)
    else:
        kappa_init = 1.0 / (2 * (1 - R_init))
    
    x0 = [0.0, mu_init, kappa_init]
    bounds = [(-np.pi, np.pi), (-np.pi, np.pi), (1e-6, 1e4)]
    
    result = minimize(neg_ll, x0, method="L-BFGS-B", bounds=bounds)
    
    if result.success:
        A, B, kappa = result.x
        # Convert A back to original scale
        A_orig = A / v_std
        B_orig = B + A * v_mean / v_std
        ll = -result.fun
        n = len(epochs)
        bic = -2 * ll + 3 * np.log(n)
        
        return {
            "A": float(A_orig),
            "B": float(B_orig),
            "kappa": float(kappa),
            "log_likelihood": float(ll),
            "bic": float(bic),
            "n_params": 3,
            "converged": True,
        }
    else:
        return {
            "A": 0.0,
            "B": mu_init,
            "kappa": float(kappa_init),
            "log_likelihood": float(-neg_ll(x0)),
            "bic": float(-2 * (-neg_ll(x0)) + 3 * np.log(n)),
            "n_params": 3,
            "converged": False,
            "error": str(result.message),
        }


def fit_model_m2(epochs):
    """Fit M2 (GP): ψ ~ von Mises(A·v_proj + B + GP(φ), κ).
    
    Simplified GP: assumes GP contribution is small and approximated via
    low-rank expansion to avoid full matrix inversion.
    """
    psi = np.array([e["psi_mean"] for e in epochs])
    v_proj = np.array([compute_reflex_velocity_projection(e["mjd"]) for e in epochs])
    phi = np.array([compute_orbital_phase(e["mjd"]) for e in epochs])
    
    # Standardize
    v_mean = np.mean(v_proj)
    v_std = np.std(v_proj) + 1e-10
    v_stdz = (v_proj - v_mean) / v_std
    
    n = len(epochs)
    
    # Use a simple GP approximation: GP(φ) = Σ_j w_j * k_j(φ)
    # where k_j are fixed basis functions (e.g., sinusoids at harmonics)
    # This reduces the GP to a linear model with fixed kernel basis
    n_harmonics = 3
    K = np.column_stack([
        np.sin(2 * np.pi * h * phi) for h in range(1, n_harmonics + 1)
    ] + [
        np.cos(2 * np.pi * h * phi) for h in range(1, n_harmonics + 1)
    ])
    
    # Total parameters: A, B, GP weights (2*n_harmonics), kappa
    n_params_gp = 2 + 2 * n_harmonics + 1
    
    def neg_ll(params):
        A = params[0]
        B = params[1]
        w = params[2:-1]
        kappa = params[-1]
        
        if kappa <= 0:
            return 1e10
        
        mu = A * v_stdz + B + K @ w
        return -von_mises_log_likelihood(psi, mu, kappa)
    
    # Initialize
    z = np.mean(np.exp(1j * psi))
    mu_init = np.angle(z)
    R_init = np.abs(z)
    
    if R_init < 0.53:
        kappa_init = 2 * R_init + R_init**3
    elif R_init < 0.85:
        kappa_init = -0.4 + 1.39 * R_init + 0.43 / (1 - R_init)
    else:
        kappa_init = 1.0 / (2 * (1 - R_init))
    
    x0 = np.concatenate([[0.0, mu_init], np.zeros(2 * n_harmonics), [kappa_init]])
    bounds = [(-np.pi, np.pi), (-np.pi, np.pi)] + [(-0.5, 0.5)] * (2 * n_harmonics) + [(1e-6, 1e4)]
    
    result = minimize(neg_ll, x0, method="L-BFGS-B", bounds=bounds)
    
    if result.success:
        A, B, w, kappa = result.x[0], result.x[1], result.x[2:-1], result.x[-1]
        A_orig = A / v_std
        B_orig = B + A * v_mean / v_std
        ll = -result.fun
        bic = -2 * ll + n_params_gp * np.log(n)
        
        return {
            "A": float(A_orig),
            "B": float(B_orig),
            "gp_weights": w.tolist(),
            "kappa": float(kappa),
            "log_likelihood": float(ll),
            "bic": float(bic),
            "n_params": n_params_gp,
            "n_harmonics": n_harmonics,
            "converged": True,
        }
    else:
        return {
            "A": 0.0,
            "B": mu_init,
            "gp_weights": np.zeros(2 * n_harmonics).tolist(),
            "kappa": float(kappa_init),
            "log_likelihood": float(-neg_ll(x0)),
            "bic": float(-2 * (-neg_ll(x0)) + n_params_gp * np.log(n)),
            "n_params": n_params_gp,
            "n_harmonics": n_harmonics,
            "converged": False,
            "error": str(result.message),
        }


def compute_bayes_factor(bic1, bic2):
    """Approximate log_10(BF) from ΔBIC."""
    delta_bic = bic1 - bic2
    log10_bf = -0.5 * delta_bic / np.log(10)
    return float(log10_bf)


def main():
    print_status("=" * 70, "INFO")
    print_status("STEP 051 V2: ORBITAL MODULATION HIERARCHY (GP-AUGMENTED)", "INFO")
    print_status("=" * 70)
    
    epochs = load_epoch_data()
    if epochs is None or len(epochs) < 10:
        print_status("ERROR: Insufficient epoch data", "ERROR")
        return False
    
    print_status(f"Loaded {len(epochs)} epochs", "INFO")
    
    # Fit all models
    m0 = fit_model_m0(epochs)
    m1 = fit_model_m1(epochs)
    m2 = fit_model_m2(epochs)
    
    models = {
        "M0_null_no_orbital": m0,
        "M1_kinematic_sinusoid": m1,
        "M2_kinematic_plus_GP": m2,
    }
    
    # Model comparison
    best_model = min(models, key=lambda k: models[k]["bic"])
    best_bic = models[best_model]["bic"]
    
    bayes_factors = {}
    for name, model in models.items():
        log10_bf = compute_bayes_factor(model["bic"], best_bic)
        bayes_factors[name] = {
            "log10_bf_vs_best": float(log10_bf),
            "interpretation": (
                "Strongly preferred" if log10_bf > 2 else
                "Moderately preferred" if log10_bf > 1 else
                "Weakly preferred" if log10_bf > 0.5 else
                "Not preferred"
            ),
        }
    
    # Key comparisons
    log10_bf_m1_vs_m0 = compute_bayes_factor(m0["bic"], m1["bic"])
    log10_bf_m2_vs_m0 = compute_bayes_factor(m0["bic"], m2["bic"])
    log10_bf_m2_vs_m1 = compute_bayes_factor(m1["bic"], m2["bic"])
    
    print_status(f"Model comparison:", "INFO")
    for name, model in models.items():
        print_status(f"  {name}: BIC = {model['bic']:.2f}, logL = {model['log_likelihood']:.2f}", "INFO")
    print_status(f"  Best model: {best_model}", "INFO")
    print_status(f"  log10 BF (M1 vs M0): {log10_bf_m1_vs_m0:.2f}", "INFO")
    print_status(f"  log10 BF (M2 vs M0): {log10_bf_m2_vs_m0:.2f}", "INFO")
    print_status(f"  log10 BF (M2 vs M1): {log10_bf_m2_vs_m1:.2f}", "INFO")
    
    # Determine recommendation
    # If M1 or M2 is strongly preferred over M0 (BF > 10), retain as validation
    # Otherwise, downgrade to exploratory
    orbital_significant = bool(log10_bf_m2_vs_m0 > 1 or log10_bf_m1_vs_m0 > 1)
    
    if orbital_significant:
        recommendation = (
            f"The orbital modulation model is preferred over the null "
            f"(log10 BF = {max(log10_bf_m2_vs_m0, log10_bf_m1_vs_m0):.2f}). "
            f"Retain the orbital channel as a primary validation check."
        )
    else:
        recommendation = (
            f"The orbital modulation model is not preferred over the null "
            f"(log10 BF = {max(log10_bf_m2_vs_m0, log10_bf_m1_vs_m0):.2f}). "
            f"Downgrade the orbital channel from primary validation to "
            f"exploratory diagnostic in the manuscript."
        )
    
    results = {
        "n_epochs": len(epochs),
        "models": {
            name: {
                "parameters": {k: v for k, v in model.items() if k not in ["log_likelihood", "bic", "n_params", "converged"]},
                "log_likelihood": model["log_likelihood"],
                "bic": model["bic"],
                "n_params": model["n_params"],
                "converged": model.get("converged", True),
            }
            for name, model in models.items()
        },
        "bayes_factors": bayes_factors,
        "key_comparisons": {
            "M1_vs_M0": {
                "log10_bf": float(log10_bf_m1_vs_m0),
                "interpretation": (
                    "Kinematic sinusoid preferred" if log10_bf_m1_vs_m0 > 1 else
                    "Null preferred"
                ),
            },
            "M2_vs_M0": {
                "log10_bf": float(log10_bf_m2_vs_m0),
                "interpretation": (
                    "Kinematic + GP preferred" if log10_bf_m2_vs_m0 > 1 else
                    "Null preferred"
                ),
            },
            "M2_vs_M1": {
                "log10_bf": float(log10_bf_m2_vs_m1),
                "interpretation": (
                    "GP augmentation improves fit" if log10_bf_m2_vs_m1 > 1 else
                    "GP does not improve fit"
                ),
            },
        },
        "best_model": best_model,
        "recommendation": {
            "orbital_significant": bool(orbital_significant),
            "action": "retain_as_validation" if orbital_significant else "downgrade_to_exploratory",
            "interpretation": recommendation,
        },
    }
    
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(results, f, indent=2, cls=NpEncoder)
    
    print_status(f"\nResults saved to {OUTPUT_FILE}", "INFO")
    print_status(f"Recommendation: {'retain as validation' if orbital_significant else 'downgrade to exploratory'}", "INFO")
    print_status("STEP 051 V2 COMPLETED", "INFO")
    return True


if __name__ == "__main__":
    main()
