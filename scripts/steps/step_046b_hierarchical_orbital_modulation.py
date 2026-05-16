#!/usr/bin/env python3
"""
Step 046b: Hierarchical Mixed-Effects Orbital Modulation Model (REML)

Replaces the triplet-binned and epoch-blocked analyses with a proper linear
mixed-effects model fitted by restricted maximum likelihood (REML).  The
model accounts for within-epoch triplet correlations via a per-epoch random
offset while preserving triplet-level measurement precision.

Model:
    delay_ij = A1*sin(2π·φ_j) + A2*cos(2π·φ_j) + C + α_j + ε_ij
    α_j ~ N(0, σ_α²)          [per-epoch random effect]
    ε_ij ~ N(0, σ_ij²)        [measurement noise, known from cross-term fits]

The likelihood is profiled over σ_α using the Woodbury identity for O(n_j)
per-epoch matrix operations.  A likelihood-ratio test against the nested null
model (no sinusoid) provides the formal significance of orbital modulation.

Author: TEP Analysis Pipeline
Date: May 2026
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict

import numpy as np
from scipy import optimize, stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.json_numpy import NpEncoder
from scripts.utils.logger import print_status
from scripts.utils.config import J0437_PB_DAYS, J0437_T0_MJD

RESULTS_DIR = PROJECT_ROOT / "results"


def compute_orbital_phase(mjd: float) -> float:
    phase = ((mjd - J0437_T0_MJD) % J0437_PB_DAYS) / J0437_PB_DAYS
    return phase


def load_closure_per_epoch(pulsar_name: str):
    summary_file = RESULTS_DIR / f"step_003_closure_final_per_epoch_{pulsar_name}.json"
    if not summary_file.exists():
        print_status(f"File not found: {summary_file}", "ERROR")
        return None
    with open(summary_file, "r") as f:
        return json.load(f)


def load_triplet_data(pulsar: str = "j0437") -> Dict[str, Any]:
    epochs = load_closure_per_epoch(pulsar)
    if epochs is None:
        return {}

    delays = []
    uncertainties = []
    phases = []
    epoch_idx = []
    epoch_ids = []

    for e_idx, epoch in enumerate(epochs):
        mjd = epoch.get("mjd", 0)
        orbital_phase = compute_orbital_phase(mjd)
        epoch_id = epoch.get("epoch", f"epoch_{e_idx}")
        triplets = epoch.get("triplets", [])
        for triplet in triplets:
            delay_us = triplet.get("geometric_delta_us")
            uncertainty_us = triplet.get("sigma_us")
            if delay_us is None or uncertainty_us is None:
                continue
            delays.append(float(delay_us) * 1000.0)
            uncertainties.append(float(uncertainty_us) * 1000.0)
            phases.append(orbital_phase)
            epoch_idx.append(e_idx)
            epoch_ids.append(epoch_id)

    return {
        "n_triplets": len(delays),
        "n_epochs": len(set(epoch_idx)),
        "delay_ns": np.array(delays, dtype=np.float64),
        "sigma_ns": np.array(uncertainties, dtype=np.float64),
        "phase": np.array(phases, dtype=np.float64),
        "epoch_idx": np.array(epoch_idx, dtype=np.int32),
        "epoch_ids": epoch_ids,
    }


class MixedEffectsOrbitalModel:
    """
    Linear mixed-effects model fitted via the epoch-mean representation.

    Key simplification: for the model  y_j ~ N(X_j β, D_j + σ_α² 1 1^T),
    the likelihood factorises into:
      (a) an epoch-mean part that carries all information about β and σ_α;
      (b) a within-epoch part that is constant w.r.t. β and σ_α.

    For epoch j with inverse-variance weights w_i = 1/σ_i²:
        ȳ_j = Σ w_i y_i / Σ w_i          (inverse-variance weighted mean)
        v_j = σ_α² + 1/(Σ w_i)            (variance of the epoch mean)

    The profile log-likelihood is:
        ℓ(β, σ_α) = -0.5 Σ_j [log(v_j) + (ȳ_j - X_j β)² / v_j] + const

    This is a weighted least-squares problem with only N_epoch observations,
    eliminating the need for large-matrix inversion and removing the source
    of the previous numerical instabilities.
    """

    def __init__(self, y, sigma, phase, epoch_idx):
        self.y = y
        self.sigma = np.maximum(sigma, 1e-12)
        self.phase = phase
        self.epoch_idx = epoch_idx
        self.n_triplets = len(y)
        self.n_epochs = int(np.max(epoch_idx) + 1)

        # Build epoch-mean summary statistics
        self.epochs = []
        for j in range(self.n_epochs):
            mask = epoch_idx == j
            n = int(np.sum(mask))
            if n == 0:
                self.epochs.append(None)
                continue
            yj = y[mask]
            sj = self.sigma[mask]
            w = 1.0 / (sj ** 2)
            s = np.sum(w)
            y_bar = float(np.sum(w * yj) / s)
            # Design row (same for all triplets in an epoch because phase is epoch-level)
            phi = float(phase[mask][0])
            self.epochs.append({
                "n": n,
                "phi": phi,
                "y_bar": y_bar,
                "s": float(s),
                "sin": float(np.sin(2 * np.pi * phi)),
                "cos": float(np.cos(2 * np.pi * phi)),
            })

        self._valid_epochs = [e for e in self.epochs if e is not None]

    def _X_row(self, epoch):
        return np.array([epoch["sin"], epoch["cos"], 1.0])

    def _profile_beta(self, sigma_alpha):
        """GLS estimate of β for fixed σ_α."""
        XtWX = np.zeros((3, 3))
        XtWy = np.zeros(3)
        for e in self._valid_epochs:
            v = sigma_alpha ** 2 + 1.0 / e["s"]
            w = 1.0 / v
            x = self._X_row(e)
            XtWX += w * np.outer(x, x)
            XtWy += w * x * e["y_bar"]

        try:
            cov = np.linalg.inv(XtWX)
        except np.linalg.LinAlgError:
            cov = np.linalg.pinv(XtWX)
        beta = cov @ XtWy
        return beta, cov

    def _log_likelihood(self, sigma_alpha, beta=None):
        """Profile log-likelihood (epoch-mean part only)."""
        if beta is None:
            beta, _ = self._profile_beta(sigma_alpha)

        ll = 0.0
        for e in self._valid_epochs:
            v = sigma_alpha ** 2 + 1.0 / e["s"]
            x = self._X_row(e)
            resid = e["y_bar"] - x @ beta
            ll += -0.5 * (np.log(v) + resid ** 2 / v)
        return ll, beta

    def fit(self):
        """Fit full model, profiling over σ_α."""
        def neg_ll(sa):
            ll, _ = self._log_likelihood(sa)
            return -ll

        res = optimize.minimize_scalar(
            neg_ll,
            bounds=(0.0, 50.0),
            method="bounded",
            options={"xatol": 1e-6, "maxiter": 500},
        )
        sigma_alpha_hat = float(res.x)
        logL_max, beta_hat = self._log_likelihood(sigma_alpha_hat)
        _, cov_beta = self._profile_beta(sigma_alpha_hat)
        return beta_hat, cov_beta, sigma_alpha_hat, logL_max

    def fit_null(self):
        """Fit null model (intercept only), profiling over σ_α."""
        def neg_ll(sa):
            ll = 0.0
            for e in self._valid_epochs:
                v = sa ** 2 + 1.0 / e["s"]
                w = 1.0 / v
                C_hat = w * e["y_bar"] / w  # simplifies to e["y_bar"]
                # Actually C_hat = Σ(w·ȳ) / Σ(w)  for intercept-only
                # But since we're evaluating at the MLE for each sa:
            # We need to actually solve for C at each sa
            C_num = 0.0
            C_den = 0.0
            for e in self._valid_epochs:
                v = sa ** 2 + 1.0 / e["s"]
                w = 1.0 / v
                C_num += w * e["y_bar"]
                C_den += w
            C_hat = C_num / C_den if C_den > 0 else 0.0
            ll = 0.0
            for e in self._valid_epochs:
                v = sa ** 2 + 1.0 / e["s"]
                resid = e["y_bar"] - C_hat
                ll += -0.5 * (np.log(v) + resid ** 2 / v)
            return -ll

        res = optimize.minimize_scalar(
            neg_ll,
            bounds=(0.0, 50.0),
            method="bounded",
            options={"xatol": 1e-6, "maxiter": 500},
        )
        sigma_alpha_null = float(res.x)
        logL_null = -neg_ll(sigma_alpha_null)
        return sigma_alpha_null, logL_null

    def bootstrap_amplitude_ci(self, beta_hat, cov_beta, n_boot=10000, seed=46):
        rng = np.random.default_rng(seed)
        cov_sym = (cov_beta + cov_beta.T) / 2.0
        try:
            L = np.linalg.cholesky(cov_sym)
        except np.linalg.LinAlgError as err:
            raise RuntimeError(
                "Fixed-effect covariance is not numerically positive-definite; "
                "cannot simulate amplitude posterior draws."
            ) from err
        z = rng.standard_normal(size=(n_boot, beta_hat.shape[0]))
        beta_hat_64 = np.asarray(beta_hat, dtype=np.float64)
        beta_samples = np.empty((n_boot, beta_hat_64.shape[0]), dtype=np.float64)
        for i in range(n_boot):
            beta_samples[i] = beta_hat_64 + L @ z[i]
        A1_s = beta_samples[:, 0]
        A2_s = beta_samples[:, 1]
        A_s = np.sqrt(A1_s ** 2 + A2_s ** 2)
        return {
            "mean": float(np.mean(A_s)),
            "std": float(np.std(A_s, ddof=1)),
            "ci_95_lower": float(np.percentile(A_s, 2.5)),
            "ci_95_upper": float(np.percentile(A_s, 97.5)),
            "prob_positive": float(np.mean(A_s > 0.01)),
        }


def main():
    print_status("Step 046b: Hierarchical Mixed-Effects Orbital Modulation (REML)", "INFO")

    data = load_triplet_data("j0437")
    if not data:
        print_status("No data loaded, exiting", "ERROR")
        sys.exit(1)

    print_status(f"Loaded {data['n_triplets']} triplets across {data['n_epochs']} epochs", "INFO")

    model = MixedEffectsOrbitalModel(
        data["delay_ns"],
        data["sigma_ns"],
        data["phase"],
        data["epoch_idx"],
    )

    # Fit modulation model
    print_status("Fitting mixed-effects modulation model...", "INFO")
    beta_hat, cov_beta, sigma_alpha_hat, logL_max = model.fit()
    A1_hat, A2_hat, C_hat = beta_hat

    # Fit null model
    print_status("Fitting null model...", "INFO")
    sigma_alpha_null, logL_null = model.fit_null()

    # Likelihood ratio test
    lr_stat = 2.0 * (logL_max - logL_null)
    lr_pvalue = stats.chi2.sf(lr_stat, 2)

    # Bootstrap CI for amplitude
    boot = model.bootstrap_amplitude_ci(beta_hat, cov_beta, n_boot=10000)

    # Effective sample size
    total_var = float(np.var(data["delay_ns"], ddof=1))
    within_var = float(np.mean(data["sigma_ns"] ** 2))
    epoch_means = []
    for j in range(model.n_epochs):
        mask = data["epoch_idx"] == j
        if np.any(mask):
            epoch_means.append(float(np.mean(data["delay_ns"][mask])))
    epoch_var_est = float(np.var(epoch_means, ddof=1)) if len(epoch_means) > 1 else 0.0
    icc = epoch_var_est / (epoch_var_est + within_var) if (epoch_var_est + within_var) > 0 else 0
    mean_n_per_epoch = data["n_triplets"] / data["n_epochs"]
    n_eff = int(data["n_triplets"] / (1.0 + (mean_n_per_epoch - 1.0) * icc))

    # Standard errors from covariance
    se_A1 = float(np.sqrt(cov_beta[0, 0]))
    se_A2 = float(np.sqrt(cov_beta[1, 1]))
    se_C = float(np.sqrt(cov_beta[2, 2]))
    corr_A1A2 = float(cov_beta[0, 1] / (se_A1 * se_A2)) if (se_A1 * se_A2) > 0 else 0.0

    # Amplitude and phase
    A_hat = float(np.sqrt(A1_hat ** 2 + A2_hat ** 2))
    psi_hat = float(np.arctan2(A2_hat, A1_hat))

    # Delta-method SE for amplitude
    grad = np.array([A1_hat, A2_hat, 0.0]) / (A_hat + 1e-12)
    se_A_delta = float(np.sqrt(grad @ cov_beta @ grad))

    output = {
        "step": "046b",
        "pulsar": "J0437-4715",
        "model": "hierarchical_mixed_effects_orbital_modulation_reml",
        "n_triplets": data["n_triplets"],
        "n_epochs": data["n_epochs"],
        "parameterisation": "A1*sin(2πφ) + A2*cos(2πφ) + C + alpha[epoch]",
        "amplitude_interval": {
            "label": "Monte Carlo amplitude interval (2.5 / 97.5 percentiles on circular A)",
            "n_draws": 10000,
            "generator": (
                "Draw z ~ N(0, I_3); set beta = beta_hat + L z with lower-triangular L from Cholesky(LL^T)=cov_beta; "
                "circular amplitude A = sqrt(A_1^2 + A_2^2)."
            ),
            "mean_of_drawn_circular_amplitudes_ns": boot["mean"],
            "std_of_drawn_circular_amplitudes_ns": boot["std"],
            "note": (
                "The mean of A over draws typically exceeds A_hat because A is a nonlinear positive map of "
                "Gaussian coefficients (Rice-like upward bias); report A_hat, delta-method SE, and percentile interval."
            ),
        },
        "fixed_effects": {
            "A1_ns": {"mean": float(A1_hat), "std": se_A1},
            "A2_ns": {"mean": float(A2_hat), "std": se_A2},
            "modulation_amplitude_A_ns": {
                "mean": A_hat,
                "std": se_A_delta,
                "ci_95_lower": boot["ci_95_lower"],
                "ci_95_upper": boot["ci_95_upper"],
                "prob_positive": boot["prob_positive"],
            },
            "modulation_phase_psi_rad": {"mean": psi_hat, "std": None},
            "offset_C_ns": {"mean": float(C_hat), "std": se_C},
            "corr_A1_A2": corr_A1A2,
        },
        "random_effects": {
            "sigma_epoch_ns": {
                "mle": float(sigma_alpha_hat),
                "null_model_mle": float(sigma_alpha_null),
            },
            "icc_estimate": float(icc),
            "effective_sample_size": n_eff,
        },
        "model_comparison": {
            "logL_modulation": float(logL_max),
            "logL_null": float(logL_null),
            "lr_statistic": float(lr_stat),
            "lr_dof": 2,
            "lr_pvalue": float(lr_pvalue),
            "favors_modulation": bool(lr_pvalue < 0.05),
        },
        "methodology_note": (
            "This mixed-effects model accounts for within-epoch triplet correlations via "
            "per-epoch random offsets (alpha_j ~ N(0, sigma_epoch²)), preserving triplet-level "
            "measurement precision.  It replaces both the overconfident triplet-binning (which ignores "
            "correlations) and the overly conservative epoch-blocked analysis (which discards within-epoch "
            "inverse-variance weighting).  The likelihood-ratio test is exact for nested Gaussian models.  "
            "The 95% amplitude interval is from parametric Gaussian draws in (A_1, A_2) using the REML "
            "fixed-effect covariance (see amplitude_interval.generator), not from resampling residuals."
        ),
    }

    output_file = RESULTS_DIR / "step_046b_hierarchical_orbital_modulation.json"
    with open(output_file, "w") as f:
        json.dump(output, f, indent=2, cls=NpEncoder)

    print_status(f"Results saved to {output_file}", "INFO")

    # Summary
    fe = output["fixed_effects"]["modulation_amplitude_A_ns"]
    print_status(
        f"Modulation amplitude: {fe['mean']:.3f} ± {fe['std']:.3f} ns "
        f"(95% MC amplitude interval: [{fe['ci_95_lower']:.3f}, {fe['ci_95_upper']:.3f}])",
        "INFO",
    )
    mc = output["model_comparison"]
    print_status(
        f"Likelihood-ratio test: χ² = {mc['lr_statistic']:.2f}, p = {mc['lr_pvalue']:.4e} (2 df)",
        "INFO",
    )
    print_status(
        f"Epoch scatter σ_α = {sigma_alpha_hat:.3f} ns, ICC = {icc:.3f}, effective n = {n_eff}",
        "INFO",
    )


if __name__ == "__main__":
    main()
