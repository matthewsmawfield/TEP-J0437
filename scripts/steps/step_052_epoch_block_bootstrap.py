#!/usr/bin/env python3
"""
Step 052: Epoch-Block Bootstrap for H_trim Excess
==================================================

Addresses the reviewer concern that the reported 79.5σ and 16.4σ H_trim excess
significances divide the excess by only the Monte Carlo noise-floor uncertainty,
ignoring the measured trimmed-amplitude uncertainty.

This step performs a proper epoch-block bootstrap that, in every resample:
  1. Resamples complete epochs (with replacement).
  2. Recomputes orientation grouping (positive/negative geom_sign).
  3. Recomputes the 10% trimmed estimator on epoch-level |H| means.
  4. Re-estimates the within-orientation noise floor.
  5. Records their difference (the excess).

The empirical significance is the fraction of bootstrap resamples whose excess
is ≤ 0, giving a one-sided p-value. The 95% CI is from the bootstrap distribution
of the excess.

Outputs:
  results/step_052_epoch_block_bootstrap.json
"""

import json
import numpy as np
from pathlib import Path
import time

RESULTS_DIR = Path(__file__).parent.parent.parent / "results"
RNG_SEED = 42 + 52  # step_052


def load_epoch_data(pulsar_file: str, min_triplets: int = 5):
    """Load per-epoch data from step_003 output."""
    data = json.load(open(pulsar_file))
    epoch_data = []
    for epoch in data:
        triplets = epoch["triplets"]
        if len(triplets) < min_triplets:
            continue
        h_vals = np.array([abs(t["geometric_delta_us"]) * 1000 for t in triplets])
        signs = np.array([t.get("geom_sign", 0) for t in triplets])
        sigmas = np.array([t.get("sigma_us", 0) * 1000 for t in triplets])

        weights = 1.0 / sigmas**2
        weights = np.where(np.isfinite(weights) & (weights > 0), weights, 0)
        if weights.sum() == 0:
            continue
        h_mean = float(np.sum(weights * h_vals) / weights.sum())

        pos_mask = signs > 0
        neg_mask = signs < 0

        def mad_sigma(vals):
            if len(vals) < 3:
                return float(np.std(vals)) if len(vals) > 0 else 0.0
            return float(np.median(np.abs(vals - np.median(vals))) * 1.4826)

        pos_sigma = mad_sigma(h_vals[pos_mask])
        neg_sigma = mad_sigma(h_vals[neg_mask])

        n_total = pos_mask.sum() + neg_mask.sum()
        if n_total > 0:
            within_orient_sigma = np.sqrt(
                (pos_sigma**2 * pos_mask.sum() + neg_sigma**2 * neg_mask.sum())
                / n_total
            )
        else:
            within_orient_sigma = 0.0

        epoch_data.append(
            {
                "h_mean": h_mean,
                "within_orient_sigma": within_orient_sigma,
                "n_triplets": len(triplets),
            }
        )
    return epoch_data


def compute_trim_excess(epoch_data):
    """Compute H_trim and within-orientation noise floor excess."""
    h_means = np.array([e["h_mean"] for e in epoch_data])
    within_sigmas = np.array([e["within_orient_sigma"] for e in epoch_data])

    n = len(h_means)
    if n < 10:
        return None

    # 10% trimmed mean
    sorted_h = np.sort(h_means)
    n_trim = int(0.1 * n)
    trimmed = sorted_h[n_trim : n - n_trim]
    H_trim = float(np.mean(trimmed))

    # Monte Carlo noise floor: E[|H|] = sigma * sqrt(2/pi) using within-orientation sigma
    mc_noise_floor = float(np.mean(within_sigmas) * np.sqrt(2.0 / np.pi))

    excess = H_trim - mc_noise_floor
    return {
        "H_trim": H_trim,
        "mc_noise_floor": mc_noise_floor,
        "excess": excess,
        "n_epochs": n,
    }


def epoch_block_bootstrap(epoch_data, n_bootstrap=10000, seed=RNG_SEED):
    """Epoch-block bootstrap: resample complete epochs, recompute excess."""
    rng = np.random.default_rng(seed)
    n = len(epoch_data)
    h_means = np.array([e["h_mean"] for e in epoch_data])
    within_sigmas = np.array([e["within_orient_sigma"] for e in epoch_data])

    excess_samples = np.zeros(n_bootstrap)
    h_trim_samples = np.zeros(n_bootstrap)
    floor_samples = np.zeros(n_bootstrap)

    for b in range(n_bootstrap):
        idx = rng.integers(0, n, size=n)
        h_resampled = h_means[idx]
        sigma_resampled = within_sigmas[idx]

        # 10% trimmed mean
        sorted_h = np.sort(h_resampled)
        n_trim = int(0.1 * n)
        trimmed = sorted_h[n_trim : n - n_trim]
        h_trim_samples[b] = np.mean(trimmed)

        # Noise floor
        floor_samples[b] = np.mean(sigma_resampled) * np.sqrt(2.0 / np.pi)

        excess_samples[b] = h_trim_samples[b] - floor_samples[b]

    return excess_samples, h_trim_samples, floor_samples


def run_bootstrap_for_pulsar(pulsar_name, pulsar_file, n_bootstrap=10000):
    """Run the full bootstrap analysis for one pulsar."""
    print(f"\n{'='*60}")
    print(f"  Epoch-Block Bootstrap: {pulsar_name}")
    print(f"{'='*60}")

    epoch_data = load_epoch_data(pulsar_file)
    n_epochs = len(epoch_data)
    print(f"  Epochs loaded: {n_epochs}")

    if n_epochs < 10:
        print(f"  SKIP: insufficient epochs ({n_epochs})")
        return None

    # Observed values
    obs = compute_trim_excess(epoch_data)
    print(f"  Observed H_trim:       {obs['H_trim']:.3f} ns")
    print(f"  Observed MC floor:     {obs['mc_noise_floor']:.3f} ns")
    print(f"  Observed excess:       {obs['excess']:.3f} ns")

    # Bootstrap
    t0 = time.time()
    excess_samples, h_trim_samples, floor_samples = epoch_block_bootstrap(
        epoch_data, n_bootstrap=n_bootstrap
    )
    t1 = time.time()
    print(f"  Bootstrap iterations:  {n_bootstrap} ({t1-t0:.1f}s)")

    # Bootstrap CI for excess
    ci_low = float(np.percentile(excess_samples, 2.5))
    ci_high = float(np.percentile(excess_samples, 97.5))
    excess_mean = float(np.mean(excess_samples))
    excess_std = float(np.std(excess_samples, ddof=1))

    # Empirical p-value: fraction of bootstrap excess <= 0
    p_value = float(np.mean(excess_samples <= 0))
    if p_value == 0:
        p_value_str = f"< {1.0/n_bootstrap:.1e}"
    else:
        p_value_str = f"{p_value:.4e}"

    # Empirical sigma: excess / std
    empirical_sigma = obs["excess"] / excess_std if excess_std > 0 else 0

    print(f"  Bootstrap excess mean: {excess_mean:.3f} ns")
    print(f"  Bootstrap excess std:  {excess_std:.3f} ns")
    print(f"  95% CI:                [{ci_low:.3f}, {ci_high:.3f}] ns")
    print(f"  Empirical p-value:     {p_value_str}")
    print(f"  Empirical significance: {empirical_sigma:.1f}σ")

    return {
        "pulsar": pulsar_name,
        "n_epochs": n_epochs,
        "n_bootstrap": n_bootstrap,
        "observed": obs,
        "bootstrap": {
            "excess_mean_ns": excess_mean,
            "excess_std_ns": excess_std,
            "ci_95_low_ns": ci_low,
            "ci_95_high_ns": ci_high,
            "p_value": p_value,
            "p_value_str": p_value_str,
            "empirical_sigma": empirical_sigma,
            "h_trim_mean_ns": float(np.mean(h_trim_samples)),
            "h_trim_std_ns": float(np.std(h_trim_samples, ddof=1)),
            "floor_mean_ns": float(np.mean(floor_samples)),
            "floor_std_ns": float(np.std(floor_samples, ddof=1)),
        },
    }


def main():
    print("=" * 60)
    print("  Step 052: Epoch-Block Bootstrap for H_trim Excess")
    print("=" * 60)

    results = {}

    # J0437-4715
    j0437 = run_bootstrap_for_pulsar(
        "J0437-4715",
        RESULTS_DIR / "step_003_closure_final_per_epoch_j0437.json",
        n_bootstrap=10000,
    )
    if j0437:
        results["J0437-4715"] = j0437

    # J1603-7202
    j1603 = run_bootstrap_for_pulsar(
        "J1603-7202",
        RESULTS_DIR / "step_003_closure_final_per_epoch_j1603.json",
        n_bootstrap=10000,
    )
    if j1603:
        results["J1603-7202"] = j1603

    # Save
    output_file = RESULTS_DIR / "step_052_epoch_block_bootstrap.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Results saved to: {output_file}")

    # Summary
    print(f"\n{'='*60}")
    print("  SUMMARY")
    print(f"{'='*60}")
    for name, r in results.items():
        b = r["bootstrap"]
        print(f"  {name}:")
        print(f"    Excess: {r['observed']['excess']:.3f} ns")
        print(f"    95% CI: [{b['ci_95_low_ns']:.3f}, {b['ci_95_high_ns']:.3f}] ns")
        print(f"    Empirical: {b['empirical_sigma']:.1f}σ (p = {b['p_value_str']})")


if __name__ == "__main__":
    main()
