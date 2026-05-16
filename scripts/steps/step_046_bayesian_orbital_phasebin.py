#!/usr/bin/env python3
"""
Step 046: Orbital Phase-Binning Analysis with Signed Delays

Performs a phase-resolved orbital modulation test using the Stokes-aligned
signed geometric delay (geometric_delta_us) rather than the sign-marginalised
magnitude |H|. TEP predicts that the signed delay varies sinusoidally with
orbital phase as the total velocity vector (proper motion + orbital) changes.

This replaces the previous sign-marginalised Bayesian approach, which fitted
a sinusoid to always-positive |H| values and was therefore methodologically
unsound for testing orbital modulation.

Author: TEP Analysis Pipeline
Date: May 2026
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
from scipy import stats

# Setup paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.json_numpy import NpEncoder
from scripts.utils.logger import print_status
from scripts.utils.config import J0437_PB_DAYS, J0437_T0_MJD

RESULTS_DIR = PROJECT_ROOT / "results"


def compute_orbital_phase(mjd: float) -> float:
    """Compute orbital phase for J0437-4715 (0-1 range)."""
    phase = ((mjd - J0437_T0_MJD) % J0437_PB_DAYS) / J0437_PB_DAYS
    return phase


def load_closure_per_epoch(pulsar_name):
    """Load per-epoch closure delay results."""
    summary_file = RESULTS_DIR / f"step_003_closure_final_per_epoch_{pulsar_name}.json"
    if not summary_file.exists():
        print_status(f"File not found: {summary_file}", "ERROR")
        return None
    with open(summary_file, "r") as f:
        return json.load(f)


def inverse_variance_weighted_mean(
    delays: np.ndarray,
    uncertainties: np.ndarray,
) -> Dict[str, float]:
    """
    Compute inverse-variance weighted mean of signed delays.

    Returns the weighted mean, its standard error, and the effective sample size.
    """
    if len(delays) == 0:
        return {
            "mean_ns": np.nan,
            "sem_ns": np.nan,
            "significance": np.nan,
            "n_eff": 0,
            "n_triplets": 0,
        }

    n = len(delays)
    sigma_safe = np.maximum(uncertainties, 1e-12)
    weights = 1.0 / (sigma_safe**2)
    weights = weights / np.sum(weights)

    weighted_mean = float(np.sum(weights * delays))
    sem = float(1.0 / np.sqrt(np.sum(1.0 / (sigma_safe**2))))
    significance = abs(weighted_mean) / sem if sem > 0 else 0.0
    n_eff = float(np.sum(weights) ** 2 / np.sum(weights**2))

    return {
        "mean_ns": weighted_mean,
        "sem_ns": sem,
        "significance": significance,
        "n_eff": n_eff,
        "n_triplets": n,
    }


def fit_sinusoid_to_bins(bin_data: List[Dict]) -> Dict[str, Any]:
    """Fit A*sin(2pi phase + phi)+offset to phase-bin means."""
    phases = np.array([b["phase"] for b in bin_data])
    means = np.array([b["mean_ns"] for b in bin_data])
    sems = np.maximum(np.array([b["sem_ns"] for b in bin_data]), 1e-12)

    def sinusoid(phase, A, phi, offset):
        return A * np.sin(2 * np.pi * phase + phi) + offset

    from scipy.optimize import curve_fit

    popt, pcov = curve_fit(
        sinusoid,
        phases,
        means,
        sigma=sems,
        p0=[np.std(means), 0, np.mean(means)],
        absolute_sigma=True,
        maxfev=10000,
    )
    A_fit, phi_fit, offset_fit = popt
    A_err = float(np.sqrt(max(pcov[0, 0], 0.0)))
    modulation_significance = abs(A_fit) / A_err if A_err > 0 else 0.0

    chi2_null = float(np.sum(((means - np.mean(means)) / sems) ** 2))
    chi2_fit = float(np.sum(((means - sinusoid(phases, *popt)) / sems) ** 2))
    delta_chi2 = max(0.0, chi2_null - chi2_fit)
    # Constant offset is nested inside the sinusoid+offset model. The modulation
    # adds two parameters (sine and cosine, equivalently amplitude and phase).
    p_val_modulation = float(stats.chi2.sf(delta_chi2, 2))

    approaching_mask = (phases >= 0) & (phases < 0.5)
    receding_mask = (phases >= 0.5) & (phases <= 1)
    if np.any(approaching_mask) and np.any(receding_mask):
        mean_approaching = float(np.mean(means[approaching_mask]))
        mean_receding = float(np.mean(means[receding_mask]))
        diff = mean_approaching - mean_receding
    else:
        mean_approaching = mean_receding = diff = np.nan

    return {
        "fitted_amplitude_ns": float(A_fit),
        "fitted_amplitude_err_ns": A_err,
        "fitted_phase_rad": float(phi_fit),
        "fitted_offset_ns": float(offset_fit),
        "modulation_significance_sigma": float(modulation_significance),
        "modulation_p_value": p_val_modulation,
        "delta_chi2": float(delta_chi2),
        "delta_chi2_dof": 2,
        "quadrant_difference_ns": float(diff),
        "mean_approaching_ns": float(mean_approaching),
        "mean_receding_ns": float(mean_receding),
        "n_bins_used": len(bin_data),
        "chi2_null": chi2_null,
        "chi2_fit": chi2_fit,
        "dof_fit": int(len(means) - 3),
    }


def aggregate_epoch_signed_delays(triplets: List[Dict]) -> List[Dict[str, Any]]:
    """Aggregate triplet delays to one independent signed-delay mean per epoch."""
    grouped: Dict[str, List[Dict]] = {}
    for t in triplets:
        grouped.setdefault(t["epoch_id"], []).append(t)

    epoch_rows = []
    for epoch_id, rows in grouped.items():
        delays = np.array([r["delay_ns"] for r in rows], dtype=float)
        uncertainties = np.array([r["uncertainty_ns"] for r in rows], dtype=float)
        if len(delays) == 0:
            continue
        sigma_safe = np.maximum(uncertainties, 1e-12)
        weights = 1.0 / (sigma_safe**2)
        weights = weights / np.sum(weights)
        mean_ns = float(np.sum(weights * delays))
        if len(delays) > 1:
            internal_sem = float(np.std(delays, ddof=1) / np.sqrt(len(delays)))
        else:
            internal_sem = float(sigma_safe[0])
        epoch_rows.append(
            {
                "epoch_id": epoch_id,
                "mjd": float(rows[0]["mjd"]),
                "orbital_phase": float(rows[0]["orbital_phase"]),
                "mean_ns": mean_ns,
                "internal_sem_ns": internal_sem,
                "n_triplets": len(delays),
            }
        )
    return epoch_rows


def epoch_block_phase_bin_analysis(
    triplets: List[Dict],
    n_phase_bins: int = 8,
    min_epochs_per_bin: int = 20,
    n_permutations: int = 2000,
    seed: int = 46,
) -> Dict[str, Any]:
    """Conservative orbital test using one independent mean per observing epoch."""
    epoch_rows = aggregate_epoch_signed_delays(triplets)
    phases = np.array([r["orbital_phase"] for r in epoch_rows])
    means = np.array([r["mean_ns"] for r in epoch_rows])
    n_triplets = np.array([r["n_triplets"] for r in epoch_rows])

    bin_edges = np.linspace(0, 1, n_phase_bins + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    phase_bins = []
    fit_bins = []

    for i in range(n_phase_bins):
        phase_min = bin_edges[i]
        phase_max = bin_edges[i + 1]
        mask = (phases >= phase_min) & (phases < phase_max)
        bin_epoch_means = means[mask]
        bin_result = {
            "bin_index": i,
            "phase_min": float(phase_min),
            "phase_max": float(phase_max),
            "phase_center": float(bin_centers[i]),
            "phase_degrees": float(bin_centers[i] * 360),
            "n_epochs": int(np.sum(mask)),
            "n_triplets": int(np.sum(n_triplets[mask])) if np.any(mask) else 0,
            "valid": False,
        }
        if len(bin_epoch_means) >= min_epochs_per_bin:
            sem = float(np.std(bin_epoch_means, ddof=1) / np.sqrt(len(bin_epoch_means)))
            mean_ns = float(np.mean(bin_epoch_means))
            bin_result.update(
                {
                    "valid": True,
                    "mean_ns": mean_ns,
                    "sem_ns": sem,
                    "epoch_scatter_ns": float(np.std(bin_epoch_means, ddof=1)),
                    "significance": float(abs(mean_ns) / sem) if sem > 0 else 0.0,
                }
            )
            fit_bins.append(
                {
                    "phase": float(bin_centers[i]),
                    "mean_ns": mean_ns,
                    "sem_ns": sem,
                    "significance": bin_result["significance"],
                }
            )
        phase_bins.append(bin_result)

    results = {
        "method": "epoch_blocked_signed_delay_phase_bins",
        "independence_unit": "observing_epoch",
        "n_epochs_total": len(epoch_rows),
        "n_bins_total": n_phase_bins,
        "n_bins_valid": len(fit_bins),
        "min_epochs_per_bin": min_epochs_per_bin,
        "phase_bins": phase_bins,
    }

    if len(fit_bins) >= 4:
        modulation = fit_sinusoid_to_bins(fit_bins)
        results["modulation_test"] = modulation

        rng = np.random.RandomState(seed)
        observed_amp = abs(modulation["fitted_amplitude_ns"])
        perm_amps = []
        for _ in range(n_permutations):
            shuffled = phases.copy()
            rng.shuffle(shuffled)
            perm_triplets = []
            for row, phase in zip(epoch_rows, shuffled):
                perm_triplets.append({**row, "orbital_phase": float(phase)})
            # Re-bin directly from permuted epoch rows to avoid re-expanding triplets.
            perm_bins = []
            perm_phases = np.array([r["orbital_phase"] for r in perm_triplets])
            perm_means = np.array([r["mean_ns"] for r in perm_triplets])
            for i in range(n_phase_bins):
                mask = (perm_phases >= bin_edges[i]) & (perm_phases < bin_edges[i + 1])
                vals = perm_means[mask]
                if len(vals) >= min_epochs_per_bin:
                    sem = float(np.std(vals, ddof=1) / np.sqrt(len(vals)))
                    perm_bins.append(
                        {
                            "phase": float(bin_centers[i]),
                            "mean_ns": float(np.mean(vals)),
                            "sem_ns": sem,
                        }
                    )
            if len(perm_bins) >= 4:
                try:
                    perm_amps.append(abs(fit_sinusoid_to_bins(perm_bins)["fitted_amplitude_ns"]))
                except Exception:
                    continue
        perm_amps = np.array(perm_amps)
        if len(perm_amps) > 0:
            results["epoch_phase_permutation"] = {
                "n_permutations": int(len(perm_amps)),
                "observed_abs_amplitude_ns": float(observed_amp),
                "null_abs_amplitude_mean_ns": float(np.mean(perm_amps)),
                "null_abs_amplitude_std_ns": float(np.std(perm_amps, ddof=1)),
                "empirical_p_value": float((np.sum(perm_amps >= observed_amp) + 1) / (len(perm_amps) + 1)),
            }

    return results


def load_triplet_data_with_orbital_phase(pulsar: str = "j0437") -> List[Dict]:
    """
    Load triplet-level closure delay data with orbital phases.

    Returns list of dicts with:
    - mjd, orbital_phase, delay_ns, uncertainty_ns, closure_snr

    Uses the Stokes-aligned signed geometric delay (geometric_delta_us)
    and its uncertainty (sigma_us) directly.
    """
    epochs = load_closure_per_epoch(pulsar)
    if epochs is None:
        return []
    triplets_with_phase = []

    for epoch in epochs:
        mjd = epoch.get("mjd", 0)
        orbital_phase = compute_orbital_phase(mjd)

        triplets = epoch.get("triplets", [])
        for triplet in triplets:
            # Use Stokes-aligned signed geometric delay directly
            delay_us = triplet.get("geometric_delta_us")
            uncertainty_us = triplet.get("sigma_us")
            if delay_us is None or uncertainty_us is None:
                continue

            # Convert from microseconds to nanoseconds
            delay_ns = float(delay_us) * 1000.0
            uncertainty_ns = float(uncertainty_us) * 1000.0

            triplets_with_phase.append(
                {
                    "mjd": mjd,
                    "orbital_phase": orbital_phase,
                    "delay_ns": delay_ns,
                    "uncertainty_ns": uncertainty_ns,
                    "epoch_id": epoch.get("epoch", "unknown"),
                }
            )

    print_status(
        f"Loaded {len(triplets_with_phase)} triplets with orbital phase info", "INFO"
    )
    return triplets_with_phase


def phase_bin_signed_analysis(
    triplets: List[Dict], n_phase_bins: int = 8, min_triplets_per_bin: int = 50
) -> Dict[str, Any]:
    """
    Perform signed-delay analysis on orbital phase-binned triplets.

    Strategy:
    1. Bin all triplets by orbital phase
    2. Compute inverse-variance weighted mean signed delay per bin
    3. Fit sinusoid to signed bin means with absolute_sigma=True
    4. Report modulation amplitude and significance
    """
    print_status(f"\n{'=' * 70}", "TITLE")
    print_status("SIGNED-DELAY ORBITAL PHASE-BINNING ANALYSIS", "TITLE")
    print_status(f"{'=' * 70}", "TITLE")
    print_status(f"Phase bins: {n_phase_bins}", "INFO")
    print_status(f"Min triplets per bin: {min_triplets_per_bin}", "INFO")

    # Extract arrays
    phases = np.array([t["orbital_phase"] for t in triplets])
    delays = np.array([t["delay_ns"] for t in triplets])
    uncertainties = np.array([t["uncertainty_ns"] for t in triplets])

    # Create phase bins (0 to 1)
    bin_edges = np.linspace(0, 1, n_phase_bins + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    results = {
        "phase_bins": [],
        "n_bins_total": n_phase_bins,
        "n_bins_valid": 0,
        "binary_period_days": J0437_PB_DAYS,
    }

    bin_means = []

    for i in range(n_phase_bins):
        phase_min = bin_edges[i]
        phase_max = bin_edges[i + 1]
        phase_center = bin_centers[i]

        # Select triplets in this phase bin
        mask = (phases >= phase_min) & (phases < phase_max)
        n_in_bin = np.sum(mask)

        bin_result = {
            "bin_index": i,
            "phase_min": float(phase_min),
            "phase_max": float(phase_max),
            "phase_center": float(phase_center),
            "phase_degrees": float(phase_center * 360),
            "n_triplets": int(n_in_bin),
            "valid": False,
        }

        if n_in_bin >= min_triplets_per_bin:
            bin_delays = delays[mask]
            bin_uncertainties = uncertainties[mask]

            ivw_result = inverse_variance_weighted_mean(
                bin_delays, bin_uncertainties
            )

            bin_result.update(
                {
                    "valid": True,
                    "mean_ns": ivw_result["mean_ns"],
                    "sem_ns": ivw_result["sem_ns"],
                    "significance": ivw_result["significance"],
                    "n_eff": ivw_result["n_eff"],
                }
            )

            bin_means.append(
                {
                    "phase": phase_center,
                    "mean_ns": ivw_result["mean_ns"],
                    "sem_ns": ivw_result["sem_ns"],
                    "significance": ivw_result["significance"],
                }
            )

            results["n_bins_valid"] += 1

            print_status(
                f"\n  Bin {i + 1}: phase {phase_min:.2f}-{phase_max:.2f} "
                f"({phase_center * 360:.0f}°)",
                "INFO",
            )
            print_status(
                f"    Triplets: {n_in_bin} (n_eff={ivw_result['n_eff']:.1f})", "INFO"
            )
            print_status(
                f"    Signed mean = {ivw_result['mean_ns']:+.3f} +/- {ivw_result['sem_ns']:.3f} ns",
                "INFO",
            )
            print_status(
                f"    Significance: {ivw_result['significance']:.1f}sigma", "INFO"
            )
        else:
            print_status(
                f"\n  Bin {i + 1}: phase {phase_min:.2f}-{phase_max:.2f} "
                f"- INSUFFICIENT DATA ({n_in_bin} triplets)",
                "WARNING",
            )

        results["phase_bins"].append(bin_result)

    # Test for orbital modulation pattern
    if len(bin_means) >= 4:
        modulation_test = test_orbital_modulation_pattern(bin_means)
        results["modulation_test"] = modulation_test

    return results


def test_orbital_modulation_pattern(bin_data: List[Dict]) -> Dict[str, Any]:
    """
    Fit sinusoidal modulation to signed mean delays per orbital phase bin.

    TEP predicts the signed delay varies as the total velocity vector
    (proper motion + orbital) changes with orbital phase.

    Uses absolute_sigma=True to ensure parameter uncertainties reflect
    the actual bin SEMs rather than being rescaled by chi-square.
    """
    phases = np.array([b["phase"] for b in bin_data])
    means = np.array([b["mean_ns"] for b in bin_data])
    sems = np.array([b["sem_ns"] for b in bin_data])

    try:
        result = fit_sinusoid_to_bins(bin_data)

        # Quadrant analysis on signed means
        approaching_mask = ((phases >= 0) & (phases < 0.5))
        receding_mask = ((phases >= 0.5) & (phases <= 1))

        if np.any(approaching_mask) and np.any(receding_mask):
            mean_approaching = float(np.mean(means[approaching_mask]))
            mean_receding = float(np.mean(means[receding_mask]))
            diff = mean_approaching - mean_receding
        else:
            mean_approaching = mean_receding = diff = np.nan

        result.update(
            {
                "quadrant_difference_ns": float(diff),
                "mean_approaching_ns": float(mean_approaching),
                "mean_receding_ns": float(mean_receding),
            }
        )
        return result
    except Exception as e:
        print_status(f"Modulation fit failed: {e}", "WARNING")
        return {"error": str(e), "n_bins_used": len(bin_data)}


def main():
    """Run signed-delay orbital phase-binned analysis."""
    print_status("\n" + "=" * 70, "TITLE")
    print_status("STEP 046: SIGNED-DELAY ORBITAL PHASE-BINNING", "TITLE")
    print_status("=" * 70, "TITLE")

    triplets = load_triplet_data_with_orbital_phase("j0437")

    if len(triplets) == 0:
        print_status("No triplet data found!", "ERROR")
        return

    results = phase_bin_signed_analysis(
        triplets, n_phase_bins=8, min_triplets_per_bin=100
    )
    epoch_blocked_results = epoch_block_phase_bin_analysis(
        triplets, n_phase_bins=8, min_epochs_per_bin=20
    )
    results["epoch_blocked_analysis"] = epoch_blocked_results
    results["recommended_modulation_test"] = "epoch_blocked_analysis.modulation_test"
    results["methodology_note"] = (
        "Triplet-binned analysis is retained as a high-precision diagnostic, but "
        "triplets within an epoch are correlated. The conservative orbital claim "
        "should use epoch_blocked_analysis, which contributes one signed-delay "
        "mean per observing epoch and estimates bin SEMs from epoch-to-epoch scatter."
    )

    # Summary
    print_status(f"\n{'=' * 70}", "TITLE")
    print_status("SUMMARY", "TITLE")
    print_status(f"{'=' * 70}", "TITLE")
    print_status(f"Total phase bins: {results['n_bins_total']}", "INFO")
    print_status(f"Valid bins: {results['n_bins_valid']}", "INFO")

    if "modulation_test" in results:
        mt = results["modulation_test"]
        print_status(f"\nModulation Test:", "INFO")
        print_status(
            f"  Fitted amplitude: {mt.get('fitted_amplitude_ns', 0):+.3f} +/- "
            f"{mt.get('fitted_amplitude_err_ns', 0):.3f} ns",
            "INFO",
        )
        print_status(
            f"  Modulation significance: {mt.get('modulation_significance_sigma', 0):.2f} sigma",
            "INFO",
        )
        print_status(
            f"  Modulation p-value: {mt.get('modulation_p_value', 1.0):.3f}",
            "INFO",
        )

        if "quadrant_difference_ns" in mt and not np.isnan(mt["quadrant_difference_ns"]):
            print_status(f"\nQuadrant Analysis:", "INFO")
            print_status(
                f"  Approaching mean: {mt['mean_approaching_ns']:+.3f} ns", "INFO"
            )
            print_status(
                f"  Receding mean: {mt['mean_receding_ns']:+.3f} ns", "INFO"
            )
            print_status(
                f"  Difference: {mt['quadrant_difference_ns']:+.3f} ns", "INFO"
            )

    if "modulation_test" in epoch_blocked_results:
        mt = epoch_blocked_results["modulation_test"]
        print_status(f"\nEpoch-blocked Modulation Test:", "INFO")
        print_status(
            f"  Fitted amplitude: {mt.get('fitted_amplitude_ns', 0):+.3f} +/- "
            f"{mt.get('fitted_amplitude_err_ns', 0):.3f} ns",
            "INFO",
        )
        print_status(
            f"  Modulation significance: {mt.get('modulation_significance_sigma', 0):.2f} sigma",
            "INFO",
        )
        print_status(
            f"  Nested-model p-value: {mt.get('modulation_p_value', 1.0):.3e}",
            "INFO",
        )
        if "epoch_phase_permutation" in epoch_blocked_results:
            print_status(
                f"  Epoch phase permutation p-value: {epoch_blocked_results['epoch_phase_permutation']['empirical_p_value']:.3e}",
                "INFO",
            )

    # Save results
    output_path = RESULTS_DIR / "step_046_bayesian_orbital_phasebin_results.json"
    with open(output_path, "w") as f:
        json.dump(results, f, indent=2, cls=NpEncoder)

    print_status(f"\nResults saved to: {output_path}", "INFO")
    print_status("=" * 70, "TITLE")


if __name__ == "__main__":
    main()
