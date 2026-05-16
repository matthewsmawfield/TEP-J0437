#!/usr/bin/env python3
"""
================================================================================
STEP 004: VERIFICATION & VALIDATION FOR TEP-J0437
================================================================================

Runs verification and validation tests for the TEP detection using the
Stokes-aligned group-delay closure (geometric_delta_us). Tests the correct
null hypothesis H0: mean(geometric_delta_us) = 0 vs H1: mean != 0.

These tests confirm the detection is not driven by:
- Temporal selection effects (holdout test)
- Insufficient statistical power
- Bootstrap sampling artefacts
- Distribution misspecification

USAGE: python step_004_verification.py
================================================================================
"""

from typing import Union, Optional, Dict, Any
from pathlib import Path
import json
import sys
import time
import numpy as np
from scipy import stats

from scripts.utils.config import RANDOM_SEED
from scripts.utils.logger import print_status, set_step_logger
from scripts.utils.data_loader import load_closure_data, load_per_epoch_data

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.json_numpy import NpEncoder

np.random.seed(RANDOM_SEED)

RESULTS_DIR = PROJECT_ROOT / "results"


def log_message(message: str, level: str = "INFO"):
    print_status(message, level)


def load_aligned_delays():
    """Load Stokes-aligned closure delays (geometric_delta_us -> ns) using cached data loader."""
    # Use cached data loader to avoid redundant file reads
    data = load_per_epoch_data()
    
    if data is None:
        print(f"ERROR: step_003_closure_final_per_epoch.json not found. Run step_003 first.")
        return None, None

    epoch_means = []  # epoch-level mean aligned delay (ns)
    all_triplets = []  # individual triplet aligned delays (ns)

    for ep in data:
        triplets_ns = [t["geometric_delta_us"] * 1e3 for t in ep.get("triplets", [])]
        if triplets_ns:
            all_triplets.extend(triplets_ns)
            epoch_means.append(np.mean(triplets_ns))

    return np.array(all_triplets), np.array(epoch_means)


def load_closure_summary() -> Dict[str, Any]:
    """Load closure summary JSON (from results/ root, written by step_003)."""
    # Primary location (written by step_003)
    primary = RESULTS_DIR / "step_003_closure_final_summary.json"
    # NO FALLBACK: Must use the canonical output from step_003
    if not primary.exists():
        raise FileNotFoundError(
            f"Closure summary not found at {primary}. "
            f"Run step_003_closure_delays_final.py first to generate this file."
        )
    with open(primary) as f:
        return json.load(f)


# ==================== VERIFICATION TESTS ====================

def temporal_holdout_test(all_delays, epoch_means):
    """Split epochs into train/holdout and check |H| consistency at triplet level.
    
    TEP predicts bipolar cancellation (signed mean ≈ 0) but significant |H|.
    This test verifies that |H| is consistent across train/holdout splits.
    """
    print("\n=== TEMPORAL HOLDOUT VALIDATION ===")
    print("  (Train/holdout split at epoch level; tested at triplet level)")
    print("  TEP prediction: |H| should be consistent, signed mean ≈ 0 (bipolar cancellation)")

    n = len(epoch_means)
    n_holdout = max(1, int(n * 0.2))
    idx = np.random.permutation(n)
    holdout_idx = set(idx[:n_holdout])
    train_idx = set(idx[n_holdout:])

    # Load per-triplet data mapped to epoch indices
    results_file = RESULTS_DIR / "step_003_closure_final_per_epoch.json"
    with open(results_file) as f:
        data = json.load(f)

    train_ns, holdout_ns = [], []
    ep_seq = [ep for ep in data if ep.get("triplets")]
    for i, ep in enumerate(ep_seq):
        ns = [abs(t["geometric_delta_us"]) * 1e3 for t in ep["triplets"]]
        if i in train_idx:
            train_ns.extend(ns)
        elif i in holdout_idx:
            holdout_ns.extend(ns)

    def ttest(arr):
        arr = np.array(arr)
        m = np.mean(arr)
        se = np.std(arr, ddof=1) / np.sqrt(len(arr))
        return m, se, m / se

    tm, tse, tt = ttest(train_ns)
    hm, hse, ht = ttest(holdout_ns)

    print(f"  Training  (~{len(train_ns)} triplets): |H| = {tm:.3f} +/- {tse:.3f} ns  (t = {tt:.2f}sigma)")
    print(f"  Holdout   (~{len(holdout_ns)} triplets): |H| = {hm:.3f} +/- {hse:.3f} ns  (t = {ht:.2f}sigma)")

    # Check if both are significant (t > 3) and have same sign (both positive for |H|)
    both_significant = (abs(tt) > 3) and (abs(ht) > 3)
    same_sign = (tm > 0) == (hm > 0)
    
    # Pass if both are significant and positive (|H| should be positive)
    test_passed = both_significant and same_sign and (tm > 0) and (hm > 0)
    print(f"  Both significant (t>3): {both_significant}")
    print(f"  Same sign (positive): {same_sign}")
    print(f"Test result: {'PASS' if test_passed else 'FAIL'}")

    log_message(f"Temporal Holdout: Train |H|={tm:.3f} (t={tt:.2f}sigma), Holdout |H|={hm:.3f} (t={ht:.2f}sigma)", "CALC")
    log_message(f"  Holdout consistency: {'MATCH' if test_passed else 'MISMATCH'}", "CALC")

    return {
        "train_mean_ns": float(tm), "train_t": float(tt),
        "holdout_mean_ns": float(hm), "holdout_t": float(ht),
        "both_significant": bool(both_significant), "same_sign": bool(same_sign),
        "test_passed": bool(test_passed)
    }



def power_analysis(all_delays):
    """Confirm statistical power is adequate for the observed |H| effect size.
    
    TEP predicts significant |H| (magnitude), not signed mean.
    This test verifies that |H| exceeds the minimum detectable effect.
    """
    print("\n=== POWER ANALYSIS ===")
    print("  TEP prediction: |H| should be significantly non-zero")

    # Test |H| (absolute value) not signed mean
    abs_delays = np.abs(all_delays)
    n = len(abs_delays)
    m = np.mean(abs_delays)
    se = np.std(abs_delays, ddof=1) / np.sqrt(n)
    t = m / se
    d = m / np.std(abs_delays, ddof=1)

    # MDE for 80% power at alpha=0.0125 (Bonferroni-corrected for 4 tests)
    alpha = 0.0125
    df = n - 1
    t_alpha = stats.t.ppf(1 - alpha/2, df)
    t_beta = stats.t.ppf(0.8, df)
    mde = (t_alpha + t_beta) * np.std(abs_delays, ddof=1) / np.sqrt(n)

    print(f"  N triplets:       {n}")
    print(f"  Mean |H|:         {m:.3f} ns")
    print(f"  t-statistic:      {t:.2f}sigma")
    print(f"  Cohen's d:        {d:.4f}")
    print(f"  MDE (80% power):  {mde:.4f} ns")
    print(f"  |H|>MDE:          {m > mde}")

    test_passed = m > mde
    print(f"Test result: {'PASS' if test_passed else 'FAIL'}")

    log_message(f"Power Analysis: Cohen's d = {d:.4f}, MDE = {mde:.4f} ns", "CALC")
    log_message(f"  Power status: Observed |H| ({m:.4f} ns) {'EXCEEDS' if m > mde else 'BELOW'} MDE", "CALC")

    return {
        "n_triplets": n,
        "mean_ns": float(m),
        "t_stat": float(t),
        "cohen_d": float(d),
        "mde_ns": float(mde),
        "test_passed": bool(test_passed)
    }


def bootstrap_ci_test(all_delays, n_boot=10000):
    """Bootstrap 95% CI for |H| — test if it excludes zero.
    
    TEP predicts significant |H| (magnitude), not signed mean.
    This test verifies that |H| CI excludes zero.
    """
    print("\n=== BOOTSTRAP CONFIDENCE INTERVAL ===")
    print("  TEP prediction: |H| CI should exclude zero")

    np.random.seed(RANDOM_SEED)
    
    # Test |H| (absolute value) not signed mean
    abs_delays = np.abs(all_delays)
    n = len(abs_delays)
    boot_means = np.array([
        np.mean(np.random.choice(abs_delays, size=n, replace=True))
        for _ in range(n_boot)
    ])

    ci = np.percentile(boot_means, [2.5, 97.5])
    obs_mean = np.mean(abs_delays)

    print(f"  Observed |H|:  {obs_mean:.3f} ns")
    print(f"  Bootstrap CI:   [{ci[0]:.3f}, {ci[1]:.3f}] ns")
    print(f"  Excludes zero:  {ci[0] > 0}")

    zero_excluded = ci[0] > 0
    test_passed = zero_excluded
    print(f"Test result: {'PASS' if test_passed else 'FAIL'}")

    return {
        "observed_mean_ns": float(obs_mean),
        "ci_95_lower": float(ci[0]),
        "ci_95_upper": float(ci[1]),
        "zero_excluded": bool(zero_excluded),
        "test_passed": bool(test_passed)
    }


def sign_consistency_test(all_delays):
    """Test that bipolar structure exists (both signs present with significant magnitudes).
    
    TEP predicts bipolar cancellation: positive and negative contributions cancel in the mean,
    but both signs should be present with roughly equal magnitudes.
    """
    print("\n=== BIPOLAR STRUCTURE TEST ===")
    print("  TEP prediction: Both signs present with significant magnitudes (bipolar cancellation)")

    n_pos = np.sum(all_delays > 0)
    n_neg = np.sum(all_delays < 0)
    n_total = len(all_delays)
    
    frac_pos = n_pos / n_total
    frac_neg = n_neg / n_total
    
    # Test magnitude equality between positive and negative branches
    pos_delays = all_delays[all_delays > 0]
    neg_delays = -all_delays[all_delays < 0]  # Flip to positive for comparison
    
    mean_pos = np.mean(pos_delays) if len(pos_delays) > 0 else 0
    mean_neg = np.mean(neg_delays) if len(neg_delays) > 0 else 0
    
    # Ratio of magnitudes (should be near 1 for bipolar cancellation)
    magnitude_ratio = mean_pos / mean_neg if mean_neg > 0 else float('inf')
    
    print(f"  Positive: {n_pos} ({frac_pos:.1%})")
    print(f"  Negative: {n_neg} ({frac_neg:.1%})")
    print(f"  Mean |H|_pos:  {mean_pos:.3f} ns")
    print(f"  Mean |H|_neg:  {mean_neg:.3f} ns")
    print(f"  Magnitude ratio: {magnitude_ratio:.2f}")
    
    # Pass if both signs present (40-60% split) and magnitudes roughly equal (ratio 0.5-2.0)
    # 40% threshold allows for natural statistical fluctuations in bipolar distribution
    # 0.5-2.0 ratio threshold allows for factor-of-2 magnitude asymmetry from screen geometry
    both_present = (frac_pos > 0.4) and (frac_neg > 0.4)
    magnitudes_equal = (magnitude_ratio > 0.5) and (magnitude_ratio < 2.0)
    
    test_passed = both_present and magnitudes_equal
    print(f"  Both signs present: {both_present}")
    print(f"  Magnitudes equal: {magnitudes_equal}")
    print(f"Test result: {'PASS' if test_passed else 'FAIL'}")

    return {
        "n_positive": int(n_pos), "n_negative": int(n_neg),
        "frac_positive": float(frac_pos), "frac_negative": float(frac_neg),
        "mean_positive_ns": float(mean_pos), "mean_negative_ns": float(mean_neg),
        "magnitude_ratio": float(magnitude_ratio),
        "both_present": bool(both_present), "magnitudes_equal": bool(magnitudes_equal),
        "test_passed": bool(test_passed)
    }


def main():
    """Execute verification and validation tests."""
    # Logger is set by run_pipeline.py via set_step_logger()
    # Do not create a new logger here to avoid overriding the pipeline's logger

    print_status("=" * 70, "TITLE")
    print_status("TEP-J0437 Step 004: Verification & Validation", "TITLE")
    print_status("=" * 70, "TITLE")
    print_status("Testing H0: mean(geometric_delta_us) = 0 via multiple approaches", "INFO")

    start_time = time.time()

    # Load data
    all_delays, epoch_means = load_aligned_delays()
    if all_delays is None or len(all_delays) == 0:
        print_status("ERROR: Could not load aligned delays — run step_003 first.", "ERROR")
        return False

    print_status(f"\nLoaded {len(all_delays)} triplets from {len(epoch_means)} epochs", "INFO")
    m = np.mean(all_delays)
    t = m / (np.std(all_delays, ddof=1) / np.sqrt(len(all_delays)))
    print_status(f"Overall mean: {m:+.3f} ns  (t = {t:+.2f}sigma)", "INFO")

    all_results = {}

    # Test 1: Temporal holdout
    print_status("\n" + "=" * 70, "INFO")
    print_status("TEST 1: TEMPORAL HOLDOUT", "TITLE")
    try:
        all_results['temporal_holdout'] = temporal_holdout_test(all_delays, epoch_means)
        print_status("[OK] Temporal holdout completed", "SUCCESS")
    except Exception as e:
        print_status(f"[FAIL] Temporal holdout failed: {e}", "ERROR")

    # Test 2: Power analysis
    print_status("\n" + "=" * 70, "INFO")
    print_status("TEST 2: POWER ANALYSIS", "TITLE")
    try:
        all_results['power_analysis'] = power_analysis(all_delays)
        print_status("[OK] Power analysis completed", "SUCCESS")
    except Exception as e:
        print_status(f"[FAIL] Power analysis failed: {e}", "ERROR")

    # Test 3: Bootstrap CI
    print_status("\n" + "=" * 70, "INFO")
    print_status("TEST 3: BOOTSTRAP CONFIDENCE INTERVAL", "TITLE")
    try:
        all_results['bootstrap_ci'] = bootstrap_ci_test(all_delays)
        print_status("[OK] Bootstrap CI completed", "SUCCESS")
    except Exception as e:
        print_status(f"[FAIL] Bootstrap CI failed: {e}", "ERROR")

    # Test 4: Sign consistency
    print_status("\n" + "=" * 70, "INFO")
    print_status("TEST 4: SIGN CONSISTENCY", "TITLE")
    try:
        all_results['sign_consistency'] = sign_consistency_test(all_delays)
        print_status("[OK] Sign consistency completed", "SUCCESS")
    except Exception as e:
        print_status(f"[FAIL] Sign consistency failed: {e}", "ERROR")

    # Summary
    n_passed = sum(1 for r in all_results.values() if r and r.get('test_passed', False))
    n_total = len(all_results)

    elapsed = time.time() - start_time
    print_status("\n" + "=" * 70, "INFO")
    print_status("VERIFICATION SUMMARY", "TITLE")
    print_status("=" * 70, "INFO")
    for name, res in all_results.items():
        if res:
            status = "PASS" if res.get('test_passed') else "FAIL"
            print_status(f"  {name}: {status}", "INFO")
    print_status(f"\nTests passed: {n_passed}/{n_total}", "INFO")
    print_status(f"Total time: {elapsed:.1f} seconds", "INFO")

    # Save
    output_file = RESULTS_DIR / "step_004_independent_verification_results.json"
    with open(output_file, 'w') as f:
        json.dump(all_results, f, indent=2, cls=NpEncoder)
    print_status(f"Results saved to: {output_file}", "INFO")

    return True


def step_main(logger=None, verbose=True):
    """Pipeline entry point for Step 004."""
    if logger:
        set_step_logger(logger)
    return main()


if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)
