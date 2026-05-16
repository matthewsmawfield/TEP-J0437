#!/usr/bin/env python3
"""
================================================================================
STEP 013: FALSIFICATION CRITERIA AUTOMATION
================================================================================
Evaluates falsification criteria for the TEP detection using Stokes-aligned
closure delays. All criteria test the signed aligned mean, not |H|.

Criteria:
  1. Phase Closure significance (primary detection metric): Rayleigh/V-test p < 1e-6
     and epoch-level circular bootstrap CI excludes zero (Step 003 summary).
  2. Phase-scramble specificity: circular concentration exceeds scramble null (Step 047).
  3. Frame invariance of unweighted Phase Closure ψ: identical under bulk-vector substitutions (Step 048).
  4. Velocity-label permutation control: reported as inconclusive at N=2 unless unique (Step 047; diagnostic).
  5. Bipolar cancellation in signed mean: |t_signed| < 3 (diagnostic; Step 003 summary).

Important: |H| magnitude diagnostics are noise-floor dominated and are not used as
primary falsification gates in this repository’s current inference logic.
================================================================================
"""

from typing import Union, Optional

import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.json_numpy import NpEncoder
from scripts.utils.config import RANDOM_SEED
from scripts.utils.logger import TEPLogger, print_status, set_step_logger

np.random.seed(RANDOM_SEED)

RAYLEIGH_P_PRIMARY = 1e-6
VTEST_P_PRIMARY = 1e-6
BOOTSTRAP_NOMINAL_CONFIDENCE = 0.95

# Logger is set by run_pipeline.py via set_step_logger()
# Do not create a new logger here to avoid overriding the pipeline's logger
RESULTS_DIR = PROJECT_ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def _load_json(path: Path) -> Optional[dict]:
    if not path.exists():
        return None
    with open(path) as fh:
        return json.load(fh)


def load_step_003_summary_j0437() -> Optional[dict]:
    return _load_json(RESULTS_DIR / "step_003_closure_final_summary_j0437.json")


def load_step_047_controls() -> Optional[dict]:
    return _load_json(RESULTS_DIR / "step_047_velocity_direction_controls_results.json")


def load_step_048_frame_analysis() -> Optional[dict]:
    return _load_json(RESULTS_DIR / "step_048_cmb_dipole_frame_analysis.json")


def evaluate_criteria():
    """Evaluate TEP falsification criteria against current data."""
    s3 = load_step_003_summary_j0437()
    s47 = load_step_047_controls()
    s48 = load_step_048_frame_analysis()

    if s3 is None:
        print_status("ERROR: step_003_closure_final_summary_j0437.json not found.", "ERROR")
        return None

    # Primary: circular-statistics significance + bootstrap CI excludes zero
    rayleigh_p = float(
        s3.get("phase_closure_rayleigh_p_unweighted", s3.get("phase_closure_rayleigh_p"))
    )
    v_p = float(s3.get("phase_closure_v_p_unweighted", s3.get("phase_closure_v_p")))
    ci_low = s3.get("phase_closure_bootstrap_ci_95_lower_rad")
    ci_high = s3.get("phase_closure_bootstrap_ci_95_upper_rad")
    ci_excludes_zero = bool(ci_low is not None and ci_high is not None and not (ci_low <= 0.0 <= ci_high))
    phase_primary_pass = bool((rayleigh_p < RAYLEIGH_P_PRIMARY or v_p < VTEST_P_PRIMARY) and ci_excludes_zero)

    # Step 047 phase-scramble specificity
    phase_scramble_p = None
    phase_scramble_pass = None
    velocity_label_status = None
    velocity_label_pass = None
    wrong_velocity_status = None
    if s47 is not None:
        phase_scramble_p = (
            s47.get("results", {})
            .get("phase_scramble", {})
            .get("p_observed_concentration_by_chance")
        )
        phase_scramble_pass = bool(
            s47.get("results", {}).get("phase_scramble", {}).get("test_passed", False)
        )
        velocity_label_status = s47.get("results", {}).get("velocity_label_permutation", {}).get("test_status")
        velocity_label_pass = s47.get("results", {}).get("velocity_label_permutation", {}).get("test_passed")
        wrong_velocity_status = s47.get("results", {}).get("wrong_velocity", {}).get("test_status")

    # Step 048 unweighted invariance
    unweighted_invariant = None
    unweighted_delta = None
    if s48 is not None:
        d = (
            s48.get("pulsars", {})
            .get("J0437-4715", {})
            .get("enhanced_tests", {})
            .get("directional_specificity", {})
        )
        unweighted_invariant = d.get("unweighted_psi_invariant_all_directions")
        unweighted_delta = d.get("unweighted_psi_max_abs_delta_rad")

    # Diagnostic: signed-mean cancellation from step_003 (already epoch-level)
    signed_t = float(s3.get("H_signed_t_statistic"))
    signed_cancel_pass = bool(abs(signed_t) < 3.0)

    # Diagnostics: store |H| quantities but do not gate falsification on them
    H_mag = float(s3.get("H_magnitude_ns"))
    H_sem = float(s3.get("H_sem_ns"))
    H_noise = float(s3.get("H_noise_bias_ns"))
    H_excess = float(s3.get("H_excess_ns"))
    H_excess_t = float(s3.get("H_excess_t_statistic"))
    H_trim = float(s3.get("H_trim_magnitude_ns"))
    H_trim_sem = float(s3.get("H_trim_sem_ns"))

    criteria = {
        "criterion_1": {
            "name": "Phase Closure significance (primary)",
            "threshold": f"(Rayleigh p < {RAYLEIGH_P_PRIMARY} OR V-test p < {VTEST_P_PRIMARY}) AND bootstrap CI excludes 0",
            "observed": {
                "phase_closure_mean_rad": s3.get("phase_closure_mean_rad"),
                "phase_closure_circ_se_rad": s3.get("phase_closure_circ_se_rad"),
                "rayleigh_p": rayleigh_p,
                "v_p": v_p,
                "bootstrap_ci_95": [ci_low, ci_high],
            },
            "passed": phase_primary_pass,
        },
        "criterion_2": {
            "name": "Phase-scramble specificity (diagnostic)",
            "threshold": "Step 047 phase_scramble p < 0.05",
            "observed": {
                "p_observed_concentration_by_chance": phase_scramble_p,
            },
            "passed": bool(phase_scramble_pass) if phase_scramble_pass is not None else False,
            "status": "missing_step_047" if s47 is None else "ok",
        },
        "criterion_3": {
            "name": "Unweighted ψ frame invariance (diagnostic)",
            "threshold": "Step 048 unweighted ψ invariant under bulk-vector substitutions",
            "observed": {
                "unweighted_psi_invariant_all_directions": unweighted_invariant,
                "unweighted_psi_max_abs_delta_rad": unweighted_delta,
            },
            "passed": bool(unweighted_invariant) if unweighted_invariant is not None else False,
            "status": "missing_step_048" if s48 is None else "ok",
        },
        "criterion_4": {
            "name": "Velocity-label permutation control (diagnostic)",
            "threshold": "Step 047 should report pass only when unique-best; N=2 typically inconclusive",
            "observed": {
                "test_status": velocity_label_status,
                "test_passed": velocity_label_pass,
            },
            "passed": None if velocity_label_status in (None, "inconclusive") else bool(velocity_label_pass),
            "status": "missing_step_047" if s47 is None else "ok",
        },
        "criterion_5": {
            "name": "Signed-mean bipolar cancellation (diagnostic)",
            "threshold": "|t_signed| < 3",
            "observed": {
                "H_signed_mean_ns": s3.get("H_signed_mean_ns"),
                "H_signed_sem_ns": s3.get("H_signed_sem_ns"),
                "H_signed_t_statistic": signed_t,
            },
            "passed": signed_cancel_pass,
        },
    }

    # Overall falsification verdict: requires criterion_1 only.
    # Other criteria are diagnostics and may be missing depending on run subset.
    not_falsified = bool(criteria["criterion_1"]["passed"])
    n_passed = sum(1 for v in criteria.values() if v.get("passed") is True)
    n_failed = sum(1 for v in criteria.values() if v.get("passed") is False)
    n_inconclusive = sum(1 for v in criteria.values() if v.get("passed") is None)

    print_status("=" * 70, "TITLE")
    print_status("FALSIFICATION CRITERIA EVALUATION (PHASE-PRIMARY)", "TITLE")
    print_status("=" * 70, "TITLE")
    for k, v in criteria.items():
        passed = v.get("passed")
        if passed is True:
            label = "PASS"
            level = "SUCCESS"
        elif passed is False:
            label = "FAIL"
            level = "WARNING"
        else:
            label = "INCONCLUSIVE"
            level = "WARNING"
        print_status(f"  {k}: {v['name']} -> {label}", level)
    print_status(f"\nPrimary criterion: {'PASS' if not_falsified else 'FAIL'}", "SUCCESS" if not_falsified else "WARNING")

    return {
        "evaluation": criteria,
        "n_passed": n_passed,
        "tep_not_falsified": bool(not_falsified),
        "summary": {
            "criteria_evaluated": 5,
            "criteria_passed": n_passed,
            "criteria_failed": n_failed,
            "criteria_inconclusive": n_inconclusive,
            "overall_status": "TEP NOT FALSIFIED" if not_falsified else "TEP FALSIFIED",
        },
        "statistics": {
            "phase_primary": {
                "phase_closure_mean_rad": s3.get("phase_closure_mean_rad"),
                "phase_closure_circ_se_rad": s3.get("phase_closure_circ_se_rad"),
                "phase_closure_rbar": s3.get("phase_closure_rbar"),
                "phase_closure_rayleigh_p": rayleigh_p,
                "phase_closure_v_p": v_p,
                "phase_closure_p_values_prefer_unweighted": True,
                "phase_closure_bootstrap_ci_95": [ci_low, ci_high],
            },
            "H_diagnostics": {
                "H_magnitude_ns": H_mag,
                "H_sem_ns": H_sem,
                "H_noise_bias_ns": H_noise,
                "H_excess_ns": H_excess,
                "H_excess_t_statistic": H_excess_t,
                "H_trim_magnitude_ns": H_trim,
                "H_trim_sem_ns": H_trim_sem,
            },
            "signed_diagnostic": {
                "H_signed_mean_ns": s3.get("H_signed_mean_ns"),
                "H_signed_sem_ns": s3.get("H_signed_sem_ns"),
                "H_signed_t_statistic": signed_t,
                "signed_threshold": 3.0,
            },
            "note_criterion_design": (
                "Criteria designed per current repo inference: Phase Closure ψ (circular statistics) is the primary "
                "falsification gate; other controls are diagnostic. |H| magnitude diagnostics are reported but are "
                "not used as primary detection gates because they are noise-floor dominated."
            ),
            "linked_outputs": {
                "step_003_summary_j0437": "results/step_003_closure_final_summary_j0437.json",
                "step_047_velocity_controls": "results/step_047_velocity_direction_controls_results.json",
                "step_048_frame_analysis": "results/step_048_cmb_dipole_frame_analysis.json",
            },
        },
    }


def step_main(logger=None, verbose=True):
    """Standard pipeline entry point for falsification criteria."""
    return main()


def main():
    from scripts.utils.logger import _active_logger

    if _active_logger is None:
        log_dir = PROJECT_ROOT / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        _logger = TEPLogger(
            "step_013",
            str(log_dir / "step_013_falsification_criteria.log"),
        )
        set_step_logger(_logger)

    print_status("=" * 70, "INFO")
    print_status("STEP 013: FALSIFICATION CRITERIA AUTOMATION", "INFO")
    print_status("=" * 70, "INFO")

    eval_results = evaluate_criteria()
    if eval_results is None:
        return False

    output_file = RESULTS_DIR / "step_013_falsification_criteria_results.json"
    with open(output_file, "w") as f:
        json.dump(eval_results, f, indent=2, cls=NpEncoder)
    print_status(f"\nResults saved to: {output_file}", "INFO")

    report_lines = [
        "=" * 70,
        "TEP FALSIFICATION CRITERIA REPORT",
        "=" * 70,
        "",
        f"Criteria passed: {eval_results['summary']['criteria_passed']}/5",
        f"Criteria failed: {eval_results['summary']['criteria_failed']}/5",
        f"Criteria inconclusive: {eval_results['summary']['criteria_inconclusive']}/5",
        f"Status: {eval_results['summary']['overall_status']}",
        "",
    ]
    for k, v in eval_results["evaluation"].items():
        passed = v["passed"]
        status = "PASS" if passed is True else "FAIL" if passed is False else "INCONCLUSIVE"
        report_lines.append(f"{k}: {v['name']}")
        report_lines.append(f"  Status: {status}")
    report = "\n".join(report_lines)

    report_file = RESULTS_DIR / "step_013_falsification_criteria_report.txt"
    with open(report_file, "w") as f:
        f.write(report)
    print_status(f"Report saved to: {report_file}", "INFO")
    print_status("STEP 013 COMPLETED SUCCESSFULLY", "INFO")
    return True


if __name__ == "__main__":
    main()
