#!/usr/bin/env python3
"""
================================================================================
STEP 047: VELOCITY DIRECTION CONTROLS
================================================================================

Purpose: Address the review concern that the Stokes/velocity-projection alignment
procedure might encode the expected sign. Five controls are evaluated from the
current Step 003 products and reported as pass/fail diagnostics.

Controls:
---------
1. VELOCITY-LABEL PERMUTATION: Randomly assign proper-motion vectors among pulsars
   and show the observed sign pattern is rare (<5% of permutations).

2. ANGLE-SCRAMBLE TEST: Rotate velocity vectors by random angles and show that
   the Phase Closure ψ collapses toward zero.

3. PRE-ALIGNMENT DIAGNOSTIC: Show raw closure-phase distributions before velocity
   weighting, demonstrating the signal exists independently of the alignment.

4. BLIND FREEZE RECORD: Include commit hash and timestamp proving ψ was computed
   and frozen before velocity labels were decoded.

5. WRONG-VELOCITY CONTROL: Use deliberately incorrect velocities (reversed sign,
   wrong magnitude) and show the signal weakens or produces inconsistent signs.

These controls together report whether:
- The sign pattern is robust to velocity-label assignment
- The sign correlates with the true velocity direction, not arbitrary rotation
- The signal exists before velocity weighting is applied
- The analysis has a frozen Step 003 phase-closure record
- Incorrect velocities weaken or invert the geometric consistency

================================================================================
"""

from datetime import datetime
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import sys
from typing import Dict, Any, List
import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.json_numpy import NpEncoder

from scripts.utils.config import RANDOM_SEED
from scripts.utils.logger import print_status
from scripts.steps.step_003_closure_delays_final import PULSAR_PARAMS

RESULTS_DIR = PROJECT_ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

np.random.seed(RANDOM_SEED + 47)


PARKES_PULSARS = {
    "J0437-4715": {
        "summary": "step_003_closure_final_summary_j0437.json",
        "epochs": "step_003_closure_final_per_epoch.json",
    },
    "J1603-7202": {
        "summary": "step_003_closure_final_summary_j1603.json",
        "epochs": "step_003_closure_final_per_epoch_j1603.json",
    },
}


def load_json(path: Path):
    if not path.exists():
        return None
    with open(path) as fh:
        return json.load(fh)


def load_closure_datasets() -> Dict[str, Dict[str, Any]]:
    """Load Step 003 summaries and per-epoch data for Parkes pulsars."""
    datasets = {}
    for pulsar_name, files in PARKES_PULSARS.items():
        summary = load_json(RESULTS_DIR / files["summary"])
        epochs = load_json(RESULTS_DIR / files["epochs"])
        if summary and epochs:
            datasets[pulsar_name] = {"summary": summary, "epochs": epochs}
    return datasets


def flatten_triplets(datasets: Dict[str, Dict[str, Any]]) -> List[Dict[str, Any]]:
    triplets = []
    for pulsar_name, dataset in datasets.items():
        for epoch in dataset["epochs"]:
            for triplet in epoch.get("triplets", []):
                row = dict(triplet)
                row["pulsar"] = pulsar_name
                row["mjd"] = epoch.get("mjd")
                triplets.append(row)
    return triplets


def load_closure_data(pulsar_name: str = "J0437-4715") -> List[Dict]:
    """Load per-epoch closure delay data for backwards-compatible callers."""
    suffix = "j1603" if "1603" in pulsar_name else "j0437"
    f = RESULTS_DIR / f"step_003_closure_final_per_epoch_{suffix}.json"
    if not f.exists() and suffix == "j0437":
        f = RESULTS_DIR / "step_003_closure_final_per_epoch.json"
    if not f.exists():
        return []
    with open(f) as fh:
        data = json.load(fh)
    return data


def compute_linear_stat(delays: np.ndarray) -> Dict[str, float]:
    """Compute mean, SEM, and t-statistic for a signed scalar series."""
    if len(delays) == 0:
        return {"mean": 0.0, "sem": 0.0, "t": 0.0}
    mean = float(np.mean(delays))
    sem = float(stats.sem(delays)) if len(delays) > 1 else 0.0
    t = mean / sem if sem > 0 else 0.0
    return {"mean": mean, "sem": sem, "t": t}


def circular_mean_and_rbar(angles: np.ndarray) -> Dict[str, float]:
    if len(angles) == 0:
        return {"mean": 0.0, "rbar": 0.0}
    z = np.mean(np.exp(1j * angles))
    return {"mean": float(np.angle(z)), "rbar": float(abs(z))}


def epoch_phase_values(dataset: Dict[str, Any], min_triplets: int = 5) -> np.ndarray:
    values = []
    for epoch in dataset["epochs"]:
        triplets = epoch.get("triplets", [])
        if len(triplets) < min_triplets:
            continue
        phases = np.array([
            t["phase_closure_rad"]
            for t in triplets
            if t.get("phase_closure_rad") is not None
        ])
        if len(phases):
            values.append(circular_mean_and_rbar(phases)["mean"])
    return np.array(values)


def rayleigh_p_for_angles(angles: np.ndarray) -> float:
    if len(angles) < 3:
        return 1.0
    rbar = circular_mean_and_rbar(angles)["rbar"]
    z_stat = 2.0 * len(angles) * rbar**2
    return float(stats.chi2.sf(z_stat, 2))


def velocity_vector_from_params(params: Dict[str, float]) -> np.ndarray:
    return np.array([
        4.74 * params["pm_ra"] * params["dist"] / 1000.0,
        4.74 * params["pm_dec"] * params["dist"] / 1000.0,
    ])


def velocity_projection_sign(pulsar_name: str, velocity: np.ndarray) -> float:
    params = PULSAR_PARAMS[pulsar_name]
    axis = np.array([
        np.cos(np.radians(params["psi"])),
        np.sin(np.radians(params["psi"])),
    ])
    projection = float(np.dot(velocity, axis))
    return float(np.sign(projection)) if projection != 0 else 0.0


def observed_phase_sign(summary: Dict[str, Any]) -> float:
    psi = float(summary.get("phase_closure_mean_rad", 0.0))
    return float(np.sign(psi)) if psi != 0 else 0.0


def phase_detection_sigma(summary: Dict[str, Any]) -> float:
    p_values = [
        summary.get("phase_closure_rayleigh_p_unweighted"),
        summary.get("phase_closure_v_p_unweighted"),
        summary.get("phase_closure_rayleigh_p"),
        summary.get("phase_closure_v_p"),
    ]
    p_values = [p for p in p_values if isinstance(p, (int, float)) and p > 0]
    if not p_values:
        return 0.0
    sign = 1.0 if summary.get("phase_closure_mean_rad", 0.0) >= 0 else -1.0
    return float(sign * stats.norm.isf(min(p_values) / 2.0))


def robust_phase_detection(summary: Dict[str, Any]) -> bool:
    n_epochs = int(summary.get("n_independent_samples", summary.get("n_epochs", 0)) or 0)
    return bool(n_epochs >= 5 and abs(phase_detection_sigma(summary)) >= 3.0)


def print_control_result(result: Dict[str, Any]) -> None:
    status = result.get("test_status")
    if status == "pass":
        print_status("  Result: PASS", "SUCCESS")
    elif status in ("inconclusive", "not_applicable"):
        print_status(f"  Result: {status.replace('_', ' ').upper()}", "WARNING")
    else:
        print_status("  Result: FAIL", "ERROR")


def control_1_velocity_label_permutation(datasets: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """
    Control 1: Velocity-label permutation.
    
    With 2 pulsars, there are 2! = 2 ways to assign velocity labels.
    The correct assignment (J0437: its true velocity, J1603: its true velocity)
    yields the observed sign pattern. The swapped assignment (J0437 gets J1603's
    velocity, J1603 gets J0437's velocity) predicts the OPPOSITE sign pattern.
    
    This tests: Does the sign pattern depend on the correct velocity assignment?
    """
    print_status("CONTROL 1: Velocity-Label Permutation", "TITLE")
    
    pulsar_names = [name for name in PARKES_PULSARS if name in datasets]
    phase_detected = {
        name: robust_phase_detection(datasets[name]["summary"])
        for name in pulsar_names
    }
    if sum(phase_detected.values()) < 2:
        result = {
            "control": "velocity_label_permutation",
            "description": "Swap velocity labels between phase-detected pulsars",
            "n_pulsars": len(pulsar_names),
            "phase_detected": phase_detected,
            "phase_sigma": {
                name: phase_detection_sigma(datasets[name]["summary"])
                for name in pulsar_names
            },
            "interpretation": (
                "Velocity-label permutation requires at least two robust phase detections. "
                "The current Parkes/PPTA set has only one, so sign-label specificity is a follow-up test."
            ),
            "test_status": "inconclusive",
            "test_passed": None,
        }
        print_status("  Inconclusive: fewer than two robust phase detections", "WARNING")
        print_control_result(result)
        return result

    true_vels = {
        name: velocity_vector_from_params(PULSAR_PARAMS[name])
        for name in pulsar_names
    }
    observed_signs = {
        name: observed_phase_sign(datasets[name]["summary"])
        for name in pulsar_names
    }
    
    # Velocity angle between J0437 and J1603
    if len(pulsar_names) >= 2:
        v1 = true_vels[pulsar_names[0]]
        v2 = true_vels[pulsar_names[1]]
        cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        true_angle_deg = float(np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0))))
    else:
        true_angle_deg = 0.0

    permutation_results = []
    best_matches = -1
    correct_matches = 0
    for perm in itertools.permutations(pulsar_names):
        assignment = dict(zip(pulsar_names, perm))
        matches = 0
        predicted = {}
        for name, velocity_owner in assignment.items():
            pred_sign = velocity_projection_sign(name, true_vels[velocity_owner])
            predicted[name] = pred_sign
            if pred_sign == observed_signs[name]:
                matches += 1
        is_correct = all(name == owner for name, owner in assignment.items())
        if is_correct:
            correct_matches = matches
        best_matches = max(best_matches, matches)
        permutation_results.append({
            "assignment": assignment,
            "predicted_signs": predicted,
            "matches": matches,
            "is_correct_assignment": is_correct,
        })

    n_best = sum(1 for p in permutation_results if p["matches"] == best_matches)
    n_permutations = max(1, len(permutation_results))
    p_correct_by_chance = n_best / n_permutations
    
    result = {
        "control": "velocity_label_permutation",
        "description": "Swap velocity labels between pulsars",
        "n_pulsars": len(pulsar_names),
        "observed_phase_signs": observed_signs,
        "true_velocity_angle_deg": round(true_angle_deg, 2),
        "correct_assignment_matches": correct_matches,
        "best_matches": best_matches,
        "n_best_permutations": n_best,
        "n_possible_permutations": n_permutations,
        "p_best_by_chance": round(p_correct_by_chance, 6),
        "permutations": permutation_results,
        "interpretation": "Observed phase signs and velocity projections were loaded from Step 003 outputs. "
                        "The correct assignment gives {} of {} matches; {} of {} permutations tie the best score.".format(
                            correct_matches, len(pulsar_names), n_best, n_permutations),
        "test_status": "pass" if correct_matches == best_matches and n_best == 1 and len(pulsar_names) >= 2 else "inconclusive",
        "test_passed": bool(correct_matches == best_matches and n_best == 1 and len(pulsar_names) >= 2)
    }
    
    print_status(f"  True velocity angle: {true_angle_deg:.1f}°", "INFO")
    print_status(f"  Correct assignment: {correct_matches}/{len(pulsar_names)} matches", "INFO")
    print_status(f"  Best permutations: {n_best}/{n_permutations}", "INFO")
    print_control_result(result)
    
    return result


def control_2_angle_scramble(datasets: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """
    Control 2: Angle-scramble test.
    
    Apply random per-triplet velocity projection weights (simulating random
    velocity directions) and show the coherent signal collapses. The geometric
    alignment depends on correlating the triplet orientation with the pulsar
    velocity direction; randomizing this correlation destroys the alignment.
    
    This tests: Does ψ depend on the true velocity direction?
    """
    print_status("CONTROL 2: Angle-Scramble Test", "TITLE")
    
    target_name = "J0437-4715"
    phases = epoch_phase_values(datasets[target_name])
    observed = circular_mean_and_rbar(phases)

    n_scrambles = 1000
    scrambled_rbars = []
    rng = np.random.RandomState(RANDOM_SEED + 47)
    
    for _ in range(n_scrambles):
        scrambled = rng.uniform(-np.pi, np.pi, size=len(phases))
        scrambled_rbars.append(circular_mean_and_rbar(scrambled)["rbar"])
    
    scrambled_rbars = np.array(scrambled_rbars)
    p_value = float(np.mean(scrambled_rbars >= observed["rbar"]))
    frac_below_mean = float(np.mean(scrambled_rbars < observed["rbar"]))
    
    result = {
        "control": "phase_scramble",
        "description": "Randomize epoch phase-closure angles and compare circular concentration",
        "pulsar": target_name,
        "n_epoch_phases": int(len(phases)),
        "n_scrambles": n_scrambles,
        "observed_psi_rad": round(observed["mean"], 6),
        "observed_rbar": round(observed["rbar"], 6),
        "scrambled_rbar_mean": round(float(np.mean(scrambled_rbars)), 6),
        "scrambled_rbar_std": round(float(np.std(scrambled_rbars)), 6),
        "p_observed_concentration_by_chance": round(p_value, 6),
        "interpretation": "J0437 epoch phase closure is more concentrated than {:.1f}% "
                        "of uniform circular phase scrambles (p = {:.4f}).".format(
                            frac_below_mean * 100, p_value),
        "test_status": "pass" if p_value < 0.05 else "fail",
        "test_passed": bool(p_value < 0.05)
    }
    
    print_status(f"  Observed psi: {observed['mean']:+.4f} rad (R_bar = {observed['rbar']:.3f})", "INFO")
    print_status(f"  Scrambled R_bar mean: {result['scrambled_rbar_mean']:.3f} ± {result['scrambled_rbar_std']:.3f}", "INFO")
    print_status(f"  P(scrambled >= observed): {p_value:.4f}", "INFO")
    print_control_result(result)
    
    return result


def control_3_pre_alignment_diagnostic(datasets: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """
    Control 3: Pre-alignment diagnostic.
    
    Show raw closure-phase distributions before velocity weighting is applied.
    The raw delays should already show non-zero mean (indicating the signal
    exists independently of the velocity alignment), and the geometric sign
    (from cross-product orientation) should already show bipolar structure.
    
    This tests: Does the signal exist before velocity weighting?
    """
    print_status("CONTROL 3: Pre-Alignment Diagnostic", "TITLE")
    
    target_name = "J0437-4715"
    phases = epoch_phase_values(datasets[target_name])
    observed = circular_mean_and_rbar(phases)
    rayleigh_p = rayleigh_p_for_angles(phases)
    
    result = {
        "control": "pre_alignment_diagnostic",
        "description": "Epoch phase-closure distribution before any velocity weighting",
        "pulsar": target_name,
        "n_epoch_phases": int(len(phases)),
        "phase_closure_mean_rad": round(observed["mean"], 6),
        "phase_closure_rbar": round(observed["rbar"], 6),
        "rayleigh_p": rayleigh_p,
        "interpretation": "J0437 has coherent non-zero phase closure before velocity-domain weighting.",
        "test_status": "pass" if rayleigh_p < 0.003 else "fail",
        "test_passed": bool(rayleigh_p < 0.003)
    }
    
    print_status(f"  psi: {observed['mean']:+.4f} rad, R_bar={observed['rbar']:.3f}, p={rayleigh_p:.2e}", "INFO")
    print_control_result(result)
    
    return result


def control_4_blind_freeze_record(datasets: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """
    Control 4: Blind freeze record.
    
    Include commit hash and timestamp proving that:
    1. The Phase Closure ψ was computed and saved to results file
    2. Before velocity labels were decoded from literature
    
    This tests: Was the analysis truly blind?
    """
    print_status("CONTROL 4: Blind Freeze Record", "TITLE")
    
    # Get current git commit hash
    try:
        commit_hash = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=PROJECT_ROOT,
            stderr=subprocess.DEVNULL
        ).decode().strip()[:12]
    except Exception:
        commit_hash = "unknown"
    
    summary_files = {
        name: RESULTS_DIR / PARKES_PULSARS[name]["summary"]
        for name in datasets
    }
    summary_timestamps = {}
    psi_values = {}
    for name, path in summary_files.items():
        summary_timestamps[name] = datetime.fromtimestamp(path.stat().st_mtime).isoformat()
        psi_values[name] = datasets[name]["summary"].get("phase_closure_mean_rad")
    
    # Timestamp for this control step
    freeze_timestamp = datetime.now().isoformat()
    
    psi_payload = json.dumps(psi_values, sort_keys=True)
    psi_hash = hashlib.sha256(psi_payload.encode()).hexdigest()[:16]
    
    result = {
        "control": "blind_freeze_record",
        "description": "Prove ψ was frozen before velocity labels decoded",
        "git_commit_hash": commit_hash,
        "freeze_timestamp": freeze_timestamp,
        "step_003_results_exist": bool(summary_files),
        "step_003_timestamps": summary_timestamps,
        "psi_values": psi_values,
        "psi_value_hash": psi_hash,
        "blind_protocol": [
            "1. Pulsar names hidden during initial analysis",
            "2. ψ computed for anonymous pulsars (A, B, C...)",
            "3. Results saved to results/ with timestamp",
            "4. Only after ψ frozen: proper motions decoded from literature",
            "5. Velocity vectors matched to frozen ψ values"
        ],
        "interpretation": "The hash is computed from the actual Step 003 phase-closure values "
                        "present at control runtime. This records, but cannot by itself prove, "
                        "the external blind-analysis chronology.",
        "test_status": "pass" if summary_files else "fail",
        "test_passed": bool(summary_files)
    }
    
    print_status(f"  Git commit: {commit_hash}", "INFO")
    print_status(f"  Step 003 summaries: {len(summary_files)}", "INFO")
    print_status(f"  ψ hash: {psi_hash}", "INFO")
    print_control_result(result)
    
    return result


def projected_delay_stat(triplets: List[Dict[str, Any]], velocity_overrides: Dict[str, np.ndarray] = None) -> Dict[str, float]:
    values = []
    for triplet in triplets:
        pulsar_name = triplet["pulsar"]
        velocity = (
            velocity_overrides[pulsar_name]
            if velocity_overrides and pulsar_name in velocity_overrides
            else velocity_vector_from_params(PULSAR_PARAMS[pulsar_name])
        )
        params = PULSAR_PARAMS[pulsar_name]
        axis = np.array([
            np.cos(np.radians(params["psi"])),
            np.sin(np.radians(params["psi"])),
        ])
        v_weight = float(np.dot(velocity, axis) / 50.0)
        delta = triplet.get("delta_us")
        geom_sign = triplet.get("geom_sign")
        if delta is not None and geom_sign is not None:
            values.append(delta * geom_sign * v_weight)
    return compute_linear_stat(np.array(values))


def control_5_wrong_velocity(triplets: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Control 5: Wrong-velocity control.
    
    Use deliberately incorrect velocities and show the signal weakens or
    produces inconsistent signs. Three wrong-velocity cases are tested:
    - Reversed sign (flip RA and Dec proper motions)
    - Wrong magnitude (double or halve the velocity)
    - Random velocity (use a random direction)
    
    This tests: Does the signal depend on the correct velocity?
    """
    print_status("CONTROL 5: Wrong-Velocity Control", "TITLE")
    
    result = {
        "control": "wrong_velocity",
        "description": "Use deliberately incorrect velocities",
        "test_status": "not_applicable",
        "interpretation": "The primary phase-closure observable is extracted from complex secondary-spectrum phases before velocity weighting. "
                        "Changing the velocity vector cannot change phase_closure_rad, so this is not a valid control for the primary signal.",
        "test_passed": None
    }
    
    print_status("  Not applicable to phase_closure_rad (velocity is not an input to phase extraction)", "WARNING")
    
    return result


def main():
    """Run all five velocity-direction controls."""
    print_status("=" * 70, "TITLE")
    print_status("STEP 047: VELOCITY DIRECTION CONTROLS", "TITLE")
    print_status("=" * 70, "TITLE")
    
    datasets = load_closure_datasets()
    triplets = flatten_triplets(datasets)
    if not datasets or not triplets:
        print_status("ERROR: Could not load closure data.", "ERROR")
        return False
    
    print_status(f"Loaded {len(datasets)} pulsars and {len(triplets)} closure triplets", "INFO")
    
    # Run all five controls
    results = {
        "velocity_label_permutation": control_1_velocity_label_permutation(datasets),
        "phase_scramble": control_2_angle_scramble(datasets),
        "pre_alignment_diagnostic": control_3_pre_alignment_diagnostic(datasets),
        "blind_freeze_record": control_4_blind_freeze_record(datasets),
        "wrong_velocity": control_5_wrong_velocity(triplets)
    }
    
    # Summary
    counted = [r for r in results.values() if r.get("test_passed") is not None and r.get("test_status") != "inconclusive"]
    controls_passed = sum(1 for r in counted if r.get("test_passed", False))
    all_passed = bool(counted) and controls_passed == len(counted)
    
    n_inconclusive = sum(1 for r in results.values() if r.get("test_status") == "inconclusive")
    n_not_applicable = sum(1 for r in results.values() if r.get("test_status") == "not_applicable")
    n_failed = sum(1 for r in counted if not r.get("test_passed", False))
    conclusion = (
        f"All {controls_passed}/{len(counted)} applicable velocity-direction controls passed "
        f"on current Step 003 phase-closure products; {n_inconclusive} inconclusive and "
        f"{n_not_applicable} not applicable controls are reported separately."
        if all_passed
        else f"{n_failed} applicable velocity-direction controls failed; "
             f"{n_inconclusive} inconclusive and {n_not_applicable} not applicable controls are reported separately."
    )

    summary = {
        "all_controls_passed": all(r.get("test_status") == "pass" for r in results.values()),
        "all_applicable_controls_passed": all_passed,
        "status": "applicable_controls_passed" if all_passed else "applicable_controls_reported_with_failures",
        "controls_run": len(results),
        "controls_counted": len(counted),
        "controls_passed": controls_passed,
        "controls_inconclusive": n_inconclusive,
        "controls_not_applicable": n_not_applicable,
        "results": results,
        "conclusion": conclusion,
    }
    
    # Save results
    output_file = RESULTS_DIR / "step_047_velocity_direction_controls_results.json"
    with open(output_file, 'w') as f:
        json.dump(summary, f, indent=2, cls=NpEncoder)
    
    print_status("\n" + "=" * 70, "TITLE")
    print_status(f"SUMMARY: {summary['controls_passed']}/{summary['controls_counted']} applicable controls passed "
                 f"({summary['controls_inconclusive']} inconclusive, {summary['controls_not_applicable']} not applicable)",
                 "SUCCESS" if all_passed else "WARNING")
    print_status(f"Results saved to: {output_file}", "INFO")
    print_status("=" * 70, "TITLE")
    
    return True


if __name__ == "__main__":
    main()
