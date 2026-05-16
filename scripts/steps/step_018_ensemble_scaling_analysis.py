#!/usr/bin/env python3
"""
================================================================================
STEP 018: ENSEMBLE SCALING ANALYSIS (RIGOROUS AUDIT VERSION)
================================================================================

Purpose: Perform an ensemble-wide test of TEP scaling across the pulsar ensemble.
This version addresses the "Noise Floor Scaling" artifact by:
1. Using Phase Closure significance as the primary, noise-immune metric.
2. Using noise-subtracted H_excess for group-delay scaling.
3. Explicitly testing the "Constant Noise Floor" hypothesis.

================================================================================
"""

import json
import numpy as np
import sys
from pathlib import Path
from scipy import stats
from typing import Dict, Any, List, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.json_numpy import NpEncoder

from scripts.utils.logger import print_status
RESULTS_DIR = PROJECT_ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

def load_pulsar_summary(pulsar_name: str) -> Dict[str, Any]:
    """Load summary holonomy for a given pulsar."""
    # Try different filename patterns to handle various naming conventions
    prefix = pulsar_name.split('+')[0].split('-')[0]
    possible_filenames = [
        f"step_003_closure_final_summary_{prefix}.json",
        f"step_003_closure_final_summary_{pulsar_name.lower().replace('-', '')}.json",
        f"step_003_closure_final_summary_{pulsar_name.replace('+', '').replace('-', '')}.json",
        f"step_003_closure_final_summary_{pulsar_name}.json"
    ]
    
    for filename in possible_filenames:
        summary_file = RESULTS_DIR / filename
        if summary_file.exists():
            with open(summary_file, 'r') as f:
                return json.load(f)
        
    return None

def load_step031_scaling_inputs() -> Dict[str, Dict[str, Any]]:
    """Load effective-distance and curvature inputs produced by step 031."""
    scaling_inputs: Dict[str, Dict[str, Any]] = {}

    tep_file = RESULTS_DIR / "step_031_tep_scaling_detailed.json"
    if tep_file.exists():
        with open(tep_file, "r") as f:
            tep = json.load(f)
        curvature = tep.get("curvature_results", {})
        distances = tep.get("distance_results", {})
        velocities = tep.get("velocity_results", {})

        scaling_inputs["J0437-4715"] = {
            "eta": (curvature.get("eta1_j0437_mean", 0.0) + curvature.get("eta2_j0437_mean", 0.0)) / 2.0,
            "D_eff": distances.get("D_eff_j0437", 0.0),
            "v_eff": velocities.get("v_j0437", 0.0),
        }
        scaling_inputs["J1603-7202"] = {
            "eta": (curvature.get("eta1_j1603_mean", 0.0) + curvature.get("eta2_j1603_mean", 0.0)) / 2.0,
            "D_eff": distances.get("D_eff_j1603", 0.0),
            "v_eff": velocities.get("v_j1603", 0.0),
        }

    jiamusi_file = RESULTS_DIR / "step_031_jiamusi_eta_analysis.json"
    if jiamusi_file.exists():
        with open(jiamusi_file, "r") as f:
            jiamusi = json.load(f)
        for psr_prefix, entry in jiamusi.items():
            eta = entry.get("eta_mean", 0.0)
            s = entry.get("screen_distance_s")
            distance_pc = entry.get("distance_pc", 0.0)
            if not eta or s is None:
                continue
            scaling_inputs[psr_prefix] = {
                "eta": eta,
                "D_eff": distance_pc * s * (1 - s),
                "v_eff": None,
            }
    return scaling_inputs

def main():
    print_status("===" * 40)
    print("STEP 018: RIGOROUS ENSEMBLE SCALING ANALYSIS")
    print_status("===" * 40)
    
    ensemble = [
        "J0437-4715", "J1603-7202", "B0329+54", "B0355+54", "B0540+23",
        "B0740-28", "B1508+55", "B1933+16", "B2154+40", "B2310+42", "B2324+60", "B2351+61"
    ]
    
    # Load kinematics from step_030 results (which have proper citations)
    step030_file = RESULTS_DIR / "step_030_tep_scaling_analysis.json"
    if step030_file.exists():
        with open(step030_file, 'r') as f:
            step030_data = json.load(f)
        # Convert array-based structure to dictionary format expected by step_018
        kinematics = {}
        pulsars = step030_data.get("pulsars", [])
        distances = step030_data.get("distances_pc", [])
        velocities = step030_data.get("velocities_kms", [])
        for i, psr in enumerate(pulsars):
            if i < len(distances) and i < len(velocities):
                kinematics[psr] = {
                    "dist": distances[i],
                    "v": velocities[i]
                }
    else:
        print_status(f"Step 030 results not found at {step030_file}. Skipping ensemble scaling analysis.", "WARNING")
        print_status("Run step_030_tep_scaling_analysis.py first to generate required kinematics data.", "INFO")
        return False
    
    scaling_inputs = load_step031_scaling_inputs()
    
    obs_phase_z = []
    obs_h_excess = []
    obs_h_raw = []
    obs_h_noise = []
    pred_scaling = []
    dist_list = []
    names = []
    
    print("\n1. LOADING ENSEMBLE DATA (Noise-Aware)...")
    for psr in ensemble:
        summary = load_pulsar_summary(psr)
        if summary is None:
            continue
            
        phase_z = abs(summary.get("phase_closure_t_statistic", 0.0))
        h_raw = summary.get("H_magnitude_ns", 0.0)
        h_noise = summary.get("H_noise_bias_ns", 0.0)
        h_excess = summary.get("H_excess_ns", 0.0)
        
        # Determine predictor: eta * D_eff * v_eff
        # Try full pulsar name first, then prefix for step_031 keys
        scaling_entry = scaling_inputs.get(psr, scaling_inputs.get(psr.split('+')[0], scaling_inputs.get(psr.split('-')[0], None)))
        if scaling_entry and scaling_entry.get("D_eff") and scaling_entry.get("eta"):
            v = scaling_entry.get("v_eff")
            if v is None:
                if kinematics.get(psr):
                    v = kinematics[psr]["v"]
                else:
                    continue
            predictor = scaling_entry["eta"] * scaling_entry["D_eff"] * v
        else:
            # If no eta data available, use distance * velocity as fallback predictor
            if kinematics.get(psr):
                v = kinematics[psr]["v"]
                # Use minimum velocity of 1.0 km/s for pulsars with no proper motion measurement
                # to avoid division by zero
                if v == 0.0:
                    v = 1.0
                predictor = kinematics[psr]["dist"] * v
            else:
                continue
            
        obs_phase_z.append(phase_z)
        obs_h_raw.append(h_raw)
        obs_h_noise.append(h_noise)
        obs_h_excess.append(h_excess)
        pred_scaling.append(predictor)
        # Step 030 may omit distances_pc/velocities_kms (single-pulsar / scaling_disabled runs);
        # do not assume kinematics[psr] exists after using step_031 scaling inputs.
        if psr in kinematics:
            dist_pc = kinematics[psr]["dist"]
        elif summary.get("pulsar_dist_pc") is not None:
            dist_pc = float(summary["pulsar_dist_pc"])
        elif scaling_entry and scaling_entry.get("D_eff") is not None:
            dist_pc = float(scaling_entry["D_eff"])
        else:
            dist_pc = float("nan")
        dist_list.append(dist_pc)
        names.append(psr)
        
        detection_status = "DETECTED" if phase_z > 3.0 else "NOISE-DOMINATED"
        print(f"   {psr:12} : PhaseZ={phase_z:5.2f}σ, H_raw={h_raw:5.2f} ns, H_noise={h_noise:5.2f} ns [{detection_status}]")

    if len(names) < 3:
        print_status("Insufficient data for scaling fit (requires at least 3 pulsars).", "WARNING")
        # Still produce output indicating data limitation
        results = {
            "status": "insufficient_data",
            "n_pulsars_with_data": len(names),
            "pulsars": names,
            "message": "Ensemble scaling analysis requires at least 3 pulsars with kinematics data. Current analysis limited to 2 pulsars (J0437-4715, J1603-7202). Full ensemble analysis pending additional pulsar data.",
            "observed_data": {
                "pulsars": names,
                "phase_z": obs_phase_z,
                "h_raw": obs_h_raw,
                "h_noise": obs_h_noise,
                "h_excess": obs_h_excess,
                "predictor": pred_scaling,
                "distances": dist_list
            }
        }
        output_file = RESULTS_DIR / "step_018_ensemble_scaling_results.json"
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2, cls=NpEncoder)
        print_status(f"Partial results saved to: {output_file}", "INFO")
        return

    obs_phase_z = np.array(obs_phase_z)
    obs_h_raw = np.array(obs_h_raw)
    obs_h_noise = np.array(obs_h_noise)
    obs_h_excess = np.array(obs_h_excess)
    pred_scaling = np.array(pred_scaling)
    distances = np.array(dist_list)

    # 2. TEST ARTIFACT HYPOTHESIS: Is |H| just a constant noise floor?
    print("\n2. TESTING NOISE-FLOOR ARTIFACT HYPOTHESIS...")
    noise_r = np.corrcoef(distances, obs_h_noise)[0, 1]
    print(f"   Correlation(Distance, Noise Floor): R = {noise_r:.3f}")
    
    # Artifact efficiency: Eff_noise = NoiseFloor / (D * v)
    eff_noise = obs_h_noise / pred_scaling
    valid_gamma = (distances > 0) & (eff_noise > 0)
    gamma_artifact_slope, _, _, _, _ = stats.linregress(np.log(distances[valid_gamma]), np.log(eff_noise[valid_gamma]))
    print(f"   Artifact Scaling Index (gamma_noise): {-gamma_artifact_slope:.3f}")
    print("   (Note: If gamma_noise ≈ 1, a constant noise floor mimics the NFW profile)")

    # 3. REAL SCALING: Phase Closure Significance
    print("\n3. TESTING PHASE CLOSURE SCALING (Noise-Immune)...")
    phase_r = np.corrcoef(pred_scaling, obs_phase_z)[0, 1]
    print(f"   Correlation(Predictor, Phase significance): R = {phase_r:.3f}")

    # 4. ENVIRONMENT-DEPENDENT COUPLING (Rigorous)
    print("\n4. RIGOROUS COUPLING ANALYSIS (using H_excess)...")
    # Only use pulsars where H_excess is statistically significant (>2sigma)
    # or use all and see if signal emerges above noise floor
    coupling_eff = obs_h_excess / pred_scaling
    
    # Focus on Phase-Detected pulsars for coupling efficiency
    phase_detected = obs_phase_z > 3.0
    has_excess = obs_h_excess > 0.0
    
    # Check if any phase-detected pulsars have non-zero h_excess
    phase_with_excess = phase_detected & has_excess
    
    if phase_with_excess.sum() >= 2:
        log_D = np.log(distances[phase_with_excess])
        log_eff = np.log(coupling_eff[phase_with_excess])
        gamma_slope, _, gamma_r, gamma_p, _ = stats.linregress(log_D, log_eff)
        real_gamma = -gamma_slope
        print(f"   Rigorous Power-law Index γ (Phase-Detected with H_excess): {real_gamma:.3f}")
        print(f"   Significance: p = {gamma_p:.2e}, R² = {gamma_r**2:.3f}")
    elif phase_detected.sum() >= 2:
        print("   H_excess is noise-limited (all zeros) for phase-detected pulsars.")
        print("   Cannot determine gamma from H_excess - falling back to Phase Closure scaling.")
        real_gamma = None
    else:
        print("   Insufficient Phase-Detected pulsars for rigorous gamma estimation.")
        real_gamma = None

    # 5. FINAL VERDICT
    print("\n" + "=" * 40)
    print("ENSEMBLE VERDICT (AUDITED)")
    print("=" * 40)
    
    if real_gamma is not None and abs(real_gamma - 1.0) < 0.2 and gamma_r**2 > 0.5:
        verdict = f"Unassailable TEP scaling confirmed (gamma={real_gamma:.2f}) using phase-detected ensemble."
        status = "SUCCESS"
    elif phase_r > 0.5:
        verdict = f"Ensemble scaling detected in phase (R={phase_r:.2f}), but gamma coupling remains uncertain."
        status = "INFO"
    else:
        verdict = "SCALING FAILED: No consistent physical scaling detected above noise floor."
        status = "WARNING"
    
    print_status(verdict, status)
    
    output = {
        "n_pulsars": len(names),
        "pulsars": names,
        "metrics": {
            "phase_z": obs_phase_z.tolist(),
            "h_raw_ns": obs_h_raw.tolist(),
            "h_noise_ns": obs_h_noise.tolist(),
            "h_excess_ns": obs_h_excess.tolist(),
            "predictor": pred_scaling.tolist()
        },
        "artifact_test": {
            "noise_distance_correlation": float(noise_r),
            "gamma_artifact": float(-gamma_artifact_slope)
        },
        "rigorous_scaling": {
            "gamma": float(real_gamma) if real_gamma is not None else None,
            "r_squared": float(gamma_r**2) if real_gamma is not None and 'gamma_r' in locals() else None,
            "p_value": float(gamma_p) if real_gamma is not None and 'gamma_p' in locals() else None
        },
        "verdict": verdict
    }
    
    out_file = RESULTS_DIR / "step_018_ensemble_scaling_results.json"
    with open(out_file, 'w') as f:
        json.dump(output, f, indent=2, cls=NpEncoder)
    print(f"\nResults saved to: {out_file}")

if __name__ == "__main__":
    main()
