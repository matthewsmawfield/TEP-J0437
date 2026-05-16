#!/usr/bin/env python3
"""
================================================================================
STEP 012: ENVIRONMENTAL DEPENDENCE ANALYSIS
================================================================================

This step analyzes environmental dependence of the holonomy magnitude to
address the independent replication concern. It establishes testable scaling
relations based on the TEP v0.8 continuous geometric screening formulation,
where Temporal Shear (∇φ) is the operative quantity driving holonomy in
low-density environments.

PURPOSE:
--------
To provide a framework for testing TEP predictions across different ISM
environments using the Temporal Topology and Temporal Shear formulation,
enabling independent replication studies to validate environment-dependent
scaling of |H| through the continuous screening paradigm.

METHODOLOGY:
-------------
1. Analyze correlation between |H| and ISM-related parameters
2. Establish scaling relations with Temporal Shear, screen distance, loop geometry
3. Provide predictions for other pulsars with known ISM properties
4. Define replication criteria for independent studies based on continuous screening

OUTPUT:
-------
- Environmental scaling relations based on Temporal Shear
- Predictions for candidate replication pulsars
- Replication framework for independent groups
- Testable hypotheses for environmental dependence in the continuous screening paradigm

AUTHOR: TEP Analysis Framework
VERSION: 1.0.0
================================================================================
"""

import json
import sys
from pathlib import Path
from typing import Optional, Union

import numpy as np

# Add parent directory to path for imports
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.json_numpy import NpEncoder

from scripts.utils.logger import TEPLogger, print_status, set_step_logger

# Setup
# Logger is set by run_pipeline.py via set_step_logger()
# Do not create a new logger here to avoid overriding the pipeline's logger
RESULTS_DIR = PROJECT_ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def analyze_density_correlation(closure_delays: dict, ism_density: float) -> dict:
    """
    Analyze correlation between holonomy magnitude and ISM properties.

    In TEP v0.8, the holonomy scales with Temporal Shear (∇φ), which in
    the low-density ISM is proportional to the density gradient. The
    correlation therefore tracks how |H| varies with the spatial variation
    of the Temporal Topology rather than with density itself.

    Parameters
    ----------
    closure_delays : dict
        Closure delay measurements
    ism_density : float
        ISM density in g/cm³

    Returns
    -------
    dict
        Correlation analysis results
    """
    # NOTE: Full correlation analysis requires multi-pulsar data
    # This would analyze how |H| varies with proxies for Temporal Shear,
    # such as arclet count, SNR, and scattering strength, which reflect
    # the spatial variation of the Temporal Topology.
    # Current limitation: Single pulsar (J0437-4715) limits environmental dependence analysis
    # Arclet count and SNR serve as proxies for ISM complexity and Temporal Shear
    # v0.8 framework: Correlation tracks Temporal Shear (∇φ), the operative quantity in continuous geometric screening

    return {
        "shear_correlation": {
            "note": "Full analysis requires multi-pulsar data",
            "current_limitation": "Single pulsar (J0437-4715) limits environmental dependence analysis",
            "proxy_analysis": "Arclet count and SNR serve as proxies for ISM complexity and Temporal Shear",
            "v0.8_framework": "Correlation tracks Temporal Shear (∇φ), the operative quantity in continuous geometric screening",
            "status": "deferred_to_multi_pulsar_analysis"
        }
    }


def predict_for_other_pulsars() -> dict:
    """
    Predict holonomy magnitudes for candidate replication pulsars.

    Based on known ISM properties of other millisecond pulsars.

    Returns
    -------
    dict
        Predictions for candidate pulsars
    """
    # Load observed H from step_003 results for J0437
    j0437_file = PROJECT_ROOT / "results" / "step_003_closure_final_summary.json"
    if not j0437_file.exists():
        raise FileNotFoundError(
            f"J0437 results not found at {j0437_file}. "
            "Run step_003_closure_delays_final.py first."
        )
    with open(j0437_file, 'r') as f:
        j0437_data = json.load(f)
        observed_h_j0437 = j0437_data.get("H_magnitude_ns")
        if observed_h_j0437 is None:
            raise ValueError("H_magnitude_ns not found in step_003 results")
    
    # Candidate pulsars with known ISM properties
    candidates = [
        {
            "name": "PSR J0437-4715",
            "distance_pc": 156,
            "flux_mjy": 150,
            "screen_distances_pc": [89.8, 124.0],
            "density_g_cm3": 1e-24,
            "observed_H_ns": observed_h_j0437,
            "note": "Current study - reference measurement. Temporal Shear unsuppressed in local ISM.",
        },
        {
            "name": "PSR J1713+0747",
            "distance_pc": 1000,
            "flux_mjy": 8,
            "screen_distances_pc": None,  # Requires actual measurement from scintillation data
            "density_g_cm3": 1e-24,
            "predicted_H_ns": None,  # Cannot predict without screen distance
            "note": "Candidate for replication - requires scintillation data to measure screen distances before prediction.",
        },
        {
            "name": "PSR J1909-3744",
            "distance_pc": 1100,
            "flux_mjy": 3,
            "screen_distances_pc": None,  # Requires actual measurement from scintillation data
            "density_g_cm3": 1e-24,
            "predicted_H_ns": None,  # Cannot predict without screen distance
            "note": "Candidate for replication - requires scintillation data to measure screen distances before prediction.",
        },
        {
            "name": "PSR J1643-1224",
            "distance_pc": 350,
            "flux_mjy": 15,
            "screen_distances_pc": None,  # Requires actual measurement from scintillation data
            "density_g_cm3": 1e-24,
            "predicted_H_ns": None,  # Cannot predict without screen distance
            "note": "Candidate for replication - requires scintillation data to measure screen distances before prediction.",
        },
    ]

    return candidates


def define_replication_framework() -> dict:
    """
    Define framework for independent replication studies.

    Returns
    -------
    dict
        Replication framework with criteria and requirements
    """
    return {
        "target_pulsar_criteria": {
            "minimum_flux_mjy": 50,
            "minimum_observation_years": 5,
            "minimum_epochs": 100,
            "minimum_arclets_per_epoch": 3,
            "well_characterized_ism": True,
        },
        "analysis_requirements": {
            "use_same_pipeline": True,
            "fixed_random_seeds": True,
            "multiple_comparison_correction": "Bonferroni",
            "significance_threshold": "5sigma after correction",
        },
        "testable_hypotheses": {
            "TEP_prediction_1": "|H| scales with Temporal Shear (∇φ), which tracks ISM density gradients in the continuous screening paradigm",
            "TEP_prediction_2": "|H| scales with screen distance through integrated Temporal Shear",
            "TEP_prediction_3": "Bipolar structure preserved across environments",
            "TEP_prediction_4": "Magnitude equality holds across environments",
            "TEP_prediction_5": "Ambient Symmetry Restoration: α_eff ∝ D^-1.6 (environment-dependent scalar coupling scaling with Galactocentric distance)",
            "alternative_prediction": "No environmental dependence expected for standard ISM effects",
        },
        "ambient_symmetry_restoration": {
            "description": "TEP v0.8 predicts environment-dependent scalar coupling α_eff ∝ D^-1.6",
            "scaling_relation": "α_eff = α_0 * (D / D_0)^-1.6",
            "physical_basis": "Non-linear superposition of field gradients (Temporal Shear) in continuous screening",
            "implications": "Coupling strength varies with Galactic environment, predicting weaker effects at larger Galactocentric radii",
            "testable": True,
        },
        "success_criteria": {
            "detection_significance": ">5sigma after multiple comparison correction",
            "bipolar_structure_confirmed": "Both signs significant, magnitudes equal",
            "environmental_scaling": "|H varies with ISM parameters as predicted",
            "reproducibility": "Independent groups obtain consistent results",
        },
    }


def step_main(logger=None, verbose=True):
    """Standard pipeline entry point for environmental dependence analysis."""
    return main()


def main():
    """Execute environmental dependence analysis."""
    print_status("=" * 80, "INFO")
    print_status("STEP 012: ENVIRONMENTAL DEPENDENCE ANALYSIS", "INFO")
    print_status("=" * 80, "INFO")

    print_status("Analyzing environmental dependence of holonomy magnitude", "INFO")
    print_status(
        "Purpose: Establish framework for independent replication studies", "INFO"
    )

    # Load closure delay results if available
    closure_file = PROJECT_ROOT / "results" / "step_003_closure_final_summary.json"
    if not closure_file.exists():
        raise FileNotFoundError(
            f"Closure delays file not found at {closure_file}. "
            "Run step_003_closure_delays_final.py first."
        )
    with open(closure_file, "r") as f:
        closure_delays = json.load(f)
    print_status(f"Loaded closure delays from {closure_file}", "INFO")

    # Analyze density correlation
    correlation_results = analyze_density_correlation(closure_delays, 1e-24)
    print_status("\nDensity Correlation Analysis:", "INFO")
    print_status(json.dumps(correlation_results, indent=2), "INFO")

    # Predict for other pulsars
    predictions = predict_for_other_pulsars()
    print_status("\nPredictions for Candidate Replication Pulsars:", "INFO")
    print_status("-" * 80, "INFO")
    for pulsar in predictions:
        print_status(f"\nPulsar: {pulsar['name']}", "INFO")
        print_status(f"  Distance: {pulsar['distance_pc']} pc", "INFO")
        print_status(f"  Flux: {pulsar['flux_mjy']} mJy", "INFO")
        print_status(f"  Screen distances: {pulsar['screen_distances_pc']} pc", "INFO")
        print_status(f"  Density: {pulsar['density_g_cm3']:.2e} g/cm³", "INFO")
        if "observed_H_ns" in pulsar:
            if pulsar['observed_H_ns'] is not None:
                print_status(f"  |Observed H|: {pulsar['observed_H_ns']:.3f} ns", "INFO")
            else:
                print_status(f"  |Observed H|: None (requires measurement)", "INFO")
        if "predicted_H_ns" in pulsar:
            if pulsar['predicted_H_ns'] is not None:
                print_status(f"  |Predicted H|: {pulsar['predicted_H_ns']:.3f} ns", "INFO")
            else:
                print_status(f"  |Predicted H|: None (requires screen distance measurement)", "INFO")
        print_status(f"  Note: {pulsar['note']}", "INFO")

    # Define replication framework
    framework = define_replication_framework()
    print_status("\nReplication Framework:", "INFO")
    print_status("-" * 80, "INFO")
    print_status(json.dumps(framework, indent=2), "INFO")

    # Compile results
    results = {
        "correlation_analysis": correlation_results,
        "candidate_pulsars": predictions,
        "replication_framework": framework,
        "current_limitations": [
            "Single pulsar analysis limits environmental dependence testing",
            "Independent replication requires multi-pulsar observational campaign",
            "Theoretical scaling relations need proper derivation from TEP field equations",
        ],
        "next_steps": [
            "Obtain scintillation data for candidate pulsars",
            "Apply same pipeline to new pulsars",
            "Test environmental scaling predictions",
            "Independent group replication using same methodology",
        ],
    }

    # Save results
    output_file = RESULTS_DIR / "step_012_environmental_dependence_results.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2, cls=NpEncoder)

    print_status(f"\nResults saved to: {output_file}", "INFO")
    print_status("=" * 80, "INFO")
    print_status("STEP 012 COMPLETED SUCCESSFULLY", "INFO")
    print_status("=" * 80, "INFO")


if __name__ == "__main__":
    main()
