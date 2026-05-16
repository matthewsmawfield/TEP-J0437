#!/usr/bin/env python3
"""
================================================================================
TEP-J0437: TEMPORAL EQUIVALENCE PRINCIPLE DETECTION PIPELINE
================================================================================

A complete, reproducible analysis pipeline for detecting temporal equivalence principle
(TEP) signals in pulsar scintillation data via closure delay triangulation.

PIPELINE ARCHITECTURE:
----------------------
    Step 000 : Data Ingestion
        Downloads J0437-4715 scintillation data from Scintools ATNF archive

    Step 001 : Parse Dynamic Spectra
        Converts ASCII dynspec files to calibrated NumPy arrays

    Step 002 : Secondary Spectrum Generation
        Computes S(τ, f_D) via 2D FFT and detects scintillation arcs

    Step 003 : Closure Delay Extraction
        Measures closure delays Δ = τ₀₁ + τ₁₂ − τ₀₂ from arclet triplets

    Step 004 : Verification & Validation
        Runs independent verification, blind analysis, and alternative ISM tests

REPRODUCIBILITY:
----------------
This pipeline is designed for full reproducibility. New users can execute:

    python scripts/run_pipeline.py

to automatically download data, process all epochs, and generate TEP detection
results. All steps are logged with timestamps for scientific traceability.

OPTIMIZATION:
-------------
Computations are optimized for Apple Silicon (M1/M2/M3/M4 Pro/Max/Ultra) via:
    • Accelerate framework vectorized FFT (Step 002)
    • Parallel epoch processing with optimal worker detection
    • Memory-mapped array operations

SCIENTIFIC OUTPUT:
------------------
The pipeline produces:
    • Secondary spectra S(τ, f_D) for each epoch
    • Detected arclet catalogs with (τᵢ, f_Dᵢ, SNR) measurements
    • Closure delay distributions with statistical significance tests
    • TEP detection verdict with confidence intervals

AUTHOR: TEP Analysis Framework
VERSION: 3.0.0 (Reorganized Structure)
================================================================================
"""

from typing import Union, Optional, Dict, Any

import argparse
import multiprocessing as mp
import os
import sys
import time

# Types updated to built-ins (Python 3.9+ compatible)
from datetime import datetime
from pathlib import Path

# Configuration
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "scripts" / "steps"))

# Import logger utility
from scripts.utils.logger import TEPLogger, set_step_logger

# Orchestration log (per-run path set in run_pipeline via pipeline_log=…)
_PIPELINE_LOG_PATH: Path = ROOT / "logs" / "pipeline.log"
# All pipeline logs live under logs/ (same directory as pipeline.log)
LOGS_DIR: Path = ROOT / "logs"

LOGS_DIR.mkdir(exist_ok=True)


def log_message(message: str, level: str = "INFO", verbose: bool = True):
    """Log message to console and orchestration log file with scientific formatting."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_line = f"[{timestamp}] [{level}] {message}"

    if verbose:
        print(log_line)

    _PIPELINE_LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(_PIPELINE_LOG_PATH, "a") as f:
        f.write(log_line + "\n")


def run_step(step_name: str, step_module: str, verbose: bool = True, force: bool = False) -> bool:
    """
    Execute a pipeline step with comprehensive logging and scientific context.

    This function orchestrates the execution of individual pipeline steps, providing
    detailed logging of the astrophysical processing, computational operations,
    and statistical analysis being performed at each stage.

    Parameters
    -----------
    step_name : str
        Human-readable step identifier describing the astrophysical operation
    step_module : str
        Python module name containing the step implementation
    verbose : bool
        Enable detailed scientific output with physical interpretations

    Returns
    --------
    success : bool
        True if step completed successfully with valid scientific outputs
    """

    # Scientific context for each step
    step_context = {
        "step_000_data_ingestion": {
            "description": "Data acquisition from CSIRO Data Access Portal",
            "physics": "Retrieving calibrated dynamic spectra for PSR J0437-4715 and PSR J1603-7202",
            "method": "DOI-resolved CSIRO DAP retrieval with SHA-256 manifest verification",
            "output": "Verified raw dynspec files for secondary spectrum computation",
        },
        "step_001_parse_dynspec": {
            "description": "Dynamic spectrum calibration and parsing",
            "physics": "Extracting intensity I(t,ν) from filterbank data",
            "method": "RFI masking, baseline subtraction, bandpass calibration",
            "output": "Calibrated dynamic spectra as NumPy arrays",
        },
        "step_002_secondary_spectra": {
            "description": "Secondary spectrum computation S(τ, f_D)",
            "physics": "2D Fourier transform of covariance C(Δt, Δν) = ⟨δI δI⟩",
            "method": "Vectorized FFT with conjugate symmetry optimization",
            "output": "Secondary spectra revealing parabolic scintillation arcs",
        },
        "step_003_closure_delays_final": {
            "description": "Closure delay extraction from arclet triplets (PARALLEL)",
            "physics": "Measuring holonomy H = ∮ dτ_prop via triangulation",
            "method": "Parabolic fitting to arclets, parallel epoch processing with ProcessPoolExecutor",
            "output": "Bipolar closure delay distribution testing TEP predictions",
        },
        "step_004_verification": {
            "description": "Verification and validation of TEP detection",
            "physics": "Testing robustness against alternative explanations and systematic errors",
            "method": "Independent verification, blind analysis, alternative ISM phenomenon tests",
            "output": "Validation report quantifying robustness and remaining alternative-explanation bounds",
        },
        "step_005_enhanced_validation": {
            "description": "Enhanced validation addressing reviewer concerns",
            "physics": "Advanced tests for systematic errors, alternative explanations, and theoretical framework",
            "method": "Improved GR null tests, advanced ISS tests, blind analysis, rigorous anisotropy derivation, weighted analysis verification",
            "output": "Comprehensive enhanced validation report addressing all reviewer concerns",
        },
        "step_006_advanced_validation": {
            "description": "Advanced validation tests",
            "physics": "Systematic error investigation for both-signs significance, aggregate mean deviation, stratified analysis",
            "method": "Both-signs systematic error test, aggregate mean deviation test, stratified analysis, sensitivity analysis, bootstrap reconciliation, orientation methodology verification",
            "output": "Advanced validation results distinguishing phase-domain robustness from delay-domain diagnostics",
        },
        "step_007_independent_validation": {
            "description": "Independent statistical validation",
            "physics": "Cross-validation using multiple independent statistical approaches",
            "method": "Effect size measures, robust statistics, permutation tests, cross-validation, Bayesian model comparison, outlier robustness analysis",
            "output": "Independent statistical validation with explicit effect-size and diagnostic caveats",
        },
        "step_008_alternative_explanations": {
            "description": "Alternative explanations testing",
            "physics": "Systematic simulation of standard ISM and instrumental effects",
            "method": "Scattering geometry simulation, screen multiplicity simulation, velocity gradient simulation, instrumental systematic simulation, bipolar structure reproduction test",
            "output": "Alternative-explanation audit identifying which standard ISM/instrumental effects are disfavored",
        },
        "step_009_parameter_sensitivity": {
            "description": "Parameter sensitivity analysis",
            "physics": "Testing robustness to parameter choices and threshold selections",
            "method": "CV threshold analysis, bootstrap sample size analysis, significance threshold analysis, edge margin analysis, SNR threshold analysis",
            "output": "Parameter sensitivity results documenting threshold robustness and weak points",
        },
        "step_010_data_quality_metrics": {
            "description": "Comprehensive data quality metrics",
            "physics": "Assessment of dataset integrity and completeness",
            "method": "Overall dataset statistics, epoch-by-epoch quality metrics, SNR distribution analysis, temporal trend analysis, outlier detection, data completeness assessment",
            "output": "Data quality metrics demonstrating excellent dataset integrity",
        },
        "step_011_ism_density_modeling": {
            "description": "ISM Temporal Topology and Temporal Shear modeling",
            "physics": "TEP v0.8 Temporal Topology: screening as a continuous gradient suppression (∇φ), replacing legacy discrete approximations",
            "method": "Compute Temporal Topology φ(r; ρ) from density profile, derive Temporal Shear ∇φ, constrain B(φ) from observation",
            "output": "TEP v0.8 framework showing holonomy driven by path-integrated Temporal Shear",
        },
        "step_012_environmental_dependence": {
            "description": "Environmental dependence and Ambient Symmetry Restoration analysis",
            "physics": "TEP v0.8 continuous gradient suppression: α_eff ∝ D^-1.6 scaling with Galactocentric distance",
            "method": "Temporal Shear scaling relations, Ambient Symmetry Restoration predictions, multi-pulsar replication framework",
            "output": "Testable framework for environment-dependent scalar coupling across Galactic environments",
        },
        "step_013_falsification_criteria": {
            "description": "Falsification criteria automation",
            "physics": "Automated quantitative assessment of TEP falsification criteria",
            "method": "Quantitative threshold evaluation, pass/fail verdict generation, comprehensive falsification report",
            "output": "Automated falsification criteria evaluation strengthening scientific rigor",
        },
        "step_014_synthetic_data_validation": {
            "description": "Synthetic data validation",
            "physics": "Pipeline integrity testing using synthetic data with known properties",
            "method": "Generate synthetic scintillation data with and without TEP effects, run through pipeline to verify detection methodology, test for false positives and false negatives",
            "output": "Synthetic validation results quantifying pipeline integrity and recovery limits",
        },
        "step_015_blind_analysis_validation": {
            "description": "Blind analysis validation and pre-registration confirmation",
            "physics": "Demonstrate that all analysis thresholds were pre-registered before examining data",
            "method": "Verify fixed random seeds, document pre-registered thresholds, confirm no post-hoc adjustment",
            "output": "Blind-analysis audit documenting fixed seeds, thresholds, and post-hoc boundaries",
        },
        "step_016_control_pulsar_analysis": {
            "description": "Control pulsar analysis",
            "physics": "Analyze J1603-7202 control pulsar results and validate opposite-sign holonomy prediction",
            "method": "Comprehensive J1603 control pulsar analysis with closure delays, verify opposite signs to J0437, validate velocity geometry predictions",
            "output": "J1603 geometric diagnostic report; not an independent phase-closure detection",
        },
        "step_017_orientation_algorithm": {
            "description": "Orientation assignment algorithm documentation",
            "physics": "Explicit definition of clockwise vs counter-clockwise triplet orientation and closure-sign conventions in (τ, f_D) space",
            "method": "Generate reproducible mathematical specification used by closure extraction and validation steps",
            "output": "Documented orientation methodology for replication and peer review",
        },
        "step_019_systematic_monte_carlo": {
            "description": "Systematic error Monte Carlo audit",
            "physics": "Quantify how large calibration bias, correlated noise, drift, and selection distortions would need to be to mimic the observed signal",
            "method": "Monte Carlo perturbations of epoch-level summaries and measurement models with fixed seeds",
            "output": "Systematic-error sensitivity bounds complementing primary inference",
        },
        "step_020_replication_readiness": {
            "description": "Replication readiness checklist",
            "physics": "Bind headline statistics to machine-readable reference values independent researchers should recover",
            "method": "Read authoritative Step 003 summaries and emit structured checklist JSON with tolerances",
            "output": "Replication checklist JSON tying code, data paths, and expected ψ / |H| checkpoints",
        },
        "step_021_epoch_level_analysis": {
            "description": "Epoch-level significance analysis",
            "physics": "Address temporal consistency falsification by testing individual epochs rather than only pooled summaries",
            "method": "Bootstrap CI for epoch-|H|, multiple-comparison correction (Bonferroni), hierarchical epoch-effect models where applicable",
            "output": "Epoch-level significance accounting for aggregation and multiplicity",
        },
        "step_022_q4_mechanism_investigation": {
            "description": "Q4 dominance mechanism investigation",
            "physics": "Explain stratification where high-triplet epochs once showed inflated unsigned |H| relative to lower quartiles",
            "method": "Correlations of |H| vs triplet count and SNR, stratified diagnostics, signed-delay checks",
            "output": "Mechanism assessment separating unsigned artifacts from signed phase-domain inference",
        },
        "step_023_snr_correlation_analysis": {
            "description": "SNR correlation analysis",
            "physics": "Test whether holonomy correlates with SNR in a pattern indicative of threshold or precision-driven bias",
            "method": "Epoch- and triplet-level SNR–delay diagnostics with robust correlation summaries",
            "output": "SNR correlation report distinguishing benign precision trends from systematic concern",
        },
        "step_024_multi_pulsar_ingestion": {
            "description": "Multi-pulsar data ingestion",
            "physics": "Ingest data for multiple pulsars (B0329, B0355, B0540)",
            "method": "Download and organize data from multiple pulsar sources",
            "output": "Multi-pulsar dataset ready for joint TEP analysis",
        },
        "step_025_random_triplet_subset_analysis": {
            "description": "Random triplet subset analysis",
            "physics": "Quantify sensitivity of summaries to which triplets enter each epoch after high-count strata are handled carefully",
            "method": "Bootstrap random subsets of triplets within epochs and compare |H| and ψ summaries",
            "output": "Selection robustness with respect to triplet subsampling",
        },
        "step_026_snr_correlation_investigation": {
            "description": "SNR correlation investigation",
            "physics": "Resolve whether prominent SNR-linked outliers reflect physics or estimator/threshold coupling",
            "method": "Regression and stratified SNR tests on signed and unsigned delays; estimator comparisons",
            "output": "Mechanism-level SNR audit supporting interpretation of primary statistics",
        },
        "step_027_alternative_selection_criteria": {
            "description": "Alternative selection criteria analysis",
            "physics": "Stress-test conclusions under alternate caps on triplets per epoch and SNR gates",
            "method": "Replay aggregation under alternate deterministic and stochastic selection rules",
            "output": "Criteria sensitivity matrix for manuscript caveats",
        },
        "step_028_bayesian_hierarchical_model": {
            "description": "Bayesian hierarchical model (alternative implementation)",
            "physics": "Alternative Bayesian model for TEP detection inference",
            "method": "Hierarchical model with conjugate priors, numerical integration",
            "output": "Bayesian posterior distributions for H magnitude and bipolar probability",
        },
        "step_029_jiamusi_analysis": {
            "description": "Jiamusi data ingestion and secondary-spectrum analysis",
            "physics": "Prepare archival Jiamusi pulsar dynamic spectra for closure-delay processing",
            "method": "Download/validate .dat files, parse to arrays, build catalogs, and compute secondary-spectrum diagnostics",
            "output": "Processed Jiamusi spectra and catalog files for downstream closure extraction",
        },
        "step_030_tep_scaling_analysis": {
            "description": "TEP scaling analysis (two-pulsar baseline plus ensemble rows)",
            "physics": "Compare holonomy magnitudes across pulsars against distance, velocity, and screen-geometry scaling expectations",
            "method": "Compute |H| ratios, velocity ratios, and scaling-motivated comparisons using regenerated summaries",
            "output": "Scaling JSON consumed by downstream ensemble and theory steps",
        },
        "step_031_tep_theoretical_predictions": {
            "description": "TEP theoretical predictions",
            "physics": "Generate theoretical TEP predictions for comparison",
            "method": "Compute |H| scaling from TEP framework parameters",
            "output": "Theoretical prediction values for model comparison",
        },
        "step_032_detailed_tep_scaling": {
            "description": "Detailed TEP scaling diagnostics",
            "physics": "Contrast multiple scaling hypotheses using J0437 vs J1603 and extended ensembles where available",
            "method": "Multi-model comparison with regenerated closure summaries and geometry covariates",
            "output": "Detailed scaling diagnostics JSON for manuscript figures",
        },
        "step_033_extract_eta_from_arclets": {
            "description": "η extraction from individual arclets",
            "physics": "Fine-grained arc curvature from arclet-level analysis",
            "method": "Individual arclet η measurement and aggregation",
            "output": "High-resolution η measurements for scaling validation",
        },
        "step_034_extract_eta": {
            "description": "Arc curvature (η) extraction",
            "physics": "Extract arc curvature parameter from secondary spectra",
            "method": "Parabolic fitting to arclet distributions, η computation",
            "output": "η values for TEP scaling analysis",
        },
        "step_035_meerkat_analysis": {
            "description": "MeerKAT telescope analysis",
            "physics": "Dedicated analysis of MeerKAT observational data",
            "method": "MeerKAT-specific calibration and processing pipeline",
            "output": "MeerKAT-optimized closure delay measurements",
        },
        "step_036_synthetic_data_injection": {
            "description": "Synthetic data injection testing",
            "physics": "Inject synthetic TEP signals for pipeline validation",
            "method": "Generate synthetic scintillation with known H, test recovery",
            "output": "Synthetic injection validation quantifying detection sensitivity",
        },
        "step_037_multi_telescope_validation": {
            "description": "Multi-telescope validation audit",
            "physics": "Audit whether current positive phase evidence is independently replicated across telescopes",
            "method": "Compare Parkes/PPTA phase detections against noise-limited Jiamusi bounding rows",
            "output": "Cross-telescope replication status and instrumental-independence caveat",
        },
        "step_038_temporal_evolution": {
            "description": "Temporal evolution analysis",
            "physics": "Analyze TEP signal stability over observation timespan",
            "method": "Time-series analysis |H| across 2008-2018 epochs",
            "output": "Temporal stability assessment and drift analysis",
        },
        "step_039_higher_order_closures": {
            "description": "Higher-order closure analysis",
            "physics": "Test quadrilateral and higher-order closure relations",
            "method": "Compute n-point closure delays beyond triplet baseline",
            "output": "Higher-order closure diagnostic for consistency with triplet evidence",
        },
        "step_040_bootstrap_resampling": {
            "description": "Bootstrap and jackknife resampling analysis",
            "physics": "Obtain robust uncertainty estimates via resampling methods",
            "method": "Bootstrap (1000 iterations) and jackknife leave-one-out resampling of epochs",
            "output": "Bootstrap confidence intervals and jackknife estimates for uncertainty auditing",
        },
        "step_041_selection_bias_analysis": {
            "description": "Selection bias comprehensive analysis",
            "physics": "Characterize selection effects in triplet formation",
            "method": "Monte Carlo tests for selection bias, completeness correction",
            "output": "Selection bias quantification and correction factors",
        },
        "step_042_systematic_errors": {
            "description": "Systematic error budget",
            "physics": "Comprehensive systematic error assessment",
            "method": "Error propagation, instrumental effects, calibration systematics",
            "output": "Complete systematic error budget for TEP measurement",
        },
        "step_043_definitive_chromatic_test": {
            "description": "Chromatic diagnostic (hierarchical model + within-source sub-band)",
            "physics": "Audit frequency dependence while keeping unsigned-|H| fits out of primary inference",
            "method": "Diagnostic model comparison with explicit boundary and calibration caveats, plus within-source J0437 sub-band test",
            "output": "Diagnostic-only chromatic report; not valid as primary evidence",
        },
        "step_018_ensemble_scaling_analysis": {
            "description": "Ensemble Scaling Analysis (Ensemble Overview)",
            "physics": "Testing the TEP scaling law |H| ∝ D × v across 8 independent pulsars",
            "method": "Multi-pulsar correlation analysis with Jiamusi offset correction",
            "output": "Multi-pulsar test of TEP environmental dependence",
        },
        "step_044_probabilistic_weighting": {
            "description": "Probabilistic weighting diagnostic (alternative to hard SNR thresholds)",
            "physics": "Characterize threshold sensitivity of unsigned-|H| diagnostics",
            "method": "Sign-marginalized diagnostic model and inverse-variance weighting, explicitly not primary inference",
            "output": "Diagnostic-only threshold sensitivity report",
        },
        "step_045_synthetic_threshold_degradation": {
            "description": "Synthetic threshold degradation testing",
            "physics": "Test detection robustness to threshold degradation",
            "method": "Synthetic data injection with varying thresholds, sensitivity analysis",
            "output": "Threshold degradation report quantifying robustness and selection sensitivity",
        },
        "step_046_bayesian_orbital_phasebin": {
            "description": "Signed-delay orbital phase-binning analysis",
            "physics": "Test orbital kinematic modulation using Stokes-aligned signed geometric delays",
            "method": "Bin triplets by orbital phase, compute inverse-variance weighted signed mean per bin, fit sinusoid with absolute_sigma=True",
            "output": "Phase-resolved signed-delay modulation testing orbital kinematic predictions",
        },
        "step_047_velocity_direction_controls": {
            "description": "Velocity-direction controls",
            "physics": "Test whether velocity-projection alignment encodes the expected sign",
            "method": "Data-driven velocity-label, scramble, pre-alignment, freeze-record, and wrong-velocity controls",
            "output": "Velocity-direction control report based on Step 003 closure products",
        },
        "step_048_cmb_dipole_frame_analysis": {
            "description": "CMB dipole bulk-vector sensitivity on Stokes-weighted closure",
            "physics": "Test invariance of phase-only Phase Closure and stability of |H| when adding the Planck 2018 kinematic dipole to Step 003 SSB velocities",
            "method": "Reload Step 003 per-epoch triplets; recompute geometric delays; Earth velocity projected onto 3D CMB dipole for annual regressions; wrong-direction and random-direction controls",
            "output": "results/step_048_cmb_dipole_frame_analysis.json",
        },
        "step_049_evidence_ledger": {
            "description": "Evidence ledger, claim hierarchy, and per-pulsar evidence tier table",
            "physics": "Separate primary phase evidence from diagnostics, follow-up targets, and non-claims",
            "method": "Read current result JSON files and write machine-readable claim-status ledger plus step_049_evidence_tier_summary.json",
            "output": "results/step_049_evidence_ledger.json",
        },
    }

    # Print header with scientific context
    context = step_context.get(step_module, {})
    log_message(f"{'=' * 80}", verbose=verbose)
    log_message(f"PIPELINE STEP: {step_name}", verbose=verbose)
    log_message(f"{'=' * 80}", verbose=verbose)

    if context and verbose:
        log_message(f"", verbose=verbose)
        log_message(f"DESCRIPTION:", verbose=verbose)
        log_message(f"    {context.get('description', 'N/A')}", verbose=verbose)
        log_message(f"", verbose=verbose)
        log_message(f"PHYSICS CONTEXT:", verbose=verbose)
        log_message(f"    {context.get('physics', 'N/A')}", verbose=verbose)
        log_message(f"", verbose=verbose)
        log_message(f"METHODOLOGY:", verbose=verbose)
        log_message(f"    {context.get('method', 'N/A')}", verbose=verbose)
        log_message(f"", verbose=verbose)
        log_message(f"EXECUTING...", verbose=verbose)
        log_message(f"", verbose=verbose)

    start_time = time.time()

    # Per-step detailed log (same directory as pipeline.log; do not duplicate via tee)
    step_log_file = LOGS_DIR / f"{step_module}.log"
    step_logger = TEPLogger(step_module, str(step_log_file), verbose=True)
    set_step_logger(step_logger)
    step_logger.info(f"{'=' * 80}")
    step_logger.info(f"STEP {step_module} STARTED")
    step_logger.info(f"Description: {context.get('description', 'N/A')}")
    step_logger.info(f"Physics: {context.get('physics', 'N/A')}")
    step_logger.info(f"Method: {context.get('method', 'N/A')}")
    step_logger.info(f"{'=' * 80}")

    # Save and clear argv to prevent argument conflicts
    original_argv = sys.argv
    argv0 = sys.argv[0] if sys.argv else "run_pipeline.py"
    skip_scintools = os.environ.get("TEP_SKIP_SCINTOOLS", "").strip().lower() in (
        "1",
        "true",
        "yes",
    )
    if step_module == "step_000_data_ingestion" and skip_scintools:
        sys.argv = [argv0, "--skip-scintools"]
    else:
        sys.argv = [argv0]

    try:
        # Dynamic import
        import importlib
        module_name = f"scripts.steps.{step_module}"
        try:
            step_module_obj = importlib.import_module(module_name)
        except ImportError:
            try:
                module_name = f"{step_module}"
                step_module_obj = importlib.import_module(module_name)
            except ImportError as e:
                log_message(f"Could not import {step_module}: {e}", "ERROR", verbose)
                return False

        # Determine entry point
        step_main = None
        if hasattr(step_module_obj, "step_main"):
            step_main = step_module_obj.step_main
        elif hasattr(step_module_obj, "main"):
            step_main = step_module_obj.main
        else:
            # Check for other known step functions
            for func_name in ["calculate_tep_scaling_ratio", "comprehensive_scaling_analysis"]:
                if hasattr(step_module_obj, func_name):
                    step_main = getattr(step_module_obj, func_name)
                    break
            
            if not step_main:
                log_message(f"Module {step_module} has no valid entry point (step_main or main)", "ERROR", verbose)
                return False
        # Execute the step
        import inspect
        sig = inspect.signature(step_main)
        
        # Handle different signature types
        kwargs: Dict[str, Any] = {}
        if "logger" in sig.parameters and "verbose" in sig.parameters:
            kwargs["logger"] = step_logger
            kwargs["verbose"] = verbose
        if "force" in sig.parameters:
            kwargs["force"] = force
        if kwargs:
            success = step_main(**kwargs)
        else:
            success = step_main()

        # If step_main returns None, assume success unless exception raised
        if success is None:
            success = True

        elapsed = time.time() - start_time

        # Success message with scientific interpretation - log to both main and step logger
        if success:
            log_message(f"{'=' * 80}", verbose=verbose)
            log_message(
                f"STEP COMPLETED SUCCESSFULLY ({elapsed:.1f} seconds)", verbose=verbose
            )
            log_message(f"{'=' * 80}", verbose=verbose)

            step_logger.info(f"{'=' * 80}")
            step_logger.info(f"STEP COMPLETED SUCCESSFULLY ({elapsed:.1f} seconds)")
            step_logger.info(f"{'=' * 80}")
        else:
            log_message(f"{'=' * 80}", verbose=verbose)
            log_message(
                f"STEP COMPLETED WITH WARNINGS ({elapsed:.1f} seconds)", "WARNING", verbose=verbose
            )
            log_message(f"{'=' * 80}", verbose=verbose)

            step_logger.warning(f"{'=' * 80}")
            step_logger.warning(f"STEP COMPLETED WITH WARNINGS ({elapsed:.1f} seconds)")
            step_logger.warning(f"{'=' * 80}")

        if step_module == "step_003_closure_delays_final" and verbose:
            log_message(f"", verbose=verbose)
            log_message(f"ACTUAL RESULTS:", verbose=verbose)
            try:
                import json
                summary_path = Path("results/step_003_closure_final_summary_j0437.json")
                if not summary_path.exists():
                    summary_path = Path("results/step_003_closure_final_summary.json")
                with open(summary_path) as f:
                    results = json.load(f)
                psi = results.get("phase_closure_mean_rad")
                psi_se = results.get("phase_closure_circ_se_rad")
                rayleigh_p = results.get("phase_closure_rayleigh_p")
                h_mag = results.get("H_magnitude_ns")
                h_trim = results.get("H_trim_magnitude_ns")
                n_epochs = results.get("n_epochs")
                n_triplets = results.get("n_total_triplets")

                log_message(f"    Summary file: {summary_path}", verbose=verbose)
                if isinstance(psi, (int, float)) and isinstance(psi_se, (int, float)):
                    log_message(
                        f"    phase_closure_mean_rad: {psi:.4f} ± {psi_se:.4f} rad",
                        verbose=verbose,
                    )
                if isinstance(rayleigh_p, (int, float)):
                    log_message(f"    phase_closure_rayleigh_p: {rayleigh_p:.2e}", verbose=verbose)
                log_message(
                    f"    H_magnitude_ns: {h_mag:.4f} ns" if isinstance(h_mag, (int, float)) else f"    H_magnitude_ns: {h_mag}",
                    verbose=verbose,
                )
                log_message(
                    f"    H_trim_magnitude_ns: {h_trim:.4f} ns" if isinstance(h_trim, (int, float)) else f"    H_trim_magnitude_ns: {h_trim}",
                    verbose=verbose,
                )
                log_message(f"    n_epochs: {n_epochs}", verbose=verbose)
                log_message(f"    n_triplets: {n_triplets}", verbose=verbose)
            except Exception as e:
                log_message(f"    Error reading results: {e}", verbose=verbose)
        elif step_module == "step_004_verification" and verbose:
            log_message(f"", verbose=verbose)
            log_message(f"VERIFICATION SUMMARY:", verbose=verbose)
            log_message(f"    All verification tests have completed.", verbose=verbose)
            log_message(f"    Results are available in:", verbose=verbose)
            log_message(
                f"        results/independent_verification_results.json",
                verbose=verbose,
            )
            log_message(f"        results/blind_analysis_results.json", verbose=verbose)
            log_message(
                f"        results/alternative_ism_tests_results.json", verbose=verbose
            )
        elif step_module == "step_005_enhanced_validation" and verbose:
            log_message(f"", verbose=verbose)
            log_message(f"ENHANCED VALIDATION SUMMARY:", verbose=verbose)
            log_message(
                f"    All enhanced validation analyses have completed.", verbose=verbose
            )
            log_message(
                f"    These address reviewer concerns and strengthen the detection.",
                verbose=verbose,
            )
            log_message(f"    Results are available in:", verbose=verbose)
            log_message(f"        results/improved_gr_null_tests.json", verbose=verbose)
            log_message(f"        results/advanced_ism_tests.json", verbose=verbose)
            log_message(
                f"        results/blind_analysis_both_signs.json", verbose=verbose
            )
            log_message(
                f"        results/rigorous_anisotropy_derivation.json", verbose=verbose
            )
            log_message(
                f"        results/weighted_analysis_verification.json", verbose=verbose
            )
            log_message(
                f"        results/tep_prediction_framework.json", verbose=verbose
            )

        return success

    except Exception as e:
        elapsed = time.time() - start_time
        log_message(f"{'=' * 80}", verbose=verbose)
        log_message(f"STEP FAILED after {elapsed:.1f}s: {e}", "ERROR", verbose)
        log_message(f"{'=' * 80}", verbose=verbose)

        # Also log error to step logger
        step_logger.error(f"{'=' * 80}")
        step_logger.error(f"STEP FAILED after {elapsed:.1f}s: {e}")
        step_logger.error(f"{'=' * 80}")

        import traceback

        tb_str = traceback.format_exc()
        log_message(tb_str, "ERROR", verbose=False)
        step_logger.error(tb_str)
        return False
    finally:
        # Restore original argv
        sys.argv = original_argv
        # Clear step logger to prevent cross-step contamination
        # Fixed: Pass valid logger or handle None properly
        logger = TEPLogger("default")
        set_step_logger(logger)


def run_pipeline(
    steps: Optional[list[str]] = None,
    verbose: bool = True,
    pipeline_log: Optional[Path] = None,
    force: bool = False,
):
    """
    Execute the complete TEP-J0437 analysis pipeline.

    This function orchestrates the full workflow from data ingestion through
    TEP detection, with comprehensive logging and error handling.

    Parameters
    -----------
    steps : list of str, optional
        Specific steps to run (e.g., ['000', '001']). Default runs all.
    verbose : bool
        Enable detailed logging output
    pipeline_log : Path, optional
        Orchestration log file (timestamps, step boundaries, summary). Defaults to
        logs/pipeline.log under the project root. Step-level detail goes only to
        logs/<step_module>.log (e.g. logs/step_000_data_ingestion.log). Avoid wrapping
        this process in tee to a second file; that duplicates stdout and is unnecessary.
    """
    global _PIPELINE_LOG_PATH

    if pipeline_log is not None:
        pl = Path(pipeline_log)
        _PIPELINE_LOG_PATH = pl.resolve() if pl.is_absolute() else (ROOT / pl).resolve()
    _PIPELINE_LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    # Fresh orchestration log for this run
    _PIPELINE_LOG_PATH.unlink(missing_ok=True)

    log_message("=" * 70, verbose=verbose)
    log_message(
        "TEP-J0437: TEMPORAL EQUIVALENCE PRINCIPLE DETECTION PIPELINE", verbose=verbose
    )
    log_message("=" * 70, verbose=verbose)
    log_message(
        f"Execution started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        verbose=verbose,
    )
    log_message(f"Working directory: {ROOT}", verbose=verbose)
    log_message(f"CPU cores available: {mp.cpu_count()}", verbose=verbose)

    # Define complete pipeline (49 steps)
    all_steps = [
        ("Step 000 - Data Ingestion (Scintools ATNF)", "step_000_data_ingestion"),
        ("Step 001 - Parse Dynamic Spectra", "step_001_parse_dynspec"),
        ("Step 002 - Secondary Spectrum Generation", "step_002_secondary_spectra"),
        ("Step 003 - Closure Delay Extraction", "step_003_closure_delays_final"),
        ("Step 004 - Verification & Validation", "step_004_verification"),
        ("Step 005 - Enhanced Validation", "step_005_enhanced_validation"),
        ("Step 006 - Advanced Validation", "step_006_advanced_validation"),
        ("Step 007 - Independent Statistical Validation", "step_007_independent_validation"),
        ("Step 008 - Alternative Explanations Testing", "step_008_alternative_explanations"),
        ("Step 009 - Parameter Sensitivity Analysis", "step_009_parameter_sensitivity"),
        ("Step 010 - Comprehensive Data Quality Metrics", "step_010_data_quality_metrics"),
        ("Step 011 - ISM Density Profile Modeling", "step_011_ism_density_modeling"),
        ("Step 012 - Environmental Dependence Analysis", "step_012_environmental_dependence"),
        ("Step 014 - Synthetic Data Validation", "step_014_synthetic_data_validation"),
        ("Step 015 - Blind Analysis Validation", "step_015_blind_analysis_validation"),
        ("Step 016 - Control Pulsar Analysis", "step_016_control_pulsar_analysis"),
        ("Step 017 - Orientation Algorithm Documentation", "step_017_orientation_algorithm"),
        ("Step 019 - Systematic Monte Carlo", "step_019_systematic_monte_carlo"),
        ("Step 020 - Replication Readiness", "step_020_replication_readiness"),
        ("Step 021 - Epoch-Level Significance Analysis", "step_021_epoch_level_analysis"),
        ("Step 022 - Q4 Dominance Mechanism Investigation", "step_022_q4_mechanism_investigation"),
        ("Step 023 - SNR Correlation Analysis", "step_023_snr_correlation_analysis"),
        ("Step 024 - Multi-Pulsar Ingestion", "step_024_multi_pulsar_ingestion"),
        ("Step 025 - Random Triplet Subset Analysis", "step_025_random_triplet_subset_analysis"),
        ("Step 026 - SNR Correlation Investigation", "step_026_snr_correlation_investigation"),
        ("Step 027 - Alternative Selection Criteria", "step_027_alternative_selection_criteria"),
        ("Step 028 - Bayesian Hierarchical Model", "step_028_bayesian_hierarchical_model"),
        ("Step 029 - Jiamusi Analysis", "step_029_jiamusi_analysis"),
        ("Step 030 - TEP Scaling Analysis", "step_030_tep_scaling_analysis"),
        ("Step 018 - Ensemble Scaling Analysis", "step_018_ensemble_scaling_analysis"),
        ("Step 031 - TEP Theoretical Predictions", "step_031_tep_theoretical_predictions"),
        ("Step 032 - Detailed TEP Scaling", "step_032_detailed_tep_scaling"),
        ("Step 033 - Extract Eta from Arclets", "step_033_extract_eta_from_arclets"),
        ("Step 034 - Extract Eta", "step_034_extract_eta"),
        ("Step 035 - MeerKAT Analysis", "step_035_meerkat_analysis"),
        ("Step 036 - Synthetic Data Injection", "step_036_synthetic_data_injection"),
        ("Step 037 - Multi-Telescope Validation", "step_037_multi_telescope_validation"),
        ("Step 038 - Temporal Evolution", "step_038_temporal_evolution"),
        ("Step 039 - Higher Order Closures", "step_039_higher_order_closures"),
        ("Step 040 - Bootstrap Resampling", "step_040_bootstrap_resampling"),
        ("Step 041 - Selection Bias Analysis", "step_041_selection_bias_analysis"),
        ("Step 042 - Systematic Errors", "step_042_systematic_errors"),
        ("Step 043 - Chromatic Diagnostic", "step_043_definitive_chromatic_test"),
        ("Step 044 - Probabilistic Weighting", "step_044_probabilistic_weighting"),
        ("Step 045 - Synthetic Threshold Degradation", "step_045_synthetic_threshold_degradation"),
        ("Step 046 - Bayesian Orbital Phase-Binning (signed delays)", "step_046_bayesian_orbital_phasebin"),
        ("Step 047 - Velocity Direction Controls", "step_047_velocity_direction_controls"),
        ("Step 048 - CMB Dipole Frame Analysis", "step_048_cmb_dipole_frame_analysis"),
        ("Step 013 - Falsification Criteria Automation", "step_013_falsification_criteria"),
        ("Step 049 - Evidence Ledger", "step_049_evidence_ledger"),
    ]

    # Filter steps if specified
    if steps:
        step_numbers = [s.zfill(3) for s in steps]
        pipeline_steps = [
            (name, module)
            for name, module in all_steps
            if any(s in name for s in step_numbers)
        ]
    else:
        pipeline_steps = all_steps

    log_message(
        f"\nPipeline configuration: {len(pipeline_steps)} steps", verbose=verbose
    )
    for name, _ in pipeline_steps:
        log_message(f"  • {name}", verbose=verbose)

    # Execute pipeline
    results = []
    total_start = time.time()

    for step_name, step_module in pipeline_steps:
        success = run_step(step_name, step_module, verbose=verbose, force=force)
        results.append((step_name, success))

        if not success:
            log_message("\n" + "=" * 70, verbose=verbose)
            log_message("PIPELINE HALTED DUE TO STEP FAILURE", "ERROR", verbose=verbose)
            log_message("=" * 70, verbose=verbose)
            break

    # Generate summary
    total_elapsed = time.time() - total_start

    log_message("\n" + "=" * 70, verbose=verbose)
    log_message("PIPELINE EXECUTION SUMMARY", verbose=verbose)
    log_message("=" * 70, verbose=verbose)

    for step_name, success in results:
        status = "[PASS]" if success else "[FAIL]"
        log_message(f"{status} - {step_name}", verbose=verbose)

    log_message(f"\nTotal execution time: {total_elapsed:.1f} seconds", verbose=verbose)
    log_message(f"Orchestration log: {_PIPELINE_LOG_PATH}", verbose=verbose)
    log_message(
        f"Per-step detail logs: {LOGS_DIR}/<step_module>.log (e.g. {LOGS_DIR}/step_001_parse_dynspec.log)",
        verbose=verbose,
    )

    # Final verdict
    all_passed = all(success for _, success in results)

    log_message("\n" + "=" * 70, verbose=verbose)
    if all_passed:
        log_message("PIPELINE COMPLETED SUCCESSFULLY", "SUCCESS", verbose=verbose)
        log_message(
            "Results available in data/secondary/ and results/", verbose=verbose
        )
        log_message(
            "Review results/step_003_closure_final_summary.json for statistical verdict.",
            verbose=verbose,
        )
    else:
        log_message("[FAIL] PIPELINE COMPLETED WITH ERRORS", "WARNING", verbose=verbose)
        log_message("Check logs for detailed error information", verbose=verbose)
    log_message("=" * 70, verbose=verbose)

    return all_passed


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="TEP-J0437: Complete TEP Detection Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
EXAMPLES:
  # Run complete pipeline (recommended)
  python scripts/run_pipeline.py

  # Run specific steps only
  python scripts/run_pipeline.py --steps 000,001,002

  # Quiet mode (minimal output)
  python scripts/run_pipeline.py --quiet

  # Skip data ingestion (use cached data)
  python scripts/run_pipeline.py --steps 001,002,003

  # Long run in background (one orchestration log; do not add tee — it duplicates)
  nohup env TEP_SKIP_SCINTOOLS=1 PYTHONPATH=. python3 scripts/run_pipeline.py &

LOGGING (by design — no third merged copy):
  Orchestration only (step list, boundaries, pass/fail, summary):
    logs/pipeline.log   (override: --pipeline-log or TEP_PIPELINE_LOG)
  Detailed stdout from each step (CALC/DATA/…):
    logs/<step_module>.log   (e.g. logs/step_002_secondary_spectra.log)

For scientific reproducibility, the complete pipeline downloads data,
processes all epochs, and generates TEP detection results automatically.

ENVIRONMENT:
  TEP_SKIP_SCINTOOLS=1   When set to 1, true, or yes, step 000 is invoked with
                         --skip-scintools (CSIRO verify still runs; Scintools
                         supplement is skipped). Use when the Scintools host
                         is unreachable; otherwise omit so failures surface.
  TEP_PIPELINE_LOG      Optional path to the orchestration log (same as --pipeline-log
                         if the CLI flag is not passed). Relative paths are under
                         the project root.
        """,
    )

    parser.add_argument(
        "--steps",
        "-s",
        type=str,
        help="Comma-separated step numbers (000-003). Default: all steps",
    )

    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        default=True,
        help="Enable verbose output with scientific logging",
    )

    parser.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="Minimal output mode (overrides --verbose)",
    )

    parser.add_argument(
        "--pipeline-log",
        type=str,
        default=None,
        metavar="PATH",
        help="Orchestration log file (default: logs/pipeline.log, or TEP_PIPELINE_LOG)",
    )

    parser.add_argument(
        "--force",
        "-f",
        action="store_true",
        help="Force re-run all steps (do not skip cached outputs)",
    )

    args = parser.parse_args()

    # Parse step selection
    step_list = None
    if args.steps:
        step_list = [s.strip() for s in args.steps.split(",")]

    # Determine verbosity
    verbose = not args.quiet and args.verbose

    raw_pl = args.pipeline_log or os.environ.get("TEP_PIPELINE_LOG", "").strip()
    if raw_pl:
        pl_path = Path(raw_pl)
        pipeline_log_path = pl_path.resolve() if pl_path.is_absolute() else (ROOT / pl_path).resolve()
    else:
        pipeline_log_path = (ROOT / "logs" / "pipeline.log").resolve()

    # Execute pipeline
    success = run_pipeline(steps=step_list, verbose=verbose, pipeline_log=pipeline_log_path, force=args.force)

    # Exit with appropriate code
    sys.exit(0 if success else 1)
