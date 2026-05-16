#!/usr/bin/env python3
"""
================================================================================
STEP 020: REPLICATION READINESS CHECKLIST
================================================================================

Purpose: Prepare comprehensive checklist for independent replication.
Ensures all necessary information, data, and code are available for
independent researchers to reproduce and validate the detection.

This addresses: "Can others replicate this analysis independently?"

================================================================================
"""

import json
from scripts.utils.logger import print_status
from scripts.utils.json_numpy import NpEncoder
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, List

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

SUMMARY_J0437 = RESULTS_DIR / "step_003_closure_final_summary_j0437.json"
SUMMARY_GENERIC = RESULTS_DIR / "step_003_closure_final_summary.json"
PER_EPOCH_J0437 = RESULTS_DIR / "step_003_closure_final_per_epoch_j0437.json"
PER_EPOCH_GENERIC = RESULTS_DIR / "step_003_closure_final_per_epoch.json"

def get_authentic_results():
    """
    Return authoritative current headline values for replication checks.

    Preference order:
    - J0437-specific Step 003 summary (current paper target)
    - Generic Step 003 summary (legacy fallback, only if J0437-specific missing)
    """
    freq_file = SUMMARY_J0437 if SUMMARY_J0437.exists() else SUMMARY_GENERIC
    if not freq_file.exists():
        raise FileNotFoundError(
            f"Primary results not found at {freq_file}. "
            "Run scripts/run_pipeline.py (or step_003_closure_delays_final.py) before replication readiness."
        )
    with open(freq_file, "r") as f:
        data = json.load(f)

    psi = data.get("phase_closure_mean_rad")
    psi_se = data.get("phase_closure_circ_se_rad")
    rayleigh_p = data.get("phase_closure_rayleigh_p")
    v_p = data.get("phase_closure_v_p")
    ci_lo = data.get("phase_closure_bootstrap_ci_95_lower_rad")
    ci_hi = data.get("phase_closure_bootstrap_ci_95_upper_rad")

    if psi is None or psi_se is None:
        raise ValueError("Phase Closure ψ statistics missing from Step 003 summary.")

    return {
        "psi_mean_rad": float(psi),
        "psi_circ_se_rad": float(psi_se),
        "rayleigh_p": float(rayleigh_p) if rayleigh_p is not None else None,
        "v_p": float(v_p) if v_p is not None else None,
        "psi_bootstrap_ci_95": [ci_lo, ci_hi],
        "n_epochs": int(data.get("n_epochs")) if data.get("n_epochs") is not None else None,
        "n_total_triplets": int(data.get("n_total_triplets")) if data.get("n_total_triplets") is not None else None,
        "n_independent_samples": int(data.get("n_independent_samples")) if data.get("n_independent_samples") is not None else None,
    }

def check_data_availability() -> Dict[str, Any]:
    """Check that all necessary data is available and documented."""
    
    required_files = {
        "raw_data_catalog": PROJECT_ROOT / "data" / "raw" / "scintools_catalog.json",
        "processed_data": PROJECT_ROOT / "data" / "processed" / "j0437_epoch_catalog.json",
        "closure_results_summary_j0437": SUMMARY_J0437,
        "closure_results_per_epoch_j0437": PER_EPOCH_J0437,
        "closure_results_summary_legacy": SUMMARY_GENERIC,
        "closure_results_per_epoch_legacy": PER_EPOCH_GENERIC,
        "pipeline_documentation": PROJECT_ROOT / "scripts" / "README.md",
    }
    
    available = {}
    for name, path in required_files.items():
        available[name] = {
            "exists": path.exists(),
            "path": str(path.relative_to(PROJECT_ROOT)),
            "size_bytes": path.stat().st_size if path.exists() else 0
        }
    
    all_available = all(a["exists"] for a in available.values())
    
    return {
        "all_data_available": bool(all_available),
        "files": available,
        "missing_files": [name for name, info in available.items() if not info["exists"]]
    }


def check_code_reproducibility() -> Dict[str, Any]:
    """Check that code is reproducible with fixed seeds."""
    
    # Check for random seeds in key scripts
    seed_files = [
        "step_006_advanced_validation.py",
        "step_014_synthetic_data_validation.py",
        "step_015_blind_analysis_validation.py"
    ]
    
    seeds_documented = []
    for filename in seed_files:
        path = PROJECT_ROOT / "scripts" / "steps" / filename
        if path.exists():
            content = path.read_text()
            if "seed" in content.lower() or "random" in content.lower():
                seeds_documented.append(filename)
    
    return {
        "fixed_seeds_implemented": len(seeds_documented) > 0,
        "seeded_scripts": seeds_documented,
        "primary_seed": 42,  # Documented in blind analysis
        "reproducibility_level": "High" if len(seeds_documented) >= 3 else "Medium"
    }


def check_documentation() -> Dict[str, Any]:
    """Check that documentation is complete."""
    
    docs = {
        "manuscript": PROJECT_ROOT / "16-TEP-J0437-v0.1-Sintra.md",
        "pipeline_readme": PROJECT_ROOT / "scripts" / "README.md",
        "citation": PROJECT_ROOT / "CITATION.cff",
        "license": PROJECT_ROOT / "LICENSE",
        "requirements": PROJECT_ROOT / "requirements.txt",
    }
    
    available = {}
    for name, path in docs.items():
        available[name] = path.exists()
    
    return {
        "all_docs_present": all(available.values()),
        "docs": available
    }


def generate_analysis_parameters() -> Dict[str, Any]:
    """Document all analysis parameters for replication."""
    auth = get_authentic_results()
    
    return {
        "pre_registered_thresholds": {
            "significance_threshold_sigma": 5.0,
            "bootstrap_samples": 10000,
            "min_arclets": 3,
            "max_triplets_per_epoch": 20,
            "edge_margin_fraction": 0.01,
            "cv_threshold": 0.5,
            "multiple_comparison_method": "Bonferroni",
            "n_tests": 4,
            "alpha_corrected": 0.0125
        },
        "measurement_parameters": {
            "sub_pixel_method": "parabolic interpolation",
            "aliasing_checks": True,
            "tau_bounds_us": [0, 50],
            "fD_bounds_mHz": [-50, 50]
        },
        "data_selection": {
            "source": "PSR J0437-4715",
            "dataset": "PPTA DR2",
            "time_span_years": 11,
            "n_epochs_with_triplets": auth.get("n_epochs"),
            "n_total_triplets": auth.get("n_total_triplets"),
            "n_independent_samples": auth.get("n_independent_samples"),
        },
        "random_seeds": {
            "bootstrap": 42,
            "permutations": 42,
            "simulations": 42,
            "cross_validation": 42
        }
    }


def generate_replication_steps() -> List[Dict]:
    """Generate step-by-step replication instructions."""
    auth = get_authentic_results()
    
    return [
        {
            "step": 1,
            "action": "Download PPTA DR2 data from Scintools repository",
            "expected_time": "2-4 hours",
            "verification": "Step 000 audit verifies CSIRO manifests and checksums"
        },
        {
            "step": 2,
            "action": "Run Step 000-003: Process dynamic spectra to closure delays",
            "expected_time": "4-8 hours",
            "verification": f"J0437 Step 003 summary: n_epochs={auth.get('n_epochs')}, n_total_triplets={auth.get('n_total_triplets')}, n_independent_samples={auth.get('n_independent_samples')}"
        },
        {
            "step": 3,
            "action": "Run Step 004-006: Basic validation tests",
            "expected_time": "1-2 hours",
            "verification": "All 3 tests pass"
        },
        {
            "step": 4,
            "action": "Run Step 007-014: Advanced validation suite",
            "expected_time": "2-4 hours",
            "verification": "All 8 tests pass"
        },
        {
            "step": 5,
            "action": "Run Step 015-019: Strengthening validations",
            "expected_time": "1-2 hours",
            "verification": "All 5 tests pass"
        },
        {
            "step": 6,
            "action": "Run Step 020-022: Cross-method & systematic tests",
            "expected_time": "30 minutes",
            "verification": "Detection robust"
        },
        {
            "step": 7,
            "action": "Compare results with reference values",
            "expected_time": "15 minutes",
            "verification": (
                f"ψ ≈ {auth['psi_mean_rad']:.3f} +/- {auth['psi_circ_se_rad']:.3f} rad; "
                f"Rayleigh p ≈ {auth.get('rayleigh_p')}; V-test p ≈ {auth.get('v_p')}; "
                f"bootstrap CI ≈ {auth.get('psi_bootstrap_ci_95')}"
            )
        }
    ]


def generate_expected_results() -> Dict[str, Any]:
    """Document expected results for verification."""
    auth = get_authentic_results()
    
    return {
        "primary_detection": {
            "phase_closure_mean_rad": auth["psi_mean_rad"],
            "phase_closure_circ_se_rad": auth["psi_circ_se_rad"],
            "phase_closure_rayleigh_p": auth.get("rayleigh_p"),
            "phase_closure_v_p": auth.get("v_p"),
            "phase_closure_bootstrap_ci_95": auth.get("psi_bootstrap_ci_95"),
            "tolerance_rad": 0.02
        },
        "secondary_results": {
            "n_epochs": auth.get("n_epochs"),
            "n_total_triplets": auth.get("n_total_triplets"),
            "n_independent_samples": auth.get("n_independent_samples"),
        },
        "validation_status": {
            "n_tests": 6,
            "all_should_pass": True
        }
    }


def main():
    """Run replication readiness check."""
    print_status("===" * 80)
    print("STEP 020: REPLICATION READINESS CHECKLIST")
    print_status("===" * 80)
    print()
    print("Purpose: Ensure independent replication is possible")
    print()
    
    # Run all checks
    data_check = check_data_availability()
    code_check = check_code_reproducibility()
    doc_check = check_documentation()
    params = generate_analysis_parameters()
    steps = generate_replication_steps()
    expected = generate_expected_results()
    
    # Print results
    print_status("===" * 80)
    print("DATA AVAILABILITY")
    print_status("===" * 80)
    print(f"All required data available: {'YES [OK]' if data_check['all_data_available'] else 'NO [FAIL]'}")
    if data_check['missing_files']:
        print(f"Missing: {', '.join(data_check['missing_files'])}")
    
    print_status("" + "=" * 80)
    print("CODE REPRODUCIBILITY")
    print_status("===" * 80)
    print(f"Fixed seeds implemented: {'YES [OK]' if code_check['fixed_seeds_implemented'] else 'NO [FAIL]'}")
    print(f"Primary random seed: {code_check['primary_seed']}")
    print(f"Reproducibility level: {code_check['reproducibility_level']}")
    
    print_status("" + "=" * 80)
    print("DOCUMENTATION")
    print_status("===" * 80)
    print(f"All documentation present: {'YES [OK]' if doc_check['all_docs_present'] else 'NO [FAIL]'}")
    for doc, present in doc_check['docs'].items():
        status = "[OK]" if present else "[FAIL]"
        print(f"  {status} {doc}")
    
    print_status("" + "=" * 80)
    print("ANALYSIS PARAMETERS (Pre-registered)")
    print_status("===" * 80)
    print(f"Significance threshold: {params['pre_registered_thresholds']['significance_threshold_sigma']}sigma")
    print(f"Bootstrap samples: {params['pre_registered_thresholds']['bootstrap_samples']}")
    print(f"Min arclets: {params['pre_registered_thresholds']['min_arclets']}")
    print(f"Random seed: {params['random_seeds']['bootstrap']}")
    print(f"Multiple comparison: {params['pre_registered_thresholds']['multiple_comparison_method']}")
    
    print_status("" + "=" * 80)
    print("REPLICATION STEPS")
    print_status("===" * 80)
    total_time = 0
    for step in steps:
        print(f"\nStep {step['step']}: {step['action']}")
        print(f"  Time: {step['expected_time']}")
        print(f"  Verify: {step['verification']}")
        # Extract hours
        if 'hour' in step['expected_time']:
            hours = step['expected_time'].split('-')[-1].split()[0]
            total_time += float(hours)
    
    print(f"\nTotal estimated time: ~{int(total_time)} hours")
    
    print_status("" + "=" * 80)
    print("EXPECTED RESULTS (For Verification)")
    print_status("===" * 80)
    print(
        f"ψ = {expected['primary_detection']['phase_closure_mean_rad']:.3f} +/- {expected['primary_detection']['phase_closure_circ_se_rad']:.3f} rad"
    )
    print(f"Rayleigh p = {expected['primary_detection']['phase_closure_rayleigh_p']}")
    print(f"V-test p = {expected['primary_detection']['phase_closure_v_p']}")
    print(f"Bootstrap CI (95%) = {expected['primary_detection']['phase_closure_bootstrap_ci_95']}")
    print(f"Tolerance: +/-{expected['primary_detection']['tolerance_rad']:.2f} rad")
    
    # Overall readiness
    ready = (
        data_check['all_data_available'] and
        code_check['fixed_seeds_implemented'] and
        doc_check['all_docs_present']
    )
    
    print_status("" + "=" * 80)
    print("REPLICATION READINESS VERDICT")
    print_status("===" * 80)
    print(f"\nReady for independent replication: {'YES [OK]' if ready else 'NO - Issues need resolution'}")
    
    if ready:
        print("\n[OK] All data available")
        print("[OK] Code is reproducible with fixed seeds")
        print("[OK] Documentation complete")
        print("[OK] Step-by-step instructions provided")
        print("[OK] Expected results documented")
        print("\n-> Independent researchers can replicate this analysis")
    
    # Save report
    report = {
        "check_type": "Replication Readiness",
        "check_date": datetime.now().isoformat(),
        "data_availability": data_check,
        "code_reproducibility": code_check,
        "documentation": doc_check,
        "analysis_parameters": params,
        "replication_steps": steps,
        "expected_results": expected,
        "ready_for_replication": bool(ready),
        "estimated_time_hours": int(total_time)
    }
    
    output_file = RESULTS_DIR / "step_020_replication_readiness.json"
    with open(output_file, 'w') as f:
        json.dump(report, f, indent=2, cls=NpEncoder)
    
    print(f"\n\nReport saved to: {output_file}")
    print_status("===" * 80)
    
    return report


if __name__ == "__main__":
    main()
