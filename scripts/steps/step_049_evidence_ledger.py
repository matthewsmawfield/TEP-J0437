#!/usr/bin/env python3
"""
STEP 049: EVIDENCE LEDGER

Create a compact claim-status ledger from the current pipeline outputs.
The goal is to make the manuscript defensible on its own terms: separate the
primary phase-closure detection from supporting diagnostics and follow-up tests.
"""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.json_numpy import NpEncoder
from scripts.utils.logger import print_status

RESULTS_DIR = PROJECT_ROOT / "results"


def load_json(name: str, *, required: bool = False) -> Dict[str, Any]:
    path = RESULTS_DIR / name
    if not path.exists():
        if required:
            raise FileNotFoundError(
                f"Required pipeline result missing: {path}. Run the upstream step first."
            )
        return {}
    with open(path) as fh:
        return json.load(fh)


def _load_json_path(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    with open(path) as fh:
        return json.load(fh)


def _discover_step003_summaries_by_pulsar() -> Dict[str, Path]:
    """Map ATNF-style pulsar name -> path of step_003 closure summary JSON."""
    mapping: Dict[str, Path] = {}
    for path in sorted(RESULTS_DIR.glob("step_003_closure_final_summary*.json")):
        if path.name == "step_003_closure_final_summary.json":
            continue
        if "_j0437_sb" in path.name:
            continue
        data = _load_json_path(path)
        pulsar = data.get("pulsar")
        if isinstance(pulsar, str) and pulsar:
            mapping[pulsar] = path
    return mapping


def _tier_for_pulsar(name: str, summary: Dict[str, Any], ingestion_role: str) -> str:
    if name == "J0437-4715":
        return "primary_phase_closure"
    if name == "J1603-7202":
        return "geometric_diagnostic_phase_noise_limited"
    if "MeerKAT" in ingestion_role or "meerkat" in ingestion_role.lower():
        return "catalog_noise_limited_telescope_auxiliary"
    if "Ensemble" in ingestion_role or "Jiamusi" in ingestion_role:
        return "catalog_noise_limited_bounding_row"
    return "catalog_noise_limited_bounding_row"


def _counts_toward_independent_phase_replication(name: str, summary: Dict[str, Any]) -> bool:
    if not summary:
        return False
    return bool(summary.get("detected_3sigma")) and name != "J0437-4715"


def build_evidence_tier_summary(
    scaling: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Per-pulsar evidence tiering for epistemic accounting (N=1 phase detection vs bounds vs geometry).
    Values are read from step_003 summaries and step_024 ingestion metadata where available.
    """
    ingestion = load_json("step_024_multi_pulsar_ingestion.json")
    summary_paths = _discover_step003_summaries_by_pulsar()
    rows: List[Dict[str, Any]] = []
    for entry in ingestion.get("pulsars_analyzed", []):
        name = entry.get("name")
        if not isinstance(name, str):
            continue
        role = entry.get("role") or ""
        path = summary_paths.get(name)
        summary = _load_json_path(path) if path else {}
        tier = _tier_for_pulsar(name, summary, role)
        row = {
            "pulsar": name,
            "evidence_tier": tier,
            "ingestion_role": role,
            "data_source": entry.get("source"),
            "step003_summary_relpath": str(path.relative_to(PROJECT_ROOT)) if path else None,
            "n_epochs": summary.get("n_epochs") if summary else entry.get("n_epochs"),
            "n_total_triplets": summary.get("n_total_triplets")
            if summary
            else entry.get("n_triplets"),
            "n_independent_samples": summary.get("n_independent_samples"),
            "phase_closure_rayleigh_p": summary.get("phase_closure_rayleigh_p"),
            "phase_closure_rbar": summary.get("phase_closure_rbar"),
            "detected_3sigma_phase_closure": summary.get("detected_3sigma"),
            "counts_toward_multi_sightline_phase_replication": _counts_toward_independent_phase_replication(
                name, summary
            ),
            "is_primary_single_sightline_phase_target": name == "J0437-4715",
        }
        rows.append(row)

    # Pulsars with Step 003 output but not listed in step_024 (e.g. auxiliary MeerKAT rows)
    listed = {r["pulsar"] for r in rows}
    for pulsar, path in summary_paths.items():
        if pulsar in listed:
            continue
        summary = _load_json_path(path)
        role_guess = "MeerKAT auxiliary" if ("J1731" in path.name or pulsar == "J1731-4744") else ""
        tier = _tier_for_pulsar(pulsar, summary, role_guess)
        rows.append(
            {
                "pulsar": pulsar,
                "evidence_tier": tier,
                "ingestion_role": role_guess or None,
                "data_source": None,
                "step003_summary_relpath": str(path.relative_to(PROJECT_ROOT)),
                "n_epochs": summary.get("n_epochs"),
                "n_total_triplets": summary.get("n_total_triplets"),
                "n_independent_samples": summary.get("n_independent_samples"),
                "phase_closure_rayleigh_p": summary.get("phase_closure_rayleigh_p"),
                "phase_closure_rbar": summary.get("phase_closure_rbar"),
                "detected_3sigma_phase_closure": summary.get("detected_3sigma"),
                "counts_toward_multi_sightline_phase_replication": _counts_toward_independent_phase_replication(
                    pulsar, summary
                ),
                "is_primary_single_sightline_phase_target": False,
            }
        )

    n_primary = sum(1 for r in rows if r["evidence_tier"] == "primary_phase_closure")
    n_geom = sum(1 for r in rows if "geometric_diagnostic" in r["evidence_tier"])
    n_bound = sum(
        1
        for r in rows
        if r["evidence_tier"] in ("catalog_noise_limited_bounding_row", "catalog_noise_limited_telescope_auxiliary")
    )
    n_multi_rep = sum(1 for r in rows if r["counts_toward_multi_sightline_phase_replication"])

    subband_paths = sorted(RESULTS_DIR.glob("step_003_closure_final_summary_j0437_sb*.json"))

    return {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "ingestion_reference": "results/step_024_multi_pulsar_ingestion.json",
        "scaling_reference_note": scaling.get("reason"),
        "epistemic_accounting": {
            "interpretation": (
                "Population-level falsification of GR or establishment of TEP as physical fact "
                "requires multiple independent high-SNR phase-closure detections on distinct sightlines. "
                "The present catalog mixes one primary phase target, one geometry-supporting Parkes pulsar, "
                "and noise-limited bounding rows."
            ),
            "n_rows_in_ingestion_catalog": len(ingestion.get("pulsars_analyzed", [])),
            "n_pulsars_with_step003_summary": sum(1 for r in rows if r.get("step003_summary_relpath")),
            "n_primary_phase_closure_targets": n_primary,
            "n_geometric_diagnostic_targets": n_geom,
            "n_noise_limited_bounding_or_auxiliary_rows": n_bound,
            "n_pulsars_counting_toward_independent_phase_replication": n_multi_rep,
        },
        "example_high_snr_replication_classes_not_in_current_catalog": [
            {
                "example": "PSR B1937+21",
                "note": (
                    "Illustrative millisecond-pulsar class often cited for scintillation SNR; "
                    "not part of the processed PPTA/Jiamusi/MeerKAT artifact set in this repository."
                ),
            }
        ],
        "j0437_subband_robustness_summaries": [str(p.relative_to(PROJECT_ROOT)) for p in subband_paths],
        "pulsars": rows,
    }


def main() -> Dict[str, Any]:
    j0437 = load_json("step_003_closure_final_summary_j0437.json", required=True)
    falsification = load_json("step_013_falsification_criteria_results.json", required=True)
    controls = load_json("step_047_velocity_direction_controls_results.json", required=True)
    frame = load_json("step_048_cmb_dipole_frame_analysis.json", required=True)
    scaling = load_json("step_030_tep_scaling_analysis.json")
    telescope = load_json("step_037_multi_telescope_results.json")
    chromatic = load_json("step_043_definitive_chromatic_test.json", required=True)
    orbital = load_json("step_046_bayesian_orbital_phasebin_results.json", required=True)
    orbital_mixed = load_json("step_046b_hierarchical_orbital_modulation.json")
    systematics = load_json("step_042_systematic_error_results.json")

    directional = (
        frame.get("pulsars", {})
        .get("J0437-4715", {})
        .get("enhanced_tests", {})
        .get("directional_specificity", {})
    )
    j0437_frame = frame.get("pulsars", {}).get("J0437-4715", {}).get("ssb_frame_summary", {})
    orbital_epoch = orbital.get("epoch_blocked_analysis", {})
    orbital_epoch_fit = orbital_epoch.get("modulation_test", {})
    orbital_epoch_perm = orbital_epoch.get("epoch_phase_permutation", {})
    orbital_mixed_fe = orbital_mixed.get("fixed_effects", {})
    orbital_mixed_amp = orbital_mixed_fe.get("modulation_amplitude_A_ns", {})
    orbital_mixed_mc = orbital_mixed.get("model_comparison", {})

    ch_mc = chromatic.get("model_comparison") or {}
    ch_tep = ch_mc.get("M_TEP") or {}
    ch_ism = ch_mc.get("M_ISM") or {}
    ch_gen = ch_mc.get("M_general") or {}
    delta_aic_ism_minus_tep = None
    try:
        if ch_tep.get("aic") is not None and ch_ism.get("aic") is not None:
            delta_aic_ism_minus_tep = float(ch_ism["aic"]) - float(ch_tep["aic"])
    except (TypeError, ValueError):
        delta_aic_ism_minus_tep = None
    ch_bf = chromatic.get("best_fit_params") or {}
    ch_delta_hat = ch_bf.get("delta")
    ch_ws = chromatic.get("within_source_test") or {}
    ch_ws_bf = ch_ws.get("best_fit") or {}

    ledger = {
        "evidence_standard": (
            "Primary claims require phase-domain circular evidence on independent epochs. "
            "Unsigned |H|, delay amplitudes, cross-pulsar scaling, chromatic fits, and orbital fits "
            "are diagnostics unless explicitly marked primary."
        ),
        "manuscript_claim_guardrails": {
            "preferred_headline": (
                "Robust single-sightline detection of non-zero Phase Closure in J0437-4715, "
                "consistent with TEP and supported by non-decisive geometric diagnostics."
            ),
            "avoid_phrases": [
                "proof of TEP",
                "decisive orbital modulation",
                "confirmed cross-telescope replication",
                "J1603 independent phase detection",
                "unsigned |H| detection",
                "raw |H| scaling law established",
            ],
            "required_caveats": [
                "J1603 is phase-noise-limited and contributes geometry, not a second phase detection.",
                "Jiamusi and MeerKAT rows are noise-limited bounds, not positive replication.",
                "Orbital modulation: mixed-effects model yields 1.14 ± 0.79 ns (LR p = 0.357, 2 df); not independently significant.",
                "Raw |H| and delay-amplitude fits are diagnostic because folded magnitudes have a noise floor.",
                "A local J0437 sightline anomaly remains a residual alternative until independent high-SNR replication.",
                "Chromatic step_043 (unsigned |H| hierarchy): valid_for_primary_inference is false in results/step_043_definitive_chromatic_test.json; cross-pulsar free exponent sits on the optimizer upper boundary; within-source J0437 sub-band comparison is directional (three bands).",
            ],
        },
        "primary_claim": {
            "claim": "J0437-4715 rejects the additive scalar path-delay null through non-zero Phase Closure psi.",
            "status": "primary_evidence",
            "pulsar": "J0437-4715",
            "phase_closure_mean_rad": j0437.get("phase_closure_mean_rad"),
            "phase_closure_circ_se_rad": j0437.get("phase_closure_circ_se_rad"),
            "unweighted_psi_rad": (
                j0437.get("phase_closure_mean_rad_unweighted")
                or j0437_frame.get("phase_closure_mean_unweighted_rad")
            ),
            "rayleigh_p_unweighted": (
                j0437.get("phase_closure_rayleigh_p_unweighted")
                or j0437_frame.get("phase_closure_rayleigh_p_unweighted")
            ),
            "rayleigh_p_weighted": j0437.get("phase_closure_rayleigh_p"),
            "v_test_p": j0437.get("phase_closure_v_p"),
            "bootstrap_ci_95_rad": [
                j0437.get("phase_closure_bootstrap_ci_95_lower_rad"),
                j0437.get("phase_closure_bootstrap_ci_95_upper_rad"),
            ],
            "falsification_status": falsification.get("summary", {}).get("overall_status"),
        },
        "supporting_diagnostics": {
            "phase_scramble": controls.get("results", {}).get("phase_scramble", {}),
            "pre_alignment_phase": controls.get("results", {}).get("pre_alignment_diagnostic", {}),
            "frame_invariance": {
                "unweighted_psi_invariant_all_directions": directional.get(
                    "unweighted_psi_invariant_all_directions"
                ),
                "unweighted_psi_max_abs_delta_rad": directional.get(
                    "unweighted_psi_max_abs_delta_rad"
                ),
            },
            "signed_bipolar_cancellation": falsification.get("statistics", {}).get(
                "signed_diagnostic", {}
            ),
            "systematics_note": systematics.get("systematic_budget", {}).get("interpretation")
            or systematics.get("interpretation"),
        },
        "explicit_non_claims": {
            "multi_pulsar_scaling": {
                "status": "not_established",
                "reason": scaling.get("reason"),
                "legacy_h_scaling_disabled": scaling.get("legacy_h_scaling_disabled"),
            },
            "cross_telescope_replication": {
                "status": telescope.get("status"),
                "instrumental_stats": telescope.get("instrumental_stats"),
            },
            "chromaticity": {
                "pipeline_step": "043",
                "results_json": "results/step_043_definitive_chromatic_test.json",
                "inference_status": chromatic.get("inference_status"),
                "valid_for_primary_inference": chromatic.get("valid_for_primary_inference"),
                "interpretation": chromatic.get("interpretation"),
                "cross_pulsar_unsigned_H_hierarchy": {
                    "preferred_model_by_aic": chromatic.get("preferred_model"),
                    "free_frequency_exponent_delta_hat": ch_delta_hat,
                    "free_delta_on_upper_boundary": ch_delta_hat == 2.0,
                    "delta_aic_M_ISM_minus_M_TEP": delta_aic_ism_minus_tep,
                    "aic_M_TEP": ch_tep.get("aic"),
                    "aic_M_ISM": ch_ism.get("aic"),
                    "aic_M_general": ch_gen.get("aic"),
                },
                "within_source_j0437_subbands": {
                    "n_points": ch_ws.get("n_points"),
                    "delta_hat": ch_ws_bf.get("delta"),
                    "delta_aic_tep_vs_free": ch_ws.get("delta_aic_tep_vs_free"),
                    "delta_aic_ism_vs_free": ch_ws.get("delta_aic_ism_vs_free"),
                },
                "limitations": chromatic.get("limitations"),
            },
            "orbital_modulation": {
                "status": "suggestive_follow_up",
                "primary_analysis": "step_046b_mixed_effects_reml",
                "mixed_effects_amplitude_ns": orbital_mixed_amp.get("mean"),
                "mixed_effects_amplitude_se_ns": orbital_mixed_amp.get("std"),
                "mixed_effects_ci_95": [
                    orbital_mixed_amp.get("ci_95_lower"),
                    orbital_mixed_amp.get("ci_95_upper"),
                ],
                "mixed_effects_lr_pvalue": orbital_mixed_mc.get("lr_pvalue"),
                "mixed_effects_favors_modulation": orbital_mixed_mc.get("favors_modulation"),
                "epoch_blocked_amplitude_ns": orbital_epoch_fit.get("fitted_amplitude_ns"),
                "epoch_blocked_amplitude_err_ns": orbital_epoch_fit.get("fitted_amplitude_err_ns"),
                "epoch_blocked_sigma": orbital_epoch_fit.get("modulation_significance_sigma"),
                "nested_model_p": orbital_epoch_fit.get("modulation_p_value"),
                "permutation_p": orbital_epoch_perm.get("empirical_p_value"),
            },
            "unsigned_H": {
                "status": "diagnostic_only",
                "reason": "Folded magnitudes have a noise floor and delay-domain systematics exceed formal SEM.",
                "H_magnitude_ns": j0437.get("H_magnitude_ns"),
                "H_noise_bias_ns": j0437.get("H_noise_bias_ns"),
                "H_excess_ns": j0437.get("H_excess_ns"),
            },
        },
        "paper_ready_summary": {
            "one_sentence": (
                "The present pipeline supports a phase-domain J0437 detection of non-zero synchronization "
                "holonomy, while classifying J1603 geometry, orbital structure, chromatic fits, scaling, "
                "and unsigned-|H| amplitudes as supporting diagnostics or follow-up tests."
            ),
            "reviewer_risk_reduction": [
                "Separates circular phase evidence from folded-magnitude noise-floor artifacts.",
                "Prevents J1603 from being over-described as an independent detection.",
                "Prevents orbital triplet-level precision from being treated as independent-epoch significance.",
                "Makes cross-telescope replication status explicit.",
                "Keeps residual local-ISM alternatives visible instead of overclaiming exclusion.",
                "Records step_043 inference_status and valid_for_primary_inference so unsigned-|H| chromatic hierarchies cannot be read as phase-primary claims.",
            ],
        },
        "replication_gates": [
            "Independent reproduction of J0437 phase closure from raw dynamic spectra.",
            "A second nearby high-SNR pulsar with robust phase closure in the predicted geometry.",
            "Cross-telescope detection of phase closure, not merely unsigned |H|.",
            "Chromatic test on same-source multi-band data with enough epochs for phase statistics.",
            "Orbital or annual modulation detected with epoch-blocked significance.",
        ],
    }

    tier_summary = build_evidence_tier_summary(scaling)
    ledger["evidence_tier_summary"] = tier_summary

    out = RESULTS_DIR / "step_049_evidence_ledger.json"
    with open(out, "w") as fh:
        json.dump(ledger, fh, indent=2, cls=NpEncoder)

    tier_out = RESULTS_DIR / "step_049_evidence_tier_summary.json"
    with open(tier_out, "w") as fh:
        json.dump(tier_summary, fh, indent=2, cls=NpEncoder)

    print_status("STEP 049: EVIDENCE LEDGER", "TITLE")
    print_status(f"Primary claim status: {ledger['primary_claim']['status']}", "SUCCESS")
    print_status(f"Falsification status: {ledger['primary_claim']['falsification_status']}", "INFO")
    print_status(f"Results saved to: {out}", "INFO")
    print_status(f"Evidence tier summary: {tier_out}", "INFO")
    return ledger


if __name__ == "__main__":
    main()
