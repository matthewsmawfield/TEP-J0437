#!/usr/bin/env python3
"""
Step 030: TEP Scaling Analysis - Multi-Pulsar Comparison

Audits whether the phase-closure detections are numerous enough to test
multi-pulsar TEP scaling.  Legacy unsigned-|H| amplitudes are retained only as
diagnostic context because mean(|delay|) is noise-floor biased.

TEP Predictions:
- Holonomy magnitude |H| should scale with pulsar distance (D)
- Should scale with proper motion velocity (v)
- Should scale with scattering strength (arc curvature eta)
- Sign depends on velocity direction relative to scattering geometry

At present this analysis should only perform a scaling fit when at least two
independent pulsars pass the phase-domain circular-statistics gate.
"""

import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.json_numpy import NpEncoder
from scripts.utils.config import (
    J0437_DIST_PC, J0437_PM_RA, J0437_PM_DEC, J0437_PB_DAYS,
    J0437_T0_MJD, J0437_A1_LC, J0437_S_SCREEN,
    J1603_DIST_PC, J1603_PM_RA, J1603_PM_DEC, J1603_K_KMS, J1603_S_SCREEN
)
from scripts.utils.logger import TEPLogger, print_status, set_step_logger

RESULTS_DIR = PROJECT_ROOT / "results"

# Physical Parameters now imported from config.py
# J0437-4715 and J1603-7202 parameters are defined in scripts/utils/config.py

# Jiamusi pulsars (from Wang et al. 2018, A&A 618, A186)
# All 10 pulsars with DM-based distances and estimated parameters

# B0329+54 - DM ~26.76 pc/cm^3, distance ~1.0 kpc (Wang et al. 2018)
B0329_DIST_PC = 1000.0  # DM distance from Wang et al. 2018
B0329_PM_RA = 7.10  # mas/yr (ATNF/Wang et al. 2018)
B0329_PM_DEC = -11.75  # mas/yr (ATNF/Wang et al. 2018)
B0329_K_KMS = 0.0  # Non-binary
B0329_S_SCREEN = 0.99  # From eta analysis (s=0.990)

# B0355+54 - DM ~57.14 pc/cm^3, distance ~1.0 kpc (Wang et al. 2018)
B0355_DIST_PC = 1000.0  # DM distance from Wang et al. 2018
B0355_PM_RA = 9.17  # mas/yr (ATNF/Wang et al. 2018)
B0355_PM_DEC = 0.70  # mas/yr (ATNF/Wang et al. 2018)
B0355_K_KMS = 0.0  # Non-binary
B0355_S_SCREEN = 0.62  # From eta analysis (s=0.619)

# B0540+23 - DM ~77.70 pc/cm^3, distance ~1.6 kpc (Wang et al. 2018)
B0540_DIST_PC = 1600.0  # DM distance from Wang et al. 2018
B0540_PM_RA = 2.50  # mas/yr (ATNF/Wang et al. 2018)
B0540_PM_DEC = -21.80  # mas/yr (ATNF/Wang et al. 2018)
B0540_K_KMS = 0.0  # Non-binary
B0540_S_SCREEN = 0.99  # From eta analysis (s=0.990)

# B0740-28 - DM ~73.78 pc/cm^3, distance ~2.0 kpc (Wang et al. 2018)
# NOTE: Screen distance not available in Wang et al. 2018
# Run step_032_extract_eta_from_arclets.py to extract eta from actual data
B0740_DIST_PC = 2000.0  # DM distance from Wang et al. 2018
B0740_PM_RA = -2.44  # mas/yr (Wang et al. 2018, Table 1)
B0740_PM_DEC = -0.09  # mas/yr (Wang et al. 2018, Table 1)
B0740_K_KMS = 0.0  # Non-binary

# B1508+55 - DM ~19.61 pc/cm^3, distance ~2.1 kpc (Wang et al. 2018)
B1508_DIST_PC = 2100.0  # DM distance from Wang et al. 2018
B1508_PM_RA = -73.70  # mas/yr (ATNF/Chatterjee et al. 2009)
B1508_PM_DEC = -62.70  # mas/yr (ATNF/Chatterjee et al. 2009)
B1508_K_KMS = 0.0  # Non-binary
B1508_S_SCREEN = 0.34  # From eta analysis (s=0.34)

# B1933+16 - DM ~158.52 pc/cm^3, distance ~3.7 kpc (Wang et al. 2018)
B1933_DIST_PC = 3700.0  # DM distance from Wang et al. 2018
B1933_PM_RA = -2.00  # mas/yr (ATNF/Chatterjee et al. 2009)
B1933_PM_DEC = -0.10  # mas/yr (ATNF/Chatterjee et al. 2009)
B1933_K_KMS = 0.0  # Non-binary
B1933_S_SCREEN = 0.94  # From eta analysis (s=0.935)

# B2154+40 - DM ~71.12 pc/cm^3, distance ~2.9 kpc (Wang et al. 2018)
B2154_DIST_PC = 2900.0  # DM distance from Wang et al. 2018
B2154_PM_RA = 14.60  # mas/yr (ATNF/Chatterjee et al. 2009)
B2154_PM_DEC = -2.60  # mas/yr (ATNF/Chatterjee et al. 2009)
B2154_K_KMS = 0.0  # Non-binary
B2154_S_SCREEN = 0.75  # From eta analysis (s=0.746)

# B2310+42 - DM ~17.27 pc/cm^3, distance ~1.06 kpc (Wang et al. 2018)
B2310_DIST_PC = 1060.0  # DM distance from Wang et al. 2018
B2310_PM_RA = -3.00  # mas/yr (ATNF/Chatterjee et al. 2009)
B2310_PM_DEC = -6.00  # mas/yr (ATNF/Chatterjee et al. 2009)
B2310_K_KMS = 0.0  # Non-binary
B2310_S_SCREEN = 0.99  # From eta analysis (s=0.990)

# B2324+60 - DM ~122.61 pc/cm^3, distance ~2.7 kpc (Wang et al. 2018)
# NOTE: Screen distance not available in Wang et al. 2018
# Run step_032_extract_eta_from_arclets.py to extract eta from actual data
B2324_DIST_PC = 2700.0  # DM distance from Wang et al. 2018
B2324_PM_RA = 0.0  # mas/yr (Wang et al. 2018, Table 1 - not measured)
B2324_PM_DEC = 0.0  # mas/yr (Wang et al. 2018, Table 1 - not measured)
B2324_K_KMS = 0.0  # Non-binary

# B2351+61 - DM ~94.66 pc/cm^3, distance ~2.4 kpc (Wang et al. 2018)
# NOTE: Screen distance not available in Wang et al. 2018
# Run step_032_extract_eta_from_arclets.py to extract eta from actual data
B2351_DIST_PC = 2400.0  # DM distance from Wang et al. 2018
B2351_PM_RA = -0.19  # mas/yr (Wang et al. 2018, Table 1)
B2351_PM_DEC = -0.01  # mas/yr (Wang et al. 2018, Table 1)
B2351_K_KMS = 0.0  # Non-binary

# Conversion factors (import from config for consistency)
from scripts.utils.config import PC_TO_KM, MAS_YR_TO_KM_S


def proper_motion_to_velocity(pm_ra, pm_dec, dist_pc):
    """Convert proper motion to transverse velocity in km/s."""
    # Convert to km/s: v = u * d * conversion factor
    # u in mas/yr, d in pc
    v_ra = pm_ra * dist_pc * MAS_YR_TO_KM_S  # km/s
    v_dec = pm_dec * dist_pc * MAS_YR_TO_KM_S  # km/s
    v_total = np.sqrt(v_ra**2 + v_dec**2)
    return v_total, v_ra, v_dec


def load_closure_summary(pulsar_name):
    """Load closure delay summary for a pulsar."""
    summary_file = RESULTS_DIR / f"step_003_closure_final_summary_{pulsar_name}.json"
    if not summary_file.exists():
        return None
    with open(summary_file, "r") as f:
        return json.load(f)


def load_secondary_catalog(pulsar_name):
    """Load secondary spectrum catalog for a pulsar."""
    # Handle Jiamusi pulsars specifically (combined catalog)
    if pulsar_name.startswith("B"):
        catalog_file = (
            PROJECT_ROOT / "data" / "secondary" / "jiamusi_secondary_catalog.json"
        )
        if not catalog_file.exists():
            return None
        with open(catalog_file, "r") as f:
            full_cat = json.load(f)

        # Filter for this specific pulsar
        prefix = pulsar_name.split("+")[0].split("-")[0]
        pulsar_epochs = [
            e for e in full_cat.get("epochs", []) if e["file"].startswith(prefix + "_")
        ]
        return {"epochs": pulsar_epochs}

    # MeerKAT pulsars (individual catalogs)
    catalog_file = (
        PROJECT_ROOT / "data" / "secondary" / f"{pulsar_name}_secondary_catalog.json"
    )
    if not catalog_file.exists():
        return None
    with open(catalog_file, "r") as f:
        return json.load(f)


def phase_detection_sigma(summary):
    """Return signed Gaussian-equivalent sigma from circular phase p-values."""
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


def is_robust_phase_detection(summary):
    """Gate scaling on phase closure, never on unsigned-|H| diagnostics."""
    n_epochs = int(summary.get("n_independent_samples", summary.get("n_epochs", 0)) or 0)
    return bool(n_epochs >= 5 and abs(phase_detection_sigma(summary)) >= 3.0)


def calculate_tep_scaling_ratio():
    """Calculate predicted TEP scaling ratio between pulsars."""

    print_status("=" * 70, "TITLE")
    print_status("TEP SCALING ANALYSIS: Multi-Pulsar Comparison", "TITLE")
    print_status("=" * 70, "TITLE")

    # Define all pulsars that might contribute to a future scaling test.
    PULSARS = {
        "j0437": {
            "name": "J0437-4715",
            "dist_pc": J0437_DIST_PC,
            "pm_ra": J0437_PM_RA,
            "pm_dec": J0437_PM_DEC,
            "s_screen": J0437_S_SCREEN,
            "telescope": "Parkes/PPTA",
        },
        "j1603": {
            "name": "J1603-7202",
            "dist_pc": J1603_DIST_PC,
            "pm_ra": J1603_PM_RA,
            "pm_dec": J1603_PM_DEC,
            "s_screen": J1603_S_SCREEN,
            "telescope": "Parkes/PPTA",
        },
        "B0329": {
            "name": "B0329+54",
            "dist_pc": B0329_DIST_PC,
            "pm_ra": B0329_PM_RA,
            "pm_dec": B0329_PM_DEC,
            "s_screen": B0329_S_SCREEN,
            "telescope": "Jiamusi",
        },
        "B0355": {
            "name": "B0355+54",
            "dist_pc": B0355_DIST_PC,
            "pm_ra": B0355_PM_RA,
            "pm_dec": B0355_PM_DEC,
            "s_screen": B0355_S_SCREEN,
            "telescope": "Jiamusi",
        },
        "B0540": {
            "name": "B0540+23",
            "dist_pc": B0540_DIST_PC,
            "pm_ra": B0540_PM_RA,
            "pm_dec": B0540_PM_DEC,
            "s_screen": B0540_S_SCREEN,
            "telescope": "Jiamusi",
        },
    }

    # Load results for all available pulsars
    available_pulsars = {}
    excluded_pulsars = []
    for pulsar_id, params in PULSARS.items():
        summary = load_closure_summary(pulsar_id.lower().replace("b", "B"))
        if summary is not None and is_robust_phase_detection(summary):
            phase_sigma = phase_detection_sigma(summary)
            available_pulsars[pulsar_id] = {"params": params, "summary": summary}
            print_status(
                f"[OK] {params['name']} ({params['telescope']}): phase closure = {phase_sigma:.2f}sigma",
                "SUCCESS",
            )
        elif summary is not None:
            phase_sigma = phase_detection_sigma(summary)
            reason = (
                "phase evidence below threshold; legacy unsigned-|H| amplitude is diagnostic only"
            )
            excluded_pulsars.append(
                {
                    "name": params["name"],
                    "telescope": params["telescope"],
                    "phase_sigma": phase_sigma,
                    "H_magnitude_ns": summary.get("H_magnitude_ns"),
                    "H_excess_ns": summary.get("H_excess_ns"),
                    "reason": reason,
                }
            )
            print_status(
                f"○ {params['name']}: not included ({reason}; phase sigma={phase_sigma:.2f})",
                "INFO",
            )
        else:
            print_status(f"[SKIP] {params['name']}: no closure summary available", "WARNING")

    n_available = len(available_pulsars)
    print_status(f"\nRobust phase detections (≥3sigma): {n_available} pulsars", "INFO")

    if n_available < 2:
        print_status(
            "Need at least 2 independent robust phase detections for scaling analysis",
            "WARNING",
        )
        print_status(
            f"Only {n_available} pulsar(s) meet the ≥3sigma threshold. Multi-pulsar scaling analysis requires at least 2.",
            "INFO",
        )
        # Save partial results for downstream steps
        partial_results = {
            "n_pulsars": n_available,
            "pulsars": [data["params"]["name"] for data in available_pulsars.values()],
            "scaling_analysis_performed": False,
            "reason": f"Insufficient robust phase detections (need ≥2, have {n_available})",
            "excluded_pulsars": excluded_pulsars,
            "legacy_h_scaling_disabled": True,
            "note": (
                "Unsigned |H| amplitudes are noise-floor-biased diagnostics and are not used to claim "
                "multi-pulsar scaling. Single-pulsar phase results are available in step_003."
            ),
        }
        output_file = RESULTS_DIR / "step_030_tep_scaling_analysis.json"
        with open(output_file, "w") as f:
            json.dump(partial_results, f, indent=2, cls=NpEncoder)
        print_status(f"Partial results saved to {output_file}", "INFO")
        return

    # Calculate velocities for available pulsars
    for pulsar_id, data in available_pulsars.items():
        params = data["params"]
        v_total, v_ra, v_dec = proper_motion_to_velocity(
            params["pm_ra"], params["pm_dec"], params["dist_pc"]
        )
        data["velocity"] = v_total
        data["v_ra"] = v_ra
        data["v_dec"] = v_dec

    # Perform multi-pulsar scaling analysis
    perform_multi_pulsar_scaling(available_pulsars)


def compute_total_sigma(H_pred_i: float, H_sem_i: float, pulsar_name: str) -> float:
    """
    Compute total uncertainty on H for model comparison.

    Combines measurement SEM with prediction uncertainty from:
    - Distance uncertainty: 2% for MeerKAT (parallax), 30% for Jiamusi (DM-based)
    - Proper-motion uncertainty: 5% for all pulsars
    - H_pred scales as D x v, so fractional uncertainty propagates in quadrature

    The prediction uncertainty sigma_pred = H_pred * sqrt(dist_frac^2 + pm_frac^2)
    is added in quadrature with H_sem.
    """
    # Distance uncertainty fractions (from literature)
    # MeerKAT (parallax): 2% from parallax error propagation (precision timing)
    # Jiamusi (DM-based): 30% from DM-distance relation uncertainty (Wang et al. 2018)
    if pulsar_name.startswith("J"):  # MeerKAT: parallax distance
        dist_frac = 0.02
    else:  # Jiamusi: DM-based distance
        dist_frac = 0.30
    pm_frac = 0.05  # 5% proper motion uncertainty (typical VLBI timing)
    sigma_pred = abs(H_pred_i) * np.sqrt(dist_frac**2 + pm_frac**2)
    return float(np.sqrt(H_sem_i**2 + sigma_pred**2))


def perform_multi_pulsar_scaling(pulsars):
    """
    Perform scaling analysis with multiple pulsars.

    Tests TEP scaling predictions across N pulsars with strong detections.
    """
    print_status("\n" + "=" * 70, "TITLE")
    print_status("MULTI-PULSAR SCALING TEST", "TITLE")
    print_status("=" * 70, "TITLE")

    n_pulsars = len(pulsars)
    print_status(f"\nAnalyzing {n_pulsars} pulsars with strong TEP detections", "INFO")

    # Extract data arrays
    names = []
    H_values = []
    H_errors = []
    distances = []
    velocities = []
    screen_factors = []
    telescopes = []

    for pulsar_id, data in pulsars.items():
        params = data["params"]
        summary = data["summary"]
        names.append(params["name"])
        # Use physical holonomy magnitude |H| (ns) for scaling analysis
        h_mag = summary.get("H_magnitude_ns", 0.0)
        h_sem = summary.get("H_sem_ns", 1.0)
        if h_sem is None or (isinstance(h_sem, float) and math.isnan(h_sem)) or h_sem <= 0:
            h_sem = 1.0
        H_values.append(h_mag)
        H_errors.append(h_sem)
        distances.append(params["dist_pc"])
        velocities.append(data["velocity"])
        s = params["s_screen"]
        screen_factors.append(s * (1 - s))
        telescopes.append(params["telescope"])

    H_values = np.array(H_values)
    H_errors = np.array(H_errors)
    distances = np.array(distances)
    velocities = np.array(velocities)
    screen_factors = np.array(screen_factors)

    # Print comparison table
    print_status("\n--- PULSAR PARAMETERS ---", "INFO")
    print_status(
        f"{'Pulsar':<12} {'H (ns)':>10} {'D (pc)':>8} {'v (km/s)':>10} {'s(1-s)':>8} {'Tel':>8}",
        "INFO",
    )
    print_status("-" * 70, "INFO")
    for i in range(n_pulsars):
        print_status(
            f"{names[i]:<12} {H_values[i]:>10.3f} {distances[i]:>8.1f} {velocities[i]:>10.1f} {screen_factors[i]:>8.3f} {telescopes[i]:>8}",
            "INFO",
        )

    # Test scaling models
    print_status("\n--- SCALING MODEL TESTS ---", "INFO")

    # Model 1: H ∝ D (distance only)
    predicted_H_d = H_values[0] * distances / distances[0]
    chi2_d = np.sum(((H_values - predicted_H_d) / H_errors) ** 2)

    # Model 2: H ∝ v (velocity only)
    predicted_H_v = H_values[0] * velocities / velocities[0]
    chi2_v = np.sum(((H_values - predicted_H_v) / H_errors) ** 2)

    # Model 3: H ∝ D x v (distance-velocity)
    dv_product = distances * velocities
    predicted_H_dv = H_values[0] * dv_product / dv_product[0]
    chi2_dv = np.sum(((H_values - predicted_H_dv) / H_errors) ** 2)

    # Model 4: H ∝ D x v x s(1-s) (full TEP)
    dvs_product = distances * velocities * screen_factors
    predicted_H_full = H_values[0] * dvs_product / dvs_product[0]
    chi2_full = np.sum(((H_values - predicted_H_full) / H_errors) ** 2)

    dof = max(n_pulsars - 1, 1)  # degrees of freedom, min 1 to avoid div by zero

    # Compute total errors for each pulsar under each model (measurement + prediction uncertainty)
    total_errors_dv = np.array(
        [
            compute_total_sigma(predicted_H_dv[i], H_errors[i], names[i])
            for i in range(n_pulsars)
        ]
    )
    total_errors_full = np.array(
        [
            compute_total_sigma(predicted_H_full[i], H_errors[i], names[i])
            for i in range(n_pulsars)
        ]
    )
    total_errors_d = np.array(
        [
            compute_total_sigma(predicted_H_d[i], H_errors[i], names[i])
            for i in range(n_pulsars)
        ]
    )
    total_errors_v = np.array(
        [
            compute_total_sigma(predicted_H_v[i], H_errors[i], names[i])
            for i in range(n_pulsars)
        ]
    )

    chi2_dv_phys = float(np.sum(((H_values - predicted_H_dv) / total_errors_dv) ** 2))
    chi2_full_phys = float(
        np.sum(((H_values - predicted_H_full) / total_errors_full) ** 2)
    )
    chi2_d_phys = float(np.sum(((H_values - predicted_H_d) / total_errors_d) ** 2))
    chi2_v_phys = float(np.sum(((H_values - predicted_H_v) / total_errors_v) ** 2))

    print_status(f"\nModel fits (χ² with {dof} dof):", "INFO")
    print_status(
        f"  H ∝ D:       χ² = {chi2_d:.2f},  reduced χ² = {chi2_d / dof:.2f}", "INFO"
    )
    print_status(
        f"  H ∝ v:       χ² = {chi2_v:.2f},  reduced χ² = {chi2_v / dof:.2f}", "INFO"
    )
    print_status(
        f"  H ∝ Dxv:     χ² = {chi2_dv:.2f},  reduced χ² = {chi2_dv / dof:.2f}", "INFO"
    )
    print_status(
        f"  H ∝ Dxvxs(1-s): χ² = {chi2_full:.2f},  reduced χ² = {chi2_full / dof:.2f}",
        "INFO",
    )

    print_status(
        f"\n  Chi2 WITH physical distance/PM uncertainties (more meaningful):", "INFO"
    )
    print_status(
        f"    Dxv model:  chi2 = {chi2_dv_phys:.1f} / {dof} DOF = {chi2_dv_phys / dof if dof > 0 else float('nan'):.2f}",
        "INFO",
    )
    comb_z = float(np.sqrt(np.sum((H_values / H_errors) ** 2)))
    n_det_5 = int(np.sum(H_values / H_errors > 5))
    print_status(
        f"\n  Combined {n_pulsars}-pulsar detection: z = {comb_z:.1f}sigma", "INFO"
    )
    print_status(f"  Pulsars detected >5sigma individually: {n_det_5}/{n_pulsars}", "INFO")

    # Determine best model
    models = {
        "Distance (D)": chi2_d,
        "Velocity (v)": chi2_v,
        "DistancexVelocity (Dxv)": chi2_dv,
        "Full TEP (Dxvxs(1-s))": chi2_full,
    }
    best_model = min(models, key=models.get)
    best_chi2 = models[best_model]

    print_status(f"\nBest model: {best_model} (χ² = {best_chi2:.2f})", "INFO")

    if best_chi2 / dof < 2.0:
        print_status(
            f"  -> Good fit (reduced χ² < 2). Scaling consistent with TEP predictions.",
            "SUCCESS",
        )
    elif best_chi2 / dof < 4.0:
        print_status(
            f"  -> Marginal fit. Some deviation from simple scaling, but within tolerance.",
            "INFO",
        )
    else:
        print_status(
            f"  -> Poor fit. Scaling may require additional physics or systematic correction.",
            "WARNING",
        )

    # Cross-telescope validation. Reaching this block only means multiple
    # phase-detected pulsars exist; it still does not prove instrument
    # replication unless both telescope groups contain phase detections.
    parkes_pulsars = [i for i in range(n_pulsars) if telescopes[i] == "Parkes/PPTA"]
    jiamusi_pulsars = [i for i in range(n_pulsars) if telescopes[i] == "Jiamusi"]

    print_status(f"\n--- CROSS-TELESCOPE VALIDATION ---", "INFO")
    print_status(f"  Parkes/PPTA: {len(parkes_pulsars)} pulsars", "INFO")
    print_status(f"  Jiamusi: {len(jiamusi_pulsars)} pulsars", "INFO")

    if len(parkes_pulsars) >= 1 and len(jiamusi_pulsars) >= 1:
        print_status(f"  -> Cross-instrument comparison possible", "SUCCESS")
        print_status(
            f"  -> Both telescope groups contain phase-detected rows; inspect step_037 before claiming replication",
            "INFO",
        )

    # Enhanced models with telescope offsets and arc curvature
    print_status(f"\n--- ENHANCED MODEL TESTS ---", "INFO")

    # Load arc curvature (eta) data from secondary catalogs
    eta_means = []
    for pulsar_id in pulsars.keys():
        catalog = load_secondary_catalog(pulsar_id.lower().replace("b", "B"))
        if catalog:
            eta1_values = [
                e["eta_screen1"] for e in catalog["epochs"] if e["eta_screen1"] > 0
            ]
            if eta1_values:
                eta_means.append(np.mean(eta1_values))
            else:
                eta_means.append(0.0)
        else:
            eta_means.append(0.0)

    eta_means = np.array(eta_means)

    # Model 5: H ∝ D x v + telescope offset (accounts for calibration differences)
    # Fit: H = A * D * v + offset_tel
    parkes_mask = np.array([t == "Parkes/PPTA" for t in telescopes])
    jiamusi_mask = np.array([t == "Jiamusi" for t in telescopes])

    # Simple 2-telescope fit: H = A * D * v + offset
    A_init = H_values[0] / (distances[0] * velocities[0])

    # Grid search for best telescope offset
    offsets = np.linspace(-2, 4, 100)
    best_offset = 0
    best_chi2_offset = 1e10

    for offset in offsets:
        H_pred = A_init * distances * velocities
        H_pred[jiamusi_mask] += offset  # Add offset to Jiamusi
        chi2 = np.sum(((H_values - H_pred) / H_errors) ** 2)
        if chi2 < best_chi2_offset:
            best_chi2_offset = chi2
            best_offset = offset

    # Calculate with best offset
    H_pred_offset = A_init * distances * velocities
    H_pred_offset[jiamusi_mask] += best_offset
    chi2_tel = np.sum(((H_values - H_pred_offset) / H_errors) ** 2)

    # Model 6: H ∝ eta (arc curvature only - scattering strength)
    if np.any(eta_means > 0):
        eta_available = eta_means > 0
        if np.sum(eta_available) >= 2:
            eta_ref = eta_means[eta_available][0]
            H_eta_ref = H_values[eta_available][0]
            H_pred_eta = np.zeros_like(H_values)
            H_pred_eta[eta_available] = H_eta_ref * eta_means[eta_available] / eta_ref
            chi2_eta = np.sum(
                (
                    (H_values[eta_available] - H_pred_eta[eta_available])
                    / H_errors[eta_available]
                )
                ** 2
            )
            dof_eta = max(np.sum(eta_available) - 1, 1)
            print_status(
                f"  H ∝ eta (arc curvature): χ² = {chi2_eta:.2f} ({dof_eta} dof)", "INFO"
            )
    else:
        chi2_eta = np.nan

    # Model 7: Combined Dxv + telescope offset
    dof_tel = max(n_pulsars - 2, 1)  # A and offset, min 1 to avoid div by zero
    print_status(
        f"  H ∝ Dxv + telescope offset: χ² = {chi2_tel:.2f}, reduced χ² = {chi2_tel / dof_tel:.2f}",
        "INFO",
    )
    print_status(f"    Best Jiamusi offset: {best_offset:+.3f} ns", "INFO")

    if chi2_tel / dof_tel < 2.0:
        print_status(f"    -> Good fit with telescope calibration offset", "SUCCESS")

    # Residual analysis by telescope
    print_status(f"\n--- RESIDUAL ANALYSIS ---", "INFO")
    residuals_dv = H_values - H_values[0] * distances * velocities / (
        distances[0] * velocities[0]
    )
    print_status(f"  Parkes/PPTA residuals (Dxv model):", "INFO")
    for i in parkes_pulsars:
        print_status(
            f"    {names[i]}: {residuals_dv[i]:+.3f} ns ({residuals_dv[i] / H_errors[i]:+.2f}sigma)",
            "INFO",
        )
    print_status(f"  Jiamusi residuals (Dxv model):", "INFO")
    for i in jiamusi_pulsars:
        print_status(
            f"    {names[i]}: {residuals_dv[i]:+.3f} ns ({residuals_dv[i] / H_errors[i]:+.2f}sigma)",
            "INFO",
        )

    parkes_mean_resid = float(np.mean(residuals_dv[parkes_mask])) if np.any(parkes_mask) else 0.0
    jiamusi_mean_resid = float(np.mean(residuals_dv[jiamusi_mask])) if np.any(jiamusi_mask) else 0.0
    print_status(f"\n  Mean Parkes/PPTA residual: {parkes_mean_resid:+.3f} ns", "INFO")
    print_status(f"  Mean Jiamusi residual: {jiamusi_mean_resid:+.3f} ns", "INFO")
    print_status(
        f"  Telescope difference: {jiamusi_mean_resid - parkes_mean_resid:+.3f} ns",
        "INFO",
    )

    # Save results
    results = {
        "n_pulsars": n_pulsars,
        "pulsars": names,
        "telescopes": telescopes,
        "H_values_ns": H_values.tolist(),
        "H_errors_ns": H_errors.tolist(),
        "distances_pc": distances.tolist(),
        "velocities_kms": velocities.tolist(),
        "screen_factors": screen_factors.tolist(),
        "scaling_models": {
            "distance_only": {
                "chi2": float(chi2_d),
                "dof": dof,
                "reduced_chi2": float(chi2_d / dof),
            },
            "velocity_only": {
                "chi2": float(chi2_v),
                "dof": dof,
                "reduced_chi2": float(chi2_v / dof),
            },
            "distance_velocity": {
                "chi2": float(chi2_dv),
                "dof": dof,
                "reduced_chi2": float(chi2_dv / dof),
            },
            "full_tep": {
                "chi2": float(chi2_full),
                "dof": dof,
                "reduced_chi2": float(chi2_full / dof),
            },
            "dv_telescope_offset": {
                "chi2": float(chi2_tel),
                "dof": dof_tel,
                "reduced_chi2": float(chi2_tel / dof_tel),
                "jiamusi_offset_ns": float(best_offset),
            },
        },
        "residual_analysis": {
            "parkes_ppta_mean_residual_ns": float(parkes_mean_resid),
            "jiamusi_mean_residual_ns": float(jiamusi_mean_resid),
            "telescope_difference_ns": float(jiamusi_mean_resid - parkes_mean_resid),
        },
        "best_model": best_model,
        "best_chi2": float(best_chi2),
        "interpretation": f"Multi-pulsar scaling with {n_pulsars} pulsars across 2 telescopes. Best simple model: {best_model}. Telescope offset of {best_offset:.2f} ns improves fit significantly, suggesting calibration differences or additional physics needed.",
        "predicted_H_by_model": {
            "distance_only": predicted_H_d.tolist(),
            "velocity_only": predicted_H_v.tolist(),
            "distance_velocity": predicted_H_dv.tolist(),
            "full_tep": predicted_H_full.tolist(),
        },
        "residuals_by_model": {
            "distance_velocity": (H_values - predicted_H_dv).tolist(),
            "full_tep": (H_values - predicted_H_full).tolist(),
        },
        "total_errors_ns": total_errors_dv.tolist(),
        "chi2_with_physical_errors": {
            "distance_only": chi2_d_phys,
            "velocity_only": chi2_v_phys,
            "distance_velocity": chi2_dv_phys,
            "full_tep": chi2_full_phys,
            "dof": int(dof),
            "reduced_chi2_dv": chi2_dv_phys / dof if dof > 0 else float('nan'),
            "note": (
                "Chi2 computed with total uncertainty = sqrt(H_sem^2 + sigma_pred^2). "
                "sigma_pred = H_pred * sqrt(dist_frac^2 + pm_frac^2); "
                "dist_frac=0.02 (MeerKAT parallax) or 0.30 (Jiamusi DM); pm_frac=0.05. "
                "This is the physically meaningful chi2 for model comparison."
            ),
        },
        "combined_detection": {
            "n_pulsars_detected_3sigma": int(
                np.sum(np.array(H_values) / np.array(H_errors) > 3)
            ),
            "n_pulsars_detected_5sigma": int(
                np.sum(np.array(H_values) / np.array(H_errors) > 5)
            ),
            "individual_significances": {
                names[i]: float(H_values[i] / H_errors[i]) for i in range(n_pulsars)
            },
            "combined_z_quadrature": float(
                np.sqrt(np.sum((np.array(H_values) / np.array(H_errors)) ** 2))
            ),
            "note": f"Combined z = sqrt(sum((H_i/SEM_i)^2)) for {n_pulsars} pulsars with strong H-magnitude detections. Note: this uses the legacy linear H/SEM significance, not circular-statistics detection status.",
        },
    }

    output_file = RESULTS_DIR / "step_030_tep_scaling_analysis.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2, cls=NpEncoder)

    print_status(f"\nResults saved to {output_file}", "SUCCESS")
    print_status("=" * 70, "TITLE")


if __name__ == "__main__":
    # Setup logging
    # Logger is set by run_pipeline.py via set_step_logger()
    # Do not create a new logger here to avoid overriding the pipeline's logger

    calculate_tep_scaling_ratio()
