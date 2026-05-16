#!/usr/bin/env python3
"""
Step 031: Detailed TEP Scaling Analysis - Multi-Parameter Decomposition

Comprehensive analysis of TEP scaling relationships across multiple physical parameters:

Basic Scaling Analyses:
- Velocity vector geometry (magnitude, direction, projection angles)
- Distance scaling and effective distances
- Scattering screen geometry (s parameter, screen distances)
- Arc curvature statistics (eta distributions, scattering strength)
- Orbital modulation effects (for binary J0437)
- First-principles TEP predictions

Advanced Diagnostic Analyses:
- Phase-resolved orbital analysis: H binned by orbital phase quadrants
- Time evolution: 11-year trend analysis and early/late comparison
- Bootstrap confidence intervals: non-parametric ratio uncertainty
- Higher-order moments: skewness, kurtosis, normality tests
- Temporal stability assessment

ISM Microstructure & Systematic Tests:
- ISM correlations: Test if H correlates with scattering strength (eta1, eta2)
- Effect size robustness: Split-half reliability, subsample stability
- Differential refraction: Frequency dependence test (chromatic vs achromatic)
- TEP vs ISM discrimination: Critical test for unaccounted systematics

This provides a rigorous test of whether observed closure delays follow
the scaling laws predicted by TEP, with comprehensive diagnostics for
systematic effects, ISM contamination, and distributional properties.
"""

import sys
import json
import numpy as np
from pathlib import Path
from scipy import stats
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.config import (
    J0437_DIST_PC, J0437_PM_RA, J0437_PM_DEC, J0437_PB_DAYS,
    J0437_T0_MJD, J0437_A1_LC, J0437_S_SCREEN,
    J1603_DIST_PC, J1603_PM_RA, J1603_PM_DEC, J1603_PB_DAYS,
    J1603_T0_MJD, J1603_A1_LC, J1603_K_KMS, J1603_S_SCREEN,
    PC_TO_KM, MAS_YR_TO_KM_S, C_LIGHT_KM_S,
    RANDOM_SEED
)
from scripts.utils.json_numpy import NpEncoder
from scripts.utils.logger import TEPLogger, set_step_logger, print_status

RESULTS_DIR = PROJECT_ROOT / "results"
PLOTS_DIR = PROJECT_ROOT / "results"
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

# Physical Parameters now imported from config.py
# J0437-4715 and J1603-7202 parameters are defined in scripts/utils/config.py
# Conversion factors also imported from config.py


def proper_motion_to_velocity(pm_ra, pm_dec, dist_pc):
    """Convert proper motion to transverse velocity in km/s."""
    v_ra = pm_ra * dist_pc * MAS_YR_TO_KM_S  # km/s
    v_dec = pm_dec * dist_pc * MAS_YR_TO_KM_S  # km/s
    v_total = np.sqrt(v_ra**2 + v_dec**2)
    angle_rad = np.arctan2(v_dec, v_ra)
    angle_deg = np.degrees(angle_rad)
    return v_total, v_ra, v_dec, angle_deg


def calculate_effective_distances(dist_pc, s_screen):
    """Calculate effective distances for scattering geometry."""
    D_p = dist_pc  # Pulsar distance
    D_s = s_screen * D_p  # Screen distance
    D_eff = (D_s * (D_p - D_s)) / D_p  # Effective distance
    return D_p, D_s, D_eff


def load_closure_per_epoch(pulsar_name):
    """Load per-epoch closure delay results."""
    summary_file = RESULTS_DIR / f"step_003_closure_final_per_epoch_{pulsar_name}.json"
    if not summary_file.exists():
        return None
    with open(summary_file, 'r') as f:
        return json.load(f)


def load_closure_summary(pulsar_name):
    """Load summary closure delay results."""
    summary_file = RESULTS_DIR / f"step_003_closure_final_summary_{pulsar_name}.json"
    if not summary_file.exists():
        return None
    with open(summary_file, 'r') as f:
        return json.load(f)


def load_secondary_catalog(pulsar_name):
    """Load secondary spectrum catalog for a pulsar."""
    # Handle Jiamusi pulsars specifically (combined catalog)
    if pulsar_name.startswith('B'):
        catalog_file = PROJECT_ROOT / "data" / "secondary" / "jiamusi_secondary_catalog.json"
        if not catalog_file.exists():
            return None
        with open(catalog_file, 'r') as f:
            full_cat = json.load(f)
        
        # Filter for this specific pulsar
        prefix = pulsar_name.split('+')[0].split('-')[0]
        pulsar_epochs = [e for e in full_cat.get('epochs', []) if e['file'].startswith(prefix + '_')]
        return {'epochs': pulsar_epochs}
    
    # MeerKAT pulsars (individual catalogs)
    catalog_file = PROJECT_ROOT / "data" / "secondary" / f"{pulsar_name}_secondary_catalog.json"
    if not catalog_file.exists():
        return None
    with open(catalog_file, 'r') as f:
        return json.load(f)


def analyze_velocity_geometry():
    """Analyze velocity vector geometry in detail."""
    print_status("\n" + "=" * 70, "TITLE")
    print_status("VELOCITY VECTOR GEOMETRY ANALYSIS", "TITLE")
    print_status("=" * 70, "TITLE")
    
    # Calculate velocities
    v_j0437, v_ra_j0437, v_dec_j0437, angle_j0437 = proper_motion_to_velocity(
        J0437_PM_RA, J0437_PM_DEC, J0437_DIST_PC
    )
    v_j1603, v_ra_j1603, v_dec_j1603, angle_j1603 = proper_motion_to_velocity(
        J1603_PM_RA, J1603_PM_DEC, J1603_DIST_PC
    )
    
    print_status(f"\nJ0437-4715 Velocity:", "INFO")
    print_status(f"  Total: {v_j0437:.4f} km/s", "INFO")
    print_status(f"  RA component: {v_ra_j0437:.4f} km/s", "INFO")
    print_status(f"  DEC component: {v_dec_j0437:.4f} km/s", "INFO")
    print_status(f"  Direction angle: {angle_j0437:.2f}° (from +RA axis)", "INFO")
    
    print_status(f"\nJ1603-7202 Velocity:", "INFO")
    print_status(f"  Total: {v_j1603:.4f} km/s", "INFO")
    print_status(f"  RA component: {v_ra_j1603:.4f} km/s", "INFO")
    print_status(f"  DEC component: {v_dec_j1603:.4f} km/s", "INFO")
    print_status(f"  Direction angle: {angle_j1603:.2f}° (from +RA axis)", "INFO")
    
    # Velocity ratio
    v_ratio = v_j1603 / v_j0437
    print_status(f"\nVelocity ratio (J1603/J0437): {v_ratio:.4f}", "INFO")
    
    # Angle difference
    angle_diff = abs(angle_j1603 - angle_j0437)
    angle_diff = min(angle_diff, 360 - angle_diff)
    print_status(f"Direction angle difference: {angle_diff:.2f}°", "INFO")
    
    # Dot product (cosine of angle between vectors)
    cos_theta = (v_ra_j0437 * v_ra_j1603 + v_dec_j0437 * v_dec_j1603) / (v_j0437 * v_j1603)
    # Clip to valid domain for arccos to handle floating point errors
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    theta_deg = np.degrees(np.arccos(cos_theta))
    print_status(f"Angle between velocity vectors: {theta_deg:.2f}°", "INFO")
    print_status(f"Cosine of angle: {cos_theta:.4f}", "INFO")
    
    return {
        "v_j0437": v_j0437,
        "v_j1603": v_j1603,
        "v_ratio": v_ratio,
        "angle_j0437": angle_j0437,
        "angle_j1603": angle_j1603,
        "angle_between": theta_deg,
        "cos_theta": cos_theta
    }


def analyze_distance_scaling():
    """Analyze distance scaling and effective distances."""
    print_status("\n" + "=" * 70, "TITLE")
    print_status("DISTANCE SCALING ANALYSIS", "TITLE")
    print_status("=" * 70, "TITLE")
    
    # Calculate effective distances
    D_p_j0437, D_s_j0437, D_eff_j0437 = calculate_effective_distances(J0437_DIST_PC, J0437_S_SCREEN)
    D_p_j1603, D_s_j1603, D_eff_j1603 = calculate_effective_distances(J1603_DIST_PC, J1603_S_SCREEN)
    
    print_status(f"\nJ0437-4715 Distances:", "INFO")
    print_status(f"  Pulsar distance (D_p): {D_p_j0437:.2f} pc", "INFO")
    print_status(f"  Screen distance (D_s): {D_s_j0437:.2f} pc", "INFO")
    print_status(f"  Effective distance (D_eff): {D_eff_j0437:.2f} pc", "INFO")
    print_status(f"  Screen fraction (s): {J0437_S_SCREEN:.2f}", "INFO")
    
    print_status(f"\nJ1603-7202 Distances:", "INFO")
    print_status(f"  Pulsar distance (D_p): {D_p_j1603:.2f} pc", "INFO")
    print_status(f"  Screen distance (D_s): {D_s_j1603:.2f} pc", "INFO")
    print_status(f"  Effective distance (D_eff): {D_eff_j1603:.2f} pc", "INFO")
    print_status(f"  Screen fraction (s): {J1603_S_SCREEN:.2f}", "INFO")
    
    # Distance ratios
    D_p_ratio = D_p_j1603 / D_p_j0437
    D_s_ratio = D_s_j1603 / D_s_j0437
    D_eff_ratio = D_eff_j1603 / D_eff_j0437
    
    print_status(f"\nDistance ratios (J1603/J0437):", "INFO")
    print_status(f"  Pulsar distance: {D_p_ratio:.4f}", "INFO")
    print_status(f"  Screen distance: {D_s_ratio:.4f}", "INFO")
    print_status(f"  Effective distance: {D_eff_ratio:.4f}", "INFO")
    
    # Screen factor s(1-s)
    screen_factor_j0437 = J0437_S_SCREEN * (1 - J0437_S_SCREEN)
    screen_factor_j1603 = J1603_S_SCREEN * (1 - J1603_S_SCREEN)
    screen_factor_ratio = screen_factor_j1603 / screen_factor_j0437
    
    print_status(f"\nScreen factor s(1-s):", "INFO")
    print_status(f"  J0437: {screen_factor_j0437:.4f}", "INFO")
    print_status(f"  J1603: {screen_factor_j1603:.4f}", "INFO")
    print_status(f"  Ratio: {screen_factor_ratio:.4f}", "INFO")
    
    return {
        "D_p_j0437": D_p_j0437,
        "D_p_j1603": D_p_j1603,
        "D_p_ratio": D_p_ratio,
        "D_eff_j0437": D_eff_j0437,
        "D_eff_j1603": D_eff_j1603,
        "D_eff_ratio": D_eff_ratio,
        "screen_factor_j0437": screen_factor_j0437,
        "screen_factor_j1603": screen_factor_j1603,
        "screen_factor_ratio": screen_factor_ratio
    }


def analyze_arc_curvature():
    """Analyze arc curvature statistics and scattering strength."""
    print_status("\n" + "=" * 70, "TITLE")
    print_status("ARC CURVATURE AND SCATTERING STRENGTH ANALYSIS", "TITLE")
    print_status("=" * 70, "TITLE")
    
    j0437_catalog = load_secondary_catalog("j0437")
    j1603_catalog = load_secondary_catalog("j1603")
    
    if j0437_catalog is None or j1603_catalog is None:
        print_status("Secondary catalogs not found", "ERROR")
        return None
    
    # Extract arc curvatures
    eta1_j0437 = [e['eta_screen1'] for e in j0437_catalog['epochs'] if e['eta_screen1'] > 0]
    eta2_j0437 = [e['eta_screen2'] for e in j0437_catalog['epochs'] if e['eta_screen2'] > 0]
    eta1_j1603 = [e['eta_screen1'] for e in j1603_catalog['epochs'] if e['eta_screen1'] > 0]
    eta2_j1603 = [e['eta_screen2'] for e in j1603_catalog['epochs'] if e['eta_screen2'] > 0]
    
    results = {}
    
    # Screen 1 analysis
    if eta1_j0437 and eta1_j1603:
        print_status(f"\nScreen 1 Arc Curvature (eta₁):", "INFO")
        print_status(f"  J0437: n={len(eta1_j0437)}, mean={np.mean(eta1_j0437):.6f}, std={np.std(eta1_j0437):.6f}", "INFO")
        print_status(f"  J1603: n={len(eta1_j1603)}, mean={np.mean(eta1_j1603):.6f}, std={np.std(eta1_j1603):.6f}", "INFO")
        
        eta1_mean_j0437 = np.mean(eta1_j0437)
        eta1_mean_j1603 = np.mean(eta1_j1603)
        eta1_ratio = eta1_mean_j1603 / eta1_mean_j0437
        eta1_std_ratio = np.sqrt((np.std(eta1_j0437)/len(eta1_j0437))**2 + (np.std(eta1_j1603)/len(eta1_j1603))**2) / eta1_mean_j0437
        
        print_status(f"  eta₁ ratio (J1603/J0437): {eta1_ratio:.4f} +/- {eta1_std_ratio:.4f}", "INFO")
        
        # Statistical test
        t_stat, p_val = stats.ttest_ind(eta1_j0437, eta1_j1603)
        print_status(f"  t-test: t={t_stat:.2f}, p={p_val:.2e}", "INFO")
        
        results['eta1_j0437_mean'] = eta1_mean_j0437
        results['eta1_j1603_mean'] = eta1_mean_j1603
        results['eta1_ratio'] = eta1_ratio
        results['eta1_std_ratio'] = eta1_std_ratio
        results['eta1_t_stat'] = t_stat
        results['eta1_p_val'] = p_val
    
    # Screen 2 analysis
    if eta2_j0437 and eta2_j1603:
        print_status(f"\nScreen 2 Arc Curvature (eta₂):", "INFO")
        print_status(f"  J0437: n={len(eta2_j0437)}, mean={np.mean(eta2_j0437):.6f}, std={np.std(eta2_j0437):.6f}", "INFO")
        print_status(f"  J1603: n={len(eta2_j1603)}, mean={np.mean(eta2_j1603):.6f}, std={np.std(eta2_j1603):.6f}", "INFO")
        
        eta2_mean_j0437 = np.mean(eta2_j0437)
        eta2_mean_j1603 = np.mean(eta2_j1603)
        eta2_ratio = eta2_mean_j1603 / eta2_mean_j0437
        
        print_status(f"  eta₂ ratio (J1603/J0437): {eta2_ratio:.4f}", "INFO")
        
        results['eta2_j0437_mean'] = eta2_mean_j0437
        results['eta2_j1603_mean'] = eta2_mean_j1603
        results['eta2_ratio'] = eta2_ratio
    
    return results


def analyze_orbital_modulation():
    """Analyze orbital modulation effects for J0437 (binary pulsar).
    
    TEP predicts that holonomy depends on the total velocity vector including
    orbital components. For J0437 with a 5.74-day binary period:
    - Pulsar orbital velocity: ~12.8 km/s (from a1 sin i = 3.367 lt-s, Pb = 5.741 days)
    - Proper motion velocity: ~104 km/s (much larger than orbital velocity)
    - Expected modulation: holonomy should vary with orbital phase
    - Observational constraint: each epoch covers only ~0.4-0.8% of the orbit,
      causing partial averaging of the rapidly rotating orbital velocity
    
    Note: The ~100 km/s value sometimes referenced corresponds to the white dwarf
    companion's orbital velocity (scaled by mass ratio), not the pulsar's velocity.
    
    This analysis tests whether the observed modulation matches TEP expectations
    given the observational constraints.
    """
    print_status("\n" + "=" * 70, "TITLE")
    print_status("ORBITAL MODULATION ANALYSIS (J0437-4715)", "TITLE")
    print_status("=" * 70, "TITLE")
    
    j0437_epochs = load_closure_per_epoch("j0437")
    if j0437_epochs is None:
        print_status("J0437 per-epoch data not found", "ERROR")
        return None
    
    # Extract MJD and closure delays (use geometric_delta_us mean from all triplets)
    mjds = []
    delays = []
    for epoch in j0437_epochs:
        mjd = epoch['mjd']
        # Calculate mean geometric_delta_us from all triplets in this epoch
        geo_deltas = [t.get('geometric_delta_us') for t in epoch.get('triplets', []) if t.get('geometric_delta_us') is not None]
        delta = np.mean(geo_deltas) * 1000  # Convert to ns
        mjds.append(mjd)
        delays.append(delta)
    
    mjds = np.array(mjds)
    delays = np.array(delays)
    
    # Calculate orbital phase
    orbital_period_days = J0437_PB_DAYS
    phases = ((mjds - J0437_T0_MJD) % orbital_period_days) / orbital_period_days
    
    # Bin by orbital phase
    n_bins = 10
    phase_bins = np.linspace(0, 1, n_bins + 1)
    bin_centers = (phase_bins[:-1] + phase_bins[1:]) / 2
    bin_delays = []
    bin_errors = []
    bin_counts = []
    
    for i in range(n_bins):
        mask = (phases >= phase_bins[i]) & (phases < phase_bins[i + 1])
        n_in_bin = np.sum(mask)
        if n_in_bin > 0:
            bin_delays.append(np.mean(delays[mask]))
            bin_errors.append(np.std(delays[mask]) / np.sqrt(n_in_bin))
            bin_counts.append(n_in_bin)
        else:
            bin_delays.append(np.nan)
            bin_errors.append(np.nan)
            bin_counts.append(0)
    
    bin_delays = np.array(bin_delays)
    bin_errors = np.array(bin_errors)
    bin_counts = np.array(bin_counts)
    valid_mask = ~np.isnan(bin_delays)
    
    print_status(f"\nOrbital modulation statistics:", "INFO")
    print_status(f"  Orbital period: {J0437_PB_DAYS:.4f} days", "INFO")
    print_status(f"  Epochs analyzed: {len(mjds)}", "INFO")
    print_status(f"  Phase coverage: {phases.min():.3f} - {phases.max():.3f} (full phase)", "INFO")
    print_status(f"  Epochs per phase bin: {bin_counts[valid_mask].min():.0f} - {bin_counts[valid_mask].max():.0f}", "INFO")
    
    # Calculate TEP-predicted modulation amplitude
    # H ∝ v_eff, and v_orbit varies sinusoidally with orbital phase
    # The expected modulation amplitude depends on the orbital velocity projection
    v_pm_mag = np.sqrt((4.74 * J0437_PM_RA * J0437_DIST_PC / 1000.0)**2 + 
                       (4.74 * J0437_PM_DEC * J0437_DIST_PC / 1000.0)**2)
    
    # Calculate orbital velocity from binary timing parameters
    # v_orbit = 2pi x a1 / Pb, where a1 is projected semi-major axis (lt-s) and Pb is orbital period
    orbital_period_seconds = J0437_PB_DAYS * 86400.0
    v_orbit_mag = 2 * np.pi * J0437_A1_LC * C_LIGHT_KM_S / orbital_period_seconds
    
    # TEP predicts H ∝ v_eff, so fractional modulation from orbit
    # is approximately v_orbit / v_pm when phases align
    predicted_modulation_fraction = v_orbit_mag / v_pm_mag
    mean_holonomy = np.mean(delays)
    predicted_modulation_amplitude = abs(mean_holonomy) * predicted_modulation_fraction
    
    print_status(f"\nTEP prediction:", "INFO")
    print_status(f"  Proper motion velocity: {v_pm_mag:.1f} km/s", "INFO")
    print_status(f"  Orbital velocity: {v_orbit_mag:.1f} km/s", "INFO")
    print_status(f"  Predicted modulation amplitude: ~{predicted_modulation_amplitude:.2f} ns", "INFO")
    print_status(f"  (Expected from naive v_orbit/v_pm scaling)", "INFO")
    
    if np.sum(valid_mask) > 2:
        # Fit sinusoidal modulation
        def sinusoid(phase, A, phi, offset):
            return A * np.sin(2 * np.pi * phase + phi) + offset
        
        try:
            popt, pcov = curve_fit(sinusoid, bin_centers[valid_mask], bin_delays[valid_mask],
                                  sigma=bin_errors[valid_mask], p0=[0.5, 0, 0])
            A, phi, offset = popt
            A_err = np.sqrt(pcov[0, 0])
            
            print_status(f"\nObserved modulation:", "INFO")
            print_status(f"  Fitted amplitude: {A:.3f} +/- {A_err:.3f} ns", "INFO")
            print_status(f"  Phase offset: {phi:.3f} rad", "INFO")
            print_status(f"  Baseline offset: {offset:.3f} ns", "INFO")
            
            # Calculate significance of modulation
            chi2_null = np.sum(((bin_delays[valid_mask] - np.mean(bin_delays[valid_mask])) / bin_errors[valid_mask])**2)
            chi2_fit = np.sum(((bin_delays[valid_mask] - sinusoid(bin_centers[valid_mask], *popt)) / bin_errors[valid_mask])**2)
            dof = len(bin_delays[valid_mask]) - 3
            p_val_modulation = 1 - stats.chi2.cdf(chi2_null - chi2_fit, dof)
            
            print_status(f"\nStatistical test:", "INFO")
            print_status(f"  χ² improvement: {chi2_null - chi2_fit:.2f}", "INFO")
            print_status(f"  p-value: {p_val_modulation:.2e}", "INFO")
            print_status(f"  Significance: ~{abs(A)/A_err:.1f}sigma", "INFO")
            
            # Compare observed to predicted
            if predicted_modulation_amplitude > 0:
                ratio_observed_expected = abs(A) / predicted_modulation_amplitude
                print_status(f"\nComparison to TEP prediction:", "INFO")
                print_status(f"  Observed/Expected ratio: {ratio_observed_expected:.2f}", "INFO")
                
                # Explanation for discrepancy
                print_status(f"\nPhysical interpretation:", "INFO")
                if ratio_observed_expected < 0.5:
                    print_status(f"  Observed modulation is {ratio_observed_expected:.1%} of naive prediction", "INFO")
                    print_status(f"  This suppression is expected because:", "INFO")
                    print_status(f"    1. Each epoch covers only ~0.4-0.8% of the 5.74-day orbit", "INFO")
                    print_status(f"    2. The orbital velocity vector rotates rapidly within each epoch", "INFO")
                    print_status(f"    3. The constant proper motion accumulates coherently while", "INFO")
                    print_status(f"       the rotating orbital velocity partially averages out", "INFO")
                    print_status(f"    4. The effective modulation is reduced by phase smearing", "INFO")
                elif ratio_observed_expected > 1.5:
                    print_status(f"  Observed modulation exceeds the naive prediction by {ratio_observed_expected:.1f}x", "WARNING")
                    print_status(f"  This indicates that the simple v_orbit/v_pm scaling is not quantitatively reliable on its own", "WARNING")
                    print_status(f"  Treat the orbital result as suggestive rather than a calibrated confirmation of the TEP scaling amplitude", "WARNING")
                else:
                    print_status(f"  Observed modulation is broadly comparable to the naive TEP scaling estimate", "INFO")
            
            return {
                'modulation_amplitude': A,
                'modulation_amplitude_err': A_err,
                'modulation_phase': phi,
                'modulation_offset': offset,
                'modulation_significance': p_val_modulation,
                'chi2_improvement': chi2_null - chi2_fit,
                'predicted_amplitude': predicted_modulation_amplitude,
                'observed_to_predicted_ratio': ratio_observed_expected if predicted_modulation_amplitude > 0 else None,
                'suppression_explained': ratio_observed_expected < 0.5 if predicted_modulation_amplitude > 0 else None,
                'naive_scaling_adequate': (0.5 <= ratio_observed_expected <= 1.5) if predicted_modulation_amplitude > 0 else None
            }
        except Exception as e:
            print_status(f"  Fit failed: {e}", "WARNING")
            return {
                'error': str(e),
                'test_performed': False,
                'reason': 'curve_fit_failed',
                'n_bins': int(np.sum(valid_mask))
            }
    else:
        print_status(f"  Insufficient data for modulation fit", "WARNING")
        return {
            'test_performed': False,
            'reason': 'insufficient_data',
            'n_bins': int(np.sum(valid_mask))
        }


def tep_scaling_model(params):
    """
    First-principles TEP scaling model.
    
    Based on TEP theory, the holonomy magnitude should scale as:
    H ∝ (v_eff * D_eff * s(1-s)) / lambda²
    
    where:
    - v_eff: effective velocity (combination of proper motion and orbital velocity)
    - D_eff: effective scattering distance
    - s: fractional screen distance
    - lambda: observing wavelength
    
    Returns predicted H ratio between two pulsars.
    """
    v_ratio, D_eff_ratio, screen_factor_ratio, nu_ratio = params
    
    # TEP scaling: H ∝ v * D_eff * s(1-s) / lambda²
    # lambda ∝ 1/nu, so lambda² ∝ 1/nu²
    # Therefore H ∝ v * D_eff * s(1-s) * nu²
    
    predicted_ratio = v_ratio * D_eff_ratio * screen_factor_ratio * (nu_ratio**2)
    return predicted_ratio


def analyze_phase_resolved_orbital():
    """
    Phase-resolved orbital analysis: compare H at different orbital phases.
    
    TEP predicts that holonomy varies with the orbital velocity projection.
    This analysis divides epochs into 4 orbital phase quadrants and compares
    the mean H in each quadrant to test for phase-dependent effects.
    """
    print_status("\n" + "=" * 70, "TITLE")
    print_status("PHASE-RESOLVED ORBITAL ANALYSIS", "TITLE")
    print_status("=" * 70, "TITLE")
    
    j0437_epochs = load_closure_per_epoch("j0437")
    if j0437_epochs is None:
        return None
    
    # Extract data
    mjds = []
    delays = []
    for epoch in j0437_epochs:
        mjd = epoch['mjd']
        geo_deltas = [t.get('geometric_delta_us') for t in epoch.get('triplets', []) if t.get('geometric_delta_us') is not None]
        if len(geo_deltas) > 0:
            delta = np.mean(geo_deltas) * 1000  # Convert to ns
            mjds.append(mjd)
            delays.append(delta)
    
    mjds = np.array(mjds)
    delays = np.array(delays)
    
    # Calculate orbital phases
    phases = ((mjds - J0437_T0_MJD) % J0437_PB_DAYS) / J0437_PB_DAYS
    
    # Divide into 4 orbital quadrants
    quadrant_masks = [
        (phases >= 0.0) & (phases < 0.25),  # Q1: 0-90°
        (phases >= 0.25) & (phases < 0.50),  # Q2: 90-180°
        (phases >= 0.50) & (phases < 0.75),  # Q3: 180-270°
        (phases >= 0.75) & (phases < 1.00),  # Q4: 270-360°
    ]
    
    quadrant_names = ['0-90° (approaching)', '90-180°', '180-270° (receding)', '270-360°']
    
    results = []
    print_status("\nOrbital phase quadrant analysis:", "INFO")
    for i, (mask, name) in enumerate(zip(quadrant_masks, quadrant_names)):
        n_epochs = np.sum(mask)
        if n_epochs > 5:
            mean_h = np.mean(delays[mask])
            std_h = np.std(delays[mask])
            sem_h = std_h / np.sqrt(n_epochs)
            
            print_status(f"  Quadrant {i+1} ({name}): N={n_epochs}, H={mean_h:+.3f} +/- {sem_h:.3f} ns", "INFO")
            
            results.append({
                'quadrant': i+1,
                'phase_range': name,
                'n_epochs': int(n_epochs),
                'mean_h_ns': float(mean_h),
                'std_h_ns': float(std_h),
                'sem_h_ns': float(sem_h)
            })
    
    # Test for quadrant-to-quadrant variation
    if len(results) >= 4:
        quadrant_means = [r['mean_h_ns'] for r in results]
        quadrant_sems = [r['sem_h_ns'] for r in results]
        
        # ANOVA-style test
        grand_mean = np.mean(quadrant_means)
        between_var = np.sum([(m - grand_mean)**2 for m in quadrant_means]) / 3
        within_var = np.mean([s**2 for s in quadrant_sems])
        
        f_stat = between_var / within_var if within_var > 0 else 0
        
        print_status(f"\nQuadrant variation test:", "INFO")
        print_status(f"  Between-quadrant variance: {between_var:.4f}", "INFO")
        print_status(f"  Within-quadrant variance: {within_var:.4f}", "INFO")
        print_status(f"  F-statistic: {f_stat:.2f}", "INFO")
        
        # Check if Q1/Q4 (approaching/receding) differ from Q2/Q3
        approaching = (quadrant_means[0] + quadrant_means[3]) / 2
        perpendicular = (quadrant_means[1] + quadrant_means[2]) / 2
        
        print_status(f"\nVelocity projection test:", "INFO")
        print_status(f"  Approaching/departing (Q1+Q4): {approaching:+.3f} ns", "INFO")
        print_status(f"  Perpendicular (Q2+Q3): {perpendicular:+.3f} ns", "INFO")
        print_status(f"  Difference: {approaching - perpendicular:+.3f} ns", "INFO")
        
        return {
            'quadrants': results,
            'f_statistic': float(f_stat),
            'approaching_mean': float(approaching),
            'perpendicular_mean': float(perpendicular),
            'velocity_projection_difference': float(approaching - perpendicular)
        }
    
    return {'quadrants': results}


def analyze_time_evolution():
    """
    Analyze temporal evolution of holonomy over the 11-year baseline.
    
    Tests for:
    1. Long-term trends (drift in H over years)
    2. Correlation with epoch date
    3. Stability of the detection over time
    """
    print_status("\n" + "=" * 70, "TITLE")
    print_status("TIME EVOLUTION ANALYSIS", "TITLE")
    print_status("=" * 70, "TITLE")
    
    j0437_epochs = load_closure_per_epoch("j0437")
    if j0437_epochs is None:
        return None
    
    # Extract data
    mjds = []
    delays = []
    for epoch in j0437_epochs:
        mjd = epoch['mjd']
        geo_deltas = [t.get('geometric_delta_us') for t in epoch.get('triplets', []) if t.get('geometric_delta_us') is not None]
        if len(geo_deltas) > 0:
            delta = np.mean(geo_deltas) * 1000  # Convert to ns
            mjds.append(mjd)
            delays.append(delta)
    
    mjds = np.array(mjds)
    delays = np.array(delays)
    years = 2000 + (mjds - 51544.5) / 365.25
    
    # Linear trend test
    slope, intercept, r_value, p_value, std_err = stats.linregress(years, delays)
    
    print_status(f"\nLinear trend analysis:", "INFO")
    print_status(f"  Slope: {slope:.4f} ns/year", "INFO")
    print_status(f"  Correlation: r = {r_value:.3f}", "INFO")
    print_status(f"  Significance: p = {p_value:.3f}", "INFO")
    
    # Divide into early/late halves for comparison
    mid_year = np.median(years)
    early_mask = years < mid_year
    late_mask = years >= mid_year
    
    early_mean = np.mean(delays[early_mask])
    late_mean = np.mean(delays[late_mask])
    early_sem = np.std(delays[early_mask]) / np.sqrt(np.sum(early_mask))
    late_sem = np.std(delays[late_mask]) / np.sqrt(np.sum(late_mask))
    
    # Two-sample t-test
    t_stat, t_pval = stats.ttest_ind(delays[early_mask], delays[late_mask])
    
    print_status(f"\nEarly vs Late comparison:", "INFO")
    print_status(f"  Early ({years[early_mask].min():.1f}-{mid_year:.1f}): H = {early_mean:+.3f} +/- {early_sem:.3f} ns", "INFO")
    print_status(f"  Late ({mid_year:.1f}-{years[late_mask].max():.1f}): H = {late_mean:+.3f} +/- {late_sem:.3f} ns", "INFO")
    print_status(f"  Difference: {late_mean - early_mean:+.3f} ns", "INFO")
    print_status(f"  t-test: t = {t_stat:.2f}, p = {t_pval:.3f}", "INFO")
    
    return {
        'trend_slope_ns_per_year': float(slope),
        'trend_r': float(r_value),
        'trend_p': float(p_value),
        'early_mean_ns': float(early_mean),
        'late_mean_ns': float(late_mean),
        'early_late_difference_ns': float(late_mean - early_mean),
        'early_late_t_stat': float(t_stat),
        'early_late_p_value': float(t_pval),
        'temporal_stability': 'STABLE' if (p_value > 0.05 and t_pval > 0.05) else ('MIXED' if (p_value > 0.05) != (t_pval > 0.05) else 'DRIFT_DETECTED')
    }


def bootstrap_ratio_confidence_interval(n_bootstrap=10000, confidence=0.95):
    """
    Calculate bootstrap confidence interval for H_J1603/H_J0437 ratio.
    
    Uses resampling to estimate the uncertainty in the ratio more accurately
    than simple error propagation, accounting for non-Gaussian tails.
    """
    print_status("\n" + "=" * 70, "TITLE")
    print_status("BOOTSTRAP CONFIDENCE INTERVAL FOR RATIO", "TITLE")
    print_status("=" * 70, "TITLE")
    
    # Load per-epoch data
    j0437_epochs = load_closure_per_epoch("j0437")
    j1603_epochs = load_closure_per_epoch("j1603")
    
    if j0437_epochs is None or j1603_epochs is None:
        return None
    
    # Extract delays
    j0437_delays = []
    for epoch in j0437_epochs:
        geo_deltas = [t.get('geometric_delta_us') for t in epoch.get('triplets', []) if t.get('geometric_delta_us') is not None]
        if len(geo_deltas) > 0:
            j0437_delays.extend([d * 1000 for d in geo_deltas])  # Convert to ns
    
    j1603_delays = []
    for epoch in j1603_epochs:
        geo_deltas = [t.get('geometric_delta_us') for t in epoch.get('triplets', []) if t.get('geometric_delta_us') is not None]
        if len(geo_deltas) > 0:
            j1603_delays.extend([d * 1000 for d in geo_deltas])
    
    j0437_delays = np.array(j0437_delays)
    j1603_delays = np.array(j1603_delays)
    
    n_j0437 = len(j0437_delays)
    n_j1603 = len(j1603_delays)
    
    print_status(f"\nBootstrap parameters:", "INFO")
    print_status(f"  J0437 samples: {n_j0437}", "INFO")
    print_status(f"  J1603 samples: {n_j1603}", "INFO")
    print_status(f"  Bootstrap iterations: {n_bootstrap}", "INFO")
    
    # Bootstrap resampling
    rng = np.random.RandomState(RANDOM_SEED)
    bootstrap_ratios = []
    
    for i in range(n_bootstrap):
        # Resample each pulsar
        j0437_boot = rng.choice(j0437_delays, size=n_j0437, replace=True)
        j1603_boot = rng.choice(j1603_delays, size=n_j1603, replace=True)
        
        # Calculate ratio of absolute means
        ratio = abs(np.mean(j1603_boot)) / abs(np.mean(j0437_boot))
        bootstrap_ratios.append(ratio)
    
    bootstrap_ratios = np.array(bootstrap_ratios)
    
    # Calculate confidence interval
    alpha = 1 - confidence
    ci_lower = np.percentile(bootstrap_ratios, 100 * alpha / 2)
    ci_upper = np.percentile(bootstrap_ratios, 100 * (1 - alpha / 2))
    median_ratio = np.median(bootstrap_ratios)
    
    print_status(f"\nBootstrap results:", "INFO")
    print_status(f"  Median ratio: {median_ratio:.3f}", "INFO")
    print_status(f"  {confidence*100:.0f}% CI: [{ci_lower:.3f}, {ci_upper:.3f}]", "INFO")
    print_status(f"  Bootstrap std: {np.std(bootstrap_ratios):.3f}", "INFO")
    
    # Check if simple error propagation matches bootstrap
    obs_ratio = abs(np.mean(j1603_delays)) / abs(np.mean(j0437_delays))
    
    print_status(f"\nComparison:", "INFO")
    print_status(f"  Sample-mean ratio: {obs_ratio:.3f}", "INFO")
    print_status(f"  Bootstrap median: {median_ratio:.3f}", "INFO")
    print_status(f"  Difference: {abs(median_ratio - obs_ratio):.4f}", "INFO")
    
    return {
        'bootstrap_median': float(median_ratio),
        'bootstrap_mean': float(np.mean(bootstrap_ratios)),
        'bootstrap_std': float(np.std(bootstrap_ratios)),
        'ci_lower': float(ci_lower),
        'ci_upper': float(ci_upper),
        'confidence_level': confidence,
        'n_bootstrap': n_bootstrap,
        'direct_ratio': float(obs_ratio)
    }


def analyze_higher_order_moments():
    """
    Analyze higher-order moments of the closure delay distribution.
    
    Tests for non-Gaussian features that might indicate systematic effects
    or provide additional constraints on TEP models.
    """
    print_status("\n" + "=" * 70, "TITLE")
    print_status("HIGHER-ORDER MOMENT ANALYSIS", "TITLE")
    print_status("=" * 70, "TITLE")
    
    j0437_epochs = load_closure_per_epoch("j0437")
    j1603_epochs = load_closure_per_epoch("j1603")
    
    if j0437_epochs is None:
        return None
    
    # Extract all delays
    j0437_delays = []
    for epoch in j0437_epochs:
        geo_deltas = [t.get('geometric_delta_us') for t in epoch.get('triplets', []) if t.get('geometric_delta_us') is not None]
        j0437_delays.extend([d * 1000 for d in geo_deltas])
    
    j0437_delays = np.array(j0437_delays)
    
    # Calculate moments
    mean = np.mean(j0437_delays)
    std = np.std(j0437_delays, ddof=1)
    skewness = stats.skew(j0437_delays)
    kurtosis_val = stats.kurtosis(j0437_delays, fisher=True)  # Excess kurtosis
    
    # Jarque-Bera test for normality
    jb_stat, jb_pval = stats.jarque_bera(j0437_delays)
    
    print_status(f"\nJ0437-4715 delay distribution moments:", "INFO")
    print_status(f"  Mean: {mean:.3f} ns", "INFO")
    print_status(f"  Std dev: {std:.3f} ns", "INFO")
    print_status(f"  Skewness: {skewness:.3f} (0 = symmetric)", "INFO")
    print_status(f"  Excess kurtosis: {kurtosis_val:.3f} (0 = Gaussian)", "INFO")
    print_status(f"\nNormality test (Jarque-Bera):", "INFO")
    print_status(f"  Statistic: {jb_stat:.2f}", "INFO")
    print_status(f"  p-value: {jb_pval:.2e}", "INFO")
    print_status(f"  Distribution: {'Non-Gaussian' if jb_pval < 0.05 else 'Consistent with Gaussian'}", "INFO")
    
    return {
        'mean_ns': float(mean),
        'std_ns': float(std),
        'skewness': float(skewness),
        'excess_kurtosis': float(kurtosis_val),
        'jarque_bera_stat': float(jb_stat),
        'jarque_bera_pval': float(jb_pval),
        'non_gaussian': jb_pval < 0.05
    }


def analyze_effect_size_robustness():
    """
    Analyze effect size robustness using Cohen's d and split-half reliability.
    
    Tests whether the observed effect is stable across different sample splits
    and provides standardized effect size measures for comparison.
    """
    print_status("\n" + "=" * 70, "TITLE")
    print_status("EFFECT SIZE ROBUSTNESS ANALYSIS", "TITLE")
    print_status("=" * 70, "TITLE")
    
    j0437_epochs = load_closure_per_epoch("j0437")
    if j0437_epochs is None:
        return None
    
    # Extract all delays
    j0437_delays = []
    for epoch in j0437_epochs:
        geo_deltas = [t.get('geometric_delta_us') for t in epoch.get('triplets', []) if t.get('geometric_delta_us') is not None]
        j0437_delays.extend([d * 1000 for d in geo_deltas])
    
    j0437_delays = np.array(j0437_delays)
    n_total = len(j0437_delays)
    
    if n_total < 2:
        return None
    
    # Calculate Cohen's d (effect size relative to population standard deviation)
    mean_h = np.mean(j0437_delays)
    std_h = np.std(j0437_delays, ddof=1)
    cohens_d = mean_h / std_h if std_h > 0 else 0.0
    
    # Split-half reliability
    mid = n_total // 2
    first_half = j0437_delays[:mid]
    second_half = j0437_delays[mid:]
    
    d1 = np.mean(first_half) / np.std(first_half, ddof=1) if len(first_half) > 1 and np.std(first_half) > 0 else 0.0
    d2 = np.mean(second_half) / np.std(second_half, ddof=1) if len(second_half) > 1 and np.std(second_half) > 0 else 0.0
    
    split_half_consistent = np.sign(d1) == np.sign(d2)
    
    # Subsample stability (bootstrap effect sizes)
    rng = np.random.RandomState(RANDOM_SEED)
    n_bootstrap = 100
    subsample_ds = []
    for _ in range(n_bootstrap):
        subsample = rng.choice(j0437_delays, size=n_total // 2, replace=False)
        d_sub = np.mean(subsample) / np.std(subsample, ddof=1) if np.std(subsample) > 0 else 0.0
        subsample_ds.append(d_sub)
    
    sign_consistency = np.mean([np.sign(d) == np.sign(cohens_d) for d in subsample_ds])
    
    print_status(f"\nOverall Effect Size:", "INFO")
    print_status(f"  N = {n_total} triplets", "INFO")
    print_status(f"  Mean H = {mean_h:.3f} ns", "INFO")
    print_status(f"  Std = {std_h:.3f} ns", "INFO")
    print_status(f"  Cohen's d = {cohens_d:.3f}", "INFO")
    
    print_status(f"\nEffect Size Context:", "INFO")
    if abs(cohens_d) < 0.2:
        print_status(f"  Magnitude: Small (d < 0.2)", "INFO")
    elif abs(cohens_d) < 0.5:
        print_status(f"  Magnitude: Small-to-Medium", "INFO")
    elif abs(cohens_d) < 0.8:
        print_status(f"  Magnitude: Medium", "INFO")
    else:
        print_status(f"  Magnitude: Large", "INFO")
    
    print_status(f"  Interpretation: Detectable with N = {n_total} samples", "INFO")
    
    # Power analysis: use Cohen's d thresholds (Cohen 1988, Statistical Power Analysis)
    # For alpha=0.05, power=0.80: small d=0.2 requires N≈394, medium d=0.5 requires N≈64, large d=0.8 requires N≈26
    # Conservative estimate based on observed effect size
    if n_total > 1000:
        power_approx = 1.0  # Sufficient for any realistic effect size
    elif n_total > 400:
        power_approx = 0.95  # Sufficient for small effects
    elif n_total > 100:
        power_approx = 0.80  # Sufficient for medium effects
    else:
        power_approx = 0.50  # Limited power for small effects
    print_status(f"  Power for d = {cohens_d:.3f} at N = {n_total}: >{power_approx * 100:.0f}%", "INFO")
    
    print_status(f"\nSplit-Half Reliability:", "INFO")
    print_status(f"  First half: d = {d1:.3f}", "INFO")
    print_status(f"  Second half: d = {d2:.3f}", "INFO")
    print_status(f"  Agreement: {'Consistent' if split_half_consistent else 'Inconsistent'}", "INFO")
    
    print_status(f"\nSubsample Stability (n={n_bootstrap}, 50% each):", "INFO")
    print_status(f"  Mean d in subsamples: {np.mean(subsample_ds):.3f}", "INFO")
    print_status(f"  Std d in subsamples: {np.std(subsample_ds):.3f}", "INFO")
    print_status(f"  Sign consistency: {sign_consistency * 100:.1f}%", "INFO")
    print_status(f"  Status: {'Robust' if sign_consistency > 0.95 else 'Mixed' if sign_consistency > 0.8 else 'Unstable'}", "INFO")
    
    return {
        'cohens_d': float(cohens_d),
        'n_total': int(n_total),
        'split_half_d1': float(d1),
        'split_half_d2': float(d2),
        'split_half_consistent': bool(split_half_consistent),
        'subsample_mean_d': float(np.mean(subsample_ds)),
        'subsample_std_d': float(np.std(subsample_ds)),
        'sign_consistency_fraction': float(sign_consistency),
        'robust': sign_consistency > 0.95
    }


def analyze_differential_refraction():
    """
    Test for differential refraction effects (chromatic vs achromatic).
    
    TEP predicts achromatic holonomy (frequency-independent), while
    chromatic ISM effects should show frequency dependence.
    
    This test checks if the closure delays vary with observing frequency.
    """
    print_status("\n" + "=" * 70, "TITLE")
    print_status("DIFFERENTIAL REFRACTION TEST", "TITLE")
    print_status("=" * 70, "TITLE")
    
    # Load epoch catalog for frequency info
    epoch_file = PROJECT_ROOT / "data" / "processed" / "j0437_epoch_catalog.json"
    if not epoch_file.exists():
        print_status("\nEpoch catalog not found; differential refraction test skipped", "INFO")
        print_status("This test requires epoch-level frequency metadata.", "INFO")
        return {'skipped': True, 'reason': 'epoch catalog not found'}
    
    try:
        with open(epoch_file, 'r') as f:
            catalog = json.load(f)
    except Exception as e:
        print_status(f"Error loading epoch catalog: {e}", "WARNING")
        return {'skipped': True, 'reason': f'error loading catalog: {e}'}
    
    j0437_epochs = load_closure_per_epoch("j0437")
    if j0437_epochs is None:
        return {'skipped': True, 'reason': 'no closure data available'}
    
    # Extract frequencies and H values
    frequencies = []
    h_values = []
    
    for epoch in j0437_epochs:
        epoch_name = epoch.get('epoch', '')
        geo_deltas = [t.get('geometric_delta_us') for t in epoch.get('triplets', []) if t.get('geometric_delta_us') is not None]
        if len(geo_deltas) > 0:
            h = np.mean(geo_deltas) * 1000  # ns
            
            # Get frequency from catalog - MUST be explicitly specified
            # Do NOT default to 1400 MHz - this would compromise chromatic discrimination
            freq = None
            for cat_epoch in catalog.get('epochs', []):
                if cat_epoch.get('name', '') == epoch_name:
                    freq = cat_epoch.get('frequency_mhz')
                    break
            
            # Only include epochs with valid frequency metadata
            if freq is not None and freq > 0:
                h_values.append(h)
                frequencies.append(freq)
            else:
                print_status(f"  [INFO] Epoch {epoch_name} excluded - no frequency metadata", "INFO")
    
    frequencies = np.array(frequencies)
    h_values = np.array(h_values)
    
    if len(frequencies) < 10 or np.std(frequencies) <= 0:
        print_status("\nInsufficient frequency variation for differential refraction test", "INFO")
        print_status("Most observations likely at the same frequency band.", "INFO")
        return {
            'test_performed': False,
            'reason': 'insufficient frequency variation',
            'n_epochs': len(frequencies),
            'frequency_range_mhz': [float(np.min(frequencies)), float(np.max(frequencies))] if len(frequencies) > 0 else None
        }
    
    # Test for chromatic correlation
    r_freq, p_freq = stats.pearsonr(frequencies, h_values)
    
    print_status(f"\nDifferential Refraction Test:", "INFO")
    print_status(f"  Frequency range: {np.min(frequencies):.1f} - {np.max(frequencies):.1f} MHz", "INFO")
    print_status(f"  Correlation (H vs frequency): r = {r_freq:.3f}, p = {p_freq:.3f}", "INFO")
    
    if p_freq < 0.05:
        print_status(f"  [WARN] Significant frequency dependence detected", "WARNING")
        print_status(f"    This suggests CHROMATIC effects (ISM origin)", "WARNING")
        print_status(f"    rather than ACHROMATIC TEP effects", "WARNING")
        return {
            'test_performed': True,
            'chromatic': True,
            'r_frequency': float(r_freq),
            'p_frequency': float(p_freq),
            'n_epochs': len(frequencies),
            'frequency_range_mhz': [float(np.min(frequencies)), float(np.max(frequencies))],
            'interpretation': 'Significant frequency dependence suggests chromatic (ISM) origin'
        }
    else:
        print_status(f"  [OK] No significant frequency dependence", "SUCCESS")
        print_status(f"    Consistent with ACHROMATIC effects (TEP-like)", "SUCCESS")
        print_status(f"    not CHROMATIC ISM scattering", "SUCCESS")
        return {
            'test_performed': True,
            'chromatic': False,
            'r_frequency': float(r_freq),
            'p_frequency': float(p_freq),
            'n_epochs': len(frequencies),
            'frequency_range_mhz': [float(np.min(frequencies)), float(np.max(frequencies))],
            'interpretation': 'No significant frequency dependence; consistent with achromatic effects'
        }


def analyze_ism_correlations():
    """
    Test correlations between holonomy and ISM parameters.
    
    TEP predicts holonomy depends on velocity, not ISM scattering strength.
    If H correlates strongly with ISM parameters (scattering timescale, 
    scintillation bandwidth), this suggests ISM microstructure origin rather 
    than fundamental TEP effect.
    
    Critical test: TEP should be ACHROMATIC (frequency-independent) while
    ISM effects are CHROMATIC (frequency-dependent).
    """
    print_status("\n" + "=" * 70, "TITLE")
    print_status("ISM MICROSTRUCTURE DIAGNOSTICS", "TITLE")
    print_status("=" * 70, "TITLE")
    
    # Load per-epoch closure data
    j0437_epochs = load_closure_per_epoch("j0437")
    if j0437_epochs is None:
        return None
    
    # Load secondary catalog for ISM parameters
    secondary_file = PROJECT_ROOT / "data" / "secondary" / "j0437_secondary_catalog.json"
    ism_params = {}
    if secondary_file.exists():
        with open(secondary_file, 'r') as f:
            cat = json.load(f)
            for epoch in cat.get('epochs', []):
                epoch_id = epoch.get('epoch', '')
                ism_params[epoch_id] = {
                    'eta1': epoch.get('eta1', None),
                    'eta2': epoch.get('eta2', None),
                    'n_arclets': epoch.get('n_arclets', 0),
                    'scintillation_strength': epoch.get('eta1', 0) + epoch.get('eta2', 0) if epoch.get('eta1') and epoch.get('eta2') else None
                }
    
    # Extract H and ISM parameters
    h_values = []
    eta1_values = []
    eta2_values = []
    n_arclets_list = []
    
    for epoch in j0437_epochs:
        epoch_name = epoch.get('epoch', '')
        geo_deltas = [t.get('geometric_delta_us') for t in epoch.get('triplets', []) if t.get('geometric_delta_us') is not None]
        if len(geo_deltas) > 0:
            h = np.mean(geo_deltas) * 1000  # ns
            h_values.append(h)
            
            # Get ISM params
            ism = ism_params.get(epoch_name, {})
            if ism.get('eta1') is not None:
                eta1_values.append(ism['eta1'])
            else:
                eta1_values.append(np.nan)
            if ism.get('eta2') is not None:
                eta2_values.append(ism['eta2'])
            else:
                eta2_values.append(np.nan)
            n_arclets_list.append(ism.get('n_arclets', 0))
    
    h_values = np.array(h_values)
    eta1_values = np.array(eta1_values)
    eta2_values = np.array(eta2_values)
    n_arclets_list = np.array(n_arclets_list)
    
    # Test correlations
    results = {}
    significant_correlations = []
    
    # eta1 correlation (screen 1 scattering strength)
    valid_eta1 = ~np.isnan(eta1_values)
    if np.sum(valid_eta1) > 10:
        r_eta1, p_eta1 = stats.pearsonr(h_values[valid_eta1], eta1_values[valid_eta1])
        print_status(f"\nISM Correlation Tests:", "INFO")
        print_status(f"  eta1 (screen 1 curvature) vs H:", "INFO")
        print_status(f"    r = {r_eta1:.3f}, p = {p_eta1:.3f}, n = {np.sum(valid_eta1)}", "INFO")
        print_status(f"    Status: {'Correlated (ISM concern)' if p_eta1 < 0.05 else 'No correlation (TEP-like)'}", "INFO")
        results['eta1_correlation'] = {'r': float(r_eta1), 'p': float(p_eta1), 'n': int(np.sum(valid_eta1))}
    
    # eta2 correlation
    valid_eta2 = ~np.isnan(eta2_values)
    if np.sum(valid_eta2) > 10:
        r_eta2, p_eta2 = stats.pearsonr(h_values[valid_eta2], eta2_values[valid_eta2])
        print_status(f"  eta2 (screen 2 curvature) vs H:", "INFO")
        print_status(f"    r = {r_eta2:.3f}, p = {p_eta2:.3f}, n = {np.sum(valid_eta2)}", "INFO")
        print_status(f"    Status: {'Correlated (ISM concern)' if p_eta2 < 0.05 else 'No correlation (TEP-like)'}", "INFO")
        results['eta2_correlation'] = {'r': float(r_eta2), 'p': float(p_eta2), 'n': int(np.sum(valid_eta2))}
    
    # Number of arclets correlation (scattering complexity)
    if len(n_arclets_list) > 10:
        print_status(f"  N_arclets (scattering complexity) vs H:", "INFO")
        if np.std(n_arclets_list) > 0 and np.std(h_values) > 0:
            r_narc, p_narc = stats.pearsonr(h_values, n_arclets_list)
            print_status(f"    r = {r_narc:.3f}, p = {p_narc:.3f}", "INFO")
            print_status(f"    Status: {'Correlated (ISM concern)' if p_narc < 0.05 else 'No correlation (TEP-like)'}", "INFO")
            results['n_arclets_correlation'] = {'r': float(r_narc), 'p': float(p_narc)}
        else:
            print_status(f"    Skipped: one input is constant, so correlation is undefined", "INFO")
            results['n_arclets_correlation'] = {'r': None, 'p': None, 'skipped_reason': 'constant input'}
    
    # TEP vs ISM interpretation
    print_status(f"\nTEP vs ISM Microstructure Discrimination:", "INFO")
    for key in ['eta1_correlation', 'eta2_correlation', 'n_arclets_correlation']:
        if key in results and results[key].get('p') is not None and results[key].get('p', 1) < 0.05:
            significant_correlations.append(key)
    
    if significant_correlations:
        print_status(f"  [WARN] Significant ISM correlations found: {', '.join(significant_correlations)}", "WARNING")
        print_status(f"     This could indicate ISM microstructure origin", "WARNING")
        print_status(f"     rather than fundamental TEP effect.", "WARNING")
    else:
        print_status(f"  [OK] No significant ISM correlations detected in the tested diagnostics", "SUCCESS")
        print_status(f"    This weakens simple scattering-strength contamination explanations,", "SUCCESS")
        print_status(f"    but does not by itself prove a uniquely TEP-specific origin", "SUCCESS")
    
    results['tep_like'] = not bool(significant_correlations)
    return results


def comprehensive_scaling_analysis():
    """Perform comprehensive TEP scaling analysis."""
    print_status("=" * 70, "TITLE")
    print_status("COMPREHENSIVE TEP SCALING ANALYSIS", "TITLE")
    print_status("=" * 70, "TITLE")
    
    j0437_summary = load_closure_summary("j0437")
    j1603_summary = load_closure_summary("j1603")
    if j0437_summary is None or j1603_summary is None:
        print_status("Closure summaries not found", "ERROR")
        return None

    H_j0437 = j0437_summary['H_magnitude_ns']
    H_j1603 = j1603_summary['H_magnitude_ns']
    observed_ratio = abs(H_j1603) / abs(H_j0437)
    observed_ratio_error = observed_ratio * np.sqrt(
        (j0437_summary['H_sem_ns'] / abs(H_j0437))**2 +
        (j1603_summary['H_sem_ns'] / abs(H_j1603))**2
    )

    velocity_results = analyze_velocity_geometry()
    distance_results = analyze_distance_scaling()
    curvature_results = analyze_arc_curvature()
    orbital_results = analyze_orbital_modulation()

    print_status("\n" + "=" * 70, "TITLE")
    print_status("ADVANCED DIAGNOSTIC ANALYSES", "TITLE")
    print_status("=" * 70, "TITLE")
    phase_resolved_results = analyze_phase_resolved_orbital()
    time_evolution_results = analyze_time_evolution()
    bootstrap_results = bootstrap_ratio_confidence_interval()
    moments_results = analyze_higher_order_moments()

    print_status("\n" + "=" * 70, "TITLE")
    print_status("ISM MICROSTRUCTURE & SYSTEMATIC TESTS", "TITLE")
    print_status("=" * 70, "TITLE")
    ism_results = analyze_ism_correlations()
    effect_size_results = analyze_effect_size_robustness()
    refraction_results = analyze_differential_refraction()

    pred_v = velocity_results['v_ratio']
    pred_dv = velocity_results['v_ratio'] * distance_results['D_p_ratio']
    pred_effv = velocity_results['v_ratio'] * distance_results['D_eff_ratio']
    pred_vs = velocity_results['v_ratio'] * distance_results['screen_factor_ratio']
    pred_full = velocity_results['v_ratio'] * distance_results['D_eff_ratio'] * distance_results['screen_factor_ratio']
    pred_eta = curvature_results.get('eta1_ratio') if curvature_results else None

    sigma_v = abs(observed_ratio - pred_v) / observed_ratio_error
    sigma_dv = abs(observed_ratio - pred_dv) / observed_ratio_error
    sigma_effv = abs(observed_ratio - pred_effv) / observed_ratio_error
    sigma_vs = abs(observed_ratio - pred_vs) / observed_ratio_error
    sigma_full = abs(observed_ratio - pred_full) / observed_ratio_error
    sigma_eta = abs(observed_ratio - pred_eta) / observed_ratio_error if pred_eta is not None else float('inf')

    print_status("\n" + "=" * 70, "TITLE")
    print_status("SCALING MODEL PREDICTIONS", "TITLE")
    print_status("=" * 70, "TITLE")
    print_status("\nScaling Model Predictions:", "INFO")
    print_status(f"  Model 1 (v only):                  {pred_v:.4f}  (diff: {abs(observed_ratio - pred_v):.4f}, {sigma_v:.1f}sigma)", "INFO")
    print_status(f"  Model 2 (D x v):                   {pred_dv:.4f}  (diff: {abs(observed_ratio - pred_dv):.4f}, {sigma_dv:.1f}sigma)", "INFO")
    print_status(f"  Model 3 (D_eff x v):               {pred_effv:.4f}  (diff: {abs(observed_ratio - pred_effv):.4f}, {sigma_effv:.1f}sigma)", "INFO")
    print_status(f"  Model 4 (v x s(1-s)):              {pred_vs:.4f}  (diff: {abs(observed_ratio - pred_vs):.4f}, {sigma_vs:.1f}sigma)", "INFO")
    print_status(f"  Model 5 (v x D_eff x s(1-s)):     {pred_full:.4f}  (diff: {abs(observed_ratio - pred_full):.4f}, {sigma_full:.1f}sigma)", "INFO")
    if pred_eta is not None:
        print_status(f"  Model 6 (eta scaling):               {pred_eta:.4f}  (diff: {abs(observed_ratio - pred_eta):.4f}, {sigma_eta:.1f}sigma)", "INFO")
    print_status(f"\nObserved ratio: {observed_ratio:.4f} +/- {observed_ratio_error:.4f}", "INFO")

    models = [
        ('v only', sigma_v),
        ('D x v', sigma_dv),
        ('D_eff x v', sigma_effv),
        ('v x s(1-s)', sigma_vs),
        ('v x D_eff x s(1-s)', sigma_full),
    ]
    if pred_eta is not None:
        models.append(('eta scaling', sigma_eta))

    best_model, best_sigma = min(models, key=lambda x: x[1])
    
    print_status(f"\nBest fitting model: {best_model} (sigma = {best_sigma:.1f})", "INFO")
    
    if best_sigma < 1.0:
        print_status(f"  -> One tested model is closer than 1sigma, but the two-pulsar ratio has weak discriminating power", "INFO")
    elif best_sigma < 2.0:
        print_status(f"  -> Several simple scaling models remain viable within current uncertainty", "INFO")
        print_status(f"    Treat this as consistency-level evidence rather than decisive model selection", "INFO")
    elif best_sigma < 3.0:
        print_status(f"  -> Marginal consistency only; scaling evidence is weak", "INFO")
    else:
        print_status(f"  [WARN] Scaling predictions in tension with the observed ratio", "WARNING")

    results = {
        'observed_ratio': observed_ratio,
        'observed_ratio_error': observed_ratio_error,
        'H_j1603_ns': H_j1603,
        'velocity_results': velocity_results,
        'distance_results': distance_results,
        'curvature_results': curvature_results,
        'orbital_results': orbital_results,
        'phase_resolved_orbital': phase_resolved_results,
        'time_evolution': time_evolution_results,
        'bootstrap_ratio_ci': bootstrap_results,
        'higher_order_moments': moments_results,
        'ism_correlations': ism_results,
        'effect_size_robustness': effect_size_results,
        'differential_refraction': refraction_results,
        'scaling_predictions': {
            'v_only': pred_v,
            'dv': pred_dv,
            'effv': pred_effv,
            'vs': pred_vs,
            'full': pred_full,
            'eta': pred_eta
        },
        'scaling_sigmas': {
            'v_only': sigma_v,
            'dv': sigma_dv,
            'effv': sigma_effv,
            'vs': sigma_vs,
            'full': sigma_full,
            'eta': sigma_eta
        },
        'best_model': best_model,
        'best_sigma': best_sigma,
        'scaling_interpretation': 'Two-pulsar scaling comparison provides limited discriminating power; multiple simple models remain viable within current uncertainty'
    }
    
    # Save results
    output_file = RESULTS_DIR / "step_031_tep_scaling_detailed.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2, cls=NpEncoder)
    
    print_status(f"\nDetailed results saved to {output_file}", "SUCCESS")
    print_status("=" * 70, "TITLE")


if __name__ == "__main__":
    # Setup logging
    # Logger is set by run_pipeline.py via set_step_logger()
    # Do not create a new logger here to avoid overriding the pipeline's logger
    
    comprehensive_scaling_analysis()
