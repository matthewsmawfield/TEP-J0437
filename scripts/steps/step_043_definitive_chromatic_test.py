#!/usr/bin/env python3
"""
Step 043: Chromatic Diagnostic Using Hierarchical Model

This is a diagnostic analysis. It models:

|H|_i = α x (eta_i)^β x (D_eff,i)^γ x (nu_i/nu_0)^δ x telescope_calibration + noise

Where:
- eta_i = arc curvature (physical scattering strength)
- D_eff,i = effective screen distance  
- nu_i = observing frequency
- δ = frequency exponent (THE KEY PARAMETER)
  * δ = 0 -> achromatic (TEP prediction)
  * δ = -2 -> chromatic ISM (nu^-2 scaling)

The telescope calibration term accounts for the +4.0 ns offset between instruments.

Model comparison:
- M_TEP: δ = 0 (fixed, no frequency dependence)
- M_ISM: δ free (estimated from data, expected ≈ -2 if chromatic)
- M_general: δ free with prior including both 0 and -2

Because the cross-pulsar fit uses unsigned |H| amplitudes, including
noise-floor-limited bounding rows, it is not valid as primary evidence for or
against chromaticity.  The within-source J0437 sub-band result is the cleaner
diagnostic, but it has few points and large errors.
"""

from typing import Union, Optional

import sys
import json
import numpy as np
from pathlib import Path
from scipy import stats
from scipy.optimize import minimize

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.json_numpy import NpEncoder

from scripts.utils.logger import print_status

RESULTS_DIR = PROJECT_ROOT / "results"
OUTPUT_FILE = RESULTS_DIR / "step_043_definitive_chromatic_test.json"


def phase_detection_sigma(summary):
    p_values = [
        summary.get('phase_closure_rayleigh_p_unweighted'),
        summary.get('phase_closure_v_p_unweighted'),
        summary.get('phase_closure_rayleigh_p'),
        summary.get('phase_closure_v_p'),
    ]
    p_values = [p for p in p_values if isinstance(p, (int, float)) and p > 0]
    if not p_values:
        return 0.0
    sign = 1.0 if summary.get('phase_closure_mean_rad', 0.0) >= 0 else -1.0
    return float(sign * stats.norm.isf(min(p_values) / 2.0))

def load_complete_data():
    """Load all pulsar data with measured eta and |H|."""
    
    # Load closure results
    # We include sub-bands for J0437 to provide a within-source achromaticity test
    closure_pulsars = [
        ('J0437', 'step_003_closure_final_summary_j0437.json', 1380, 'Parkes'),
        ('J0437_sb0', 'step_003_closure_final_summary_j0437_sb0.json', 1089, 'Parkes'),
        ('J0437_sb1', 'step_003_closure_final_summary_j0437_sb1.json', 1476, 'Parkes'),
        ('J1603', 'step_003_closure_final_summary_j1603.json', 1380, 'Parkes'),
        ('B0329', 'step_003_closure_final_summary_B0329.json', 2250, 'Jiamusi'),
        ('B0355', 'step_003_closure_final_summary_B0355.json', 2250, 'Jiamusi'),
        ('B0540', 'step_003_closure_final_summary_B0540.json', 2250, 'Jiamusi'),
        ('B0740', 'step_003_closure_final_summary_B0740.json', 2250, 'Jiamusi'),
        ('B1933', 'step_003_closure_final_summary_B1933.json', 2250, 'Jiamusi'),
        ('B2154', 'step_003_closure_final_summary_B2154.json', 2250, 'Jiamusi'),
        ('B2310', 'step_003_closure_final_summary_B2310.json', 2250, 'Jiamusi'),
        ('J0908', 'step_003_closure_final_summary_J0908.json', 1300, 'MeerKAT'),
        ('J0922', 'step_003_closure_final_summary_J0922.json', 1300, 'MeerKAT'),
        ('J1731', 'step_003_closure_final_summary_J1731.json', 1300, 'MeerKAT'),
    ]
    
    # Load eta and D_eff from detailed TEP scaling.  The current producer
    # (step_032_detailed_tep_scaling.py) writes step_031_tep_scaling_detailed.json;
    # keep the older name as a fallback for archived runs.
    tep_file = RESULTS_DIR / "step_031_tep_scaling_detailed.json"
    if not tep_file.exists():
        tep_file = RESULTS_DIR / "step_031_detailed_tep_scaling.json"
    tep_data = {}
    if tep_file.exists():
        with open(tep_file) as f:
            tep = json.load(f)
            # Map results to pulsar names
            if 'curvature_results' in tep:
                c = tep['curvature_results']
                d = tep['distance_results']
                tep_data['J0437'] = {
                    'eta': (c['eta1_j0437_mean'] + c['eta2_j0437_mean']) / 2,
                    'D_eff': d['D_eff_j0437'],
                    'v': 104.38
                }
                tep_data['J1603'] = {
                    'eta': (c['eta1_j1603_mean'] + c['eta2_j1603_mean']) / 2,
                    'D_eff': d['D_eff_j1603'],
                    'v': 31.67
                }
    
    # Load Jiamusi specific geometry
    jiamusi_file = RESULTS_DIR / "step_031_jiamusi_eta_analysis.json"
    jiamusi_data = {}
    if jiamusi_file.exists():
        with open(jiamusi_file) as f:
            jiamusi = json.load(f)
            for p, jd in jiamusi.items():
                if isinstance(jd, dict) and jd.get('eta_mean', 0) > 0:
                    s = jd.get('screen_distance_s', 0.5) or 0.5
                    D_p = jd.get('distance_pc', 1000)
                    D_eff = D_p * s * (1 - s)
                    jiamusi_data[p] = {
                        'eta': jd['eta_mean'],
                        'D_eff': D_eff,
                        'v': 50.0
                    }

    # Load MeerKAT geometry from secondary catalog
    meerkat_file = PROJECT_ROOT / "data" / "secondary" / "meerkat_secondary_catalog.json"
    meerkat_data = {}
    if meerkat_file.exists():
        with open(meerkat_file) as f:
            meerkat = json.load(f)
        for epoch in meerkat.get('epochs', []):
            fname = epoch['file']
            # Extract pulsar name from filename (e.g., J1731-4744_secondary.npz)
            pname = fname.split('_')[0]
            if pname not in meerkat_data:
                meerkat_data[pname] = {'etas': [], 'n_arclets': []}
            if epoch.get('eta_screen1', 0) > 0:
                meerkat_data[pname]['etas'].append(epoch['eta_screen1'])
            if epoch.get('eta_screen2', 0) > 0:
                meerkat_data[pname]['etas'].append(epoch['eta_screen2'])
            meerkat_data[pname]['n_arclets'].append(epoch.get('n_arclets', 0))
        # Compute means
        for pname, data in list(meerkat_data.items()):
            if data['etas']:
                eta_mean = np.mean(data['etas'])
                # Approximate screen distance from eta using standard thin-screen model
                # For a rough estimate, use s = 0.5 for unknown geometries
                s = 0.5
                # Get pulsar distance from PULSAR_PARAMS if available
                from scripts.steps.step_003_closure_delays_final import PULSAR_PARAMS
                full_name = pname
                if full_name not in PULSAR_PARAMS:
                    # Try with hyphen if needed
                    for k in PULSAR_PARAMS:
                        if k.startswith(pname):
                            full_name = k
                            break
                dist_pc = PULSAR_PARAMS.get(full_name, {}).get('dist', 400.0)
                D_eff = dist_pc * s * (1 - s)
                meerkat_data[pname] = {
                    'eta': eta_mean,
                    'D_eff': D_eff,
                    'v': 50.0
                }
            else:
                meerkat_data[pname] = None

    full_data = []
    for pulsar_id, filename, freq, telescope in closure_pulsars:
        filepath = RESULTS_DIR / filename
        if not filepath.exists():
            continue

        with open(filepath) as f:
            cd = json.load(f)

        # Quality assessment: significant vs bounding constraint
        # Bounding pulsars are noise-limited but still provide valid |H|
        # measurements that constrain amplitude models
        rayleigh_p = cd.get('phase_closure_rayleigh_p', 1.0)
        n_epochs = cd.get('n_epochs', 0)
        phase_sigma = phase_detection_sigma(cd)
        quality = 'phase_significant' if abs(phase_sigma) >= 3 and n_epochs >= 5 else 'bounding'
        if quality == 'bounding':
            print_status(f"  Including {pulsar_id} as bounding constraint (Rayleigh p={rayleigh_p:.3f}, n_epochs={n_epochs})", "INFO")

        h_err = cd.get('H_sem_ns')
        # For single-epoch pulsars, approximate H_err from triplet count
        if h_err is None or not np.isfinite(h_err) or h_err <= 0:
            n_triplets = cd.get('n_total_triplets', 0)
            h_mag = cd.get('H_magnitude_ns', 0.0)
            if n_triplets > 0 and h_mag > 0:
                h_err = h_mag / np.sqrt(n_triplets)
                print_status(f"  {pulsar_id}: approximated H_err={h_err:.3f} from {n_triplets} triplets", "INFO")
            else:
                print_status(f"  Skipping {pulsar_id}: invalid H_err={h_err} and no triplets for approximation", "WARNING")
                continue

        h_mag = cd.get('H_magnitude_ns')
        if h_mag is None or not np.isfinite(h_mag) or h_mag <= 0:
            print_status(f"  Skipping {pulsar_id}: invalid |H|={h_mag} (need finite > 0)", "WARNING")
            continue

        # Get geometry for this pulsar
        base_name = pulsar_id.split('_')[0]
        geom = tep_data.get(base_name) or jiamusi_data.get(base_name) or meerkat_data.get(base_name)

        if geom:
            full_data.append({
                'pulsar': pulsar_id,
                'H': h_mag,
                'H_excess': cd.get('H_excess_ns', 0.0),
                'H_err': h_err,
                'eta': geom['eta'],
                'D_eff': geom['D_eff'],
                'v': geom['v'],
                'freq': freq,
                'telescope': telescope,
                'quality': quality,
                'phase_sigma': phase_sigma,
                'valid_for_primary_inference': False,
                'log_eta': np.log(geom['eta']),
                'log_D_eff': np.log(geom['D_eff']),
                'log_v': np.log(geom['v']),
                'log_freq': np.log(freq),
                'log_H': np.log(h_mag)
            })

    return full_data

def log_likelihood(params, data, model_type='general', verbose=False):
    """
    Compute log-likelihood for the model.
    
    params = [log_α, β, γ, δ, cal_Jiamusi]
    Where:
    - log_α: overall amplitude (log scale)
    - β: eta exponent  
    - γ: D_eff exponent
    - δ: FREQUENCY EXPONENT (key parameter!)
    - cal_Jiamusi: telescope calibration offset relative to Parkes
    """
    log_alpha, beta, gamma, delta, cal_jiamusi = params
    
    if verbose:
        print(f"\n    [LL CALCULATION] Parameters: log_α={log_alpha:.4f}, β={beta:.4f}, γ={gamma:.4f}, δ={delta:.4f}")
        print(f"    [LL CALCULATION] Calibration: Jiamusi={cal_jiamusi:.4f} (Parkes anchored to 0)")
    
    ll = 0
    if verbose:
        print(f"    [LL CALCULATION] Computing likelihood for {len(data)} pulsars:")
        print(f"    {'Pulsar':<10} {'log_H_obs':>10} {'log_H_pred':>10} {'sigma':>10} {'contrib':>12}")
        print("    " + "-" * 60)
    
    for i, d in enumerate(data):
        # Model prediction (in log space)
        log_H_pred = (log_alpha + 
                     beta * d['log_eta'] + 
                     gamma * d['log_D_eff'] + 
                     delta * d['log_freq'])
        
        # Add telescope calibration
        if d['telescope'] != 'Parkes':
            log_H_pred += cal_jiamusi
        
        # Gaussian likelihood
        sigma = d['H_err'] / d['H']  # Convert to log-scale uncertainty
        residual = d['log_H'] - log_H_pred
        chi2_term = residual**2 / sigma**2
        log_term = np.log(2*np.pi*sigma**2)
        contrib = -0.5 * (chi2_term + log_term)
        ll += contrib
        
        if verbose and i < 10:  # Show first 10 for brevity
            print(f"    {d['pulsar']:<10} {d['log_H']:>10.4f} {log_H_pred:>10.4f} {sigma:>10.4f} {residual:>10.4f} {contrib:>12.4f}")
        elif verbose and i == 10:
            print(f"    ... ({len(data)-10} more pulsars)")
    
    if verbose:
        print(f"    [LL CALCULATION] Total log-likelihood: {ll:.4f}")
    
    return ll

def fit_model(data, model_type='general', delta_fixed=None, verbose=False):
    """
    Fit model to data using maximum likelihood.
    
    If delta_fixed is provided, fix δ to that value (for model comparison).
    """
    
    if verbose:
        delta_str = f"δ={delta_fixed:.2f}" if delta_fixed is not None else "δ=free"
        print(f"\n  [FIT] Starting model fit: {delta_str}, {len(data)} data points")
    
    # Initial guess
    if delta_fixed is not None:
        # 4 params: log_α, β, γ, cal_J (δ fixed)
        # Full order: log_α, β, γ, δ, cal_J
        x0 = [0.0, 1.0, 1.0, 0.0]  # log_α, β, γ, cal_J
        bounds = [(-5, 5), (-2, 3), (-2, 3), (-1, 1)]
        
        def neg_ll_fixed(params):
            # params: [log_α, β, γ, cal_J]
            # full:   [log_α, β, γ, δ, cal_J]
            full_params = [params[0], params[1], params[2], delta_fixed, params[3]]
            return -log_likelihood(full_params, data, model_type, verbose=False)
        
        result = minimize(neg_ll_fixed, x0, bounds=bounds, method='L-BFGS-B')
        
        # Reconstruct full param vector
        opt_params = [result.x[0], result.x[1], result.x[2], delta_fixed, result.x[3]]
        
    else:
        # 5 params: log_α, β, γ, δ, cal_J
        x0 = [0.0, 1.0, 1.0, -1.0, 0.0]
        bounds = [(-5, 5), (-2, 3), (-2, 3), (-4, 2), (-1, 1)]
        
        def neg_ll_free(params):
            return -log_likelihood(params, data, model_type, verbose=False)
        
        result = minimize(neg_ll_free, x0, bounds=bounds, method='L-BFGS-B')
        opt_params = result.x
    
    max_ll = -result.fun
    
    return opt_params, max_ll, result.success

def compute_aic_bic(params, max_ll, data, n_params):
    """Compute AIC and BIC for model comparison."""
    n = len(data)
    if n == 0:
        return float('inf'), float('inf')
    aic = 2 * n_params - 2 * max_ll
    bic = n_params * np.log(n) - 2 * max_ll
    return aic, bic

def within_source_chromatic_test(all_data):
    """
    Within-source chromatic test using J0437 sub-bands.

    Same pulsar (J0437) at three frequencies eliminates:
      - Geometry differences (same η, D_eff, v)
      - Telescope calibration differences (same Parkes backend)
      - ISM path differences (same line of sight)

    Model: log H_i = log A + δ · log(ν_i/ν_ref) + noise
    Where A absorbs all common pulsar-dependent factors.
    """
    j0437_points = [d for d in all_data if d['pulsar'].startswith('J0437')]
    if len(j0437_points) < 3:
        return None

    j0437_points.sort(key=lambda d: d['freq'])

    # Reference frequency = full band
    nu_ref = 1380.0
    log_nus = [np.log(d['freq'] / nu_ref) for d in j0437_points]
    log_Hs = [d['log_H'] for d in j0437_points]
    sigmas = [d['H_err'] / d['H'] for d in j0437_points]

    def neg_ll(params):
        log_A, delta = params
        ll = 0.0
        for ln, lH, sig in zip(log_nus, log_Hs, sigmas):
            pred = log_A + delta * ln
            resid = lH - pred
            ll += -0.5 * (resid**2 / sig**2 + np.log(2 * np.pi * sig**2))
        return -ll

    from scipy.optimize import minimize

    # TEP model: δ = 0
    res_tep = minimize(lambda p: neg_ll([p[0], 0.0]), [np.mean(log_Hs)], method='L-BFGS-B')
    ll_tep = -res_tep.fun

    # Free δ model
    res_free = minimize(neg_ll, [np.mean(log_Hs), 0.0], bounds=[(None, None), (-4, 2)], method='L-BFGS-B')
    log_A, delta = res_free.x
    ll_free = -res_free.fun

    # ISM model: δ = -2
    res_ism = minimize(lambda p: neg_ll([p[0], -2.0]), [np.mean(log_Hs)], method='L-BFGS-B')
    ll_ism = -res_ism.fun

    aic_tep = 2 * 1 - 2 * ll_tep
    aic_free = 2 * 2 - 2 * ll_free
    aic_ism = 2 * 1 - 2 * ll_ism

    # Predictions
    preds = {d['pulsar']: float(np.exp(log_A + delta * np.log(d['freq'] / nu_ref))) for d in j0437_points}

    return {
        'n_points': len(j0437_points),
        'frequencies_mhz': [d['freq'] for d in j0437_points],
        'H_observed_ns': [d['H'] for d in j0437_points],
        'H_err_ns': [d['H_err'] for d in j0437_points],
        'H_predicted_ns': [preds[d['pulsar']] for d in j0437_points],
        'best_fit': {
            'log_A': float(log_A),
            'A': float(np.exp(log_A)),
            'delta': float(delta),
        },
        'model_comparison': {
            'M_TEP_delta0': {'log_likelihood': float(ll_tep), 'aic': float(aic_tep)},
            'M_ISM_delta-2': {'log_likelihood': float(ll_ism), 'aic': float(aic_ism)},
            'M_free': {'log_likelihood': float(ll_free), 'aic': float(aic_free)},
        },
        'delta_aic_tep_vs_free': float(aic_tep - aic_free),
        'delta_aic_ism_vs_free': float(aic_ism - aic_free),
        'conclusion': (
            'Achromatic TEP preferred' if aic_tep < aic_ism and aic_tep - aic_free < 2
            else 'Chromatic ISM preferred' if aic_ism < aic_tep and aic_ism - aic_free < 2
            else 'Inconclusive'
        ),
        'note': 'Within-source test eliminates geometry, telescope calibration, and ISM path differences. Limitation: sub-band errors are large (few epochs), so constraints on δ are weak.'
    }

def main():
    print("=" * 70)
    print("STEP 043: CHROMATIC DIAGNOSTIC (Hierarchical Model)")
    print("=" * 70)
    print("\nModel: log |H| = logα + β·logeta + γ·logD_eff + δ·lognu + cal + noise")
    print("Status: diagnostic only; unsigned |H| is noise-floor biased and not primary evidence")
    print("\nKey parameter δ (frequency exponent):")
    print("  δ = 0    -> Achromatic (TEP prediction)")
    print("  δ = -2   -> Chromatic ISM (nu^-2 scaling)")
    
    # Load data
    print("\n[1] Loading data...")
    data = load_complete_data()
    
    if len(data) == 0:
        print("  ERROR: No data loaded! Cannot perform chromatic test.")
        print("  Missing required data files:")
        print("    - step_031_tep_scaling_detailed.json")
        print("    - step_031_jiamusi_eta_analysis.json")
        print("    - step_003_closure_final_summary_*.json files")
        
        # Save empty but valid result
        output = {
            'status': 'failed',
            'error': 'No data available - missing prerequisite result files',
            'missing_files': [
                'step_031_tep_scaling_detailed.json',
                'step_031_jiamusi_eta_analysis.json'
            ],
            'best_fit_params': {'delta': None},
            'model_comparison': {},
            'data': []
        }
        with open(OUTPUT_FILE, 'w') as f:
            json.dump(output, f, indent=2, cls=NpEncoder)
        print(f"\n  Saved failure status to {OUTPUT_FILE}")
        return 1
    
    # Minimum data requirement relaxed: with quality-weighted bounding pulsars,
    # a stable 5-parameter fit can run with >= 3 data points (though constraints
    # are weak). Below 3 the model is underdetermined.
    if len(data) < 3:
        print_status(f"  [SKIP] Insufficient data: {len(data)} pulsars available, need >=3 for any fit.", "WARNING")
        output = {
            'status': 'skipped',
            'reason': f'Insufficient pulsars with valid data ({len(data)} < 3 required)',
            'n_pulsars_available': len(data),
            'best_fit_params': {'delta': None},
            'model_comparison': {},
            'data': data,
            'note': 'Chromatic test requires at least 3 pulsars with valid H, H_err, eta, D_eff, and frequency measurements.'
        }
        with open(OUTPUT_FILE, 'w') as f:
            json.dump(output, f, indent=2, cls=NpEncoder)
        print_status(f"  Saved skip status to {OUTPUT_FILE}", "INFO")
        return 1

    print_status(f"  Loaded {len(data)} pulsars:")
    print_status(f"\n    {'Pulsar':<10} {'Telescope':<10} {'Freq':<8} {'eta':<10} {'|H| (ns)':<12} {'log(H)':<10}")
    print_status("    " + "-" * 70)
    for d in data:
        log_h = np.log(d['H'])
        print_status(f"    {d['pulsar']:<10} {d['telescope']:<10} {d['freq']:<8} "
              f"{d['eta']:<10.4f} {d['H']:>12.2f} {log_h:>10.4f}")
    
    # Show detailed data transformations
    print_status("\n  [DATA PREPARATION] Computing log-transformed variables:")
    print_status(f"    {'Pulsar':<10} {'log(eta)':<10} {'log(D_eff)':<12} {'log(v)':<10} {'log(freq)':<10}")
    print_status("    " + "-" * 60)
    for d in data:
        print_status(f"    {d['pulsar']:<10} {d['log_eta']:<10.4f} {d['log_D_eff']:<12.4f} "
              f"{d['log_v']:<10.4f} {d['log_freq']:<10.4f}")
    
    # Fit three models
    print_status("\n[2] Fitting models with maximum likelihood...")
    
    # Model 1: TEP (δ = 0 fixed)
    print_status("\n  Model M_TEP: δ = 0 (achromatic, fixed)")
    print_status("    [FIT] Initial parameters: log_α=0.0, β=1.0, γ=1.0, δ=0.0, cal_J=0.0")
    print_status("    [FIT] Running L-BFGS-B optimization...")
    params_tep, ll_tep, success_tep = fit_model(data, delta_fixed=0.0)
    print_status(f"    [FIT] Optimization converged: {success_tep}")
    print_status(f"    [FIT] Optimal parameters: log_α={params_tep[0]:.4f}, β={params_tep[1]:.4f}, γ={params_tep[2]:.4f}")
    print_status(f"    [FIT] Fixed δ=0.0, cal_J={params_tep[4]:.4f}")
    aic_tep, bic_tep = compute_aic_bic(params_tep, ll_tep, data, 4)
    print_status(f"    [RESULT] Maximum log-likelihood: {ll_tep:.4f}")
    print_status(f"    [RESULT] AIC: {aic_tep:.4f}, BIC: {bic_tep:.4f}")
    
    # Model 2: General (δ free)
    print_status("\n  Model M_general: δ free (estimated from data)")
    print_status("    [FIT] Initial parameters: log_α=0.0, β=1.0, γ=1.0, δ=-1.0, cal_J=0.0")
    print_status("    [FIT] Running L-BFGS-B optimization with δ free...")
    params_gen, ll_gen, success_gen = fit_model(data, delta_fixed=None)
    print_status(f"    [FIT] Optimization converged: {success_gen}")
    print_status(f"    [FIT] Optimal parameters: log_α={params_gen[0]:.4f}, β={params_gen[1]:.4f}, γ={params_gen[2]:.4f}")
    print_status(f"    [FIT] Fitted δ={params_gen[3]:.4f} (KEY PARAMETER), cal_J={params_gen[4]:.4f}")
    aic_gen, bic_gen = compute_aic_bic(params_gen, ll_gen, data, 5)
    print_status(f"    [RESULT] Maximum log-likelihood: {ll_gen:.4f}")
    print_status(f"    [RESULT] AIC: {aic_gen:.4f}, BIC: {bic_gen:.4f}")
    
    # Model 3: ISM-constrained (δ = -2 fixed)
    print_status("\n  Model M_ISM: δ = -2 (chromatic, fixed)")
    print_status("    [FIT] Initial parameters: log_α=0.0, β=1.0, γ=1.0, δ=-2.0, cal_J=0.0")
    print_status("    [FIT] Running L-BFGS-B optimization...")
    params_ism, ll_ism, success_ism = fit_model(data, delta_fixed=-2.0)
    print_status(f"    [FIT] Optimization converged: {success_ism}")
    print_status(f"    [FIT] Optimal parameters: log_α={params_ism[0]:.4f}, β={params_ism[1]:.4f}, γ={params_ism[2]:.4f}")
    print_status(f"    [FIT] Fixed δ=-2.0, cal_J={params_ism[4]:.4f}")
    aic_ism, bic_ism = compute_aic_bic(params_ism, ll_ism, data, 4)
    print_status(f"    [RESULT] Maximum log-likelihood: {ll_ism:.4f}")
    print_status(f"    [RESULT] AIC: {aic_ism:.4f}, BIC: {bic_ism:.4f}")
    
    # Extract results for general model
    log_alpha, beta, gamma, delta, cal_j = params_gen
    
    print_status("\n[3] General model best-fit parameters:")
    print_status(f"  log α = {log_alpha:.3f}  (α = {np.exp(log_alpha):.3f})")
    print_status(f"  β (eta exponent)     = {beta:.3f}")
    print_status(f"  γ (D_eff exponent) = {gamma:.3f}")
    print_status(f"  δ (FREQUENCY)      = {delta:.3f}  ← KEY PARAMETER")
    print_status(f"  cal_Parkes         = 0.000  (1.000 multiplicative, reference telescope)")
    print_status(f"  cal_Jiamusi        = {cal_j:.3f}  ({np.exp(cal_j):.3f} multiplicative)")
    
    # Model comparison
    print_status("\n[4] Model comparison:")
    delta_aic_tep = aic_tep - aic_gen
    delta_aic_ism = aic_ism - aic_gen
    
    print_status(f"  DeltaAIC (M_TEP vs M_general):   {delta_aic_tep:+.1f}")
    print_status(f"  DeltaAIC (M_ISM vs M_general):   {delta_aic_ism:+.1f}")
    
    if delta_aic_tep < 2 and delta_aic_ism > 2:
        winner = "M_TEP"
        interpretation = "TEP achromatic model is preferred"
    elif delta_aic_ism < 2 and delta_aic_tep > 2:
        winner = "M_ISM"
        interpretation = "ISM chromatic model is preferred"
    elif abs(delta_aic_tep) < 2 and abs(delta_aic_ism) < 2:
        winner = "AMBIGUOUS"
        interpretation = "Data cannot distinguish between models"
    else:
        winner = "M_general"
        interpretation = "General model with free δ is preferred"
    
    print_status(f"\n  Preferred model: {winner}")
    print_status(f"  Interpretation: {interpretation}")
    
    # Evidence for/against specific δ values
    print_status("\n[5] Evidence for specific δ values:")
    print_status(f"  δ = 0 (TEP):     log-likelihood = {ll_tep:.2f}")
    print_status(f"  δ = -2 (ISM):    log-likelihood = {ll_ism:.2f}")
    print_status(f"  δ = {delta:.2f} (best-fit): log-likelihood = {ll_gen:.2f}")
    
    # Likelihood ratio test for δ = 0 vs δ = -2
    if ll_tep > ll_ism:
        lr_winner = "TEP (δ=0)"
        lr_diff = ll_tep - ll_ism
    else:
        lr_winner = "ISM (δ=-2)"
        lr_diff = ll_ism - ll_tep
    
    print_status(f"\n  Likelihood ratio: {lr_diff:.2f} in favor of {lr_winner}")
    
    # Interpret δ value
    print_status(f"\n[6] Interpretation of best-fit δ = {delta:.3f}:")
    if abs(delta) < 0.5:
        freq_conclusion = "Consistent with achromatic TEP (δ = 0)"
    elif delta < -1.5:
        freq_conclusion = "Consistent with chromatic ISM (δ ≈ -2)"
    elif delta < -0.5:
        freq_conclusion = "Intermediate: weak frequency dependence"
    else:
        freq_conclusion = "Unexpected positive frequency dependence"
    
    if np.isclose(delta, 2.0, atol=1e-3) or np.isclose(delta, -4.0, atol=1e-3):
        freq_conclusion = (
            f"{freq_conclusion}; fitted exponent is on an optimizer boundary, "
            "so the cross-pulsar frequency exponent is not interpretable as evidence."
        )
    print_status(f"  {freq_conclusion}")
    
    # Model predictions
    print_status("\n[7] Model predictions vs observations:")
    print_status(f"  {'Pulsar':<10} {'Freq':>6} {'|H _obs':>10} {'|H _pred':>10} {'Residual':>10}")
    print_status("  " + "-" * 60)
    
    for d in data:
        log_H_pred = (log_alpha + 
                     beta * d['log_eta'] + 
                     gamma * d['log_D_eff'] + 
                     delta * d['log_freq'])
        if d['telescope'] != 'Parkes':
            log_H_pred += cal_j

        H_pred = np.exp(log_H_pred)
        residual = d['H'] - H_pred
        
        print_status(f"  {d['pulsar']:<10} {d['freq']:>6} {d['H']:>10.2f} {H_pred:>10.2f} {residual:>10.2f}")
    
    # Within-source J0437 sub-band test
    print_status("\n[8] Within-source chromatic test (J0437 sub-bands)...")
    ws_result = within_source_chromatic_test(data)
    if ws_result:
        print_status(f"  J0437 sub-band test ({ws_result['n_points']} frequency points):")
        print_status(f"  Best-fit δ = {ws_result['best_fit']['delta']:.3f}")
        print_status(f"  ΔAIC (TEP vs free):   {ws_result['delta_aic_tep_vs_free']:+.1f}")
        print_status(f"  ΔAIC (ISM vs free):   {ws_result['delta_aic_ism_vs_free']:+.1f}")
        print_status(f"  Conclusion: {ws_result['conclusion']}")
        for d, obs, pred in zip(ws_result['frequencies_mhz'], ws_result['H_observed_ns'], ws_result['H_predicted_ns']):
            print_status(f"    {d} MHz: |H|_obs = {obs:.2f} ns, |H|_pred = {pred:.2f} ns")
    else:
        print_status("  Insufficient J0437 sub-band data for within-source test")

    # Final conclusion
    print_status("\n" + "=" * 70)
    print_status("CONCLUSION")
    print_status("=" * 70)
    ws_conclusion = ws_result['conclusion'] if ws_result else 'Not available (insufficient sub-band data)'
    ws_delta = f"{ws_result['best_fit']['delta']:.3f}" if ws_result else 'N/A'

    conclusion_text = f"""
HIERARCHICAL DIAGNOSTIC MODEL RESULTS:

Best-fit frequency exponent: δ = {delta:.3f}

Interpretation:
- δ = 0 corresponds to achromatic TEP (no frequency dependence)
- δ = -2 corresponds to chromatic ISM (nu^-2 scaling)
- Best-fit δ = {delta:.3f} indicates: {freq_conclusion}

Model comparison (AIC):
- M_TEP (δ=0):    AIC = {aic_tep:.1f}
- M_ISM (δ=-2):   AIC = {aic_ism:.1f}
- M_general:      AIC = {aic_gen:.1f}

Preferred model: {winner}
{interpretation}

Telescope calibration difference (Parkes/Jiamusi): {np.exp(-cal_j):.3f}
(Calibration uncertainty is a major limiting factor)

WITHIN-SOURCE J0437 SUB-BAND TEST:
{ws_conclusion}
Best-fit δ = {ws_delta}

OVERALL ASSESSMENT:
The cross-pulsar chromatic fit is diagnostic only because it uses unsigned |H|
amplitudes, includes noise-floor-limited bounding rows, and is strongly limited
by telescope calibration and geometry differences. The fitted free exponent
should not be used as primary evidence when it lands on an optimizer boundary.
The cleaner within-source J0437 sub-band test removes telescope and geometry
differences, but with only a few frequency points and large sub-band errors it
is best described as weakly achromatic/inconclusive rather than decisive.
"""

    print_status(conclusion_text)

    # Save results
    output = {
        'status': 'complete',
        'inference_status': 'diagnostic_only',
        'valid_for_primary_inference': False,
        'method': 'Hierarchical diagnostic model with frequency exponent',
        'best_fit_params': {
            'log_alpha': float(log_alpha),
            'alpha': float(np.exp(log_alpha)),
            'beta': float(beta),
            'gamma': float(gamma),
            'delta': float(delta),
            'cal_parkes': 0.0,
            'cal_jiamusi': float(cal_j)
        },
        'model_comparison': {
            'M_TEP': {'aic': float(aic_tep), 'bic': float(bic_tep), 'log_likelihood': float(ll_tep)},
            'M_ISM': {'aic': float(aic_ism), 'bic': float(bic_ism), 'log_likelihood': float(ll_ism)},
            'M_general': {'aic': float(aic_gen), 'bic': float(bic_gen), 'log_likelihood': float(ll_gen)}
        },
        'preferred_model': winner,
        'interpretation': freq_conclusion,
        'limitations': [
            'Cross-pulsar fit uses unsigned |H| amplitudes, which are folded-noise-floor biased.',
            'Bounding rows constrain amplitudes but are not positive phase detections.',
            'Telescope calibration and geometry differences dominate the cross-pulsar model.',
            'A free exponent on an optimizer boundary is not interpretable as strong chromatic evidence.',
            'Within-source J0437 sub-band test has few points and large errors.'
        ],
        'within_source_test': ws_result,
        'data': data
    }

    with open(OUTPUT_FILE, 'w') as f:
        json.dump(output, f, indent=2, cls=NpEncoder)

    print(f"\nResults saved: {OUTPUT_FILE}")
    print("=" * 70)

if __name__ == "__main__":
    main()
