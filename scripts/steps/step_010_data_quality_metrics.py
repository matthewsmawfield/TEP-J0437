#!/usr/bin/env python3
"""
STEP 010: COMPREHENSIVE DATA QUALITY METRICS

This step provides comprehensive data quality metrics and epoch-by-epoch
reporting to ensure transparency and enable detailed inspection of the dataset.

Key metrics:
1. Overall dataset statistics
2. Epoch-by-epoch quality metrics
3. SNR distribution analysis
4. Temporal quality trends
5. Outlier detection and reporting
6. Data completeness assessment
"""

import sys
import json
from pathlib import Path
import numpy as np
from scipy import stats
from typing import Dict, Any, Optional
from datetime import datetime

# Project configuration
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.json_numpy import NpEncoder
from scripts.utils.config import RANDOM_SEED
from scripts.utils.logger import print_status
RESULTS_DIR = PROJECT_ROOT / "results"

np.random.seed(RANDOM_SEED)


def load_closure_results():
    """Load closure delay results."""
    results_file = RESULTS_DIR / "step_003_closure_final_per_epoch.json"
    
    if not results_file.exists():
        print("ERROR: step_003_closure_final_per_epoch.json not found. Run step_003 first.")
        return None
    
    with open(results_file, 'r') as f:
        data = json.load(f)
    
    results = {}
    for epoch_data in data:
        epoch_name = epoch_data.get("epoch", "unknown")
        
        # Track data quality - warn if using defaults
        triplets = epoch_data.get("triplets", [])
        n_arclets = epoch_data.get("n_arclets")
        mjd = epoch_data.get("mjd")
        
        if n_arclets is None:
            print(f"[WARN] Epoch {epoch_name}: missing n_arclets field")
            n_arclets = 0
        if mjd is None:
            print(f"[WARN] Epoch {epoch_name}: missing mjd field")
            mjd = np.nan  # Use NaN instead of 0 (which is year 1858!)
        
        results[epoch_name] = {
            "triplet_delays": [t.get("geometric_delta_us", t.get("delta_us")) for t in triplets],
            "triplet_snrs": [t.get("snr") for t in triplets],  # Use None if missing, not 0
            "n_triplets": len(triplets),
            "n_arclets": n_arclets,
            "mjd": mjd
        }
    
    return results


def generate_overall_statistics(results):
    """Generate overall dataset statistics."""
    print_status("" + "=" * 70)
    print("OVERALL DATASET STATISTICS")
    print_status("===" * 70)
    
    # Collect all data
    all_delays = []
    all_snrs = []
    epoch_stats = []
    
    for epoch_name, epoch_data in results.items():
        # Filter out None values from missing data
        delays_raw = epoch_data.get('triplet_delays', [])
        snrs_raw = epoch_data.get('triplet_snrs', [])
        
        delays = np.array([d for d in delays_raw if d is not None])
        snrs = np.array([s for s in snrs_raw if s is not None])
        
        if len(delays) > 0:
            all_delays.extend(delays)
            all_snrs.extend(snrs)
            
            # Handle std calculation for single-triplet case (ddof=1 returns NaN)
            if len(delays) > 1:
                std_delay = float(np.std(delays, ddof=1))
            else:
                std_delay = 0.0  # Single measurement has no variance
            
            mjd = epoch_data.get('mjd')
            # Handle NaN MJDs gracefully
            if np.isnan(mjd) if isinstance(mjd, float) else mjd is None:
                mjd = None
            
            epoch_stats.append({
                'epoch': epoch_name,
                'n_triplets': len(delays),
                'n_arclets': epoch_data.get('n_arclets', 0),
                'mean_delay': float(np.mean(delays)),
                'std_delay': std_delay,
                'mean_snr': float(np.mean(snrs)) if len(snrs) > 0 else None,  # Use None if no SNR data
                'mjd': mjd
            })
    
    all_delays = np.array(all_delays)
    all_snrs = np.array(all_snrs)
    
    print(f"\nDataset overview:")
    print_status(f"Total epochs: {len(results)}")
    print_status(f"Total triplets: {len(all_delays)}")
    print_status(f"Mean triplets per epoch: {np.mean([e['n_triplets'] for e in epoch_stats]):.1f}")
    print_status(f"Median triplets per epoch: {np.median([e['n_triplets'] for e in epoch_stats]):.1f}")
    
    print(f"\nClosure delay statistics:")
    print_status(f"Mean: {np.mean(all_delays)*1e3:.3f} ns")
    print_status(f"Std: {np.std(all_delays, ddof=1)*1e3:.3f} ns")
    print_status(f"Median: {np.median(all_delays)*1e3:.3f} ns")
    print_status(f"Min: {np.min(all_delays)*1e3:.3f} ns")
    print_status(f"Max: {np.max(all_delays)*1e3:.3f} ns")
    print_status(f"Skewness: {stats.skew(all_delays):.3f}")
    print_status(f"Kurtosis: {stats.kurtosis(all_delays):.3f}")
    
    print(f"\nSNR statistics:")
    print_status(f"Mean SNR: {np.mean(all_snrs):.2f}")
    print_status(f"Median SNR: {np.median(all_snrs):.2f}")
    print_status(f"Min SNR: {np.min(all_snrs):.2f}")
    print_status(f"Max SNR: {np.max(all_snrs):.2f}")
    
    return {
        'n_epochs': len(results),
        'n_triplets_total': len(all_delays),
        'mean_triplets_per_epoch': float(np.mean([e['n_triplets'] for e in epoch_stats])),
        'median_triplets_per_epoch': float(np.median([e['n_triplets'] for e in epoch_stats])),
        'delay_mean_ns': float(np.mean(all_delays) * 1e3),
        'delay_std_ns': float(np.std(all_delays, ddof=1) * 1e3),
        'delay_median_ns': float(np.median(all_delays) * 1e3),
        'delay_min_ns': float(np.min(all_delays) * 1e3),
        'delay_max_ns': float(np.max(all_delays) * 1e3),
        'delay_skewness': float(stats.skew(all_delays)),
        'delay_kurtosis': float(stats.kurtosis(all_delays)),
        'snr_mean': float(np.mean(all_snrs)),
        'snr_median': float(np.median(all_snrs)),
        'snr_min': float(np.min(all_snrs)),
        'snr_max': float(np.max(all_snrs)),
        'epoch_stats': epoch_stats
    }


def generate_epoch_by_epoch_report(results):
    """Generate detailed epoch-by-epoch quality report."""
    print_status("" + "=" * 70)
    print("EPOCH-BY-EPOCH QUALITY REPORT")
    print_status("===" * 70)
    
    epoch_report = []
    
    for epoch_name, epoch_data in sorted(results.items(), key=lambda x: x[1].get('mjd', 0)):
        delays = np.array(epoch_data.get('triplet_delays', []))
        snrs = np.array(epoch_data.get('triplet_snrs', []))
        
        if len(delays) > 0:
            # Calculate epoch-specific metrics
            mean_delay = np.mean(delays)
            # Handle std calculation for single-triplet case (ddof=1 returns NaN)
            if len(delays) > 1:
                std_delay = np.std(delays, ddof=1)
            else:
                std_delay = 0.0  # Single measurement has no variance
            abs_mean_delay = abs(np.mean(delays))
            
            # Sign distribution
            n_neg = np.sum(delays < 0)
            n_pos = np.sum(delays > 0)
            neg_frac = n_neg / len(delays)
            
            # SNR metrics
            mean_snr = np.mean(snrs) if len(snrs) > 0 else 0.0
            min_snr = np.min(snrs) if len(snrs) > 0 else 0.0
            
            # Quality flags
            low_snr_flag = mean_snr < 5.0 if len(snrs) > 0 else False
            outlier_flag = abs(mean_delay) > 3 * std_delay if std_delay > 0 else False
            low_triplet_flag = len(delays) < 5
            
            epoch_report.append({
                'epoch': epoch_name,
                'mjd': epoch_data.get('mjd', 0),
                'n_triplets': len(delays),
                'n_arclets': epoch_data.get('n_arclets', 0),
                'mean_delay_ns': float(mean_delay * 1e3),
                'std_delay_ns': float(std_delay * 1e3),
                'abs_mean_delay_ns': float(abs_mean_delay * 1e3),
                'n_negative': int(n_neg),
                'n_positive': int(n_pos),
                'negative_fraction': float(neg_frac),
                'mean_snr': float(mean_snr),
                'min_snr': float(min_snr),
                'quality_flags': {
                    'low_snr': bool(low_snr_flag),
                    'outlier': bool(outlier_flag),
                    'low_triplet_count': bool(low_triplet_flag)
                }
            })
    
    # Print summary table
    print(f"\n{'Epoch':<20} {'MJD':>10} {'N_trip':>7} {'N_arc':>6} {'Mean(ns)':>10} {'Std(ns)':>9} {'Neg%':>6} {'SNR':>6} {'Flags':>10}")
    print("-" * 90)
    
    for ep in epoch_report[:20]:  # Print first 20
        flags = []
        if ep['quality_flags']['low_snr']:
            flags.append('L')
        if ep['quality_flags']['outlier']:
            flags.append('O')
        if ep['quality_flags']['low_triplet_count']:
            flags.append('T')
        flag_str = ''.join(flags) if flags else '-'
        
        print(f"{ep['epoch']:<20} {ep['mjd']:>10.1f} {ep['n_triplets']:>7} {ep['n_arclets']:>6} "
              f"{ep['mean_delay_ns']:>10.3f} {ep['std_delay_ns']:>9.3f} {ep['negative_fraction']*100:>5.1f}% "
              f"{ep['mean_snr']:>6.1f} {flag_str:>10}")
    
    if len(epoch_report) > 20:
        print(f"... ({len(epoch_report) - 20} more epochs)")
    
    # Quality summary
    low_snr_count = sum(1 for ep in epoch_report if ep['quality_flags']['low_snr'])
    outlier_count = sum(1 for ep in epoch_report if ep['quality_flags']['outlier'])
    low_triplet_count = sum(1 for ep in epoch_report if ep['quality_flags']['low_triplet_count'])
    
    print(f"\nQuality flags summary:")
    print_status(f"Low SNR epochs: {low_snr_count}/{len(epoch_report)} ({100*low_snr_count/len(epoch_report):.1f}%)")
    print_status(f"Outlier epochs: {outlier_count}/{len(epoch_report)} ({100*outlier_count/len(epoch_report):.1f}%)")
    print_status(f"Low triplet count epochs: {low_triplet_count}/{len(epoch_report)} ({100*low_triplet_count/len(epoch_report):.1f}%)")
    
    return {
        'epoch_report': epoch_report,
        'quality_summary': {
            'low_snr_count': low_snr_count,
            'outlier_count': outlier_count,
            'low_triplet_count': low_triplet_count,
            'total_epochs': len(epoch_report)
        }
    }


def analyze_snr_distribution(results):
    """Analyze SNR distribution across the dataset."""
    print_status("" + "=" * 70)
    print("SNR DISTRIBUTION ANALYSIS")
    print_status("===" * 70)
    
    all_snrs = []
    for epoch_data in results.values():
        all_snrs.extend(epoch_data.get('triplet_snrs', []))
    
    all_snrs = np.array(all_snrs)
    
    # SNR bins
    bins = [0, 3, 5, 10, 20, 50, 100, float('inf')]
    bin_labels = ['0-3', '3-5', '5-10', '10-20', '20-50', '50-100', '>100']
    
    hist, _ = np.histogram(all_snrs, bins=bins)
    
    print(f"\nSNR distribution:")
    print("-" * 70)
    for i, (label, count) in enumerate(zip(bin_labels, hist)):
        frac = 100 * count / len(all_snrs)
        print(f"  {label:>6}: {count:>6} ({frac:>5.1f}%)")
    
    print(f"\nSNR statistics:")
    print_status(f"Mean: {np.mean(all_snrs):.2f}")
    print_status(f"Median: {np.median(all_snrs):.2f}")
    print_status(f"Std: {np.std(all_snrs):.2f}")
    
    # SNR vs delay correlation
    all_delays = []
    for epoch_data in results.values():
        all_delays.extend(epoch_data.get('triplet_delays', []))
    
    all_delays = np.array(all_delays)
    
    if len(all_snrs) == len(all_delays):
        r, p = stats.pearsonr(all_snrs, np.abs(all_delays))
        print(f"\nSNR |delay| correlation:")
        print(f"  r = {r:.3f}, p = {p:.2e}")
        
        if abs(r) < 0.1:
            print(f"  No significant correlation (good - SNR independent of effect)")
        else:
            print(f"  Correlation detected - may indicate SNR bias")
    
    return {
        'snr_distribution': {
            'bins': bin_labels,
            'counts': [int(c) for c in hist]
        },
        'snr_mean': float(np.mean(all_snrs)),
        'snr_median': float(np.median(all_snrs)),
        'snr_std': float(np.std(all_snrs))
    }


def analyze_temporal_trends(results):
    """Analyze temporal quality trends."""
    print_status("" + "=" * 70)
    print("TEMPORAL QUALITY TRENDS")
    print_status("===" * 70)
    
    # Collect temporal data
    temporal_data = []
    for epoch_name, epoch_data in sorted(results.items(), key=lambda x: x[1].get('mjd', 0)):
        mjd = epoch_data.get('mjd', 0)
        if mjd > 0:
            delays = np.array(epoch_data.get('triplet_delays', []))
            if len(delays) > 0:
                # Handle std calculation for single-triplet case
                if len(delays) > 1:
                    std_delay = np.std(delays, ddof=1)
                else:
                    std_delay = 0.0
                temporal_data.append({
                    'mjd': mjd,
                    'n_triplets': len(delays),
                    'mean_abs_delay': abs(np.mean(delays)),
                    'std_delay': std_delay
                })
    
    if len(temporal_data) < 10:
        print("Insufficient temporal data for trend analysis")
        return None
    
    temporal_data = np.array([(d['mjd'], d['n_triplets'], d['mean_abs_delay'], d['std_delay']) 
                             for d in temporal_data])
    
    mjds = temporal_data[:, 0]
    n_triplets = temporal_data[:, 1]
    mean_abs_delays = temporal_data[:, 2]
    std_delays = temporal_data[:, 3]
    
    print(f"\nTemporal range:")
    print_status(f"Start MJD: {np.min(mjds):.1f}")
    print_status(f"End MJD: {np.max(mjds):.1f}")
    print_status(f"Span: {np.max(mjds) - np.min(mjds):.1f} days ({(np.max(mjds) - np.min(mjds))/365.25:.1f} years)")
    
    # Trend analysis
    r_n, p_n = stats.pearsonr(mjds, n_triplets)
    r_h, p_h = stats.pearsonr(mjds, mean_abs_delays)
    # Handle case where std_delays has no variance (all zeros from single-triplet epochs)
    if np.std(std_delays) > 0:
        r_s, p_s = stats.pearsonr(mjds, std_delays)
    else:
        r_s, p_s = 0.0, 1.0  # No correlation when no variance
    
    print(f"\nTemporal correlations:")
    print_status(f"MJD vs n_triplets: r = {r_n:.3f}, p = {p_n:.2e}")
    print_status(f"MJD |H|: r = {r_h:.3f}, p = {p_h:.2e}")
    print_status(f"MJD vs std: r = {r_s:.3f}, p = {p_s:.2e}")
    
    print(f"\nInterpretation:")
    if p_n > 0.05 and p_h > 0.05:
        print(f"  No significant temporal trends detected")
        print(f"  Data quality stable over time (good)")
    else:
        print(f"  Significant temporal trends detected")
        if p_n < 0.05:
            print(f"  - Triplet count varies with time")
        if p_h < 0.05:
            print(f"  - |H| varies with time")
    
    return {
        'temporal_span_days': float(np.max(mjds) - np.min(mjds)),
        'correlation_n_triplets': float(r_n),
        'correlation_h_magnitude': float(r_h),
        'correlation_std': float(r_s),
        'p_n_triplets': float(p_n),
        'p_h_magnitude': float(p_h),
        'p_std': float(p_s)
    }


def detect_outliers(results):
    """Detect and report outliers in the dataset."""
    print_status("" + "=" * 70)
    print("OUTLIER DETECTION")
    print_status("===" * 70)
    
    all_delays = []
    epoch_delays = {}
    
    for epoch_name, epoch_data in results.items():
        delays = np.array(epoch_data['triplet_delays'])
        if len(delays) > 0:
            all_delays.extend(delays)
            epoch_delays[epoch_name] = delays
    
    all_delays = np.array(all_delays)
    
    # Detect outliers using IQR method
    q25, q75 = np.percentile(all_delays, [25, 75])
    iqr = q75 - q25
    lower_bound = q25 - 3 * iqr
    upper_bound = q75 + 3 * iqr
    
    outliers = all_delays[(all_delays < lower_bound) | (all_delays > upper_bound)]
    outlier_fraction = len(outliers) / len(all_delays)
    
    print(f"\nOutlier statistics (IQR method, 3xIQR):")
    print_status(f"Lower bound: {lower_bound*1e3:.3f} ns")
    print_status(f"Upper bound: {upper_bound*1e3:.3f} ns")
    print_status(f"Outliers detected: {len(outliers)}/{len(all_delays)} ({100*outlier_fraction:.2f}%)")
    
    # Detect epoch-level outliers
    epoch_means = [np.mean(delays) for delays in epoch_delays.values()]
    epoch_stds = [np.std(delays, ddof=1) if len(delays) > 1 else 0.0 for delays in epoch_delays.values()]
    
    q25_mean, q75_mean = np.percentile(epoch_means, [25, 75])
    iqr_mean = q75_mean - q25_mean
    lower_mean = q25_mean - 1.5 * iqr_mean
    upper_mean = q75_mean + 1.5 * iqr_mean
    
    outlier_epochs = []
    for epoch_name, delays in epoch_delays.items():
        mean_delay = np.mean(delays)
        if mean_delay < lower_mean or mean_delay > upper_mean:
            outlier_epochs.append({
                'epoch': epoch_name,
                'mean_delay_ns': float(mean_delay * 1e3),
                'n_triplets': len(delays)
            })
    
    print(f"\nEpoch-level outliers (1.5xIQR):")
    print_status(f"Outlier epochs: {len(outlier_epochs)}/{len(epoch_delays)}")
    
    if len(outlier_epochs) > 0:
        print(f"\nOutlier epoch details:")
        for ep in outlier_epochs[:10]:  # Print first 10
            print(f"  {ep['epoch']}: mean = {ep['mean_delay_ns']:.3f} ns, n_triplets = {ep['n_triplets']}")
        if len(outlier_epochs) > 10:
            print(f"  ... ({len(outlier_epochs) - 10} more)")
    
    return {
        'outlier_count': len(outliers),
        'outlier_fraction': float(outlier_fraction),
        'outlier_epochs': outlier_epochs
    }


def assess_data_completeness(results):
    """Assess data completeness and missing data."""
    print_status("" + "=" * 70)
    print("DATA COMPLETENESS ASSESSMENT")
    print_status("===" * 70)
    
    # Check for missing MJD values
    epochs_with_mjd = sum(1 for e in results.values() if e.get('mjd', 0) > 0)
    epochs_without_mjd = len(results) - epochs_with_mjd
    
    print(f"\nMJD completeness:")
    print_status(f"Epochs with MJD: {epochs_with_mjd}/{len(results)} ({100*epochs_with_mjd/len(results):.1f}%)")
    print_status(f"Epochs without MJD: {epochs_without_mjd}/{len(results)} ({100*epochs_without_mjd/len(results):.1f}%)")
    
    # Check for empty epochs
    empty_epochs = sum(1 for e in results.values() if len(e['triplet_delays']) == 0)
    non_empty_epochs = len(results) - empty_epochs
    
    print(f"\nEpoch data completeness:")
    print_status(f"Non-empty epochs: {non_empty_epochs}/{len(results)} ({100*non_empty_epochs/len(results):.1f}%)")
    print_status(f"Empty epochs: {empty_epochs}/{len(results)} ({100*empty_epochs/len(results):.1f}%)")
    
    # Check SNR availability
    epochs_with_snr = sum(1 for e in results.values() if len(e['triplet_snrs']) > 0 and any(s > 0 for s in e['triplet_snrs']))
    epochs_without_snr = len(results) - epochs_with_snr
    
    print(f"\nSNR data completeness:")
    print_status(f"Epochs with SNR: {epochs_with_snr}/{len(results)} ({100*epochs_with_snr/len(results):.1f}%)")
    print_status(f"Epochs without SNR: {epochs_without_snr}/{len(results)} ({100*epochs_without_snr/len(results):.1f}%)")
    
    # Overall completeness score
    completeness_score = (epochs_with_mjd + non_empty_epochs + epochs_with_snr) / (3 * len(results))
    
    print(f"\nOverall completeness score: {completeness_score:.1%}")
    
    if completeness_score > 0.95:
        print(f"  Excellent data completeness")
    elif completeness_score > 0.90:
        print(f"  Good data completeness")
    elif completeness_score > 0.80:
        print(f"  Acceptable data completeness")
    else:
        print(f"  Poor data completeness - review needed")
    
    return {
        'mjd_completeness': float(epochs_with_mjd / len(results)),
        'data_completeness': float(non_empty_epochs / len(results)),
        'snr_completeness': float(epochs_with_snr / len(results)),
        'overall_completeness': float(completeness_score)
    }


def step_main(logger=None, verbose=True):
    """Standard pipeline entry point for data quality metrics."""
    return main()


def main():
    """Run comprehensive data quality metrics analysis."""
    print_status("===" * 70)
    print("STEP 010: COMPREHENSIVE DATA QUALITY METRICS")
    print_status("===" * 70)
    print()
    print("This step provides comprehensive data quality metrics and")
    print("epoch-by-epoch reporting for transparency and inspection.")
    print()
    
    results = load_closure_results()
    if results is None:
        return
    
    all_results = {}
    
    # Test 1: Overall statistics
    try:
        all_results['overall_statistics'] = generate_overall_statistics(results)
    except Exception as e:
        print(f"[FAIL] Overall statistics failed: {e}")
        all_results['overall_statistics'] = None
    
    # Test 2: Epoch-by-epoch report
    try:
        all_results['epoch_report'] = generate_epoch_by_epoch_report(results)
    except Exception as e:
        print(f"[FAIL] Epoch report failed: {e}")
        all_results['epoch_report'] = None
    
    # Test 3: SNR distribution
    try:
        all_results['snr_distribution'] = analyze_snr_distribution(results)
    except Exception as e:
        print(f"[FAIL] SNR distribution failed: {e}")
        all_results['snr_distribution'] = None
    
    # Test 4: Temporal trends
    try:
        all_results['temporal_trends'] = analyze_temporal_trends(results)
    except Exception as e:
        print(f"[FAIL] Temporal trends failed: {e}")
        all_results['temporal_trends'] = None
    
    # Test 5: Outlier detection
    try:
        all_results['outlier_detection'] = detect_outliers(results)
    except Exception as e:
        print(f"[FAIL] Outlier detection failed: {e}")
        all_results['outlier_detection'] = None
    
    # Test 6: Data completeness
    try:
        all_results['data_completeness'] = assess_data_completeness(results)
    except Exception as e:
        print(f"[FAIL] Data completeness failed: {e}")
        all_results['data_completeness'] = None
    
    # Save results
    output_file = RESULTS_DIR / "step_010_data_quality_metrics.json"
    with open(output_file, 'w') as f:
        json.dump(all_results, f, indent=2, cls=NpEncoder)
    print_status(f"Results saved to: {output_file}")
    
    # Save epoch report separately for easier inspection
    if all_results.get('epoch_report'):
        epoch_report_file = RESULTS_DIR / "step_010_epoch_by_epoch_report.json"
        with open(epoch_report_file, 'w') as f:
            json.dump(all_results['epoch_report'], f, indent=2, cls=NpEncoder)
        print(f"Epoch report saved to: {epoch_report_file}")
    
    # Summary
    print_status("" + "=" * 70)
    print("DATA QUALITY SUMMARY")
    print_status("===" * 70)
    
    if all_results.get('overall_statistics'):
        print(f"\nDataset: {all_results['overall_statistics']['n_epochs']} epochs, "
              f"{all_results['overall_statistics']['n_triplets_total']} triplets")
    
    if all_results.get('data_completeness'):
        completeness = all_results['data_completeness']['overall_completeness']
        print(f"Completeness: {completeness:.1%}")
    
    if all_results.get('outlier_detection'):
        outlier_frac = all_results['outlier_detection']['outlier_fraction']
        print(f"Outlier fraction: {outlier_frac:.2%}")
    
    print_status("Quality metrics generated:")
    print("  1. Overall dataset statistics")
    print("  2. Epoch-by-epoch quality report")
    print("  3. SNR distribution analysis")
    print("  4. Temporal quality trends")
    print("  5. Outlier detection and reporting")
    print("  6. Data completeness assessment")
    
    print_status("" + "=" * 70)


if __name__ == "__main__":
    main()
