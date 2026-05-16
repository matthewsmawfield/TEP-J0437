#!/usr/bin/env python3
"""
Step 032: Jiamusi Pulsar Data Analysis

Downloads and analyzes Jiamusi 66m telescope dynamic spectra data.
This is a MAJOR data source with MULTIPLE EPOCHS for 10 pulsars.

Source: http://zmtt.bao.ac.cn/psr-jms/
Paper: Wang et al. (2018), A&A 618, A186
"""

import sys
import json
import numpy as np
from pathlib import Path
from collections import defaultdict
from typing import List

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.json_numpy import NpEncoder

from scripts.utils.logger import TEPLogger, set_step_logger, print_status

# Directories
JIAMUSI_RAW_DIR = PROJECT_ROOT / "data" / "raw" / "jiamusi"
JIAMUSI_PROC_DIR = PROJECT_ROOT / "data" / "processed" / "jiamusi"
JIAMUSI_RESULTS_DIR = PROJECT_ROOT / "results"

for d in [JIAMUSI_RAW_DIR, JIAMUSI_PROC_DIR, JIAMUSI_RESULTS_DIR]:
    d.mkdir(parents=True, exist_ok=True)


# Jiamusi data URLs - ALL 10 pulsars from Wang et al. (2018)
# Total: 33 epochs across 10 pulsars (2015-2017)
JIAMUSI_PULSARS = {
    'B0329+54': [
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B0329_201602211028d.dat',
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B0329_201602241326d.dat',
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B0329_201711081846d.dat',
    ],
    'B0355+54': [
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B0355_201508202106d.dat',
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B0355_201601290850d.dat',
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B0355_201601292120d.dat',
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B0355_201711060010d.dat',
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B0355_201711091756d.dat',
    ],
    'B0540+23': [
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B0540_201506251100d.dat',
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B0540_201608090445d.dat',
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B0540_201711012106d.dat',
    ],
    'B0740-28': [
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B0740_201512120030d.dat',
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B0740_201601262105d.dat',
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B0740_201711020255d.dat',
    ],
    'B1508+55': [
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B1508_201512140350d.dat',
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B1508_201710310630d.dat',
    ],
    'B1933+16': [
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B1933_201506160205d.dat',
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B1933_201602190440d.dat',
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B1933_201602260440d.dat',
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B1933_201605210150d.dat',
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B1933_201605212310d.dat',
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B1933_201605222340d.dat',
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B1933_201605240125d.dat',
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B1933_201605260210d.dat',
    ],
    'B2154+40': [
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B2154_201601251455d.dat',
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B2154_201710312035d.dat',
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B2154_201711101140d.dat',
    ],
    'B2310+42': [
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B2310_201507202135d.dat',
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B2310_201710252101d.dat',
    ],
    'B2324+60': [
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B2324_201703291426d.dat',
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B2324_201711041040d.dat',
    ],
    'B2351+61': [
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B2351_201506190015d.dat',
        'http://zmtt.bao.ac.cn/psr-jms/paper2/B2351_201608110515d.dat',
    ],
}


def parse_jiamusi_file(filepath: Path) -> dict:
    """Parse Jiamusi format dynamic spectrum.
    
    Format: 3 columns - frequency(MHz), MJD, intensity
    Data is in long format (one row per time-freq pixel).
    """
    print_status(f"  Parsing {filepath.name}...", "INFO")
    
    # Read data
    data = np.loadtxt(filepath)
    
    # Extract columns
    freq = data[:, 0]  # MHz
    mjd = data[:, 1]
    intensity = data[:, 2]
    
    # Get unique values
    unique_freqs = np.unique(freq)
    unique_mjds = np.unique(mjd)
    
    n_freq = len(unique_freqs)
    n_time = len(unique_mjds)
    
    print_status(f"    Reshaping to ({n_time}, {n_freq})...", "INFO")
    
    # Validate array size before reshaping
    expected_size = n_time * n_freq
    actual_size = len(intensity)
    if actual_size != expected_size:
        print_status(f"    Array size mismatch: expected {expected_size}, got {actual_size}. Skipping file.", "WARNING")
        return None
    
    # Reshape to 2D array (time, freq)
    dynspec = intensity.reshape((n_time, n_freq))
    
    # Calculate time resolution
    if n_time > 1:
        dt_s = (unique_mjds[1] - unique_mjds[0]) * 86400  # Convert MJD diff to seconds
    else:
        raise ValueError(
            f"Cannot determine time sampling: only {len(unique_mjds)} unique MJD in data. "
            "Jiamusi data should have multiple time samples per observation."
        )
    
    return {
        'dynspec': dynspec,
        'freq_MHz': unique_freqs,
        'mjd_start': float(unique_mjds[0]),
        'dt_s': float(dt_s),
        'n_time': n_time,
        'n_freq': n_freq,
    }


def download_jiamusi_data():
    """Download Jiamusi dynamic spectra data."""
    import urllib.request

    print_status("Downloading Jiamusi pulsar data...", "INFO")

    expected_files: List[str] = []
    failures: List[str] = []
    for pulsar, urls in JIAMUSI_PULSARS.items():
        print_status(f"  {pulsar}: {len(urls)} epochs", "INFO")

        for url in urls:
            filename = Path(url).name
            expected_files.append(filename)
            filepath = JIAMUSI_RAW_DIR / filename

            if filepath.exists():
                continue

            print_status(f"    Downloading {filename}...", "INFO")
            try:
                urllib.request.urlretrieve(url, filepath)
            except Exception as exc:
                failures.append(f"{filename}: {exc}")
                print_status(f"    Failed to download {filename}: {exc}", "ERROR")

    present = sorted(path.name for path in JIAMUSI_RAW_DIR.glob("*.dat"))
    missing = sorted(set(expected_files) - set(present))
    if failures or missing:
        raise RuntimeError(
            "Jiamusi archive download incomplete. "
            f"missing={missing[:10]}{'...' if len(missing) > 10 else ''}; "
            f"failures={failures[:3]}{'...' if len(failures) > 3 else ''}"
        )

    print_status(f"Downloaded/verified {len(present)} Jiamusi epochs", "SUCCESS")
    return len(present)


def process_jiamusi_data():
    """Process all Jiamusi dynamic spectra."""
    print_status("Processing Jiamusi dynamic spectra...", "INFO")
    
    catalog = {'n_epochs': 0, 'epochs': []}
    
    for dat_file in sorted(JIAMUSI_RAW_DIR.glob('*.dat')):
        try:
            # Parse file
            data = parse_jiamusi_file(dat_file)
            
            # Skip if file could not be parsed (e.g., array size mismatch)
            if data is None:
                print_status(f"  Skipping {dat_file.name} (parse returned None)", "WARNING")
                continue
            
            # Save as NPZ
            out_name = dat_file.stem + '.npz'
            out_path = JIAMUSI_PROC_DIR / out_name
            
            np.savez(out_path,
                     dynspec=data['dynspec'],
                     freq_MHz=data['freq_MHz'],
                     mjd_start=data['mjd_start'],
                     dt_s=data['dt_s'],
                     n_time=data['n_time'],
                     n_freq=data['n_freq'])
            
            # Add to catalog
            catalog['epochs'].append({
                'file': out_name,
                'source': dat_file.name.split('_')[0],
                'mjd_start': data['mjd_start'],
                'n_time': data['n_time'],
                'n_freq': data['n_freq'],
                'dt_s': data['dt_s'],
                'freq_min': float(data['freq_MHz'].min()),
                'freq_max': float(data['freq_MHz'].max()),
                'frequency_mhz': float((data['freq_MHz'].min() + data['freq_MHz'].max()) / 2),
            })
            catalog['n_epochs'] += 1
            
            print_status(f"  Saved {out_name}", "SUCCESS")
            
        except Exception as e:
            print_status(f"  Failed to process {dat_file.name}: {e}", "ERROR")
    
    # Save catalog
    catalog_path = JIAMUSI_PROC_DIR / 'jiamusi_epoch_catalog.json'
    with open(catalog_path, 'w') as f:
        json.dump(catalog, f, indent=2, cls=NpEncoder)
    
    print_status(f"Catalog saved: {catalog_path}", "SUCCESS")
    return catalog


def compute_secondary_spectrum(dynspec: np.ndarray, dt_s: float, freq_MHz: np.ndarray) -> tuple:
    """Compute secondary spectrum from dynamic spectrum.

    Correct definition: S = |FFT2(delta I)|^2 where delta I is the
    mean-subtracted (normalized) dynamic spectrum.
    """
    # Normalize: subtract mean to get intensity fluctuations
    delta_I = dynspec - np.mean(dynspec)

    # 2D FFT and shift
    F = np.fft.fft2(delta_I)
    F = np.fft.fftshift(F)

    # Secondary spectrum is squared magnitude
    secondary = np.abs(F) ** 2

    # Convert to delay-Doppler coordinates
    n_time, n_freq = dynspec.shape
    df = freq_MHz[1] - freq_MHz[0]  # MHz
    dt = dt_s  # seconds

    # Delay axis (tau) in microseconds
    tau = np.fft.fftshift(np.fft.fftfreq(n_freq, df * 1e6)) * 1e6
    # Doppler axis (f_D) in mHz
    f_d = np.fft.fftshift(np.fft.fftfreq(n_time, dt)) * 1000

    return secondary, tau, f_d


def analyze_jiamusi_data():
    """Analyze processed Jiamusi data."""
    print_status("Analyzing Jiamusi data...", "INFO")
    
    # Load catalog
    catalog_path = JIAMUSI_PROC_DIR / 'jiamusi_epoch_catalog.json'
    with open(catalog_path, 'r') as f:
        catalog = json.load(f)
    
    results = {}
    
    for epoch in catalog['epochs']:
        file_path = JIAMUSI_PROC_DIR / epoch['file']
        
        print_status(f"  Processing {epoch['file']}...", "INFO")
        
        # Load data
        arr = np.load(file_path)
        dynspec = arr['dynspec']
        freq_MHz = arr['freq_MHz']
        dt_s = float(arr['dt_s'])
        
        # Compute secondary spectrum
        secondary, tau, f_d = compute_secondary_spectrum(dynspec, dt_s, freq_MHz)
        
        # Store results
        results[epoch['file']] = {
            'secondary_max': float(np.max(secondary)),
            'n_time': int(arr['n_time']),
            'n_freq': int(arr['n_freq']),
            'mjd_start': float(arr['mjd_start']),
            'source': epoch['source'],
        }
        
        print_status(f"    Max secondary: {results[epoch['file']]['secondary_max']:.2f}", "INFO")
    
    # Save results
    results_path = JIAMUSI_RESULTS_DIR / 'step_029_jiamusi_secondary_spectra.json'
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2, cls=NpEncoder)
    
    print_status(f"Results saved: {results_path}", "SUCCESS")
    
    # Summary by pulsar
    by_pulsar = defaultdict(list)
    for fname, res in results.items():
        by_pulsar[res['source']].append(res)
    
    print_status("\n  Jiamusi Analysis Summary:", "INFO")
    print_status("  " + "="*50, "INFO")
    for pulsar, epochs in sorted(by_pulsar.items()):
        print_status(f"  {pulsar}: {len(epochs)} epochs", "INFO")
    print_status("  " + "="*50, "INFO")
    
    return results


def main():
    """Main function for Jiamusi analysis."""
    # Logger is set by run_pipeline.py via set_step_logger()
    # Do not create a new logger here to avoid overriding the pipeline's logger
    
    print_status("="*70, "INFO")
    print_status("Jiamusi Pulsar Data Analysis", "INFO")
    print_status("="*70, "INFO")
    print_status("Source: Wang et al. (2018), A&A 618, A186", "INFO")
    print_status("Data: http://zmtt.bao.ac.cn/psr-jms/", "INFO")
    print_status("="*70, "INFO")
    
    # Download data
    download_jiamusi_data()
    
    # Process data
    catalog = process_jiamusi_data()
    
    # Analyze data
    results = analyze_jiamusi_data()
    
    # Final summary
    print_status("="*70, "INFO")
    print_status("Jiamusi Analysis Complete", "INFO")
    print_status(f"  Total epochs: {catalog['n_epochs']}", "INFO")
    print_status(f"  Raw data: {JIAMUSI_RAW_DIR}", "INFO")
    print_status(f"  Processed: {JIAMUSI_PROC_DIR}", "INFO")
    print_status(f"  Results: {JIAMUSI_RESULTS_DIR}", "INFO")
    print_status("="*70, "INFO")


if __name__ == "__main__":
    main()
