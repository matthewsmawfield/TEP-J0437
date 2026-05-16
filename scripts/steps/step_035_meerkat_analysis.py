#!/usr/bin/env python3
"""
STEP 035: MeerKAT Data Availability Audit

This step inventories the MeerKAT data present in the repository and
produces an honest report on what is available vs what is missing for
a J0437-4715 TEP analysis.

DATA SITUATION:
- The repository contains 5 MeerKAT .dynspec files in data/raw/meerkat/data/pdfb4/
  (J0908-1739, J0922-0638, J1731-4744, sample, test1).
- NONE of these are J0437-4715.
- J0437-4715 MeerKAT data from the Thousand-Pulsar-Array has not been sourced.
- The meerkat_all.tar.gz file is an HTML error page, not a real archive.

Consequently, a dedicated MeerKAT TEP analysis for J0437 cannot be performed.
The primary J0437 analysis uses PPTA DR2 Parkes/CASPSR data, which is the
currently operative data source.
"""

import json
from pathlib import Path

import numpy as np
from scripts.utils.logger import print_status
from scripts.utils.json_numpy import NpEncoder

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

MEERKAT_RAW_DIR = PROJECT_ROOT / "data" / "raw" / "meerkat" / "data" / "pdfb4"
MEERKAT_PROC_DIR = PROJECT_ROOT / "data" / "processed" / "meerkat"
CATALOG_PATH = PROJECT_ROOT / "data" / "processed" / "meerkat_epoch_catalog.json"


def inventory_meerkat_data():
    """Scan raw and processed MeerKAT directories and return findings."""
    raw_files = []
    if MEERKAT_RAW_DIR.exists():
        raw_files = sorted([f.name for f in MEERKAT_RAW_DIR.glob("*.dynspec")])

    proc_files = []
    if MEERKAT_PROC_DIR.exists():
        proc_files = sorted([f.name for f in MEERKAT_PROC_DIR.glob("*.npz")])

    catalog_exists = CATALOG_PATH.exists()

    return {
        "raw_dynspec_files": raw_files,
        "processed_npz_files": proc_files,
        "catalog_exists": catalog_exists,
    }


def compute_secondary_spectrum(dynspec, dt_s, freq_MHz):
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

    return secondary


def main():
    print_status("=" * 70, "INFO")
    print_status("STEP 035: MeerKAT Data Availability Audit", "TITLE")
    print_status("=" * 70, "INFO")

    inventory = inventory_meerkat_data()
    raw_files = inventory["raw_dynspec_files"]
    proc_files = inventory["processed_npz_files"]
    catalog_exists = inventory["catalog_exists"]

    print_status(f"Raw .dynspec files found: {len(raw_files)}", "INFO")
    for f in raw_files:
        print_status(f"  - {f}", "INFO")

    print_status(f"Processed .npz files found: {len(proc_files)}", "INFO")
    if proc_files:
        for f in proc_files:
            print_status(f"  - {f}", "INFO")
    else:
        print_status("  (none)", "INFO")

    print_status(f"Epoch catalog exists: {catalog_exists}", "INFO")

    # Determine J0437 availability
    j0437_present = any("J0437" in f for f in raw_files)
    if not j0437_present:
        print_status("", "INFO")
        print_status("J0437-4715 MeerKAT data: NOT AVAILABLE (expected for this repository)", "WARNING")
        print_status("The raw MeerKAT files present are for other pulsars", "INFO")
        print_status("(J0908-1739, J0922-0638, J1731-4744, sample, test1).", "INFO")
        print_status("A dedicated MeerKAT TEP analysis for J0437 cannot be performed.", "INFO")
        print_status("The primary J0437 analysis uses PPTA DR2 Parkes/CASPSR data.", "INFO")

    print_status("", "INFO")
    print_status("=" * 70, "INFO")
    print_status("STEP 035 SUMMARY", "INFO")
    print_status("=" * 70, "INFO")
    if j0437_present:
        print_status("J0437 MeerKAT data found. Analysis would proceed.", "SUCCESS")
    else:
        print_status("J0437 MeerKAT data missing. Step returns no-op.", "WARNING")
    print_status("=" * 70, "INFO")

    # Save audit report
    results = {
        "step": "035_meerkat_analysis",
        "status": "no_op" if not j0437_present else "completed",
        "j0437_data_available": j0437_present,
        "raw_files_present": raw_files,
        "processed_files_present": proc_files,
        "catalog_exists": catalog_exists,
        "note": (
            "J0437-4715 MeerKAT data is not present in this repository. "
            "Available MeerKAT files are for other pulsars. "
            "Primary J0437 TEP analysis uses PPTA DR2 Parkes/CASPSR data."
        ),
    }

    results_path = RESULTS_DIR / "step_035_meerkat_audit.json"
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2, cls=NpEncoder)
    print_status(f"Audit report saved to {results_path}", "INFO")

    return True  # Audit completed, whether or not J0437 MeerKAT data are present.


if __name__ == "__main__":
    main()
