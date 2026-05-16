#!/usr/bin/env python3
"""
Step 001: Parse & Calibrate Dynamic Spectra

Reads .dynspec files produced by psrflux (PPTA DR2) and outputs
calibrated, RFI-masked dynamic spectra as NumPy .npz files ready
for secondary-spectrum analysis.

A .dynspec file is a binary format written by the psrflux command
from PSRCHIVE. Its structure is:
  - ASCII header section terminated by "# start_mjd ..."
  - Binary float32 data block: shape (n_time, n_freq) for intensity,
    with an optional parallel variance/weight array.

If psrchive / psrflux Python bindings are NOT installed, we fall back
to reading the binary block directly using numpy.

Outputs (per epoch):
  data/processed/j0437/<epoch>.npz
    keys: dynspec    (n_time x n_freq float32)
          freq_MHz   (n_freq,)
          dt_s       (time resolution in seconds)
          mjd_start  (float MJD)
          filename   (str)

Usage:
    python step_001_parse_dynspec.py [--workers N]
"""

import argparse
import json
import re
import struct
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Union

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.json_numpy import NpEncoder
from scripts.utils.logger import TEPLogger, print_status, set_step_logger
from scripts.utils.parallel_workers import configure_blas_thread_env, worker_count

RAW_DIR_SCINTOOLS = PROJECT_ROOT / "data" / "raw" / "scintools"
RAW_DIR_J0437 = PROJECT_ROOT / "data" / "raw" / "j0437"
RAW_DIR_CASPSR = PROJECT_ROOT / "data" / "raw" / "j0437" / "caspsr"
RAW_DIR_MPTA = (
    PROJECT_ROOT / "data" / "raw" / "mpta" / "data_august23_32ch" / "J0437-4715"
)
RAW_DIR_PPTA_DR2 = PROJECT_ROOT / "data" / "raw" / "j0437"
RAW_DIR_J1603 = PROJECT_ROOT / "data" / "raw" / "j1603"
RAW_DIR_MEERKAT = PROJECT_ROOT / "data" / "raw" / "meerkat" / "data" / "pdfb4"
RAW_DIR_JIAMUSI = PROJECT_ROOT / "data" / "raw" / "jiamusi"
PROC_DIR = PROJECT_ROOT / "data" / "processed" / "j0437"
PROC_DIR_J1603 = PROJECT_ROOT / "data" / "processed" / "j1603"
PROC_DIR_MEERKAT = PROJECT_ROOT / "data" / "processed" / "meerkat"
PROC_DIR_JIAMUSI = PROJECT_ROOT / "data" / "processed" / "jiamusi"

# --- ASCII dynspec reader (Scintools format) ---------------------------------─


def parse_ascii_dynspec(
    path: Path, exclusion_log: Optional[list] = None
) -> Optional[dict]:
    """Parse ASCII format dynamic spectrum from Scintools.

    Format: isub ichan time(min) freq(MHz) flux flux_err

    Parameters
    ----------
    path : Path
        Path to the dynspec file
    exclusion_log : list, optional
        List to append exclusion reasons for audit trail
    """

    def log_exclusion(reason: str):
        if exclusion_log is not None:
            exclusion_log.append(
                {"file": path.name, "reason": reason, "stage": "ascii_parse"}
            )

    fname = path.name

    try:
        fname = path.name

        # Read MJD from header
        mjd = None
        with open(path, "r") as f:
            for line in f:
                if line.startswith("# MJD0:"):
                    mjd = float(line.split(":")[1].strip())
                    break
        
        if mjd is None:
            log_exclusion("missing_mjd_in_header")
            return None

        # Load data (skip header lines starting with #)
        data = np.loadtxt(path, comments="#")

        if data.ndim != 2 or data.shape[1] < 6:
            log_exclusion(
                f"invalid_data_shape: ndim={data.ndim}, shape={data.shape if hasattr(data, 'shape') else 'N/A'}"
            )
            return None

        # Extract columns
        isub = data[:, 0].astype(int)
        ichan = data[:, 1].astype(int)
        time_min = data[:, 2]
        freq_MHz_data = data[:, 3]
        flux = data[:, 4]

        # Determine dimensions
        n_time = len(np.unique(isub))
        n_freq = len(np.unique(ichan))

        # Need at least 10 time samples
        if n_time < 10:
            log_exclusion(f"insufficient_time_samples: n_time={n_time} < 10")
            return None

        # Reshape into 2D dynamic spectrum
        dynspec = flux.reshape(n_time, n_freq)

        # Get frequency axis
        freq_MHz = freq_MHz_data.reshape(n_time, n_freq)[0, :]

        # Calculate time resolution (convert minutes to seconds)
        time_axis = time_min.reshape(n_time, n_freq)[:, 0]
        if n_time > 1:
            dt_s = np.mean(np.diff(time_axis)) * 60.0
        else:
            log_exclusion("insufficient_time_samples_for_dt")
            return None

        # Clean data
        dynspec = np.nan_to_num(dynspec, nan=0.0, posinf=0.0, neginf=0.0)

        # Check for valid data (reject overflow/corrupted values)
        if not np.all(np.isfinite(dynspec)):
            log_exclusion("non_finite_values_in_dynspec")
            return None

        max_val = np.max(np.abs(dynspec))
        if max_val > 1e10 or max_val == 0:
            log_exclusion(f"invalid_max_value: max_val={max_val:.2e}")
            return None

        if np.std(dynspec) < 1e-10:
            log_exclusion(f"insufficient_variation: std={np.std(dynspec):.2e}")
            return None

        # Bandpass normalization
        col_med = np.median(dynspec, axis=0)
        col_med[col_med == 0] = 1
        dynspec = dynspec / col_med[np.newaxis, :]

        # RFI masking
        dynspec = np.where(dynspec <= 0, np.nan, dynspec)

        chan_med = np.nanmedian(dynspec, axis=0)
        global_med = np.nanmedian(chan_med)
        if global_med > 0:
            rfi_chans = (chan_med > 10 * global_med) | (chan_med < global_med / 10)
            dynspec[:, rfi_chans] = np.nan

        row_means = np.nanmean(dynspec, axis=1)
        mu, sigma = np.nanmean(row_means), np.nanstd(row_means)
        if sigma > 0:
            bad_times = np.abs(row_means - mu) > 5 * sigma
            dynspec[bad_times, :] = np.nan

        # Final normalization
        col_med = np.nanmedian(dynspec, axis=0)
        col_med[col_med == 0] = np.nan
        with np.errstate(divide="ignore", invalid="ignore"):
            dynspec = dynspec / col_med[np.newaxis, :]
        dynspec = np.nan_to_num(dynspec, nan=0.0)

        print_status(f"Parsed ASCII dynspec: {n_time}x{n_freq}, dt={dt_s:.2f}s, mjd={mjd:.5f}", "DATA")
        print_status(f"  Intensity stats: mean={np.mean(dynspec):.3f}, std={np.std(dynspec):.3f}, max={np.max(dynspec):.3f}", "CALC")
        print_status(f"  RFI Masking: {np.sum(np.isnan(dynspec)) / dynspec.size * 100:.1f}% pixels masked", "CALC")

        return {
            "dynspec": dynspec.astype(np.float32),
            "freq_MHz": freq_MHz.astype(np.float64),
            "dt_s": float(dt_s),
            "mjd_start": float(mjd),
            "n_time": int(n_time),
            "n_freq": int(n_freq),
            "filename": fname,
        }

    except Exception as e:
        log_exclusion(f"parse_error: {str(e)}")
        return None


# --- PSRFITS dynspec reader (legacy .dly format) ---------------------------──


def parse_psrfits_dynspec(
    path: Path, exclusion_log: Optional[list] = None
) -> Optional[dict]:
    """Parse PSRFITS .dly dynamic spectrum files.

    Parameters
    ----------
    path : Path
        Path to the .dly file
    exclusion_log : list, optional
        List to append exclusion reasons for audit trail
    """
    from astropy.io import fits

    def log_exclusion(reason: str):
        if exclusion_log is not None:
            exclusion_log.append(
                {"file": path.name, "reason": reason, "stage": "psrfits_parse"}
            )

    try:
        fname = path.name
        mjd = None  # Must be extracted from header

        with fits.open(path, memmap=False) as hdul:
            primary = hdul[0].header
            imjd = primary.get("STT_IMJD")
            smjd = primary.get("STT_SMJD", 0)
            offs = primary.get("STT_OFFS", 0)
            
            if imjd is not None:
                mjd = float(imjd) + (float(smjd) + float(offs)) / 86400.0
            else:
                log_exclusion("missing_stt_imjd_in_primary_header")
                return None

            if "SUBINT" not in hdul:
                return None

            subint = hdul["SUBINT"]
            raw_data = subint.data["DATA"]

            if raw_data.ndim != 4:
                log_exclusion(f"invalid_data_ndim: ndim={raw_data.ndim} != 4")
                return None

            nsub, npol, nchan, nbin = raw_data.shape

            # Need at least 10 time samples
            if nsub < 10:
                log_exclusion(f"insufficient_time_samples: nsub={nsub} < 10")
                return None

            # Get scaling and extract intensity
            dat_scl = subint.data["DAT_SCL"]
            dat_offs = subint.data["DAT_OFFS"]

            scl_expanded = dat_scl[:, :, :, np.newaxis]
            offs_expanded = dat_offs[:, :, :, np.newaxis]
            scaled_data = raw_data.astype(np.float32) * scl_expanded + offs_expanded

            if npol >= 2:
                intensity = scaled_data[:, 0, :, :] + scaled_data[:, 1, :, :]
            else:
                intensity = scaled_data[:, 0, :, :]

            dynspec = np.mean(intensity, axis=-1).astype(np.float32)

            dat_freq = subint.data["DAT_FREQ"][0]
            freq_MHz = dat_freq.astype(np.float64)

            tsubint = (
                subint.data["TSUBINT"][0] if "TSUBINT" in subint.columns.names else None
            )
            if tsubint is None or nsub == 0:
                log_exclusion("missing_tsubint_or_nsub")
                return None
            dt_s = float(tsubint / nsub)

            # Check for valid data (reject overflow/corrupted values)
            if not np.all(np.isfinite(dynspec)):
                log_exclusion("non_finite_values_in_dynspec")
                return None

            max_val = np.max(np.abs(dynspec))
            if max_val > 1e10 or max_val == 0:
                # Reject files with overflow values or all zeros
                log_exclusion(f"invalid_max_value: max_val={max_val:.2e}")
                return None

            if np.std(dynspec) < 1e-10:
                # No variation - not useful for scintillation analysis
                log_exclusion(f"insufficient_variation: std={np.std(dynspec):.2e}")
                return None

            dynspec = np.nan_to_num(dynspec, nan=0.0, posinf=0.0, neginf=0.0)

            col_med = np.median(dynspec, axis=0)
            col_med[col_med == 0] = 1
            dynspec = dynspec / col_med[np.newaxis, :]

            dynspec = np.where(dynspec <= 0, np.nan, dynspec)

            chan_med = np.nanmedian(dynspec, axis=0)
            global_med = np.nanmedian(chan_med)
            if global_med > 0:
                rfi_chans = (chan_med > 10 * global_med) | (chan_med < global_med / 10)
                dynspec[:, rfi_chans] = np.nan

            row_means = np.nanmean(dynspec, axis=1)
            mu, sigma = np.nanmean(row_means), np.nanstd(row_means)
            if sigma > 0:
                bad_times = np.abs(row_means - mu) > 5 * sigma
                dynspec[bad_times, :] = np.nan

            col_med = np.nanmedian(dynspec, axis=0)
            col_med[col_med == 0] = np.nan
            with np.errstate(divide="ignore", invalid="ignore"):
                dynspec = dynspec / col_med[np.newaxis, :]
            dynspec = np.nan_to_num(dynspec, nan=0.0)

            print_status(f"Parsed PSRFITS dynspec: {nsub}x{nchan}, dt={dt_s:.2f}s, mjd={mjd:.5f}", "DATA")
            print_status(f"  Intensity stats: mean={np.mean(dynspec):.3f}, std={np.std(dynspec):.3f}, max={np.max(dynspec):.3f}", "CALC")
            print_status(f"  RFI Masking: {np.sum(np.isnan(dynspec)) / dynspec.size * 100:.1f}% pixels masked", "CALC")

            return {
                "dynspec": dynspec.astype(np.float32),
                "freq_MHz": freq_MHz.astype(np.float64),
                "dt_s": float(dt_s),
                "mjd_start": float(mjd),
                "n_time": int(nsub),
                "n_freq": int(nchan),
                "filename": fname,
            }

    except Exception as e:
        log_exclusion(f"parse_error: {str(e)}")
        return None


# --- Jiamusi format parser (ASCII 3-column) ---------------------------------


def parse_jiamusi_dynspec(
    path: Path, exclusion_log: Optional[list] = None
) -> Optional[dict]:
    """Parse Jiamusi 66m telescope dynamic spectrum format.

    Format: 3 columns - frequency(MHz), MJD, normalized_intensity
    Data is in long format (one row per time-freq pixel).

    Parameters
    ----------
    path : Path
        Path to the .dat file
    exclusion_log : list, optional
        List to append exclusion reasons for audit trail
    """

    def log_exclusion(reason: str):
        if exclusion_log is not None:
            exclusion_log.append(
                {"file": path.name, "reason": reason, "stage": "jiamusi_parse"}
            )

    fname = path.name

    try:
        # Read data
        data = np.loadtxt(path)

        # Extract columns
        freq = data[:, 0]  # MHz
        mjd = data[:, 1]
        intensity = data[:, 2]

        # Get unique values
        unique_freqs = np.unique(freq)
        unique_mjds = np.unique(mjd)

        n_freq = len(unique_freqs)
        n_time = len(unique_mjds)

        if n_time < 10 or n_freq < 10:
            log_exclusion(f"insufficient_data: n_time={n_time}, n_freq={n_freq}")
            return None

        # Reshape to 2D array (time, freq)
        dynspec = intensity.reshape((n_time, n_freq))

        # Calculate time resolution (convert MJD diff to seconds)
        if n_time > 1:
            dt_s = (unique_mjds[1] - unique_mjds[0]) * 86400.0
        else:
            dt_s = 30.0

        # Clean data
        dynspec = np.nan_to_num(dynspec, nan=0.0, posinf=0.0, neginf=0.0)

        # Check for valid data
        if not np.all(np.isfinite(dynspec)):
            log_exclusion("non_finite_values_in_dynspec")
            return None

        max_val = np.max(np.abs(dynspec))
        if max_val > 1e10 or max_val == 0:
            log_exclusion(f"invalid_max_value: max_val={max_val:.2e}")
            return None

        if np.std(dynspec) < 1e-10:
            log_exclusion(f"insufficient_variation: std={np.std(dynspec):.2e}")
            return None

        # Bandpass normalization
        col_med = np.median(dynspec, axis=0)
        col_med[col_med == 0] = 1
        dynspec = dynspec / col_med[np.newaxis, :]

        # RFI masking
        dynspec = np.where(dynspec <= 0, np.nan, dynspec)

        chan_med = np.nanmedian(dynspec, axis=0)
        global_med = np.nanmedian(chan_med)
        if global_med > 0:
            rfi_chans = (chan_med > 10 * global_med) | (chan_med < global_med / 10)
            dynspec[:, rfi_chans] = np.nan

        row_means = np.nanmean(dynspec, axis=1)
        mu, sigma = np.nanmean(row_means), np.nanstd(row_means)
        if sigma > 0:
            bad_times = np.abs(row_means - mu) > 5 * sigma
            dynspec[bad_times, :] = np.nan

        # Final normalization
        col_med = np.nanmedian(dynspec, axis=0)
        col_med[col_med == 0] = np.nan
        with np.errstate(divide="ignore", invalid="ignore"):
            dynspec = dynspec / col_med[np.newaxis, :]
        dynspec = np.nan_to_num(dynspec, nan=0.0)

        print_status(f"Parsed Jiamusi dynspec: {n_time}x{n_freq}, dt={dt_s:.2f}s, mjd={unique_mjds[0]:.5f}", "DATA")
        print_status(f"  Intensity stats: mean={np.mean(dynspec):.3f}, std={np.std(dynspec):.3f}, max={np.max(dynspec):.3f}", "CALC")
        print_status(f"  RFI Masking: {np.sum(np.isnan(dynspec)) / dynspec.size * 100:.1f}% pixels masked", "CALC")

        return {
            "dynspec": dynspec.astype(np.float32),
            "freq_MHz": unique_freqs.astype(np.float64),
            "dt_s": float(dt_s),
            "mjd_start": float(unique_mjds[0]),
            "n_time": int(n_time),
            "n_freq": int(n_freq),
            "filename": fname,
        }

    except Exception as e:
        log_exclusion(f"parse_error: {str(e)}")
        return None


def parse_dynspec_file(
    path: Path, exclusion_log: Optional[list] = None
) -> Optional[dict]:
    """Parse dynamic spectrum file (auto-detects format).

    Parameters
    ----------
    path : Path
        Path to the data file
    exclusion_log : list, optional
        List to append exclusion reasons for audit trail
    """
    # Check file extension and content to determine format
    if path.suffix == ".dynspec" or path.name.endswith(".dynspec"):
        # ASCII format from Scintools
        return parse_ascii_dynspec(path, exclusion_log=exclusion_log)
    elif path.suffix == ".dly" or path.name.endswith(".dly"):
        # PSRFITS format
        return parse_psrfits_dynspec(path, exclusion_log=exclusion_log)
    elif path.suffix == ".dat" and "jiamusi" in str(path).lower():
        # Jiamusi 66m telescope format
        return parse_jiamusi_dynspec(path, exclusion_log=exclusion_log)
    elif path.suffix == ".dat":
        # Try Jiamusi format first for .dat files
        result = parse_jiamusi_dynspec(path, exclusion_log=exclusion_log)
        if result is not None:
            return result
        # Fall back to other formats
        result = parse_ascii_dynspec(path, exclusion_log=exclusion_log)
        if result is not None:
            return result
        return parse_psrfits_dynspec(path, exclusion_log=exclusion_log)
    else:
        # Try ASCII first, then FITS
        result = parse_ascii_dynspec(path, exclusion_log=exclusion_log)
        if result is not None:
            return result
        return parse_psrfits_dynspec(path, exclusion_log=exclusion_log)


def process_one(src_path: Path, out_dir: Path, force: bool = False) -> str:
    """Worker: parse one .dynspec, save .npz, return status string."""
    stem = src_path.stem
    out_path = out_dir / f"{stem}.npz"
    if out_path.exists() and not force:
        return f"SKIP {src_path.name}"

    excl: list = []
    result = parse_dynspec_file(src_path, exclusion_log=excl)
    if result is None:
        reason = excl[-1].get("reason", "unknown") if excl else "unknown"
        safe = str(reason).replace("\t", " ").replace("\n", " ")
        return f"FAIL\t{src_path.name}\t{safe}"

    np.savez_compressed(out_path, **result)
    return f"OK   {src_path.name}  ({result['n_time']}x{result['n_freq']})"


def _record_parse_failure(exclusions: list, msg: str) -> None:
    """Parse process_one return lines into structured exclusion records."""
    if msg.startswith("FAIL\t"):
        parts = msg.split("\t", 2)
        if len(parts) >= 3:
            exclusions.append({"file": parts[1], "reason": parts[2]})
        elif len(parts) >= 2:
            exclusions.append({"file": parts[1], "reason": "unknown"})
    elif msg.startswith("FAIL "):
        exclusions.append({"file": msg[5:].strip(), "reason": "unknown"})


# --- Main ---------------------------------------------------------------------


def main(workers: int = None, force: bool = False):
    # Logger is set by run_pipeline.py via set_step_logger()
    # Do not create a new logger here to avoid overriding the pipeline's logger

    configure_blas_thread_env()
    if workers is None:
        workers = worker_count(role="io_bound", reserve=2)

    PROC_DIR.mkdir(parents=True, exist_ok=True)
    PROC_DIR_J1603.mkdir(parents=True, exist_ok=True)

    print_status("=" * 60, "TITLE")
    print_status("TEP-J0437 Step 001: Parse Dynamic Spectra", "TITLE")
    print_status("=" * 60, "TITLE")
    print_status(f"Workers: {workers}", "INFO")

    # Collect files from all directories
    j0437_files = []
    j1603_files = []

    # Scintools ASCII dynspec files (preferred - proper time sampling)
    if RAW_DIR_SCINTOOLS.exists():
        scintools_files = list(RAW_DIR_SCINTOOLS.glob("*.dynspec"))
        j0437_files.extend(scintools_files)
        print_status(
            f"Found {len(scintools_files)} .dynspec files in scintools/", "INFO"
        )

    # CASPSR ASCII dynspec files (additional epochs from 2010)
    if RAW_DIR_CASPSR.exists():
        caspsr_files = list(RAW_DIR_CASPSR.glob("*.dynspec"))
        j0437_files.extend(caspsr_files)
        print_status(
            f"Found {len(caspsr_files)} .dynspec files in j0437/caspsr/", "INFO"
        )

    # MPTA PSRFITS .dly files (fallback - may have insufficient time samples)
    if RAW_DIR_MPTA.exists():
        mpta_files = list(RAW_DIR_MPTA.rglob("*.dly"))
        j0437_files.extend(mpta_files)
        print_status(f"Found {len(mpta_files)} .dly files in mpta/", "INFO")

    # PPTA DR2 pdfb4 ASCII dynspec files (2296 epochs from 2008-2018)
    if RAW_DIR_PPTA_DR2.exists():
        ppta_files = list(RAW_DIR_PPTA_DR2.rglob("*.dynspec"))
        j0437_files.extend(ppta_files)
        print_status(
            f"Found {len(ppta_files)} .dynspec files in PPTA DR2 pdfb4/", "INFO"
        )

    # J1603-7202 ASCII dynspec files (control pulsar, 2002-2016)
    if RAW_DIR_J1603.exists():
        j1603_files = list(RAW_DIR_J1603.rglob("*.dynspec"))
        print_status(
            f"Found {len(j1603_files)} .dynspec files in j1603/ (control pulsar)",
            "INFO",
        )

    j0437_files = sorted(set(j0437_files))  # Remove duplicates
    j1603_files = sorted(set(j1603_files))

    print_status(f"Total J0437 files: {len(j0437_files)}", "INFO")
    print_status(f"Total J1603 files: {len(j1603_files)}", "INFO")

    if not j0437_files and not j1603_files:
        print_status(
            "No dynamic spectrum files found\n"
            "  -> Check data/raw/scintools/, data/raw/j0437/, and data/raw/j1603/",
            "WARNING",
        )
        return

    # Process J0437 files
    if j0437_files:
        print_status(f"Processing {len(j0437_files)} J0437 files...", "INFO")
        ok_count = fail_count = skip_count = 0
        exclusions: list = []

        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(process_one, f, PROC_DIR, force): f for f in j0437_files}
            for fut in as_completed(futures):
                msg = fut.result()
                if msg.startswith("OK"):
                    ok_count += 1
                elif msg.startswith("SKIP"):
                    skip_count += 1
                else:
                    fail_count += 1
                    _record_parse_failure(exclusions, msg)
                print_status(msg, "INFO")

        print_status(f"\nJ0437 Processing Summary:", "INFO")
        print_status(
            f"  OK: {ok_count}, Failed: {fail_count}, Skipped: {skip_count}", "INFO"
        )

        # Write J0437 epoch catalogue
        epochs = []
        for npz_path in sorted(PROC_DIR.glob("*.npz")):
            arr = np.load(npz_path, allow_pickle=False)
            epochs.append(
                {
                    "file": npz_path.name,
                    "mjd_start": float(arr["mjd_start"]),
                    "n_time": int(arr["n_time"]),
                    "n_freq": int(arr["n_freq"]),
                    "dt_s": float(arr["dt_s"]),
                    "freq_min": float(arr["freq_MHz"].min()),
                    "freq_max": float(arr["freq_MHz"].max()),
                }
            )
        epochs.sort(key=lambda x: x["mjd_start"])

        catalog_path = PROJECT_ROOT / "data" / "processed" / "j0437_epoch_catalog.json"
        with open(catalog_path, "w") as cf:
            json.dump({"n_epochs": len(epochs), "epochs": epochs}, cf, indent=2, cls=NpEncoder)

        excl_path = PROJECT_ROOT / "data" / "processed" / "j0437_parse_exclusions.json"
        with open(excl_path, "w", encoding="utf-8") as xf:
            json.dump(
                {"n_excluded": len(exclusions), "exclusions": exclusions},
                xf,
                indent=2,
                cls=NpEncoder,
            )

        print_status(
            f"J0437: {ok_count} new, {skip_count} cached, {fail_count} failed.\n"
            f"  Epoch catalogue -> {catalog_path}\n"
            f"  Parse exclusion audit -> {excl_path}",
            "SUCCESS",
        )

    # Process J1603 files
    if j1603_files:
        print_status(
            f"Processing {len(j1603_files)} J1603 files (control pulsar)...", "INFO"
        )
        ok_count = fail_count = skip_count = 0
        exclusions: list = []

        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(process_one, f, PROC_DIR_J1603, force): f for f in j1603_files
            }
            for fut in as_completed(futures):
                msg = fut.result()
                if msg.startswith("OK"):
                    ok_count += 1
                elif msg.startswith("SKIP"):
                    skip_count += 1
                else:
                    fail_count += 1
                    _record_parse_failure(exclusions, msg)
                print_status(msg, "INFO")

        print_status(f"\nJ1603 Processing Summary:", "INFO")
        print_status(
            f"  OK: {ok_count}, Failed: {fail_count}, Skipped: {skip_count}", "INFO"
        )

        # Write J1603 epoch catalogue
        epochs = []
        for npz_path in sorted(PROC_DIR_J1603.glob("*.npz")):
            arr = np.load(npz_path, allow_pickle=False)
            epochs.append(
                {
                    "file": npz_path.name,
                    "mjd_start": float(arr["mjd_start"]),
                    "n_time": int(arr["n_time"]),
                    "n_freq": int(arr["n_freq"]),
                    "dt_s": float(arr["dt_s"]),
                    "freq_min": float(arr["freq_MHz"].min()),
                    "freq_max": float(arr["freq_MHz"].max()),
                }
            )
        epochs.sort(key=lambda x: x["mjd_start"])

        catalog_path = PROJECT_ROOT / "data" / "processed" / "j1603_epoch_catalog.json"
        with open(catalog_path, "w") as cf:
            json.dump({"n_epochs": len(epochs), "epochs": epochs}, cf, indent=2, cls=NpEncoder)

        excl_path = PROJECT_ROOT / "data" / "processed" / "j1603_parse_exclusions.json"
        with open(excl_path, "w", encoding="utf-8") as xf:
            json.dump(
                {"n_excluded": len(exclusions), "exclusions": exclusions},
                xf,
                indent=2,
                cls=NpEncoder,
            )

        print_status(
            f"J1603: {ok_count} new, {skip_count} cached, {fail_count} failed.\n"
            f"  Epoch catalogue -> {catalog_path}\n"
            f"  Parse exclusion audit -> {excl_path}",
            "SUCCESS",
        )

    print_status("\n" + "=" * 70, "TITLE")
    print_status("STEP 001 COMPLETED SUCCESSFULLY", "SUCCESS")
    print_status("=" * 70, "TITLE")

    # Process MeerKAT files
    meerkat_files = (
        list(RAW_DIR_MEERKAT.glob("*.dynspec")) if RAW_DIR_MEERKAT.exists() else []
    )
    if meerkat_files:
        print_status(f"Processing {len(meerkat_files)} MeerKAT files...", "INFO")
        ok_count = fail_count = skip_count = 0
        exclusions: list = []

        PROC_DIR_MEERKAT.mkdir(parents=True, exist_ok=True)

        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(process_one, f, PROC_DIR_MEERKAT, force): f for f in meerkat_files
            }
            for fut in as_completed(futures):
                msg = fut.result()
                if msg.startswith("OK"):
                    ok_count += 1
                elif msg.startswith("SKIP"):
                    skip_count += 1
                else:
                    fail_count += 1
                    _record_parse_failure(exclusions, msg)
                print_status(msg, "INFO")

        # Write MeerKAT epoch catalogue
        epochs = []
        for npz_path in sorted(PROC_DIR_MEERKAT.glob("*.npz")):
            arr = np.load(npz_path, allow_pickle=False)
            epochs.append(
                {
                    "file": npz_path.name,
                    "mjd_start": float(arr["mjd_start"]),
                    "n_time": int(arr["n_time"]),
                    "n_freq": int(arr["n_freq"]),
                    "dt_s": float(arr["dt_s"]),
                    "freq_min": float(arr["freq_MHz"].min()),
                    "freq_max": float(arr["freq_MHz"].max()),
                }
            )
        epochs.sort(key=lambda x: x["mjd_start"])

        catalog_path = (
            PROJECT_ROOT / "data" / "processed" / "meerkat_epoch_catalog.json"
        )
        with open(catalog_path, "w") as cf:
            json.dump({"n_epochs": len(epochs), "epochs": epochs}, cf, indent=2, cls=NpEncoder)

        excl_path = PROJECT_ROOT / "data" / "processed" / "meerkat_parse_exclusions.json"
        with open(excl_path, "w", encoding="utf-8") as xf:
            json.dump(
                {"n_excluded": len(exclusions), "exclusions": exclusions},
                xf,
                indent=2,
                cls=NpEncoder,
            )

        print_status(
            f"MeerKAT: {ok_count} new, {skip_count} cached, {fail_count} failed.\n"
            f"  Epoch catalogue -> {catalog_path}\n"
            f"  Parse exclusion audit -> {excl_path}",
            "SUCCESS",
        )

    # Process Jiamusi files
    jiamusi_files = (
        list(RAW_DIR_JIAMUSI.glob("*.dat")) if RAW_DIR_JIAMUSI.exists() else []
    )
    if jiamusi_files:
        print_status(f"Processing {len(jiamusi_files)} Jiamusi files...", "INFO")
        ok_count = fail_count = skip_count = 0
        exclusions: list = []

        PROC_DIR_JIAMUSI.mkdir(parents=True, exist_ok=True)

        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(process_one, f, PROC_DIR_JIAMUSI, force): f for f in jiamusi_files
            }
            for fut in as_completed(futures):
                msg = fut.result()
                if msg.startswith("OK"):
                    ok_count += 1
                elif msg.startswith("SKIP"):
                    skip_count += 1
                else:
                    fail_count += 1
                    _record_parse_failure(exclusions, msg)
                print_status(msg, "INFO")

        # Write Jiamusi epoch catalogue
        epochs = []
        for npz_path in sorted(PROC_DIR_JIAMUSI.glob("*.npz")):
            arr = np.load(npz_path, allow_pickle=False)
            epochs.append(
                {
                    "file": npz_path.name,
                    "mjd_start": float(arr["mjd_start"]),
                    "n_time": int(arr["n_time"]),
                    "n_freq": int(arr["n_freq"]),
                    "dt_s": float(arr["dt_s"]),
                    "freq_min": float(arr["freq_MHz"].min()),
                    "freq_max": float(arr["freq_MHz"].max()),
                }
            )
        epochs.sort(key=lambda x: x["mjd_start"])

        catalog_path = (
            PROJECT_ROOT / "data" / "processed" / "jiamusi_epoch_catalog.json"
        )
        with open(catalog_path, "w") as cf:
            json.dump({"n_epochs": len(epochs), "epochs": epochs}, cf, indent=2, cls=NpEncoder)

        excl_path = PROJECT_ROOT / "data" / "processed" / "jiamusi_parse_exclusions.json"
        with open(excl_path, "w", encoding="utf-8") as xf:
            json.dump(
                {"n_excluded": len(exclusions), "exclusions": exclusions},
                xf,
                indent=2,
                cls=NpEncoder,
            )

        print_status(
            f"Jiamusi: {ok_count} new, {skip_count} cached, {fail_count} failed.\n"
            f"  Epoch catalogue -> {catalog_path}\n"
            f"  Parse exclusion audit -> {excl_path}",
            "SUCCESS",
        )


def step_main(logger=None, verbose=True, force=False):
    """Pipeline entry point for Step 001."""
    if logger:
        set_step_logger(logger)
    return main(force=force)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse J0437 .dynspec files")
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Process pool size (default: auto from CPU count; override via TEP_WORKERS / TEP_STEP001_MAX_WORKERS)",
    )
    args = parser.parse_args()

    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger("step_001", str(log_dir / "step_001_parse_dynspec.log"))
    set_step_logger(logger)

    main(workers=args.workers)
