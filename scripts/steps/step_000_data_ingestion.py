#!/usr/bin/env python3
"""
================================================================================
STEP 000: DATA INGESTION FOR TEP-J0437
================================================================================

Acquires and verifies Parkes/PPTA DR2 dynamic spectra from the CSIRO Data
Access Portal (DAP) for PSR J0437-4715 and PSR J1603-7202, then optionally
fetches supplementary J0437 epochs from the Scintools ATNF archive.

Primary PPTA collections are resolved by DOI, downloaded through the public
/ws/v2 API, and verified against the CSIRO import SHA-256 manifests shipped
under data/raw/{j0437,j1603}/metadata/. Any checksum mismatch or missing
manifest entry fails the step.

DATA SOURCES:
-------------
1. CSIRO DAP PPTA DR2 dynamic spectra
   - J0437-4715: DOI 10.25919/5f3cd2bc1c213 -> data/raw/j0437/
   - J1603-7202: DOI 10.25919/82f5-mh79 -> data/raw/j1603/

2. Scintools ATNF archive (supplementary J0437 epochs)
   - URL: https://scintools.atnf.csiro.au/data/J0437-4715/
   - Files are stored under data/raw/scintools/ for step_001 discovery.

3. Jiamusi: handled separately by step_029_jiamusi_analysis.py

AUTHOR: TEP Analysis Framework
VERSION: 0.4
================================================================================
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.json_numpy import NpEncoder
from scripts.utils.logger import set_step_logger, print_status

RAW_DIR_J0437 = PROJECT_ROOT / "data" / "raw" / "j0437"
RAW_DIR_J1603 = PROJECT_ROOT / "data" / "raw" / "j1603"
# Scintools supplement files must live under data/raw/scintools/ so step_001 discovers them.
RAW_DIR_SCINTOOLS = PROJECT_ROOT / "data" / "raw" / "scintools"
CATALOG_PATH = PROJECT_ROOT / "data" / "raw" / "ingestion_catalog.json"
INGESTION_AUDIT_PATH = PROJECT_ROOT / "data" / "raw" / "ingestion_audit.json"
SCINTOOLS_BASE = "https://scintools.atnf.csiro.au/data/J0437-4715/"
DAP_API_BASE = "https://data.csiro.au/dap/ws/v2/collections"
USER_AGENT = "TEP-J0437-Data-Ingestion/1.0"
FILE_PAGE_SIZE = 1000
DOWNLOAD_TIMEOUT_S = 300
DEFAULT_WORKERS = max(4, min(16, (os.cpu_count() or 4)))
PROGRESS_INTERVAL = 25

SCINTOOLS_TARGET_FILES = [
    "p111220_074112.rf.pcm.dynspec",
    "p111220_084944.rf.pcm.dynspec",
    "p111220_095816.rf.pcm.dynspec",
    "p111220_110656.rf.pcm.dynspec",
    "p111220_121536.rf.pcm.dynspec",
    "p111221_074112.rf.pcm.dynspec",
    "p111222_074112.rf.pcm.dynspec",
    "p111223_074112.rf.pcm.dynspec",
    "p120119_074112.rf.pcm.dynspec",
    "p120218_074112.rf.pcm.dynspec",
]


@dataclass(frozen=True)
class PPTACollection:
    pulsar: str
    doi: str
    raw_dir: Path
    checksum_manifest: Path


PPTA_COLLECTIONS: Tuple[PPTACollection, ...] = (
    PPTACollection(
        pulsar="J0437-4715",
        doi="10.25919/5f3cd2bc1c213",
        raw_dir=RAW_DIR_J0437,
        checksum_manifest=RAW_DIR_J0437 / "metadata" / "collection_import_sha256sum.txt",
    ),
    PPTACollection(
        pulsar="J1603-7202",
        doi="10.25919/82f5-mh79",
        raw_dir=RAW_DIR_J1603,
        checksum_manifest=RAW_DIR_J1603 / "metadata" / "collection_import_sha256sum.txt",
    ),
)


def log_message(message: str, level: str = "INFO") -> None:
    print_status(message, level)


def compute_sha256(file_path: Path) -> str:
    sha256_hash = hashlib.sha256()
    with open(file_path, "rb") as handle:
        for byte_block in iter(lambda: handle.read(1024 * 1024), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()


def fetch_json(url: str, timeout: int = 120) -> dict:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.load(response)


def canonical_manifest_path(raw_dir: Path, manifest_rel: str) -> Path:
    relative = manifest_rel.strip()
    if relative.startswith("./"):
        relative = relative[2:]
    return raw_dir / relative


def resolve_local_path(raw_dir: Path, manifest_rel: str) -> Optional[Path]:
    primary = canonical_manifest_path(raw_dir, manifest_rel)
    if primary.exists():
        return primary

    relative = manifest_rel.strip()
    if relative.startswith("./"):
        relative = relative[2:]
    parts = Path(relative).parts
    if len(parts) == 2 and parts[0] == "data":
        alternate = raw_dir / parts[1]
        if alternate.exists():
            return alternate
    return None


def parse_checksum_manifest(manifest_path: Path) -> List[Tuple[str, str]]:
    if not manifest_path.exists():
        raise FileNotFoundError(
            f"CSIRO checksum manifest not found: {manifest_path}. "
            "Import the collection metadata from CSIRO DAP before running step 000."
        )

    entries: List[Tuple[str, str]] = []
    with open(manifest_path, "r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            parts = stripped.split(maxsplit=1)
            if len(parts) != 2:
                raise ValueError(
                    f"Invalid checksum manifest line {line_number} in {manifest_path}: {line!r}"
                )
            entries.append((parts[0], parts[1]))
    if not entries:
        raise ValueError(f"Checksum manifest is empty: {manifest_path}")
    return entries


def filter_csiro_data_manifest_entries(
    entries: List[Tuple[str, str]],
) -> Tuple[List[Tuple[str, str]], int]:
    """Keep only CSIRO DAP data members; metadata sidecars are not API-downloadable."""
    data_entries = [
        (checksum, manifest_rel)
        for checksum, manifest_rel in entries
        if manifest_rel.startswith("./data/")
    ]
    skipped = len(entries) - len(data_entries)
    if not data_entries:
        raise ValueError("Checksum manifest contains no ./data/ members for CSIRO DAP retrieval.")
    return data_entries, skipped


def resolve_doi_to_data_endpoint(doi: str) -> Tuple[int, str]:
    metadata = fetch_json(f"{DAP_API_BASE}/{doi}.json")
    data_endpoint = metadata.get("data")
    collection_id = metadata.get("dataCollectionId")
    if not data_endpoint or collection_id is None:
        raise RuntimeError(f"CSIRO DAP metadata for DOI {doi} did not expose a data endpoint.")
    return int(collection_id), str(data_endpoint)


def manifest_rel_to_api_key(manifest_rel: str) -> str:
    relative = manifest_rel.strip()
    if relative.startswith("./"):
        relative = relative[2:]
    return relative


def fetch_file_links(collection_id: int, data_endpoint: str, refresh_index: bool) -> Dict[str, dict]:
    index_path = (
        PROJECT_ROOT
        / "data"
        / "raw"
        / f"csiro_dap_file_index_{collection_id}.json"
    )
    if index_path.exists() and not refresh_index:
        with open(index_path, "r", encoding="utf-8") as handle:
            cached = json.load(handle)
        if cached.get("collection_id") == collection_id and isinstance(cached.get("files"), dict):
            log_message(
                f"Using cached CSIRO file index for collection {collection_id}: "
                f"{cached.get('n_files', len(cached['files']))} files",
                "INFO",
            )
            return cached["files"]

    log_message(
        f"Fetching CSIRO DAP file index for collection {collection_id} (single API request)...",
        "INFO",
    )
    page_url = f"{data_endpoint}.json?offset=0&limit={FILE_PAGE_SIZE}"
    page = fetch_json(page_url, timeout=DOWNLOAD_TIMEOUT_S)
    page_files = page.get("file") or []
    if not page_files:
        raise RuntimeError(f"CSIRO DAP returned no files for collection {collection_id}.")

    files: Dict[str, dict] = {}
    for item in page_files:
        filename = item.get("filename")
        if not filename:
            raise RuntimeError(
                f"CSIRO DAP file listing for collection {collection_id} "
                "returned an entry without filename."
            )
        files[manifest_rel_to_api_key(f"./data/{filename}")] = item

    index_path.parent.mkdir(parents=True, exist_ok=True)
    with open(index_path, "w", encoding="utf-8") as handle:
        json.dump(
            {
                "collection_id": collection_id,
                "generated": datetime.now().isoformat(),
                "n_files": len(files),
                "files": files,
            },
            handle,
            indent=2,
            cls=NpEncoder,
        )
    log_message(
        f"Indexed {len(files)} CSIRO DAP files for collection {collection_id}.",
        "DATA",
    )
    return files


def download_bytes(
    url: str,
    output_path: Path,
    timeout: int = DOWNLOAD_TIMEOUT_S,
    log_progress: bool = True,
) -> int:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        total_size = int(response.headers.get("Content-Length", 0))
        downloaded = 0
        chunk_size = 1024 * 1024
        with open(output_path, "wb") as handle:
            while True:
                chunk = response.read(chunk_size)
                if not chunk:
                    break
                handle.write(chunk)
                downloaded += len(chunk)
                if (
                    log_progress
                    and total_size > 0
                    and downloaded % (10 * chunk_size) == 0
                ):
                    percent = 100 * downloaded / total_size
                    log_message(
                        f"  {output_path.name}: {percent:.1f}% "
                        f"({downloaded / (1024 * 1024):.1f} / "
                        f"{total_size / (1024 * 1024):.1f} MB)",
                        "INFO",
                    )

    if output_path.stat().st_size == 0:
        output_path.unlink(missing_ok=True)
        raise RuntimeError(f"Downloaded file is empty: {output_path}")
    return downloaded


def classify_manifest_file(
    collection: PPTACollection,
    expected_sha256: str,
    manifest_rel: str,
) -> Dict[str, object]:
    canonical = canonical_manifest_path(collection.raw_dir, manifest_rel)
    local_path = resolve_local_path(collection.raw_dir, manifest_rel)
    is_dynspec = manifest_rel.endswith(".dynspec")
    if local_path is None:
        return {
            "manifest_rel": manifest_rel,
            "state": "missing",
            "canonical_path": str(canonical.relative_to(PROJECT_ROOT)),
            "is_dynspec": is_dynspec,
        }

    if local_path.stat().st_size == 0:
        return {
            "manifest_rel": manifest_rel,
            "state": "corrupted",
            "reason": "empty_file",
            "expected_sha256": expected_sha256,
            "local_path": str(local_path.relative_to(PROJECT_ROOT)),
            "canonical_path": str(canonical.relative_to(PROJECT_ROOT)),
            "is_dynspec": is_dynspec,
        }

    observed = compute_sha256(local_path)
    if observed != expected_sha256:
        return {
            "manifest_rel": manifest_rel,
            "state": "corrupted",
            "expected_sha256": expected_sha256,
            "observed_sha256": observed,
            "local_path": str(local_path.relative_to(PROJECT_ROOT)),
            "canonical_path": str(canonical.relative_to(PROJECT_ROOT)),
            "is_dynspec": is_dynspec,
        }

    if local_path != canonical:
        canonical.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(local_path), str(canonical))
    return {
        "manifest_rel": manifest_rel,
        "state": "verified",
        "canonical_path": str(canonical.relative_to(PROJECT_ROOT)),
        "is_dynspec": is_dynspec,
    }


def classify_manifest_entries(
    collection: PPTACollection,
    manifest_entries: List[Tuple[str, str]],
    workers: int,
) -> List[Dict[str, object]]:
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [
            executor.submit(classify_manifest_file, collection, expected_sha256, manifest_rel)
            for expected_sha256, manifest_rel in manifest_entries
        ]
        return [future.result() for future in futures]


def summarize_manifest_states(classifications: List[Dict[str, object]]) -> Dict[str, int]:
    counts = {
        "verified": 0,
        "missing": 0,
        "corrupted": 0,
        "downloaded": 0,
        "dynspec": 0,
    }
    for item in classifications:
        state = str(item["state"])
        if state in counts:
            counts[state] += 1
        if item.get("downloaded_in_pass"):
            counts["downloaded"] += 1
        if item.get("is_dynspec"):
            counts["dynspec"] += 1
    return counts


def write_collection_ingestion_state(
    collection: PPTACollection,
    collection_id: int,
    classifications: List[Dict[str, object]],
    counts: Dict[str, int],
) -> Path:
    state_path = collection.raw_dir / "metadata" / "ingestion_state.json"
    state = {
        "generated": datetime.now().isoformat(),
        "pulsar": collection.pulsar,
        "doi": collection.doi,
        "collection_id": collection_id,
        "checksum_manifest": str(collection.checksum_manifest.relative_to(PROJECT_ROOT)),
        "manifest_fingerprint_sha256": compute_sha256(collection.checksum_manifest),
        "counts": counts,
        "files": classifications,
    }
    state_path.parent.mkdir(parents=True, exist_ok=True)
    with open(state_path, "w", encoding="utf-8") as handle:
        json.dump(state, handle, indent=2, cls=NpEncoder)
    return state_path


def remove_corrupted_paths(classification: Dict[str, object]) -> None:
    for key in ("local_path", "canonical_path"):
        path_value = classification.get(key)
        if not path_value:
            continue
        path = PROJECT_ROOT / str(path_value)
        if path.exists():
            path.unlink()


def download_manifest_file(
    collection: PPTACollection,
    expected_sha256: str,
    manifest_rel: str,
    file_links: Dict[str, dict],
) -> Dict[str, object]:
    canonical = canonical_manifest_path(collection.raw_dir, manifest_rel)
    api_key = manifest_rel_to_api_key(manifest_rel)
    file_entry = file_links.get(api_key)
    if file_entry is None:
        raise RuntimeError(
            f"{collection.pulsar}: CSIRO DAP file index has no entry for manifest path {manifest_rel}."
        )

    download_url = file_entry.get("presignedLink", {}).get("href") or file_entry.get("link", {}).get("href")
    if not download_url:
        raise RuntimeError(
            f"{collection.pulsar}: CSIRO DAP entry for {manifest_rel} did not expose a download URL."
        )

    expected_size = file_entry.get("fileSize")
    log_message(
        f"Downloading {collection.pulsar}: {canonical.name} "
        f"({expected_size / (1024 * 1024):.2f} MB expected)"
        if isinstance(expected_size, int)
        else f"Downloading {collection.pulsar}: {canonical.name}",
        "INFO",
    )
    bytes_written = download_bytes(download_url, canonical, log_progress=False)
    observed = compute_sha256(canonical)
    if observed != expected_sha256:
        canonical.unlink(missing_ok=True)
        raise RuntimeError(
            f"{collection.pulsar}: checksum mismatch after download for {manifest_rel}. "
            f"expected {expected_sha256}, observed {observed}."
        )
    return {
        "manifest_rel": manifest_rel,
        "state": "verified",
        "bytes": bytes_written,
        "canonical_path": str(canonical.relative_to(PROJECT_ROOT)),
        "is_dynspec": manifest_rel.endswith(".dynspec"),
        "downloaded_in_pass": True,
    }


def report_manifest_progress(
    collection: PPTACollection,
    completed: int,
    total: int,
    counts: Dict[str, int],
    phase: str,
) -> None:
    log_message(
        f"{collection.pulsar}: {phase} {completed}/{total} "
        f"(verified={counts['verified']}, missing={counts['missing']}, "
        f"corrupted={counts['corrupted']}, downloaded={counts['downloaded']}, "
        f"dynspec={counts['dynspec']})",
        "INFO",
    )


def ingest_ppta_collection(
    collection: PPTACollection,
    force: bool,
    skip_download: bool,
    refresh_index: bool,
    workers: int,
) -> Dict[str, object]:
    log_message("=" * 70, "INFO")
    log_message(f"CSIRO DAP ingestion: {collection.pulsar} ({collection.doi})", "TITLE")
    log_message("=" * 70, "INFO")

    manifest_entries, skipped_metadata = filter_csiro_data_manifest_entries(
        parse_checksum_manifest(collection.checksum_manifest)
    )
    total_entries = len(manifest_entries)
    if skipped_metadata:
        log_message(
            f"{collection.pulsar}: manifest includes {skipped_metadata} non-API metadata "
            "sidecar entries; enforcing CSIRO DAP ./data/ members only.",
            "INFO",
        )

    collection_id, data_endpoint = resolve_doi_to_data_endpoint(collection.doi)
    log_message(
        f"{collection.pulsar}: classifying {total_entries} manifest files with {workers} workers.",
        "INFO",
    )
    classifications = classify_manifest_entries(collection, manifest_entries, workers)
    counts = summarize_manifest_states(classifications)
    log_message(
        f"{collection.pulsar}: inventory verified={counts['verified']}, "
        f"missing={counts['missing']}, corrupted={counts['corrupted']}.",
        "INFO",
    )

    if skip_download:
        if counts["missing"] or counts["corrupted"]:
            raise RuntimeError(
                f"{collection.pulsar}: verify-only mode found "
                f"{counts['missing']} missing and {counts['corrupted']} corrupted manifest files."
            )
        state_path = write_collection_ingestion_state(
            collection, collection_id, classifications, counts
        )
        log_message(f"{collection.pulsar}: ingestion state saved: {state_path}", "INFO")
        return {
            "pulsar": collection.pulsar,
            "doi": collection.doi,
            "collection_id": collection_id,
            "raw_dir": str(collection.raw_dir.relative_to(PROJECT_ROOT)),
            "checksum_manifest": str(collection.checksum_manifest.relative_to(PROJECT_ROOT)),
            "manifest_entries": total_entries,
            "dynspec_files": counts["dynspec"],
            "verified_files": counts["verified"],
            "missing_files": counts["missing"],
            "corrupted_files": counts["corrupted"],
            "downloaded_files": 0,
            "cached_files": counts["verified"],
            "status": "verified",
        }

    if force:
        download_targets = classifications
        log_message(
            f"{collection.pulsar}: --force enabled; queueing all {len(download_targets)} manifest files.",
            "WARNING",
        )
    else:
        download_targets = [
            item for item in classifications if item["state"] in {"missing", "corrupted"}
        ]
        if not download_targets:
            log_message(
                f"{collection.pulsar}: all manifest files already verified; no downloads required.",
                "SUCCESS",
            )
        else:
            log_message(
                f"{collection.pulsar}: queueing {len(download_targets)} manifest files "
                f"({counts['missing']} missing, {counts['corrupted']} corrupted).",
                "INFO",
            )

    file_links = fetch_file_links(collection_id, data_endpoint, refresh_index=refresh_index)
    manifest_lookup = {
        manifest_rel: expected_sha256 for expected_sha256, manifest_rel in manifest_entries
    }
    completed = 0
    total_targets = len(download_targets)
    if download_targets:
        for item in download_targets:
            if item["state"] == "corrupted":
                remove_corrupted_paths(item)

        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(
                    download_manifest_file,
                    collection,
                    manifest_lookup[str(item["manifest_rel"])],
                    str(item["manifest_rel"]),
                    file_links,
                ): str(item["manifest_rel"])
                for item in download_targets
            }
            for future in as_completed(futures):
                result = future.result()
                manifest_rel = str(result["manifest_rel"])
                for index, item in enumerate(classifications):
                    if item["manifest_rel"] == manifest_rel:
                        classifications[index] = result
                        break
                counts = summarize_manifest_states(classifications)
                completed += 1
                if completed == 1 or completed == total_targets or completed % PROGRESS_INTERVAL == 0:
                    report_manifest_progress(
                        collection, completed, total_targets, counts, phase="download"
                    )

    classifications = classify_manifest_entries(collection, manifest_entries, workers)
    counts = summarize_manifest_states(classifications)
    if counts["missing"] or counts["corrupted"]:
        raise RuntimeError(
            f"{collection.pulsar}: ingestion incomplete after download pass "
            f"(missing={counts['missing']}, corrupted={counts['corrupted']})."
        )

    state_path = write_collection_ingestion_state(
        collection, collection_id, classifications, counts
    )
    summary = {
        "pulsar": collection.pulsar,
        "doi": collection.doi,
        "collection_id": collection_id,
        "raw_dir": str(collection.raw_dir.relative_to(PROJECT_ROOT)),
        "checksum_manifest": str(collection.checksum_manifest.relative_to(PROJECT_ROOT)),
        "manifest_entries": total_entries,
        "dynspec_files": counts["dynspec"],
        "verified_files": counts["verified"],
        "missing_files": counts["missing"],
        "corrupted_files": counts["corrupted"],
        "downloaded_files": counts["downloaded"],
        "cached_files": counts["verified"] - counts["downloaded"],
        "status": "verified",
        "ingestion_state": str(state_path.relative_to(PROJECT_ROOT)),
    }
    log_message(
        f"{collection.pulsar}: verified {summary['verified_files']} manifest files "
        f"({counts['dynspec']} dynamic spectra; downloaded={counts['downloaded']}, "
        f"cached={counts['verified']}, missing={counts['missing']}, "
        f"corrupted={counts['corrupted']}).",
        "SUCCESS",
    )
    return summary


def download_scintools_file(filename: str, force: bool) -> Dict[str, object]:
    RAW_DIR_SCINTOOLS.mkdir(parents=True, exist_ok=True)
    output_path = RAW_DIR_SCINTOOLS / filename
    if output_path.exists() and not force:
        return {
            "filename": filename,
            "path": str(output_path.relative_to(PROJECT_ROOT)),
            "size_bytes": output_path.stat().st_size,
            "status": "cached",
            "sha256": compute_sha256(output_path),
            "source": "scintools",
        }

    url = f"{SCINTOOLS_BASE}{filename}"
    try:
        download_bytes(url, output_path)
    except (urllib.error.URLError, urllib.error.HTTPError, RuntimeError) as exc:
        return {
            "filename": filename,
            "status": "failed",
            "source": "scintools",
            "error": str(exc),
        }

    return {
        "filename": filename,
        "path": str(output_path.relative_to(PROJECT_ROOT)),
        "size_bytes": output_path.stat().st_size,
        "status": "downloaded",
        "sha256": compute_sha256(output_path),
        "source": "scintools",
        "url": url,
    }


def ingest_scintools_supplement(force: bool, workers: int) -> List[Dict[str, object]]:
    log_message("Supplementary Scintools ATNF epochs for J0437-4715", "INFO")
    RAW_DIR_SCINTOOLS.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=min(workers, len(SCINTOOLS_TARGET_FILES))) as executor:
        futures = [
            executor.submit(download_scintools_file, filename, force)
            for filename in SCINTOOLS_TARGET_FILES
        ]
        return [future.result() for future in as_completed(futures)]


def audit_csiro_collection(
    collection: PPTACollection,
    verify_checksums: bool,
    workers: int,
) -> Dict[str, object]:
    manifest_entries, skipped_metadata = filter_csiro_data_manifest_entries(
        parse_checksum_manifest(collection.checksum_manifest)
    )
    dynspec_entries = [
        (checksum, manifest_rel)
        for checksum, manifest_rel in manifest_entries
        if manifest_rel.endswith(".dynspec")
    ]

    classifications = classify_manifest_entries(collection, dynspec_entries, workers)
    counts = summarize_manifest_states(classifications)
    sample_missing = [
        str(item["manifest_rel"])
        for item in classifications
        if item["state"] == "missing"
    ][:10]
    sample_mismatch = [
        str(item["manifest_rel"])
        for item in classifications
        if item["state"] == "corrupted"
    ][:10]

    if not verify_checksums:
        present_dynspec = counts["verified"] + counts["corrupted"]
        complete = counts["missing"] == 0 and present_dynspec == len(dynspec_entries)
    else:
        complete = counts["missing"] == 0 and counts["corrupted"] == 0 and counts["verified"] == len(dynspec_entries)

    return {
        "pulsar": collection.pulsar,
        "doi": collection.doi,
        "source": "CSIRO DAP",
        "raw_dir": str(collection.raw_dir.relative_to(PROJECT_ROOT)),
        "manifest_data_members": len(manifest_entries),
        "skipped_metadata_sidecars": skipped_metadata,
        "expected_dynspec": len(dynspec_entries),
        "present_dynspec": counts["verified"] + counts["corrupted"],
        "verified_dynspec": counts["verified"],
        "missing_dynspec": counts["missing"],
        "checksum_mismatch_dynspec": counts["corrupted"],
        "complete": complete,
        "sample_missing": sample_missing,
        "sample_checksum_mismatch": sample_mismatch,
    }


def audit_jiamusi_archive() -> Dict[str, object]:
    from scripts.steps.step_029_jiamusi_analysis import JIAMUSI_PULSARS

    expected_files: List[str] = []
    for urls in JIAMUSI_PULSARS.values():
        expected_files.extend(Path(url).name for url in urls)

    raw_dir = PROJECT_ROOT / "data" / "raw" / "jiamusi"
    present = sorted(path.name for path in raw_dir.glob("*.dat")) if raw_dir.exists() else []
    missing = sorted(set(expected_files) - set(present))
    return {
        "source": "Jiamusi 66m archive",
        "raw_dir": str(raw_dir.relative_to(PROJECT_ROOT)),
        "expected_epochs": len(expected_files),
        "present_epochs": len(present),
        "missing_epochs": len(missing),
        "complete": not missing,
        "sample_missing": missing[:10],
    }


def audit_optional_archives() -> List[Dict[str, object]]:
    audits: List[Dict[str, object]] = []
    j0613_dir = PROJECT_ROOT / "data" / "raw" / "j0613"
    j0613_files = list(j0613_dir.rglob("*.dynspec")) if j0613_dir.exists() else []
    audits.append(
        {
            "pulsar": "J0613-0200",
            "source": "CSIRO DAP (manual acquisition)",
            "raw_dir": str(j0613_dir.relative_to(PROJECT_ROOT)),
            "present_dynspec": len(j0613_files),
            "complete": bool(j0613_files),
            "note": "Control pulsar raw data is optional until acquired into data/raw/j0613/.",
        }
    )

    meerkat_dir = PROJECT_ROOT / "data" / "raw" / "meerkat" / "data" / "pdfb4"
    meerkat_files = list(meerkat_dir.glob("*.dynspec")) if meerkat_dir.exists() else []
    audits.append(
        {
            "pulsar": "MeerKAT archive (non-J0437 inventory)",
            "source": "Repository inventory",
            "raw_dir": str(meerkat_dir.relative_to(PROJECT_ROOT)),
            "present_dynspec": len(meerkat_files),
            "complete": True,
            "note": "MeerKAT step inventories non-J0437 files; J0437 MeerKAT data is not in this repository.",
        }
    )
    return audits


def write_ingestion_audit(
    verify_checksums: bool,
    workers: int,
) -> Dict[str, object]:
    ppta_audits = [
        audit_csiro_collection(collection, verify_checksums=verify_checksums, workers=workers)
        for collection in PPTA_COLLECTIONS
    ]
    audit = {
        "generated": datetime.now().isoformat(),
        "verify_checksums": verify_checksums,
        "ppta_collections": ppta_audits,
        "jiamusi": audit_jiamusi_archive(),
        "optional_archives": audit_optional_archives(),
    }
    audit["all_required_complete"] = all(
        item["complete"] for item in ppta_audits
    ) and audit["jiamusi"]["complete"]

    INGESTION_AUDIT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(INGESTION_AUDIT_PATH, "w", encoding="utf-8") as handle:
        json.dump(audit, handle, indent=2, cls=NpEncoder)

    for item in ppta_audits:
        log_message(
            f"AUDIT {item['pulsar']}: {item['present_dynspec']}/{item['expected_dynspec']} "
            f"dynspec present; verified={item['verified_dynspec']}; "
            f"missing={item['missing_dynspec']}; mismatch={item['checksum_mismatch_dynspec']}",
            "SUCCESS" if item["complete"] else "WARNING",
        )
    jiamusi = audit["jiamusi"]
    log_message(
        f"AUDIT Jiamusi: {jiamusi['present_epochs']}/{jiamusi['expected_epochs']} epochs present",
        "SUCCESS" if jiamusi["complete"] else "WARNING",
    )
    log_message(f"Audit saved: {INGESTION_AUDIT_PATH}", "INFO")
    return audit


def write_ingestion_catalog(
    ppta_summaries: Iterable[Dict[str, object]],
    scintools_files: Iterable[Dict[str, object]],
) -> Dict[str, object]:
    scintools_files = list(scintools_files)
    catalog = {
        "catalog_metadata": {
            "generated": datetime.now().isoformat(),
            "pipeline_version": "0.3",
            "checksum_algorithm": "SHA256",
        },
        "ppta_collections": list(ppta_summaries),
        "scintools_supplement": {
            "n_target_files": len(SCINTOOLS_TARGET_FILES),
            "n_available": len(
                [item for item in scintools_files if item.get("status") in {"cached", "downloaded"}]
            ),
            "files": scintools_files,
        },
    }
    CATALOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(CATALOG_PATH, "w", encoding="utf-8") as handle:
        json.dump(catalog, handle, indent=2, cls=NpEncoder)
    return catalog


def main() -> bool:
    parser = argparse.ArgumentParser(
        description="Acquire and verify PPTA DR2 dynamic spectra for TEP-J0437."
    )
    parser.add_argument(
        "--force",
        "-f",
        action="store_true",
        help="Re-download PPTA and Scintools files even when local copies exist.",
    )
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="Verify existing PPTA files against CSIRO manifests without downloading.",
    )
    parser.add_argument(
        "--skip-scintools",
        action="store_true",
        help="Skip supplementary Scintools ATNF downloads.",
    )
    parser.add_argument(
        "--refresh-index",
        action="store_true",
        help="Refresh cached CSIRO DAP file index JSON.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=DEFAULT_WORKERS,
        help=f"Parallel workers for manifest verification/downloads (default: {DEFAULT_WORKERS}).",
    )
    parser.add_argument(
        "--audit-only",
        action="store_true",
        help="Audit pulsar raw-data coverage without downloading or verifying the full manifest.",
    )
    parser.add_argument(
        "--audit-full",
        action="store_true",
        help="During --audit-only, verify SHA-256 checksums for all present dynspec files.",
    )
    args = parser.parse_args()
    if args.workers < 1:
        raise ValueError("--workers must be at least 1.")

    if args.audit_only:
        write_ingestion_audit(verify_checksums=args.audit_full, workers=args.workers)
        return True

    print_status("=" * 70, "TITLE")
    print_status("STEP 000: DATA INGESTION", "TITLE")
    print_status("=" * 70, "TITLE")

    ppta_summaries: List[Dict[str, object]] = []
    for collection in PPTA_COLLECTIONS:
        ppta_summaries.append(
            ingest_ppta_collection(
                collection,
                force=args.force,
                skip_download=args.skip_download,
                refresh_index=args.refresh_index,
                workers=args.workers,
            )
        )

    scintools_files: List[Dict[str, object]] = []
    if not args.skip_scintools:
        scintools_files = ingest_scintools_supplement(force=args.force, workers=args.workers)
        failed = [item for item in scintools_files if item.get("status") == "failed"]
        if failed:
            print_status(
                f"WARNING: Scintools supplementary download failed for {len(failed)} files; "
                "continuing with primary PPTA data only.",
                "WARNING",
            )

    catalog = write_ingestion_catalog(ppta_summaries, scintools_files)
    write_ingestion_audit(verify_checksums=True, workers=args.workers)

    print_status("\n" + "=" * 70, "TITLE")
    print_status("DATA INGESTION SUMMARY", "TITLE")
    print_status("=" * 70, "TITLE")
    for summary in ppta_summaries:
        print_status(
            f"{summary['pulsar']}: {summary['verified_files']} manifest files verified "
            f"({summary['dynspec_files']} dynamic spectra)",
            "SUCCESS",
        )
    if not args.skip_scintools:
        available = catalog["scintools_supplement"]["n_available"]
        print_status(
            f"Scintools supplement: {available}/{len(SCINTOOLS_TARGET_FILES)} files available",
            "SUCCESS" if available > 0 else "WARNING",
        )
    print_status(f"Catalog saved: {CATALOG_PATH}", "INFO")
    return True


def step_main(logger=None, verbose=True):
    if logger:
        set_step_logger(logger)
    return main()


if __name__ == "__main__":
    try:
        success = main()
    except Exception as exc:
        print_status(f"STEP 000 failed: {exc}", "ERROR")
        success = False
    sys.exit(0 if success else 1)
