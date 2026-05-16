#!/usr/bin/env python3
"""
Shared Data Loading Utility with Caching

This module provides cached data loading functions to eliminate redundant
file reads across pipeline steps, improving I/O efficiency and reducing
memory pressure.

Usage:
    from scripts.utils.data_loader import load_closure_data, load_per_epoch_data
    
    # First call loads from file
    data1 = load_closure_data()
    
    # Subsequent calls return cached data
    data2 = load_closure_data()  # Returns cached data
"""

import json
import functools
from pathlib import Path
from typing import Dict, Any, Optional, List
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"


def _triplet_delay_us(triplet: Dict[str, Any]) -> Optional[float]:
    """Stokes-aligned delay in microseconds; None if neither field is present."""
    g = triplet.get("geometric_delta_us")
    if g is not None:
        return float(g)
    d = triplet.get("delta_us")
    if d is not None:
        return float(d)
    return None


def _j0437_summary_path() -> Path:
    j0437 = RESULTS_DIR / "step_003_closure_final_summary_j0437.json"
    if j0437.exists():
        return j0437
    return RESULTS_DIR / "step_003_closure_final_summary.json"


def _j0437_per_epoch_path() -> Path:
    j0437 = RESULTS_DIR / "step_003_closure_final_per_epoch_j0437.json"
    if j0437.exists():
        return j0437
    return RESULTS_DIR / "step_003_closure_final_per_epoch.json"


@functools.lru_cache(maxsize=1)
def load_closure_data() -> Optional[Dict[str, Any]]:
    """
    Load J0437 closure delay summary data with caching.

    Prefers step_003_closure_final_summary_j0437.json (canonical for this paper).
    """
    results_file = _j0437_summary_path()

    if not results_file.exists():
        return None

    with open(results_file, "r") as f:
        return json.load(f)


@functools.lru_cache(maxsize=1)
def load_per_epoch_data() -> Optional[list]:
    """
    Load J0437 per-epoch closure delay data with caching.

    Prefers step_003_closure_final_per_epoch_j0437.json.
    """
    results_file = _j0437_per_epoch_path()

    if not results_file.exists():
        return None

    with open(results_file, "r") as f:
        return json.load(f)


def clear_cache():
    """Clear all cached data from memory."""
    load_closure_data.cache_clear()
    load_per_epoch_data.cache_clear()


def get_all_triplet_delays() -> Optional[np.ndarray]:
    """
    Extract all triplet delays from per-epoch data.
    
    Returns:
        NumPy array of all triplet delays in microseconds, or None if data not available
    """
    per_epoch_data = load_per_epoch_data()
    
    if per_epoch_data is None:
        return None
    
    all_delays: List[float] = []
    for epoch in per_epoch_data:
        for triplet in epoch.get("triplets", []):
            v = _triplet_delay_us(triplet)
            if v is not None:
                all_delays.append(v)
    
    return np.array(all_delays, dtype=float)


def get_epoch_means() -> Optional[np.ndarray]:
    """
    Extract epoch-level mean delays from per-epoch data.
    
    Returns:
        NumPy array of epoch mean delays in microseconds, or None if data not available
    """
    per_epoch_data = load_per_epoch_data()
    
    if per_epoch_data is None:
        return None
    
    epoch_means = []
    for epoch in per_epoch_data:
        delays = []
        for t in epoch.get("triplets", []):
            v = _triplet_delay_us(t)
            if v is not None:
                delays.append(v)
        if delays:
            epoch_means.append(np.mean(delays))
    
    return np.array(epoch_means, dtype=float)


def get_epoch_statistics() -> Optional[Dict[str, Any]]:
    """
    Get comprehensive epoch-level statistics.
    
    Returns:
        Dictionary with epoch statistics or None if data not available
    """
    per_epoch_data = load_per_epoch_data()
    
    if per_epoch_data is None:
        return None
    
    epoch_stats = []
    for epoch in per_epoch_data:
        delays = []
        for t in epoch.get("triplets", []):
            v = _triplet_delay_us(t)
            if v is not None:
                delays.append(v)
        if delays:
            delays_arr = np.array(delays)
            epoch_stats.append({
                "epoch": epoch.get("epoch", "unknown"),
                "n_triplets": len(delays),
                "n_arclets": epoch.get("n_arclets", 0),
                "mean_delay_us": float(np.mean(delays_arr)),
                "std_delay_us": float(np.std(delays_arr, ddof=1)),
                "median_delay_us": float(np.median(delays_arr)),
                "min_delay_us": float(np.min(delays_arr)),
                "max_delay_us": float(np.max(delays_arr))
            })
    
    return {
        "n_epochs": len(epoch_stats),
        "epochs": epoch_stats
    }
