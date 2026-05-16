"""
Parallel worker counts and BLAS/OpenMP threading for multiprocess pipeline steps.

When using ProcessPoolExecutor, each child process must not spawn a full BLAS thread
pool, or machines oversubscribe (slow on Apple Silicon and many-core Linux).

Environment overrides (optional):
  TEP_WORKERS             — positive int: same cap for all steps that read it
  TEP_STEP001_MAX_WORKERS — cap for Step 001 (I/O + parse)
  TEP_STEP002_MAX_WORKERS — cap for Step 002 (large FFTs, memory-sensitive)
  TEP_STEP003_MAX_WORKERS — cap for Step 003 (epoch parallel closure)
  TEP_BLAS_THREADS        — BLAS threads per process (default 1 under multiprocessing)
"""

from __future__ import annotations

import os
from typing import Literal, Optional

Role = Literal["io_bound", "fft_heavy", "cpu_bound"]

_BLAS_ENV_KEYS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
)


def _parse_positive_int(name: str, raw: str) -> int:
    v = raw.strip()
    if not v:
        raise ValueError(f"{name} is empty")
    n = int(v)
    if n < 1:
        raise ValueError(f"{name} must be >= 1, got {n}")
    return n


def _env_optional_positive(name: str) -> Optional[int]:
    raw = os.environ.get(name)
    if raw is None:
        return None
    s = str(raw).strip()
    if not s:
        return None
    return _parse_positive_int(name, s)


def logical_cpu_count() -> int:
    n = os.cpu_count()
    if n is None or n < 1:
        raise RuntimeError("os.cpu_count() returned invalid value; set TEP_WORKERS explicitly")
    return n


def configure_blas_thread_env() -> None:
    """Set BLAS/OpenMP thread counts for the current process (call before spawning pools)."""
    raw = os.environ.get("TEP_BLAS_THREADS")
    if raw is None or not str(raw).strip():
        threads = 1
    else:
        threads = _parse_positive_int("TEP_BLAS_THREADS", str(raw))
    ts = str(threads)
    for key in _BLAS_ENV_KEYS:
        os.environ[key] = ts


def worker_count(
    *,
    role: Role,
    reserve: int = 2,
) -> int:
    """
    Choose ProcessPoolExecutor max_workers.

    reserve: leave this many logical CPUs unclaimed for the OS / main process.
    """
    override = _env_optional_positive("TEP_WORKERS")
    cpu = logical_cpu_count()
    spare = max(1, cpu - max(0, reserve))

    if override is not None:
        return min(override, spare)

    if role == "fft_heavy":
        cap_env = _env_optional_positive("TEP_STEP002_MAX_WORKERS")
        default_cap = min(spare, 8)
        cap = cap_env if cap_env is not None else default_cap
        return max(1, min(cap, spare))

    if role == "io_bound":
        cap_env = _env_optional_positive("TEP_STEP001_MAX_WORKERS")
        cap = cap_env if cap_env is not None else spare
        return max(1, min(cap, spare))

    # cpu_bound (Step 003)
    cap_env = _env_optional_positive("TEP_STEP003_MAX_WORKERS")
    cap = cap_env if cap_env is not None else spare
    return max(1, min(cap, spare))
