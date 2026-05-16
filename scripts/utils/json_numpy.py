"""JSON helpers for pipeline outputs that contain NumPy scalars and arrays."""

from __future__ import annotations

import json
from typing import Any

import numpy as np


class NpEncoder(json.JSONEncoder):
    """JSON encoder that handles NumPy scalars (including numpy.bool_) and arrays."""

    def default(self, obj: Any) -> Any:
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.generic):
            return obj.item()
        if isinstance(obj, (bytes, memoryview)):
            return obj.decode("utf-8", errors="replace")
        return super().default(obj)
