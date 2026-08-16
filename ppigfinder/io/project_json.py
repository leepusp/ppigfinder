#!/usr/bin/env python3
"""
Project JSON reader/writer for ppigFinder.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any
import json

from ppigfinder.domain.project import ProjectState


def make_json_safe(value: Any) -> Any:
    """
    Convert common non-JSON objects into JSON-safe values.
    """
    if value is None or isinstance(value, (str, int, float, bool)):
        return value

    if isinstance(value, dict):
        return {
            str(make_json_safe(k)): make_json_safe(v)
            for k, v in value.items()
        }

    if isinstance(value, (list, tuple, set)):
        return [make_json_safe(item) for item in value]

    # NumPy scalar support without importing NumPy.
    if hasattr(value, "item") and callable(value.item):
        try:
            return make_json_safe(value.item())
        except Exception:
            pass

    # pathlib.Path and other simple objects.
    if isinstance(value, Path):
        return str(value)

    return str(value)


def read_project_json(path: str | Path) -> ProjectState:
    path = Path(path)

    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)

    if not isinstance(data, dict):
        raise ValueError("Invalid ppigFinder project: top-level JSON must be an object.")

    return ProjectState.from_dict(data)


def write_project_json(path: str | Path, project: ProjectState, indent: int = 2) -> None:
    path = Path(path)
    payload = make_json_safe(project.to_dict())

    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=indent, ensure_ascii=False)
        handle.write("\n")


def validate_project_json(path: str | Path) -> None:
    read_project_json(path)
