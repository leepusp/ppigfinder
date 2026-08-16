#!/usr/bin/env python3
"""
Small JSON cache utilities for ppigFinder.

Used to avoid repeated parsing of expensive result folders.
"""

from __future__ import annotations

from pathlib import Path
import hashlib
import json
import os
import time


def file_signature(path: str | Path) -> dict:
    """
    Return a stable signature for one file based on size and mtime.
    """
    path = Path(path)
    stat = path.stat()

    return {
        "path": str(path),
        "size": stat.st_size,
        "mtime": stat.st_mtime,
    }


def directory_signature(
    path: str | Path,
    suffixes: tuple[str, ...] = (".json", ".cif", ".csv", ".tsv"),
) -> str:
    """
    Compute a lightweight signature for relevant files in a directory.
    """
    path = Path(path)
    payload = []

    for root, _, files in os.walk(path):
        for name in sorted(files):
            file_path = Path(root) / name

            if suffixes and file_path.suffix.lower() not in suffixes:
                continue

            try:
                rel = file_path.relative_to(path)
                stat = file_path.stat()
            except OSError:
                continue

            payload.append((str(rel), stat.st_size, stat.st_mtime))

    raw = json.dumps(payload, sort_keys=True).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


class JsonCache:
    """
    Simple JSON cache stored in a project/cache directory.
    """

    def __init__(self, cache_dir: str | Path):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def _path_for_key(self, key: str) -> Path:
        safe = hashlib.sha256(key.encode("utf-8")).hexdigest()
        return self.cache_dir / f"{safe}.json"

    def get(self, key: str):
        path = self._path_for_key(key)

        if not path.exists():
            return None

        try:
            with path.open("r", encoding="utf-8") as handle:
                return json.load(handle)
        except Exception:
            return None

    def set(self, key: str, value) -> None:
        path = self._path_for_key(key)
        payload = {
            "created_at": time.time(),
            "value": value,
        }

        with path.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)

    def get_value(self, key: str):
        payload = self.get(key)

        if not payload:
            return None

        return payload.get("value")
