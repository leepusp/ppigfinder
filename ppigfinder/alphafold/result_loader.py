#!/usr/bin/env python3
"""
AlphaFold result loading helpers.

This module will progressively receive AF3 result parsing logic currently
stored in legacy_v20.py.
"""

from __future__ import annotations

from pathlib import Path

from ppigfinder.infrastructure.cache import JsonCache, directory_signature
from ppigfinder.infrastructure.parallel import parallel_map


def find_af3_job_dirs(root: str | Path) -> list[Path]:
    """
    Return candidate AF3 job directories.

    A directory is considered a candidate when it contains JSON, CIF or CSV
    files commonly produced by AlphaFold 3.
    """
    root = Path(root)

    if not root.exists():
        return []

    candidates: list[Path] = []

    for path in root.rglob("*"):
        if not path.is_dir():
            continue

        try:
            names = {item.suffix.lower() for item in path.iterdir() if item.is_file()}
        except OSError:
            continue

        if names & {".json", ".cif", ".csv"}:
            candidates.append(path)

    return candidates


def cached_directory_signature(path: str | Path) -> str:
    """
    Return a directory signature suitable for cache keys.
    """
    return directory_signature(path, suffixes=(".json", ".cif", ".csv", ".tsv"))


def load_directories_parallel(
    directories: list[str | Path],
    parser,
    workers: int | None = None,
    cache_dir: str | Path | None = None,
) -> list:
    """
    Parse result directories in parallel.

    parser must be a callable accepting one Path and returning a serializable
    result or domain object.
    """
    dirs = [Path(item) for item in directories]

    if cache_dir is None:
        return parallel_map(parser, dirs, workers=workers, mode="thread")

    cache = JsonCache(cache_dir)

    def parse_with_cache(path: Path):
        signature = cached_directory_signature(path)
        key = f"af3-result:{path}:{signature}"

        cached = cache.get_value(key)
        if cached is not None:
            return cached

        result = parser(path)
        cache.set(key, result)
        return result

    return parallel_map(parse_with_cache, dirs, workers=workers, mode="thread")
