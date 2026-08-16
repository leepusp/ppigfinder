#!/usr/bin/env python3
"""
AlphaFold result parsing service.

Provides cached and parallel parsing of AF3 result directories.
"""

from __future__ import annotations

from pathlib import Path

from ppigfinder.alphafold.result_loader import find_af3_job_dirs, load_directories_parallel
from ppigfinder.alphafold.result_parser import parse_af3_result_directory, parse_af3_results_root


class AlphaFoldResultsService:
    """
    High-level service for AF3 result parsing.
    """

    def parse_directory(self, path: str | Path) -> dict:
        return parse_af3_result_directory(path)

    def parse_root_sequential(self, root: str | Path) -> list[dict]:
        return parse_af3_results_root(root)

    def parse_root_parallel(
        self,
        root: str | Path,
        workers: int | None = None,
        cache_dir: str | Path | None = None,
    ) -> list[dict]:
        root = Path(root)

        candidates = find_af3_job_dirs(root)

        if not candidates:
            return parse_af3_results_root(root)

        return load_directories_parallel(
            candidates,
            parser=parse_af3_result_directory,
            workers=workers,
            cache_dir=cache_dir,
        )
