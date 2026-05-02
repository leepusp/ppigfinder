#!/usr/bin/env python3
"""
Performance defaults for ppigFinder.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(slots=True)
class PerformanceConfig:
    """
    Runtime performance configuration.
    """

    workers: int | None = None
    use_processes_for_cpu_tasks: bool = False
    enable_result_cache: bool = True
    max_cache_age_days: int = 30


DEFAULT_PERFORMANCE = PerformanceConfig()
