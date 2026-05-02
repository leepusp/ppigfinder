#!/usr/bin/env python3
"""
Parallel execution helpers for ppigFinder.

This module provides conservative defaults suitable for desktops and HPC
login/interactive sessions. Heavy jobs submitted to schedulers should still
respect the resources requested from SLURM/PBS/LSF.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from os import cpu_count
from typing import Callable, Iterable, Literal, TypeVar


T = TypeVar("T")
R = TypeVar("R")


def recommended_workers(
    requested: int | None = None,
    max_default: int = 8,
    reserve_cpus: int = 1,
) -> int:
    """
    Return a safe worker count.

    If requested is provided, it is respected but clamped to at least 1.
    Otherwise use available CPUs minus reserve_cpus, limited by max_default.
    """
    if requested is not None:
        return max(1, int(requested))

    detected = cpu_count() or 1
    usable = max(1, detected - reserve_cpus)
    return max(1, min(max_default, usable))


def chunked(items: list[T], chunk_size: int) -> list[list[T]]:
    """
    Split a list into fixed-size chunks.
    """
    chunk_size = max(1, int(chunk_size))
    return [items[i:i + chunk_size] for i in range(0, len(items), chunk_size)]


def parallel_map(
    function: Callable[[T], R],
    items: Iterable[T],
    workers: int | None = None,
    mode: Literal["thread", "process"] = "thread",
) -> list[R]:
    """
    Run function over items in parallel while preserving result order.

    mode='thread' is best for I/O or external command orchestration.
    mode='process' is best for pure Python CPU-heavy algorithms.
    """
    item_list = list(items)

    if not item_list:
        return []

    n_workers = recommended_workers(workers)

    if n_workers <= 1 or len(item_list) == 1:
        return [function(item) for item in item_list]

    executor_cls = ThreadPoolExecutor if mode == "thread" else ProcessPoolExecutor
    results: list[R | None] = [None] * len(item_list)

    with executor_cls(max_workers=n_workers) as executor:
        future_to_index = {
            executor.submit(function, item): index
            for index, item in enumerate(item_list)
        }

        for future in as_completed(future_to_index):
            index = future_to_index[future]
            results[index] = future.result()

    return results  # type: ignore[return-value]
