#!/usr/bin/env python3
"""
Built-in PSSM/domain scanner infrastructure.

This module is intentionally lightweight for now. It will receive the domain
fallback logic currently embedded in legacy_v20.py.
"""

from __future__ import annotations

from dataclasses import dataclass
from ppigfinder.infrastructure.parallel import parallel_map


@dataclass(slots=True)
class PSSMScanTask:
    """
    One protein sequence to scan.
    """

    orf_index: int
    protein: str


@dataclass(slots=True)
class PSSMScanHit:
    """
    One built-in scanner hit.
    """

    orf_index: int
    domain: str
    start: int | None = None
    end: int | None = None
    score: float | None = None


def scan_one_protein(task: PSSMScanTask, scanner) -> list[PSSMScanHit]:
    """
    Run a scanner callable against one protein.
    """
    result = scanner(task.protein)
    hits: list[PSSMScanHit] = []

    for item in result or []:
        hits.append(
            PSSMScanHit(
                orf_index=task.orf_index,
                domain=item.get("domain", item.get("name", "unknown")),
                start=item.get("start"),
                end=item.get("end"),
                score=item.get("score"),
            )
        )

    return hits


def scan_proteins_parallel(
    proteins: list[str],
    scanner,
    workers: int | None = None,
) -> list[PSSMScanHit]:
    """
    Scan many proteins in parallel using the provided scanner callable.
    """
    tasks = [
        PSSMScanTask(orf_index=index, protein=protein)
        for index, protein in enumerate(proteins)
    ]

    nested = parallel_map(
        lambda task: scan_one_protein(task, scanner),
        tasks,
        workers=workers,
        mode="thread",
    )

    hits: list[PSSMScanHit] = []
    for group in nested:
        hits.extend(group)

    return hits
