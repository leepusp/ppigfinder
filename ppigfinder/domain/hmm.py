#!/usr/bin/env python3
"""
Domain models for HMM/domain search results.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(slots=True)
class HMMHit:
    """
    Conserved domain hit detected by HMMER or fallback scanner.
    """

    orf_index: int
    domain: str
    start: int | None = None
    end: int | None = None
    score: float | None = None
    evalue: float | None = None
    source: str = "unknown"


def hmm_hit_from_legacy_dict(data: dict, orf_index: int = 0) -> HMMHit:
    """
    Convert legacy domain dictionary into HMMHit.
    """
    return HMMHit(
        orf_index=int(data.get("orf_index", orf_index)),
        domain=data.get("domain") or data.get("name") or "unknown",
        start=data.get("start"),
        end=data.get("end"),
        score=data.get("score"),
        evalue=data.get("evalue"),
        source=data.get("source", "unknown"),
    )
