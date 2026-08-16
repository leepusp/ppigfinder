#!/usr/bin/env python3
"""
Domain models for BLAST-like search results.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class BlastHit:
    """
    Protein-vs-ORF similarity hit.
    """

    orf_index: int
    score: float
    identity: float
    positives: float | None = None
    alignment_length: int | None = None
    evalue: float | None = None
    query_start: int | None = None
    query_end: int | None = None
    subject_start: int | None = None
    subject_end: int | None = None
    method: str = "unknown"


def blast_hit_from_legacy_dict(data: dict) -> BlastHit:
    """
    Convert a legacy BLAST hit dictionary into a BlastHit dataclass.

    This adapter is intentionally permissive because older ppigFinder
    versions used slightly different key names.
    """
    return BlastHit(
        orf_index=int(data.get("orf_index", data.get("index", data.get("subject_index", 0)))),
        score=float(data.get("score", data.get("bitscore", 0.0))),
        identity=float(data.get("identity", data.get("pident", 0.0))),
        positives=data.get("positives"),
        alignment_length=data.get("length"),
        evalue=data.get("evalue"),
        query_start=data.get("qstart"),
        query_end=data.get("qend"),
        subject_start=data.get("sstart"),
        subject_end=data.get("send"),
        method=data.get("method", "unknown"),
    )
