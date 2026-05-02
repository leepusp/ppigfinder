#!/usr/bin/env python3
"""
Domain models related to open reading frames.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(slots=True)
class DomainHit:
    """
    Conserved domain annotation assigned to an ORF.
    """

    name: str
    start: int | None = None
    end: int | None = None
    score: float | None = None
    evalue: float | None = None
    source: str = "unknown"


@dataclass(slots=True)
class ORF:
    """
    Open reading frame represented in genomic coordinates.
    """

    id: str
    start: int
    end: int
    strand: str
    frame: int
    dna: str
    protein: str
    length: int
    gc: float
    source: str = "unknown"
    candidate_score: float = 0.0
    domains: list[DomainHit] = field(default_factory=list)
    neighborhood: list[str] = field(default_factory=list)
    gene_name: str | None = None
    putative_function: str | None = None

    @property
    def protein_length(self) -> int:
        return len(self.protein.rstrip("*"))

    @property
    def genomic_length(self) -> int:
        return self.end - self.start


@dataclass(slots=True)
class ORFSet:
    """
    Collection of ORFs from the same genome or sequence record.
    """

    sequence_name: str
    orfs: list[ORF] = field(default_factory=list)

    def sort_by_position(self) -> None:
        self.orfs.sort(key=lambda item: item.start)

    def __len__(self) -> int:
        return len(self.orfs)

    def __iter__(self):
        return iter(self.orfs)
