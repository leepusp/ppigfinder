#!/usr/bin/env python3
"""
BLAST search service for ppigFinder.

This service separates BLAST algorithm selection from the GUI. During the
transition, it can still call methods from the legacy AdvancedORFAnalyzer,
but the GUI no longer needs to know how each backend is selected.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from ppigfinder.infrastructure.backends import BACKENDS


@dataclass(slots=True)
class BlastSearchParams:
    """
    Parameters used for protein-vs-ORF similarity search.
    """

    threshold: float = 30.0
    gap_open: int = -10
    gap_extend: int = -1
    evalue: float = 0.05
    word_size: int = 3
    max_targets: int = 50
    matrix: str = "BLOSUM62"
    low_complexity: bool = True

    def to_legacy_dict(self) -> dict[str, Any]:
        return {
            "threshold": self.threshold,
            "gap_open": self.gap_open,
            "gap_extend": self.gap_extend,
            "evalue": self.evalue,
            "word_size": self.word_size,
            "max_targets": self.max_targets,
            "matrix": self.matrix,
            "low_complexity": self.low_complexity,
        }


@dataclass(slots=True)
class BlastSearchResult:
    """
    Result returned by BlastSearchService.
    """

    hits: list[dict] = field(default_factory=list)
    algorithm_used: str = ""
    query_sequence: str = ""
    backend_error: str = ""


class BlastSearchService:
    """
    Select and run BLAST-like search backends.

    The service currently delegates the actual algorithms to the legacy
    analyzer object. In the next refactor pass, run_ncbi_blast, kmer_blast
    and sw_blast should move into ppigfinder.annotation.
    """

    def __init__(self, analyzer):
        self.analyzer = analyzer

    def search(
        self,
        query_sequence: str,
        orfs: list[dict],
        algorithm: str = "Auto",
        params: BlastSearchParams | None = None,
    ) -> BlastSearchResult:
        params = params or BlastSearchParams()
        legacy_params = params.to_legacy_dict()

        query_sequence = query_sequence.strip().upper()
        proteins = [orf.get("protein", "") for orf in orfs]

        hits = None
        algorithm_used = algorithm
        backend_error = ""

        self.analyzer._last_blast_error = ""

        if algorithm.startswith("Auto"):
            if BACKENDS.get("blast+", {}).get("available"):
                hits = self.analyzer.run_ncbi_blast(
                    query_sequence,
                    orfs,
                    legacy_params,
                )
                algorithm_used = "NCBI BLAST+"

            if hits is None:
                backend_error = self.analyzer._last_blast_error
                hits = self.analyzer.kmer_blast(
                    query_sequence,
                    proteins,
                    legacy_params,
                )
                algorithm_used = "K-mer Filter"

        elif algorithm.startswith("NCBI"):
            hits = self.analyzer.run_ncbi_blast(
                query_sequence,
                orfs,
                legacy_params,
            )
            algorithm_used = "NCBI BLAST+"

            if hits is None:
                backend_error = self.analyzer._last_blast_error
                hits = self.analyzer.kmer_blast(
                    query_sequence,
                    proteins,
                    legacy_params,
                )
                algorithm_used = "K-mer fallback"

        elif algorithm.startswith("K-mer"):
            hits = self.analyzer.kmer_blast(
                query_sequence,
                proteins,
                legacy_params,
            )
            algorithm_used = "K-mer Filter"

        else:
            hits = self.analyzer.sw_blast(
                query_sequence,
                proteins,
                legacy_params,
            )
            algorithm_used = "Smith-Waterman"

        return BlastSearchResult(
            hits=hits or [],
            algorithm_used=algorithm_used,
            query_sequence=query_sequence,
            backend_error=backend_error,
        )
