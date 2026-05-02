#!/usr/bin/env python3
"""
ORF prediction service.

This layer orchestrates ORF prediction and can be used by the GUI, CLI,
or batch workflows.
"""

from __future__ import annotations

from ppigfinder.bioseq.orf_finder import (
    find_orfs,
    find_orfs_hybrid,
    find_orfs_pyrodigal,
)


class ORFPredictionService:
    """
    High-level service for ORF prediction.
    """

    def predict_six_frame(
        self,
        dna_sequence: str,
        min_aa: int = 30,
        start_codons: set[str] | None = None,
    ) -> list[dict]:
        return find_orfs(
            dna_sequence,
            min_aa=min_aa,
            start_codons=start_codons,
        )

    def predict_pyrodigal(
        self,
        dna_sequence: str,
        meta: bool = True,
        min_aa: int = 30,
        closed_ends: bool = False,
        translation_table: int = 11,
        mask: bool = False,
    ) -> list[dict]:
        return find_orfs_pyrodigal(
            dna_sequence,
            meta=meta,
            min_aa=min_aa,
            closed_ends=closed_ends,
            translation_table=translation_table,
            mask=mask,
        )

    def predict_hybrid(
        self,
        dna_sequence: str,
        min_aa: int = 30,
        start_codons: set[str] | None = None,
        pyro_meta: bool = True,
        pyro_min_aa: int = 30,
        pyro_closed: bool = False,
        pyro_translation_table: int = 11,
        pyro_mask: bool = False,
        pyro_start_filter: dict | None = None,
    ) -> list[dict]:
        return find_orfs_hybrid(
            dna_sequence,
            min_aa=min_aa,
            start_codons=start_codons,
            pyro_meta=pyro_meta,
            pyro_min_aa=pyro_min_aa,
            pyro_closed=pyro_closed,
            pyro_translation_table=pyro_translation_table,
            pyro_mask=pyro_mask,
            pyro_start_filter=pyro_start_filter,
        )
