#!/usr/bin/env python3
"""
Genome loading service.

This module centralizes genome/sequence loading independently from the GUI.
"""

from __future__ import annotations

from pathlib import Path

from ppigfinder.domain import GenomeRecord
from ppigfinder.io import (
    choose_longest_record,
    parse_genbank,
    parse_snapgene_dna,
    read_fasta,
)


class GenomeService:
    """
    High-level service for loading genome files.
    """

    FASTA_EXTENSIONS = {".fasta", ".fa", ".fna", ".faa"}
    GENBANK_EXTENSIONS = {".gb", ".gbk", ".genbank"}
    SNAPGENE_EXTENSIONS = {".dna"}

    def load_fasta(self, path: str | Path) -> GenomeRecord:
        path = Path(path)
        record = choose_longest_record(read_fasta(path))

        return GenomeRecord(
            name=record.name or path.stem,
            sequence=record.sequence.upper(),
            source_path=path,
        )

    def load_snapgene(self, path: str | Path) -> GenomeRecord:
        path = Path(path)
        result = parse_snapgene_dna(str(path))

        return GenomeRecord(
            name=result.get("name") or path.stem,
            sequence=result.get("sequence", "").upper(),
            source_path=path,
            topology=result.get("topology", "linear"),
            features=result.get("features", []),
            primers=result.get("primers", []),
            notes=result.get("notes", {}),
        )

    def load_genbank(self, path: str | Path) -> GenomeRecord:
        path = Path(path)
        result = parse_genbank(str(path))

        return GenomeRecord(
            name=result.get("name") or path.stem,
            sequence=result.get("sequence", "").upper(),
            source_path=path,
            topology=result.get("topology", "linear"),
            features=result.get("features", []),
            primers=result.get("primers", []),
            notes=result.get("notes", {}),
        )

    def load_by_extension(self, path: str | Path) -> GenomeRecord:
        path = Path(path)
        suffix = path.suffix.lower()

        if suffix in self.FASTA_EXTENSIONS:
            return self.load_fasta(path)

        if suffix in self.GENBANK_EXTENSIONS:
            return self.load_genbank(path)

        if suffix in self.SNAPGENE_EXTENSIONS:
            return self.load_snapgene(path)

        # Default fallback: try FASTA/plain sequence.
        return self.load_fasta(path)
