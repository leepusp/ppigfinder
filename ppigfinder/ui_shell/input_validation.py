#!/usr/bin/env python3
"""
Input validation and lightweight metadata extraction for the guided UI shell.

This module is intentionally small and independent from the legacy GUI.
It gives the guided shell useful feedback before the full backend state
binding is connected.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re


@dataclass(slots=True)
class GenomeInputSummary:
    path: str
    file_type: str
    name: str
    sequence_count: int = 0
    total_length: int = 0
    longest_length: int = 0
    gc_percent: float | None = None
    valid: bool = False
    message: str = ""


def _gc_percent(sequence: str) -> float | None:
    sequence = sequence.upper()
    bases = [b for b in sequence if b in "ACGT"]
    if not bases:
        return None

    gc = sum(1 for b in bases if b in "GC")
    return round((gc / len(bases)) * 100.0, 2)


def _read_text(path: Path) -> str:
    for encoding in ("utf-8", "latin-1", "cp1252"):
        try:
            return path.read_text(encoding=encoding)
        except UnicodeDecodeError:
            continue
    return path.read_text(errors="ignore")


def _parse_fasta(path: Path) -> GenomeInputSummary:
    text = _read_text(path)

    if not text.lstrip().startswith(">"):
        return GenomeInputSummary(
            path=str(path),
            file_type="FASTA",
            name=path.name,
            valid=False,
            message="FASTA file does not start with a header line beginning with '>'.",
        )

    sequences = []
    current = []

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            if current:
                sequences.append("".join(current))
                current = []
        else:
            current.append(re.sub(r"[^A-Za-z]", "", line))

    if current:
        sequences.append("".join(current))

    lengths = [len(seq) for seq in sequences if seq]

    if not lengths:
        return GenomeInputSummary(
            path=str(path),
            file_type="FASTA",
            name=path.name,
            valid=False,
            message="No nucleotide sequence could be parsed from the FASTA file.",
        )

    longest = max(sequences, key=len)

    return GenomeInputSummary(
        path=str(path),
        file_type="FASTA",
        name=path.name,
        sequence_count=len(lengths),
        total_length=sum(lengths),
        longest_length=max(lengths),
        gc_percent=_gc_percent(longest),
        valid=True,
        message="FASTA input parsed successfully.",
    )


def _parse_genbank(path: Path) -> GenomeInputSummary:
    text = _read_text(path)

    locus = path.name
    total_length = 0

    for line in text.splitlines():
        if line.startswith("LOCUS"):
            parts = line.split()
            if len(parts) >= 2:
                locus = parts[1]
            for part in parts:
                if part.isdigit():
                    total_length = int(part)
                    break
            break

    origin_index = text.find("ORIGIN")
    sequence = ""

    if origin_index != -1:
        sequence_text = text[origin_index:]
        sequence = re.sub(r"[^A-Za-z]", "", sequence_text.replace("ORIGIN", ""))
        sequence = sequence.upper()

        # Remove common end marker residue if it entered as letters.
        sequence = sequence.replace("END", "")

    if sequence:
        total_length = len(sequence)

    if total_length <= 0:
        return GenomeInputSummary(
            path=str(path),
            file_type="GenBank",
            name=locus,
            valid=False,
            message="Could not detect sequence length in GenBank file.",
        )

    return GenomeInputSummary(
        path=str(path),
        file_type="GenBank",
        name=locus,
        sequence_count=1,
        total_length=total_length,
        longest_length=total_length,
        gc_percent=_gc_percent(sequence) if sequence else None,
        valid=True,
        message="GenBank input parsed successfully.",
    )


def validate_genome_input(path: str) -> GenomeInputSummary:
    file_path = Path(path)

    if not file_path.exists():
        return GenomeInputSummary(
            path=path,
            file_type="Unknown",
            name=file_path.name,
            valid=False,
            message="File does not exist.",
        )

    suffix = file_path.suffix.lower()

    if suffix in {".fa", ".fasta", ".fna", ".ffn"}:
        return _parse_fasta(file_path)

    if suffix in {".gb", ".gbk", ".genbank"}:
        return _parse_genbank(file_path)

    if suffix == ".dna":
        return GenomeInputSummary(
            path=str(file_path),
            file_type="SnapGene",
            name=file_path.name,
            valid=True,
            message=(
                "SnapGene file selected. Full SnapGene parsing is available through "
                "the ppigFinder backend and will be connected to the guided shell."
            ),
        )

    return GenomeInputSummary(
        path=str(file_path),
        file_type="Unknown",
        name=file_path.name,
        valid=False,
        message="Unsupported file extension for guided genome validation.",
    )


def summary_to_state(summary: GenomeInputSummary) -> dict:
    return {
        "genome_name": summary.name,
        "genome_file_type": summary.file_type,
        "genome_sequence_count": summary.sequence_count,
        "genome_total_length": summary.total_length,
        "genome_longest_length": summary.longest_length,
        "genome_gc_percent": summary.gc_percent,
        "genome_valid": summary.valid,
        "genome_message": summary.message,
    }
