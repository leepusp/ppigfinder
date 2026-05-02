#!/usr/bin/env python3
"""
FASTA input/output utilities for ppigFinder.

This module contains FASTA parsing and writing logic independent from the
PyQt interface.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re


@dataclass
class FastaRecord:
    """
    A single FASTA record.
    """

    name: str
    sequence: str
    description: str = ""


def clean_sequence(sequence: str, alphabet: str = "letters") -> str:
    """
    Remove non-sequence characters.

    alphabet='letters' keeps alphabetic symbols only.
    """
    if alphabet == "letters":
        return re.sub(r"[^A-Za-z]", "", sequence).upper()

    return re.sub(r"\s+", "", sequence).upper()


def parse_fasta_content(content: str, allow_plain_sequence: bool = True) -> list[FastaRecord]:
    """
    Parse FASTA content into records.

    If no FASTA header is found and allow_plain_sequence=True, the whole
    content is treated as a plain sequence.
    """
    records: list[FastaRecord] = []
    current_header: str | None = None
    current_chunks: list[str] = []

    def flush_record() -> None:
        nonlocal current_header, current_chunks

        if current_header is None:
            return

        header = current_header.strip()
        parts = header.split(maxsplit=1)
        name = parts[0] if parts else "sequence"
        description = parts[1] if len(parts) > 1 else ""
        sequence = clean_sequence("".join(current_chunks))

        if sequence:
            records.append(FastaRecord(name=name, description=description, sequence=sequence))

    for raw_line in content.splitlines():
        line = raw_line.strip()

        if not line:
            continue

        if line.startswith(">"):
            flush_record()
            current_header = line[1:].strip()
            current_chunks = []
        elif current_header is not None:
            current_chunks.append(line)

    flush_record()

    if not records and allow_plain_sequence:
        sequence = clean_sequence(content)
        if sequence:
            records.append(FastaRecord(name="sequence", sequence=sequence))

    return records


def read_fasta(path: str | Path, allow_plain_sequence: bool = True) -> list[FastaRecord]:
    """
    Read a FASTA file.
    """
    path = Path(path)

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        content = handle.read()

    records = parse_fasta_content(
        content,
        allow_plain_sequence=allow_plain_sequence,
    )

    if len(records) == 1 and records[0].name == "sequence":
        records[0].name = path.stem

    return records


def choose_longest_record(records: list[FastaRecord]) -> FastaRecord:
    """
    Return the longest sequence from a list of FASTA records.
    """
    if not records:
        raise ValueError("No FASTA sequence found.")

    return max(records, key=lambda record: len(record.sequence))


def wrap_sequence(sequence: str, width: int = 80) -> str:
    """
    Wrap a biological sequence to fixed-width FASTA lines.
    """
    return "\n".join(sequence[i:i + width] for i in range(0, len(sequence), width))


def write_fasta_records(
    path: str | Path,
    records: list[FastaRecord],
    width: int = 80,
) -> None:
    """
    Write FASTA records to disk.
    """
    path = Path(path)

    with path.open("w", encoding="utf-8") as handle:
        for record in records:
            header = record.name
            if record.description:
                header += f" {record.description}"

            handle.write(f">{header}\n")
            handle.write(wrap_sequence(record.sequence, width=width))
            handle.write("\n")


def make_orf_protein_records(orfs: list[dict]) -> list[FastaRecord]:
    """
    Convert legacy ORF dictionaries into protein FASTA records.
    """
    records: list[FastaRecord] = []

    for index, orf in enumerate(orfs, start=1):
        domains = "|".join(
            domain.get("domain", domain.get("name", ""))
            for domain in orf.get("domains", [])
            if isinstance(domain, dict)
        )
        domains = domains.strip("|")

        header = (
            f"ORF{index}"
            f"|F{orf.get('frame', '')}{orf.get('strand', '')}"
            f"|{orf.get('start', '')}-{orf.get('end', '')}"
        )

        if domains:
            header += f"|{domains}"

        protein = str(orf.get("protein", "")).rstrip("*")

        records.append(
            FastaRecord(
                name=header,
                sequence=protein,
            )
        )

    return records


def write_orf_protein_fasta(
    path: str | Path,
    orfs: list[dict],
    width: int = 80,
) -> None:
    """
    Write ORF protein sequences in FASTA format.
    """
    records = make_orf_protein_records(orfs)
    write_fasta_records(path, records, width=width)
