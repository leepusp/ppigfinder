#!/usr/bin/env python3
"""
Basic biological sequence utilities used by ppigFinder.
"""

from __future__ import annotations

from .genetic_code import STANDARD_CODON_TABLE


_DNA_COMPLEMENT = str.maketrans({
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "a": "t",
    "t": "a",
    "g": "c",
    "c": "g",
    "N": "N",
    "n": "n",
})


def normalize_dna(sequence: str) -> str:
    """
    Return an uppercase DNA sequence containing only canonical symbols
    and N for unknown bases.
    """
    sequence = sequence.upper()
    return "".join(base if base in {"A", "T", "G", "C", "N"} else "N" for base in sequence)


def reverse_complement(sequence: str) -> str:
    """
    Return the reverse complement of a DNA sequence.
    """
    return sequence.translate(_DNA_COMPLEMENT)[::-1]


def gc_content(sequence: str) -> float:
    """
    Return GC content as percentage.
    """
    sequence = sequence.upper()
    if not sequence:
        return 0.0

    gc = sequence.count("G") + sequence.count("C")
    return gc / len(sequence) * 100.0


def translate_dna(
    dna_sequence: str,
    codon_table: dict[str, str] | None = None,
    stop_at_stop: bool = True,
) -> str:
    """
    Translate a DNA sequence into a protein sequence.
    """
    codon_table = codon_table or STANDARD_CODON_TABLE
    dna_sequence = dna_sequence.upper()

    protein: list[str] = []

    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i + 3]
        if len(codon) != 3:
            break

        amino_acid = codon_table.get(codon, "X")
        protein.append(amino_acid)

        if stop_at_stop and amino_acid == "*":
            break

    return "".join(protein)
