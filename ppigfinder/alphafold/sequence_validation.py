#!/usr/bin/env python3
"""
Protein sequence validation for AlphaFold workflows.

AlphaFold Server JSON export should not silently submit malformed protein
sequences. This module validates and normalizes protein sequences before
job creation.
"""

from __future__ import annotations

from dataclasses import dataclass
import re


STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


@dataclass(slots=True)
class ProteinValidationResult:
    """
    Validation result for one protein sequence.
    """

    original_sequence: str
    cleaned_sequence: str
    valid: bool
    invalid_residues: list[str]
    original_length: int
    cleaned_length: int
    message: str = ""


def normalize_protein_sequence(sequence: str) -> str:
    """
    Remove whitespace, gaps and terminal stop codons, then uppercase.
    """
    sequence = sequence or ""
    sequence = re.sub(r"\s+", "", sequence)
    sequence = sequence.replace("-", "")
    sequence = sequence.upper()

    # Terminal stop codon from translated ORFs is acceptable to remove.
    sequence = sequence.rstrip("*")

    return sequence


def validate_protein_sequence(sequence: str) -> ProteinValidationResult:
    """
    Validate a protein sequence for AlphaFold Server export.
    """
    original = sequence or ""
    cleaned = normalize_protein_sequence(original)

    invalid = sorted(set(cleaned) - STANDARD_AA)

    if not cleaned:
        return ProteinValidationResult(
            original_sequence=original,
            cleaned_sequence=cleaned,
            valid=False,
            invalid_residues=[],
            original_length=len(original),
            cleaned_length=0,
            message="Protein sequence is empty after normalization.",
        )

    if invalid:
        return ProteinValidationResult(
            original_sequence=original,
            cleaned_sequence=cleaned,
            valid=False,
            invalid_residues=invalid,
            original_length=len(original),
            cleaned_length=len(cleaned),
            message=(
                "Protein sequence contains unsupported residues: "
                + ", ".join(invalid)
            ),
        )

    return ProteinValidationResult(
        original_sequence=original,
        cleaned_sequence=cleaned,
        valid=True,
        invalid_residues=[],
        original_length=len(original),
        cleaned_length=len(cleaned),
        message="OK",
    )


def clean_protein_sequence(sequence: str, strict: bool = True) -> str:
    """
    Return a normalized valid protein sequence.

    strict=True raises ValueError on unsupported residues.
    """
    result = validate_protein_sequence(sequence)

    if not result.valid and strict:
        raise ValueError(result.message)

    if not result.valid:
        return "".join(aa for aa in result.cleaned_sequence if aa in STANDARD_AA)

    return result.cleaned_sequence
