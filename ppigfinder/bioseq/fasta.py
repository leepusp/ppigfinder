#!/usr/bin/env python3
"""
FASTA utilities exposed from the bioseq namespace.

The implementation lives in ppigfinder.io.fasta because FASTA parsing and
writing are file I/O concerns.
"""

from ppigfinder.io.fasta import (
    FastaRecord,
    clean_sequence,
    parse_fasta_content,
    read_fasta,
    choose_longest_record,
    wrap_sequence,
    write_fasta_records,
    make_orf_protein_records,
    write_orf_protein_fasta,
)

__all__ = [
    "FastaRecord",
    "clean_sequence",
    "parse_fasta_content",
    "read_fasta",
    "choose_longest_record",
    "wrap_sequence",
    "write_fasta_records",
    "make_orf_protein_records",
    "write_orf_protein_fasta",
]
