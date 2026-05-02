"""
Input/output modules for ppigFinder.
"""

from .fasta import (
    FastaRecord,
    read_fasta,
    parse_fasta_content,
    choose_longest_record,
    write_fasta_records,
    write_orf_protein_fasta,
)
from .snapgene import parse_snapgene_dna, write_snapgene_dna
from .genbank import parse_genbank, write_genbank

__all__ = [
    "FastaRecord",
    "read_fasta",
    "parse_fasta_content",
    "choose_longest_record",
    "write_fasta_records",
    "write_orf_protein_fasta",
    "parse_snapgene_dna",
    "write_snapgene_dna",
    "parse_genbank",
    "write_genbank",
]
