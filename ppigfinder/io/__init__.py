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

from .html_report import render_basic_report, write_basic_report

__all__ += ["render_basic_report", "write_basic_report"]

from .project_json import read_project_json, write_project_json, validate_project_json

__all__ += ["read_project_json", "write_project_json", "validate_project_json"]

from .html_report import (
    render_project_report,
    write_project_report,
    write_report_from_project_json,
)
__all__ += [
    "render_project_report",
    "write_project_report",
    "write_report_from_project_json",
]

from .af3_table_export import write_af3_results_table, AF3_RESULT_COLUMNS

__all__ += ["write_af3_results_table", "AF3_RESULT_COLUMNS"]
