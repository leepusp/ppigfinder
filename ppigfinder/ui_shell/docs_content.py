#!/usr/bin/env python3
"""
Documentation/help content for workspace modules.
"""

from __future__ import annotations


MODULE_DOCS = {
    "overview": {
        "purpose": "General overview of the analysis flow.",
        "input": "Project state",
        "output": "Guided navigation across modules",
        "next": "Start with Data / Project if beginning a new analysis.",
    },
    "data": {
        "purpose": "Create, open or restore the analysis dataset.",
        "input": "Genome files, project files or snapshot files",
        "output": "Loaded project or genome workspace",
        "next": "After loading data, proceed to DNA / Genome and Protein / ORFs.",
    },
    "genome": {
        "purpose": "Load and inspect genomic DNA sequence data.",
        "input": "FASTA, GenBank or SnapGene DNA files",
        "output": "Genome loaded into the analysis workspace",
        "next": "Proceed to Protein / ORFs after loading data.",
    },
    "orfs": {
        "purpose": "Predict and inspect protein-coding ORFs.",
        "input": "Loaded genome sequence",
        "output": "ORF table and protein sequences",
        "next": "Continue to Annotation or export protein FASTA.",
    },
    "annotation": {
        "purpose": "Assign similarity hits, domains and genomic context.",
        "input": "Predicted protein sequences and ORF coordinates",
        "output": "BLAST/HMM/neighborhood annotations",
        "next": "Select candidates for AlphaFold / PPI analysis.",
    },
    "alphafold": {
        "purpose": "Prepare and analyze structural interaction predictions.",
        "input": "Protein pairs or AF3 result folders",
        "output": "Interaction metrics and classified results",
        "next": "Export AF3 tables or move to Reports.",
    },
    "reports": {
        "purpose": "Generate final project outputs.",
        "input": "Current project state",
        "output": "HTML, JSON and TSV/CSV exports",
        "next": "Archive or share results.",
    },
}


def docs_for(module_id: str) -> dict:
    return MODULE_DOCS.get(
        module_id,
        {
            "purpose": "Module documentation not defined yet.",
            "input": "N/A",
            "output": "N/A",
            "next": "Continue through the workflow.",
        },
    )
