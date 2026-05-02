#!/usr/bin/env python3
"""
Documentation/help content for workspace modules.
"""

from __future__ import annotations


MODULE_DOCS = {
    "overview": {
        "purpose": "General overview of the analysis flow.",
        "input": "Project state",
        "output": "Guided navigation across Data, Genome, ORFs, Annotation, AlphaFold/PPI and Reports",
        "next": "Start with Data / Project, then proceed through Genome, ORFs, Annotation, AlphaFold/PPI and Reports.",
    },
    "data": {
        "purpose": "Create a new analysis session, validate input data, or restore an existing project. This is the entry point for inserting genome data, opening saved projects, or importing a reproducible project snapshot.",
        "input": "FASTA, multi-FASTA, GenBank, SnapGene, ppigFinder project files, or Project Snapshot v3 JSON files",
        "output": "Loaded genome/project state ready for ORF prediction, annotation, AlphaFold/PPI analysis, and reporting",
        "next": "After selecting and validating data, continue to DNA / Genome for inspection or Protein / ORFs for ORF prediction.",
    },
    "genome": {
        "purpose": "Load and inspect genomic DNA sequence data.",
        "input": "FASTA, GenBank or SnapGene DNA files",
        "output": "Genome loaded into the analysis workspace",
        "next": "Proceed to Protein / ORFs after loading data.",
    },
    "orfs": {
        "purpose": "Predict, inspect and export protein-coding ORFs. The guided shell currently provides a lightweight six-frame ORF scan and will progressively connect Pyrodigal/hybrid backends.",
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
