#!/usr/bin/env python3
"""
Embedded guidance content for ppigFinder analysis steps.
"""

from __future__ import annotations


HELP_CONTENT = {
    "open_genome": {
        "purpose": "Load the DNA/genome sequence that will be used as the basis for ORF prediction and downstream analyses.",
        "when": "Use this as the first step for a new analysis.",
        "data": "FASTA, GenBank or SnapGene nucleotide files.",
        "next": "After loading the genome, inspect the genome overview and run ORF prediction.",
    },
    "open_project": {
        "purpose": "Restore a previously saved ppigFinder session.",
        "when": "Use this when continuing a previous analysis.",
        "data": "ppigFinder project files.",
        "next": "Review loaded ORFs, annotations or AlphaFold results.",
    },
    "predict_orfs": {
        "purpose": "Detect protein-coding regions in the loaded genome.",
        "when": "Run this after loading a DNA/genome sequence.",
        "data": "DNA/genome sequence.",
        "next": "Review the ORF table, export proteins or run annotation.",
    },
    "annotation": {
        "purpose": "Associate predicted proteins with similarity hits, conserved domains and genomic context.",
        "when": "Use this after ORFs have been predicted.",
        "data": "Predicted proteins and ORF genomic coordinates.",
        "next": "Use BLAST/HMM/domain/neighborhood results to select candidates for AlphaFold/PPI analysis.",
    },
    "alphafold": {
        "purpose": "Prepare protein pairs for AlphaFold 3 and import predicted interaction results.",
        "when": "Use this after selecting candidate ORFs or protein pairs.",
        "data": "Protein sequences, AlphaFold Server JSON or AF3 result folders.",
        "next": "Score interactions using ipTM, cp_ipTM, PAE_min, PAE_inter and contact percentage.",
    },
    "reports": {
        "purpose": "Generate final outputs for documentation, interpretation and reproducibility.",
        "when": "Use this after ORF prediction, annotation or AF3 result import.",
        "data": "Current project state.",
        "next": "Export HTML, TSV/CSV tables or Project Snapshot v3.",
    },
}


def help_for(action_id: str) -> dict:
    return HELP_CONTENT.get(
        action_id,
        {
            "purpose": "This step is part of the ppigFinder guided workflow.",
            "when": "Use it when appropriate for the current analysis stage.",
            "data": "Depends on the selected workflow step.",
            "next": "Continue with the next analysis step.",
        },
    )
