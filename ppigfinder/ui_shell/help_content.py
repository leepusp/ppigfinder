#!/usr/bin/env python3
"""
Embedded guidance content for ppigFinder analysis steps.
"""

from __future__ import annotations


HELP_CONTENT = {
    "open_current_interface": {
        "purpose": "Open the current complete ppigFinder interface while the guided interface is still under development.",
        "when": "Use this when you need access to all existing tools in the current Qt interface.",
        "data": "Current ppigFinder application.",
        "next": "Use the current interface for complete analysis execution.",
    },
    "open_genome": {
        "purpose": "Start a new project by loading DNA/genome data.",
        "when": "Use this as the first step for a new genome-based analysis.",
        "data": "FASTA, GenBank or SnapGene nucleotide files.",
        "next": "Continue to the Data / Project module to load files, then proceed to ORF prediction.",
    },
    "open_project": {
        "purpose": "Restore a previously saved ppigFinder analysis.",
        "when": "Use this when continuing a previous project.",
        "data": "ppigFinder project files or project snapshots.",
        "next": "Review loaded data, ORFs, annotations or AlphaFold results.",
    },
    "predict_orfs": {
        "purpose": "Move to the Protein / ORFs module, where ORF prediction and protein export are organized.",
        "when": "Use this after loading a DNA/genome sequence.",
        "data": "Loaded DNA/genome sequence.",
        "next": "Predict ORFs, inspect the ORF table and export protein sequences if needed.",
    },
    "annotation": {
        "purpose": "Move to the Annotation module for BLAST, HMM/domain and neighborhood analyses.",
        "when": "Use this after ORFs have been predicted.",
        "data": "Predicted proteins and ORF genomic coordinates.",
        "next": "Run similarity/domain/context analyses to select candidates for downstream interpretation.",
    },
    "alphafold": {
        "purpose": "Move to the AlphaFold / PPI module for interaction-oriented structural analysis.",
        "when": "Use this after selecting candidate ORFs or protein pairs.",
        "data": "Protein sequences, AlphaFold Server JSON or AF3 result folders.",
        "next": "Prepare AF3 input, import AF3 outputs and interpret interaction metrics.",
    },
    "guided_workspace": {
        "purpose": "Open the modular workspace preview organized by data type and analysis stage.",
        "when": "Use this to explore the future guided interface.",
        "data": "Current or future project state.",
        "next": "Choose a module such as Data, Genome, ORFs, Annotation, AlphaFold or Reports.",
    },
    "reports": {
        "purpose": "Move to the Reports module for final documentation and exports.",
        "when": "Use this after loading data, predicting ORFs, annotating proteins or importing AF3 results.",
        "data": "Current project state.",
        "next": "Export HTML reports, project snapshots and tabular results.",
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
