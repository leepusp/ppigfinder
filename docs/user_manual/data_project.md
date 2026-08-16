# Data / Project

## Purpose

This module is the entry point for a ppigFinder analysis. Use it to create a new analysis, open genome data, restore a previous project, or import a reproducible project snapshot.

## Input data

Supported inputs include:

- FASTA nucleotide files
- multi-FASTA genome files
- GenBank files
- SnapGene files
- ppigFinder project files
- Project Snapshot v3 JSON files

## Main actions

### Open genome file

Load the main nucleotide dataset for a new bacterial genome analysis.

### Open project

Resume a previous ppigFinder session with genome data, ORFs, annotations and saved analysis state.

### Import Project Snapshot v3

Load a portable JSON snapshot for reproducible analysis, reporting or continuation of previous work.

### Open full current interface

Advanced option for accessing the complete current interface while the guided workflow is still under development.

## Expected output

After this step, ppigFinder should have a loaded genome or restored project state ready for:

- DNA / Genome inspection
- ORF prediction
- annotation
- AlphaFold / PPI analysis
- reports and exports

## Recommended next step

After loading data, continue to **DNA / Genome** or **Protein / ORFs**.
