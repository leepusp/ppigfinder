# Protein / ORFs Module

## Purpose

The Protein / ORFs module predicts protein-coding open reading frames from the loaded genome.

## Input data

- Loaded DNA/genome sequence

## Main tasks

- Predict ORFs with Pyrodigal, six-frame scanning or hybrid logic
- Inspect ORF coordinates, strand, frame, size and GC content
- Export predicted protein sequences

## Expected output

- ORF table
- Protein sequence set
- FASTA export for downstream tools

## Recommended next step

After predicting ORFs, continue to Annotation for BLAST, HMM/domain and neighborhood analyses.
