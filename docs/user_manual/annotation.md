# Annotation Module

## Purpose

The Annotation module links predicted ORFs to functional evidence using similarity search, domain detection and genomic context.

## Input data

- Predicted ORF proteins
- Protein query sequence
- HMM/domain profiles
- ORF genomic coordinates

## Main tasks

- Run BLAST or fallback similarity search
- Annotate conserved domains with HMMER or internal scanner
- Inspect genomic neighborhood context

## Expected output

- Similarity hits
- HMM/domain annotations
- Neighborhood/operon context

## Recommended next step

Use annotation evidence to select candidate ORFs or protein pairs for AlphaFold / PPI analysis.
