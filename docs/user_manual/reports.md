# Reports Module

## Purpose

The Reports module generates final outputs for interpretation, sharing and reproducibility.

## Input data

- loaded genome
- predicted ORFs
- BLAST results
- HMM/domain annotations
- AlphaFold jobs
- parsed AF3 results

## Main outputs

ppigFinder can generate:

- HTML reports
- Project Snapshot v3 JSON files
- ORF tables
- FASTA exports
- AF3 JSON files
- TSV/CSV result tables

## Project Snapshot v3

The Project Snapshot is a versioned JSON representation of the analysis state. It is intended to support reproducible, interruptible and shareable workflows.

## Recommended use

Generate reports after completing the major analysis steps or whenever you need to archive the current state of a project.
