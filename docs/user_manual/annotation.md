# Annotation Module

## Purpose

The Annotation module links predicted ORFs to functional evidence using similarity search, conserved domain detection and genomic neighbourhood context.

## Input data

- predicted ORF proteins
- protein query sequence
- HMM/domain profiles
- ORF genomic coordinates

## BLAST / similarity search

ppigFinder can use local BLAST+ when available. It can also use built-in fallback approaches such as k-mer filtering and Smith-Waterman alignment.

Typical output includes:

- percent identity
- alignment score
- coverage
- E-value when available
- mapped ORF candidates

## HMM/domain annotation

HMM profile scanning can use HMMER3 when available. Profiles may come from Pfam, TIGRFAM or custom databases.

Output includes domain assignments that can be propagated to ORF tables and genome visualizations.

## Neighbourhood analysis

The neighbourhood step inspects ORFs around a selected candidate within a configurable genomic window. This is central to the ppigFinder strategy because co-localized bacterial genes often encode functionally or physically associated proteins.

## Recommended next step

Use annotation and neighbourhood evidence to select candidate protein pairs for AlphaFold / PPI analysis.


## Embedded candidate table

The normal guided workflow now embeds candidate ORFs directly inside the Annotation page. This is more intuitive than opening a separate table window because candidate review is part of the annotation process.

The separate candidate window is kept as an optional full-screen inspection mode for large projects.
