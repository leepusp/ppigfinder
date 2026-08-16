# Protein / ORFs Module

## Purpose

The Protein / ORFs module predicts protein-coding regions from a loaded genome sequence.

## Input data

- DNA/genome sequence

## ORF prediction modes

ppigFinder supports multiple strategies:

- Pyrodigal / Prodigal-like gene calling
- exhaustive six-frame translation
- hybrid ORF prediction

## Output

The ORF table can include:

- ORF identifier
- start and end coordinates
- strand
- reading frame
- size
- GC content
- predicted protein sequence
- prediction source
- annotation fields

## Recommended next step

After predicting ORFs, continue to Annotation to identify homologs, domains and genomic neighbourhood context.
