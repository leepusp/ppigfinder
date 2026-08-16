# ppigFinder

**ppigFinder** — **Protein–Protein Interaction Genomic Finder** — is a computational tool developed to support the investigation of protein–protein interactions in bacterial genomes.

The software integrates genome loading, ORF prediction, sequence-based annotation, genomic neighbourhood analysis, AlphaFold 3 job preparation, AlphaFold result parsing, and report generation into a desktop graphical workflow.

## Scientific motivation

Protein–protein interactions are fundamental to cellular function, but systematic discovery of candidate interactions remains challenging in bacterial genomes. Exhaustive all-versus-all structural prediction across complete proteomes is often computationally expensive and difficult to interpret.

ppigFinder follows a targeted, neighbourhood-guided strategy. Instead of testing every possible protein pair, it helps the user identify biologically plausible candidates by combining:

- ORF prediction from nucleotide sequence data
- protein homology search
- HMM/domain annotation
- genomic neighbourhood context
- AlphaFold 3-compatible job generation
- interaction metric interpretation

This design reduces the structural search space while preserving biologically motivated candidates derived from genomic context.

## Main features

- Load bacterial genome data from FASTA, GenBank, and SnapGene-compatible files.
- Predict ORFs using Pyrodigal, exhaustive six-frame scanning, or hybrid approaches.
- Translate nucleotide sequences into predicted proteins.
- Search predicted ORFs using BLAST+ or built-in fallback approaches.
- Annotate conserved domains using HMMER3 or internal fallback scanners.
- Inspect genomic neighbourhoods around ORFs of interest.
- Build AlphaFold 3 prediction jobs for selected protein pairs or groups.
- Export AlphaFold Server-compatible JSON files.
- Import AlphaFold 3 result folders and extract interaction metrics.
- Interpret ipTM, pTM, cp_ipTM, PAE_inter, PAE_min, contact percentage, and interaction confidence.
- Export ORF tables, FASTA files, TSV/CSV results, project snapshots, and HTML reports.
- Run as a local desktop application without requiring internet access for core analyses.

## Typical workflow

A typical ppigFinder session follows these steps:

1. Load a bacterial genome sequence.
2. Predict ORFs using Pyrodigal, six-frame translation, or hybrid mode.
3. Annotate ORFs using BLAST and/or HMM profile scanning.
4. Inspect the genomic neighbourhood of candidate ORFs.
5. Generate AlphaFold 3 jobs for neighbourhood-selected protein pairs.
6. Submit or run AlphaFold 3 predictions externally or through configured HPC workflows.
7. Import and analyse AlphaFold 3 outputs.
8. Export tables, reports, and project snapshots.

## Current interface

The current stable interface is launched with:

```bash
python main.py
