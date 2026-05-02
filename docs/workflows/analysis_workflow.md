# ppigFinder Analysis Workflow

ppigFinder is designed as a neighbourhood-guided workflow for investigating candidate protein-protein interactions directly from bacterial genome sequence data.

The central idea is to reduce the AlphaFold 3 search space before structural prediction. Instead of testing all possible protein pairs in a bacterial proteome, ppigFinder combines ORF prediction, functional annotation and genomic neighbourhood context to prioritize biologically plausible candidate interactions.

## Workflow summary

A typical ppigFinder session follows seven major steps:

1. Load a genome sequence.
2. Predict ORFs using Pyrodigal, six-frame translation or hybrid mode.
3. Annotate ORFs using HMM profile search and/or BLASTp.
4. Inspect genomic neighbourhoods around candidate ORFs.
5. Generate AlphaFold 3 jobs for selected protein pairs or groups.
6. Submit or run AlphaFold 3 predictions externally or through configured HPC workflows.
7. Import and analyse AlphaFold 3 outputs using confidence metrics and visual reports.

```text
Genome input
  -> ORF prediction
  -> BLAST / HMM annotation
  -> neighbourhood inspection
  -> AF3 job construction
  -> AF3 execution / submission
  -> AF3 result analysis
  -> reports and exports
```

## 1. Genome input

### Purpose

Load the bacterial genome or nucleotide sequence that will be used as the basis for the analysis.

### Accepted input

- FASTA nucleotide files
- multi-FASTA files
- GenBank files
- SnapGene files

For multi-FASTA input, the longest contig can be selected automatically for the main analysis workflow.

### Output

- loaded DNA sequence
- genome metadata
- coordinate system for ORF prediction and visualization

## 2. ORF prediction

### Purpose

Identify protein-coding open reading frames from the loaded genome sequence.

### Supported strategies

ppigFinder supports complementary ORF prediction strategies:

- **Pyrodigal / Prodigal-like gene calling**: gene prediction based on GC-content modelling, coding potential and ribosome binding site scoring.
- **Exhaustive six-frame translation**: searches all six reading frames for ORFs defined by configurable start and stop codons.
- **Hybrid mode**: combines external gene-calling logic with built-in ORF scanning.

### Typical codons

Default translated ORF search uses bacterial start codons such as:

- ATG
- GTG
- TTG

and standard stop codons.

### Output

- ORF table
- genomic coordinates
- strand and reading frame
- predicted protein sequences
- GC content
- prediction source

## 3. Functional annotation

### Purpose

Annotate predicted ORFs using sequence similarity, conserved domain detection and genomic context.

### BLAST / similarity search

Protein homology search can be performed using:

- local NCBI BLAST+ when available
- k-mer pre-filter plus Smith-Waterman fallback
- full Smith-Waterman alignment fallback

Typical output includes:

- ranked hits
- percent identity
- coverage
- score
- E-value when available

### HMM/domain annotation

HMM profile scanning can use HMMER3 when available. Profiles may come from:

- Pfam
- TIGRFAM
- custom HMM profile databases

When HMMER3 is unavailable, internal fallback scanners can be used for selected workflows.

### Output

- annotated ORFs
- HMM/domain hits
- similarity hits
- colour-coded annotation layers for the genome map

## 4. Genomic neighbourhood analysis

### Purpose

Use prokaryotic genome organization to prioritize candidate protein-protein interactions.

In bacterial genomes, functionally related and physically interacting proteins are often encoded close to one another, frequently in operons or conserved genomic neighbourhoods. ppigFinder uses this principle to reduce the candidate interaction search space.

### Typical analysis

For a selected ORF, the neighbourhood module can inspect coding sequences within a configurable genomic window, for example:

- upstream ORFs
- downstream ORFs
- ORFs sharing domain annotations
- ORFs co-localized around a query locus

### Output

- neighbourhood table
- candidate partner list
- exported neighbourhood FASTA
- candidates for AlphaFold 3 prediction

## 5. AlphaFold 3 job construction

### Purpose

Prepare protein-protein interaction prediction jobs from biologically prioritized candidates.

### Supported job concepts

ppigFinder is designed to support multiple AlphaFold 3 job construction modes, including:

- query ORF versus neighbouring ORFs
- query ORF plus homodimer jobs
- three-chain complexes involving neighbouring ORFs
- selected ORFs all-versus-all
- ORFs sharing positive HMM annotation
- selected annotation-table ORFs
- symmetric homodimers
- genome-wide interactome screens when computationally feasible
- custom stoichiometries with multiple protein chains

### AlphaFold Server JSON

Jobs are exported using an AlphaFold Server-compatible JSON structure. Each job includes:

- unique job name
- protein chain entries
- amino acid sequences
- chain copy number
- model seeds

### Output

- individual JSON job files
- batch JSON arrays
- partitioned job chunks for large sessions

## 6. HPC or external execution

### Purpose

Run AlphaFold 3 predictions externally or through configured infrastructure.

### Supported directions

ppigFinder can support workflows involving:

- AlphaFold Server upload
- local AlphaFold 3 installations
- HPC job submission through SSH/SFTP
- SLURM, PBS/Torque or LSF scheduler integrations

Large sessions can be partitioned into sequential chunks to reduce memory and submission problems.

### Output

- submitted jobs
- result folders
- downloadable AF3 outputs

## 7. AF3 result analysis

### Purpose

Interpret AlphaFold 3 predictions using interaction confidence metrics.

### Metrics

ppigFinder extracts and reports metrics such as:

- **ipTM**: interface predicted TM-score
- **pTM**: global predicted TM-score
- **ranking_score**: model ranking score
- **PAE_inter**: mean inter-chain predicted aligned error
- **PAE_min**: minimum focal inter-chain PAE
- **cp_ipTM**: chain-pair ipTM when available
- **contact percentage**: percentage of inter-chain cells below a PAE threshold
- **interaction classification**: confidence class derived from selected metrics

### Interpretation

Low inter-chain PAE values and high interface confidence support the presence of a stable predicted interaction interface. These scores should be interpreted as structural evidence, not as definitive proof of in vivo interaction.

### Output

- AF3 result table
- interaction confidence summaries
- PAE heatmaps
- pLDDT plots
- TSV/CSV exports
- HTML reports

## Reports and reproducibility

ppigFinder can export:

- ORF tables
- FASTA files
- AF3 Server JSON
- AF3 result tables
- Project Snapshot v3 JSON
- HTML reports

The Project Snapshot format is intended to support reproducible and interruptible workflows.
