# ppigFinder v2.9.14 Quick Start Guide

This page is a cleaned, Git-ready quick-start guide derived from the initial incomplete slide guide and updated for the current ppigFinder v2.9.14 / v29.14 release.

## Purpose

ppigFinder is open-source software for scientific and educational use in molecular biology, microbiology, and bioinformatics. It helps users start from a bacterial genome sequence and organize a workflow for protein-protein interaction discovery without requiring advanced command-line or high-performance-computing experience.

The workflow integrates:

- bacterial genome loading from nucleotide FASTA files;
- open reading frame (ORF) prediction;
- gene and protein annotation;
- BLAST-based homolog search;
- HMM/profile-based domain search;
- genomic-neighborhood inspection;
- AlphaFold 3 job generation and export;
- remote HPC submission workflows;
- downstream organization and interpretation of structural-prediction results.

## Requirements

Core runtime requirements:

- Python >= 3.8
- PyQt6 >= 6.4, or PyQt5 >= 5.15
- matplotlib >= 3.5
- NumPy >= 1.21

Recommended analysis tools:

- Pyrodigal >= 2.0 for prokaryotic ORF prediction
- NCBI BLAST+ >= 2.12 for local sequence similarity searches
- HMMER3 >= 3.0 for profile-HMM searches
- Paramiko >= 2.9 for SSH/SFTP-based HPC submission

## Basic workflow

### 1. Launch ppigFinder

Start the graphical interface using the installed entry point or the repository launcher:

```bash
ppigfinder.ui
```

or, from the repository root:

```bash
python main.py
```

After launch, verify the backend status indicators. Green check marks indicate that optional tools such as BLAST+, HMMER3, and Pyrodigal were detected.

### 2. Load a genome FASTA file

Use:

```text
File > Open FASTA
```

Select a bacterial genome sequence file. The Genome tab should update with information such as genome length, GC content, ORF count, and annotation status.

Example used in the initial guide:

- Organism: *Xanthomonas citri* strain 306
- NCBI accession: AE008923.1

### 3. Save the project

Before running longer analyses, save the project:

```text
File > Save Project
```

ppigFinder creates a structured project directory containing subdirectories for genome files, BLAST results, HMM searches, AlphaFold predictions, and final results.

### 4. Predict ORFs

Use:

```text
Translate genome
```

Available ORF-prediction strategies include:

- Pyrodigal
- Hybrid
- Automatic

For complete analyses, compare different ORF-prediction strategies when appropriate. For a simple tutorial, Pyrodigal is a good first choice.

After prediction, the genome map becomes visible and the ORF table is populated with predicted coding sequences.

### 5. Explore ORFs

The ORF table and genome map are synchronized:

- clicking an ORF in the table centers it on the genome map;
- clicking an ORF on the map highlights it in the table;
- right-clicking ORFs opens context-menu actions for annotation, sequence extraction, neighborhood inspection, and structural-prediction setup.

Use zoom controls gradually for large genomes.

### 6. Identify proteins of interest with BLAST

Use the BLAST query panel to paste a protein sequence and search for putative homologs among predicted ORFs. Parameters such as E-value and scoring options can be adjusted in the interface.

BLAST hits are linked back to the ORF table and genome map, allowing spatial inspection of candidate homologs.

### 7. Identify domains or protein families with HMM profiles

HMM profiles can be custom-built from aligned sequences or obtained from public resources such as Pfam/EBI.

Basic HMM workflow:

1. Load one or more `.hmm` profiles.
2. Run the search across all predicted ORFs.
3. Review significant hits in the HMM results panel.
4. Use the Parameters menu to adjust HMM reporting thresholds.

HMM hits can be visualized on the genome map and used to select proteins for downstream AlphaFold 3 job generation.

### 8. Build AlphaFold 3 jobs

ORFs can be sent to the AlphaFold job list from the ORF table context menu. Common job-building options include:

- add selected ORFs to the AF3 list;
- predict selected proteins against genomic neighbors;
- predict all pairwise combinations among selected proteins;
- define custom multi-chain complexes.

After generating the job list, export AF3-compatible JSON files for local, server, or HPC-based prediction workflows.

### 9. Export and submit jobs

Use the AlphaFold export tools to generate JSON batches for AlphaFold 3 prediction. Depending on the target execution environment, jobs may be submitted manually to an AlphaFold server or prepared for private/HPC execution workflows such as DaVinci.

For DaVinci-specific workflows, refer to:

```text
docs/user_manual/davinci_hpc.md
docs/developer/structure_prediction/davinci_af3_validation.md
```

## Notes for documentation maintenance

The initial slide guide contains useful workflow content but includes screenshots from older ppigFinder versions. This Markdown guide should be kept synchronized with the current v2.9.14 / v29.14 interface and updated whenever menus, job-generation modes, or output columns change.
