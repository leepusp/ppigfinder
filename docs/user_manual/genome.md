# DNA / Genome Module

## Purpose

The DNA / Genome module loads and inspects nucleotide sequence data before downstream protein-level analyses.

## Input data

- FASTA nucleotide files
- multi-FASTA files
- GenBank files
- SnapGene files

## Main tasks

- load genome sequence
- inspect genome length and metadata
- provide coordinates for ORF prediction
- translate genome regions when needed
- support genome map visualization and export

## Expected output

- loaded genome sequence
- genome metadata
- coordinate system for ORFs
- genome map data

## Recommended next step

After loading and inspecting the genome, proceed to Protein / ORFs.


## Automatic guided loading

In the guided shell, selecting a genome file triggers immediate validation and metadata extraction. A genome inspection window opens with the loaded file, organism/header information, sequence preview, length, GC content and validation status. If the file is valid, the workflow automatically advances to Protein / ORFs.


## Numbered sequence preview

The guided genome inspector displays sequence previews with coordinate numbering and fixed-width formatting. This makes the loaded nucleotide sequence easier to inspect, similar to genome database views. The inspector also exposes the next workflow actions, including moving to Protein / ORFs or starting ORF prediction directly.
