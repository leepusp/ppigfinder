# AlphaFold / PPI Module

## Purpose

The AlphaFold / PPI module supports preparation and interpretation of protein-protein interaction predictions.

## Input data

- selected ORF protein sequences
- selected ORF pairs or groups
- AlphaFold Server JSON files
- AlphaFold 3 output folders

## Job construction

ppigFinder can generate AlphaFold Server-compatible JSON jobs. Each job contains:

- a unique job name
- protein chain definitions
- amino acid sequences
- copy numbers for chains
- model seed configuration

Supported prediction concepts include:

- pairwise heterodimers
- homodimers
- selected ORF combinations
- neighbourhood-guided pairs
- HMM-positive ORF sets
- custom multichain complexes
- broader all-versus-all screens when computationally feasible

## Output files

ppigFinder can export:

- individual AlphaFold Server JSON files
- consolidated batch JSON arrays
- partitioned JSON batches for larger sessions
- TSV/CSV result tables after parsing AF3 output

## Result metrics

Parsed AF3 result folders may include:

- ipTM
- pTM
- ranking_score
- cp_ipTM
- PAE_inter
- PAE_min
- contact percentage
- interaction classification

## Interpretation

High ipTM or cp_ipTM and low PAE_inter / PAE_min values support a more confident predicted interaction interface. These values should be interpreted as structural evidence and should ideally be combined with genomic context, annotation evidence and experimental validation.

## Recommended next step

After importing AF3 results, export the AF3 results table and generate an HTML report.


## Guided candidate pair builder

The guided shell includes an initial AlphaFold/PPI pair builder. The first strategy generates adjacent ORF pairs from the guided ORF prediction table. This supports the neighbourhood-guided concept of ppigFinder and prepares the interface for more advanced strategies such as query-vs-neighbours, HMM-positive ORFs, selected ORFs all-vs-all, homodimers and custom stoichiometry.

The generated pairs can be exported as AlphaFold Server-compatible JSON.
