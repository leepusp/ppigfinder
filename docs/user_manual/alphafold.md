# AlphaFold / PPI Module

## Purpose

The AlphaFold / PPI module supports preparation and interpretation of protein-protein interaction predictions.

## Input data

- Selected ORF protein pairs
- AlphaFold Server JSON
- AlphaFold 3 output folders

## Main tasks

- Export AlphaFold Server-compatible JSON
- Import AF3 result folders
- Parse interaction metrics
- Classify interaction confidence

## Metrics

Important metrics include:

- ipTM
- pTM
- cp_ipTM
- PAE_inter
- PAE_min
- contact percentage
- interaction classification

## Expected output

- Parsed AF3 result table
- Interaction confidence classification
- TSV/CSV export
- Data for reports

## Recommended next step

After importing AF3 results, continue to Reports to export tables and HTML summaries.
