# Protein / ORFs



## ORF Discovery results window

After ORF prediction, the guided shell opens an ORF Discovery results window. This window summarizes the predicted protein-coding regions with:

- total predicted ORFs;
- mean and maximum ORF length;
- plus/minus strand counts;
- a coordinate table containing ORF ID, start, end, strand, frame and sequence preview;
- a compact ORF map preview;
- next workflow options for Annotation, FASTA export and AlphaFold/PPI analysis.

This step produces the protein sequence set used by BLAST, HMM/domain annotation, neighbourhood candidate selection and AlphaFold 3 job construction.


## Interactive ORF browser

The guided ORF results view includes an interactive browser inspired by locus visualization tools.

Current features:

- filter ORFs by ID or protein sequence fragment;
- filter by strand;
- filter by minimum amino-acid length;
- inspect a selected ORF;
- focus the selected ORF on the map;
- zoom and pan across the genomic coordinate map;
- view protein and DNA sequence previews;
- move directly to annotation, candidate review, FASTA export or AlphaFold/PPI.

This is the first step toward a lovis4u-like neighbourhood visualization layer inside ppigFinder.
