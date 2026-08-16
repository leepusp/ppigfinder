# Troubleshooting

## Backend not detected

Check whether BLAST+, HMMER3, Pyrodigal and Python dependencies are installed and available in the active environment. Restart ppigFinder after changing PATH, modules or conda environments.

## Genome does not load

Confirm that the input file is a valid nucleotide FASTA, multi-FASTA, GenBank or supported project file. Very large files may take longer to render in the genome map.

## ORFs do not appear

Run an ORF prediction mode after loading the genome. If the map does not update, wait for rendering to finish or reopen the saved project.

## BLAST or HMM search fails

Confirm that the query sequence or HMM file is valid and that the external tool is available. Review E-value thresholds when no hits are detected.

## AlphaFold 3 JSON is rejected

Check whether the exported JSON matches the intended submission target. Public AlphaFold servers and private/HPC wrappers may enforce different limits, naming rules or batch sizes.

## AF3 result has missing metrics

Some AF3 folders may be incomplete or may lack files such as `model.cif`, `summary.json`, `confidences.json` or `ranking_scores.csv`. ppigFinder can still parse partial results, but plotting and hotspot analysis require the heavier confidence files.

## PPI result interpretation

Treat ipTM, PAE_min, PAE_inter, HotSpotPAE and Contact% as structural evidence, not as experimental proof. Candidate interactions should be interpreted together with genomic context, annotation evidence and biological knowledge.
