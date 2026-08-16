# FoldCP input writer

The FoldCP input writer creates initial backend-specific input files for planned FoldCP jobs.

Current output per FoldCP job:

- input/foldcp_input.fasta
- input/foldcp_job_spec.yaml

This is intentionally conservative. FoldCP may require a structure-comparison workflow rather than a pure sequence prediction workflow depending on the local installation and use case.

The FASTA and YAML files provide traceability for the biological target set even if later execution uses predicted structures from AF3 or Boltz-2 as FoldCP inputs.

Future work:

- validate the local FoldCP command and accepted input types;
- define whether FoldCP consumes sequences, structures, or both in this workflow;
- connect AF3/Boltz-2 outputs to FoldCP comparison jobs;
- record FoldCP result paths back into the manifest.
