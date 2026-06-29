# Structural prediction input writers

The generic input writer creates common input files for each planned structural prediction job.

For each job layout, it writes:

- input/targets.fasta
- input/job_spec.json

These files are backend-agnostic and are useful for traceability before converting into backend-specific inputs.

Backend-specific writers can later derive from this layer:

- AF3 JSON writer
- Boltz-2 input writer
- FoldCP input writer

The intended workflow is:

1. build PredictionJobSpec objects
2. plan resources
3. build output layouts
4. write manifests
5. write generic inputs
6. write backend-specific inputs
7. render Slurm scripts
8. submit or review manually
