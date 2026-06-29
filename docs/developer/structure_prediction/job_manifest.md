# Structural prediction job manifests

The job manifest records one planned structural prediction per line.

It is designed to support:

- AF3 jobs
- Boltz-2 jobs
- FoldCP jobs
- token-aware Slurm arrays
- failed-job reprocessing
- reproducible output directories
- downstream structural analysis

Manifest columns include:

- job_id
- backend_id
- model_mode
- target_count
- estimated_tokens
- token_class
- partition_hint
- gpu_hint
- cpus_per_task
- memory_gb
- time_limit
- array_chunk_size
- job_dir
- input_dir
- result_dir
- log_dir
- retry_dir
- status

The manifest is intentionally backend-agnostic. Backend-specific input files can live under input/, while backend outputs should be collected under results/.
