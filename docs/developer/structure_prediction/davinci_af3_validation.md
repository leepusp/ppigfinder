# DaVinci AF3 wrapper validation

This document records the validation of the ppigFinder structural prediction AF3 path on DaVinci.

Validated chain:

    PredictionJobSpec
    -> input/targets.fasta
    -> input/job_spec.json
    -> input/af3_input.json
    -> DaVinci alphafold3/1.0 module
    -> /home/public/davinci/bin/af3 wrapper
    -> local dry-run preparation

Validated command pattern:

    module purge
    module load alphafold3/1.0

    af3 \
      --json-path <input/af3_input.json> \
      --job-name <job_id> \
      --workdir <workdir> \
      --stage all \
      --executor local \
      --dry-run \
      --force \
      --image "$AF3_IMAGE" \
      --model-dir "$AF3_MODEL_DIR" \
      --db-dir "$AF3_DB_DIR"

Observed dry-run result:

    Submitted/prepared: 1
    Failed: 0
    tokens: 280
    sequences: 3
    multimer: True
    executor: local

The wrapper generated:

    <workdir>/<job_id>/<job_id>_input.json
    <workdir>/<job_id>/run_<job_id>_all.sh

Implementation consequence:

The generated ppigFinder Slurm scripts should call the DaVinci AF3 wrapper with
`--executor local`, because the outer ppigFinder script is already an sbatch job.
This avoids nested Slurm submission from inside the allocated job.

Current status:

- AF3 JSON generation is compatible with the DaVinci wrapper in dry-run mode.
- Real inference was not executed during this validation.
- Boltz-2 and FoldCP remain conservative placeholders until their local commands are validated.
