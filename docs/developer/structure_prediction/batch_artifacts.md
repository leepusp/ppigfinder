# Structural prediction batch artifacts

The batch artifacts writer runs the complete local preparation step for a planned structural prediction batch.

It creates:

- output layouts
- job_manifest.tsv
- manifest validation report
- generic inputs
- AF3-specific inputs
- Boltz-2-specific inputs
- FoldCP-specific inputs

This layer still does not submit jobs. It prepares a reviewable batch directory that can later be connected to Slurm rendering and submission.

Conceptual workflow:

    build jobs
    plan batch
    write batch artifacts
    review manifest and inputs
    render Slurm scripts
    submit manually or through a controlled workflow

The artifact writer is intended to be the main backend entry point for future GUI actions.
