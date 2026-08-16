# Structural prediction Slurm batch writer

The Slurm batch writer generates reviewable Slurm scripts and task tables from a PredictionBatchPlan.

It writes files under:

    <batch_root>/slurm/

For each backend and token class group, it creates:

    tasks_<backend>_<token_class>.tsv
    submit_<backend>_<token_class>.slurm

The scripts are intentionally conservative. They do not yet execute final AF3, Boltz-2 or FoldCP commands. Instead, they load or mark the backend environment and print the exact input paths that should be consumed.

This separates three phases:

1. planning and artifact generation;
2. Slurm script rendering;
3. final backend command integration and submission.

The DaVinci profile resolves generic resource hints into local Slurm directives such as partition and GRES.
