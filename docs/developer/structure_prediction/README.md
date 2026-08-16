# Structural prediction backends

This layer prepares ppigFinder to support multiple structural prediction backends:

- AlphaFold 3 / AF3
- Boltz-2
- FoldCP
- token-aware job partitioning
- Slurm array planning
- CPU/GPU resource planning
- failed-job reprocessing
- reproducible result organization

The GUI should not contain backend-specific hardcoded logic. The GUI should call small services that produce prediction job specifications, token estimates, HPC resource plans, backend-specific input files and output manifests.

Initial module:

    ppigfinder/structure_prediction/

Current files:

    models.py
    token_budget.py
    hpc_planner.py

The resource planner separates CPU preprocessing, GPU inference and CPU postprocessing. This does not assume CPU + GPU always accelerate a single model simultaneously. Instead, it supports staged and pipelined resource use across many prediction jobs.

Planned backend roles:

- AF3: high-confidence structural prediction and multicomponent complex modeling.
- Boltz-2: complementary prediction backend for comparison of structural hypotheses.
- FoldCP: optional structural comparison or prediction-support backend, depending on local deployment.

HPC planning goals:

- choose CPU or GPU partitions according to backend and token load;
- keep Slurm arrays reproducible;
- route large-token jobs to stronger GPU-memory profiles;
- retry failed jobs with stronger memory/GPU settings;
- organize outputs into auditable per-job directories.
