# Slurm rendering and output layout

This note documents the lightweight helpers added for structural prediction planning.

Files:

- ppigfinder/structure_prediction/output_layout.py
- ppigfinder/structure_prediction/slurm_renderer.py

The output layout standardizes prediction results by backend and job:

    <root>/<backend>/<job_id>/
      input/
      results/
      logs/
      manifests/
      retries/

The Slurm renderer converts an HPCResourcePlan into a reproducible sbatch header. It does not submit jobs. Submission should remain explicit and reviewable.

Planned use cases:

- generate AF3 JSON batches and Slurm arrays;
- generate Boltz-2 input batches and Slurm arrays;
- keep FoldCP as backend-specific CPU or support workflow;
- route large-token jobs to stronger GPU-memory profiles;
- keep failed-job retries organized under retries/.
