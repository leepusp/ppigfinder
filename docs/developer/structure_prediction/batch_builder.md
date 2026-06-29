# Structural prediction batch builder

The batch builder creates one planned job per backend for the same target set.

Example conceptual flow:

    targets -> AF3 job
            -> Boltz-2 job
            -> FoldCP job

Each job receives:

- backend_id
- job_id
- estimated token count
- HPC resource plan
- output layout
- manifest row

This layer is still independent from the GUI. It prepares the structure needed for future buttons or workflows that submit AF3, Boltz-2 and FoldCP jobs from the same biological interaction candidate.

Recommended future use:

- build one multibackend batch for each ORF pair or multicomponent complex
- write one manifest per batch
- render one Slurm script per backend or per token class
- keep retries linked to the original manifest row
