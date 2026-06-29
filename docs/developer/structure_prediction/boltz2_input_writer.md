# Boltz-2 input writer

The Boltz-2 input writer creates initial backend-specific input files for planned Boltz-2 jobs.

Current output per Boltz-2 job:

- input/boltz2_input.fasta
- input/boltz2_job_spec.yaml

This is intentionally conservative. The writer does not assume a final Boltz-2 CLI syntax yet. It creates reproducible intermediate inputs that can later be connected to the installed DaVinci Boltz-2 module or command.

The generic input files are still written separately:

- input/targets.fasta
- input/job_spec.json

Future work:

- detect the local Boltz-2 command/module;
- render a real Boltz-2 execution script;
- capture output locations into the manifest;
- compare Boltz-2 results with AF3 outputs.
