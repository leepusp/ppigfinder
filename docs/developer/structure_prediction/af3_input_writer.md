# AF3 input writer

The AF3 input writer converts a backend-agnostic PredictionJobSpec into an initial AlphaFold 3 style JSON input.

Current output per AF3 job:

- input/af3_input.json

Current scope:

- protein targets only
- one AF3 job per PredictionJobSpec
- configurable model seeds
- optional top-level list wrapping

The writer intentionally rejects non-protein targets for now. DNA, RNA, ligands and covalent modifications should be added later with explicit schema handling.

This layer should be used after generic inputs are written. The generic files remain useful for traceability:

- input/targets.fasta
- input/job_spec.json
- input/af3_input.json

Before production submission, the generated AF3 JSON should be compared against the local AlphaFold 3 runner or server input requirements used on DaVinci.
