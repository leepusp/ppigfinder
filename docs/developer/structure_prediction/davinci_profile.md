# DaVinci Slurm profile

This profile converts generic ppigFinder resource hints into local DaVinci Slurm directives.

Internal hint mapping:

- gpu_shared_or_standard -> max50
- gpu_standard -> max50
- gpu_high_memory_preferred -> max90
- gpu_high_memory_required -> unrestricted
- cpu_or_backend_specific -> basic

GPU/GRES mapping:

- standard_gpu_or_slice -> shard:10
- one_standard_gpu -> shard:32
- one_high_memory_gpu -> shard:50
- largest_available_gpu_memory -> shard:50
- none_by_default -> no GRES directive

These values are defaults for script generation and previews. They should be reviewed against current sinfo and scontrol output before production runs.
