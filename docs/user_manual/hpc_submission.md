# Generic HPC / Server Submission

## Purpose

ppigFinder is not tied to a specific HPC environment. It can prepare structural-prediction jobs for different execution targets, including generic SSH/Slurm servers, institution-specific HPC profiles, and manual AlphaFold 3 Server submission.

## Supported execution targets

ppigFinder supports three broad workflows:

1. **AlphaFold 3 Server export**
   - Generate AlphaFold-compatible JSON files.
   - Split large job lists into batches when needed.
   - Upload JSON files manually to a compatible AlphaFold 3 server.
   - Import completed result folders back into ppigFinder for analysis.

2. **Generic HPC / SSH / Slurm**
   - Generate job manifests and Slurm submission scripts.
   - Upload files to a remote Linux/HPC server when SSH/SFTP is configured.
   - Adapt partition, memory, CPU, wall-time and GPU/GRES directives to the target cluster.

3. **Preconfigured cluster profiles**
   - Use a named profile when available.
   - The DaVinci profile is provided as an institutional example.
   - Additional clusters should be represented as separate profiles rather than hard-coded into the core workflow.

## Generic Slurm profile

The default profile is `generic_slurm`. It does not assume a specific partition name, GPU model, GRES syntax or module system. Users should review the generated `.slurm` files before production submission.

Typical items to edit for a new cluster include:

- partition or queue name
- CPU count
- memory request
- wall-time limit
- GPU/GRES syntax
- module or conda activation commands
- AlphaFold 3 runner command
- database/model paths
- scratch/work directory layout

## DaVinci profile

`davinci` is an optional preconfigured profile for the DaVinci cluster. It maps ppigFinder resource hints to DaVinci Slurm partitions and GRES values. It should be treated as an example profile, not as a required runtime environment.

## AlphaFold 3 Server export

Users without HPC access can still use ppigFinder to generate AlphaFold 3 JSON files. These files can be submitted manually to a compatible AlphaFold 3 server, subject to that server's current limits and accepted JSON dialect.

## Recommended practice

Always inspect generated JSON, manifest and Slurm files before submission. Cluster policies vary, and resource directives that work on one system may need adjustment on another.
