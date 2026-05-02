# DaVinci / HPC Module

## Purpose

The DaVinci / HPC module provides an optional connection point for running or preparing project workflows on the DaVinci server or compatible HPC environments.

## Why this module exists

ppigFinder can generate data and jobs that may become computationally demanding, especially AlphaFold 3 pairwise or multichain prediction batches. The HPC module is designed to make this part of the workflow explicit and graphical.

## Supported concepts

- detect whether ppigFinder is already running inside a DaVinci/GN node
- configure SSH host, user, port and key path
- test SSH connectivity without freezing the GUI with password prompts
- prepare Slurm templates for AF3 workflows
- keep HPC execution optional in the guided workflow

## Recommended workflow

1. Prepare ORFs and candidates.
2. Build AF3 jobs.
3. Open DaVinci / HPC.
4. Test connection or detect local cluster mode.
5. Prepare Slurm submission template.
6. Run or submit jobs according to the configured environment.
7. Import AF3 results back into ppigFinder.

## Notes

The guided shell currently starts with safe connection testing and template generation. Full submission and monitoring will be connected progressively through the service layer.
