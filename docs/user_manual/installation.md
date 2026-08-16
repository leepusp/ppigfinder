# Installation and Dependencies

## Purpose

This section summarizes the software requirements needed to run ppigFinder.

## Python environment

ppigFinder is a Python desktop application. A supported environment should include:

- Python 3.8 or newer
- PyQt6 6.4 or newer, or PyQt5 5.15 or newer
- NumPy
- Matplotlib
- Pyrodigal
- Paramiko, when SSH/HPC submission is used

## External tools

Some analysis modules require external command-line tools:

- BLAST+ for similarity searches
- HMMER3 for profile-HMM searches
- AlphaFold 3 access, either through exported JSON files or a configured server/HPC workflow

## Backend check

After launching ppigFinder, inspect the Backends panel. Green check marks indicate that required tools were detected. Missing tools should be installed or added to the system PATH before running the corresponding analysis.

## Windows notes

On Windows, ppigFinder can be launched from a Python IDE such as Spyder or from a configured Python environment. BLAST+ and HMMER3 must still be available locally or through the configured workflow.

## Linux / HPC notes

On Linux or HPC systems, prefer a reproducible environment with explicit module or conda activation. For DaVinci-specific execution, see `davinci_hpc.md`.
