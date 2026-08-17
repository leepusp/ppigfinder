# ppigFinder Software Architecture

ppigFinder is currently being refactored from a monolithic desktop application into a modular architecture while preserving the functional Qt interface.

## Current stable entry point

```text
main.py
  -> ppigfinder.app
      -> ppigfinder.legacy_v29_14.ppigFinderApp
```

The file `ppigfinder/legacy_v29_14.py` contains the current complete working interface. It is preserved during the refactor and should not be modified unless necessary. Older standalone and packaged legacy sources are retained for provenance and compatibility checks.

## Experimental guided interface

The experimental interface shell can be launched with:

```bash
python -m ppigfinder.ui_shell.launcher
```

It currently provides:

- splash screen
- home/start screen
- guided workspace preview
- module-oriented navigation
- embedded documentation panels

This interface does not replace the current stable interface yet.

## Major modules

| Package | Responsibility |
|---|---|
| `domain/` | Data models for genome, ORF, BLAST, HMM, AlphaFold and project state. |
| `services/` | Orchestration layer connecting domain logic, IO, infrastructure and GUI bridges. |
| `io/` | File readers, writers and exporters. |
| `bioseq/` | Sequence manipulation, ORF prediction, translation and Pyrodigal integration. |
| `annotation/` | BLAST, HMMER, fallback search and domain annotation logic. |
| `alphafold/` | AF3 JSON generation, job building, result parsing, metrics and classification. |
| `hpc/` | SSH, transfer and scheduler abstractions for SLURM, PBS and LSF workflows. |
| `infrastructure/` | Backend detection, caching, subprocess execution, logging and parallelism. |
| `visualization/` | Layout and data preparation for genome maps, PAE, pLDDT and PPI visualization. |
| `ui/` | Current Qt add-ons, widgets, dialogs and modular adapters. |
| `ui_shell/` | Experimental guided interface and future modular workspace. |

## Backend detection

At startup, ppigFinder can detect optional tools such as:

- BLAST+
- HMMER3
- Pyrodigal
- Paramiko / SSH support

The application should remain functional when optional tools are unavailable by using internal fallback implementations where feasible.

## Background execution

Computationally intensive tasks should run outside the main GUI thread. The current architecture uses Qt worker threads and modular helper services for long-running tasks such as:

- ORF prediction
- BLAST or fallback search
- HMM scanning
- AF3 result parsing
- file export

## Refactor principle

The refactor should progressively move logic out of `legacy_v29_14.py` into smaller modules:

```text
legacy GUI event
  -> ui bridge
      -> service
          -> domain / io / backend module
```

New code should avoid importing Qt unless it belongs in `ui/` or `ui_shell/`.

## Future direction

The long-term interface direction is a modular, guided workflow:

```text
Start
  -> Data / Project
  -> DNA / Genome
  -> Protein / ORFs
  -> Annotation
  -> AlphaFold / PPI
  -> Reports
```

Future visualization modules may use HTML/D3/Observable-style components embedded in Qt through web views, while keeping Python as the scientific backend.
