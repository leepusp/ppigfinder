# ppigFinder Directory Layout

## Root files

| Path | Purpose |
|---|---|
| `main.py` | Stable entry point for the current application. |
| `pyproject.toml` | Python packaging/build configuration and future command entry-points. |
| `environment.yml` | Conda environment definition. |
| `README.md` | Short project overview. |
| `ARCHITECTURE.md` | High-level technical architecture notes. |
| `docs/` | User manuals, workflow documentation and developer notes. |
| `legacy_sources/` | Historical standalone scripts and archived upstream legacy sources preserved for provenance. |

## Main package layout

| Package | Purpose |
|---|---|
| `ppigfinder/legacy_v29_14.py` | Current active monolithic GUI kept during the refactor. |
| `ppigfinder/legacy_v20.py` | Historical packaged legacy GUI retained for provenance and compatibility checks. |
| `ppigfinder/app.py` | Current desktop launcher/bootstrap. |
| `ppigfinder/ui/` | Current Qt adapters, modular actions, dialogs and widgets. |
| `ppigfinder/ui_shell/` | Experimental guided interface: splash, home and modular workspace. |
| `ppigfinder/services/` | Application orchestration layer. |
| `ppigfinder/domain/` | Serializable data models. |
| `ppigfinder/io/` | File readers, writers and exporters. |
| `ppigfinder/bioseq/` | Sequence, ORF and translation logic. |
| `ppigfinder/annotation/` | BLAST, HMMER, fallback and domain annotation logic. |
| `ppigfinder/alphafold/` | AF3 JSON, parser, metrics and classifier logic. |
| `ppigfinder/hpc/` | HPC, SSH and scheduler integrations. |
| `ppigfinder/infrastructure/` | Cache, subprocess, logging, backend detection and parallelism. |
| `ppigfinder/visualization/` | Layout/data preparation for future visualizations. |
| `ppigfinder/resources/` | Static resources, translations and help content. |
| `ppigfinder/config/` | Defaults, thresholds and runtime constants. |

## Interface separation

Current stable GUI:

```text
main.py
  -> ppigfinder.app
      -> ppigfinder.legacy_v29_14.ppigFinderApp
```

## Legacy source archive

Historical standalone scripts and upstream source snapshots are preserved under
`legacy_sources/`. They are not the active application entry point, but they are
kept to support reproducibility, provenance tracking and comparison against
earlier ppigFinder versions.
