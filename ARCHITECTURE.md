# ppigFinder Architecture

ppigFinder is organized by scientific and operational responsibility.

## Main layers

### apps/
Application entry points:
- desktop.py: graphical PyQt application
- cli.py: command-line interface
- batch.py: non-interactive batch execution

### domain/
Core data models shared by the entire application:
- Genome
- ORF
- BlastHit
- HMMHit
- DomainHit
- AF3Job
- InteractionResult
- ProjectState

### bioseq/
Biological sequence logic:
- DNA/protein sequence utilities
- genetic code
- FASTA parsing/writing
- ORF prediction
- Pyrodigal integration
- hybrid ORF prediction

### annotation/
Functional annotation:
- BLAST+
- built-in BLAST fallback
- HMMER3
- built-in PSSM fallback
- domain profiles
- annotation orchestration

### neighborhood/
Genomic context analysis:
- neighboring ORF windows
- operon-like context
- neighborhood exports

### alphafold/
AlphaFold 3 workflows:
- pair generation
- AF3 JSON generation
- ColabFold FASTA generation
- batch partitioning
- AF3 result parsing
- PAE/ipTM/contact metrics
- motif detection
- interaction classification

### hpc/
Remote cluster integration:
- SSH/SFTP
- file transfer
- SLURM/PBS/LSF schedulers
- job monitoring
- cluster profiles such as DaVinci

### io/
File input/output:
- FASTA
- GenBank
- SnapGene
- project JSON
- TSV/CSV export
- reports

### visualization/
Plot and layout logic independent from the main PyQt window:
- genome map layout
- PPI arc layout
- PAE heatmaps
- pLDDT plots
- color palettes

### ui/
PyQt interface:
- main window
- menus/actions
- dialogs
- tabs
- reusable widgets

### services/
Workflow orchestration:
- genome loading
- ORF prediction
- annotation
- AlphaFold job creation and analysis
- HPC submission
- project save/load

### infrastructure/
System-level integration:
- external backend detection
- subprocess calls
- platform-specific paths
- logging

### resources/
Static resources:
- help content
- translations
- icons

### config/
Constants and defaults:
- default parameters
- thresholds
- application constants

## Refactor rule

New code should not be added to legacy_v20.py unless it is a temporary bridge.
Every extracted feature should move toward one of the modules above.
