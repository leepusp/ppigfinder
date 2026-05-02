# DaVinci / HPC Execution Center

The DaVinci / HPC module is optional.

Users with access to DaVinci or another server can use this module to:

- configure SSH/local cluster access;
- detect local cluster mode when running directly on a DaVinci/GN node;
- prepare editable Slurm templates;
- choose a workflow target such as AlphaFold 3, GROMACS, BLAST+, HMMER3, CryoSPARC context or custom Python scripts;
- export scheduler scripts for manual review and execution.

Users without server access can still use ppigFinder normally and export:

- AlphaFold 3 Server JSON files;
- HTML reports;
- TSV/CSV tables;
- SVG/PDF vector figures;
- project snapshots.

Future report exports should favor interactive HTML for result exploration and vector graphics for publication-ready editing.
