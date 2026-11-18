# ppigFinder — User Guide

**Protein–Protein Interaction Genomic Finder · Documentation draft for v2.9.14**

> **Documentation status:** This guide is based on the 27-page preliminary guide by Gabriel U. Oka, together with the repository's v2.9.14 code and technical documentation reviewed on 18 November 2026. Screenshots in the preliminary guide show an earlier interface. Menu names may differ in newer builds; verify against the installed application. The generated file is a proposed documentation update, not evidence of a new software release.

## 1. Purpose, availability and scientific scope

ppigFinder is an open-source desktop application designed to help researchers investigate candidate protein–protein interactions in bacterial genomes. Its main workflow combines ORF prediction, homology and domain annotation, genomic-neighbourhood evidence, selection of plausible interaction partners, AlphaFold 3-compatible job preparation, and subsequent result interpretation.

The key scientific principle is **candidate prioritization**, not assuming that genomic adjacency proves physical interaction. Structural predictions and scores are computational evidence, not experimental confirmation. ppigFinder is distributed under the repository's MIT license. See the root `CITATION.cff` for the software authors and release metadata.

Repository: https://github.com/leepusp/ppigfinder

## 2. Install and launch

### 2.1 Supported setup

Use a Python environment. The packaging metadata specifies Python >=3.8, NumPy >=1.21 and Matplotlib >=3.5; PyQt6 >=6.4 is supplied as the `gui` extra. Pyrodigal and Paramiko are recommended optional packages. BLAST+ and HMMER3 are **external executables** and are not installed by pip with ppigFinder.

```bash
# From a checked-out repository:
python -m venv .venv
source .venv/bin/activate            # Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e '.[gui,recommended]'
python main.py
```

On Windows, the source can also be opened in Spyder, as shown in the original guide. Launching from `main.py` remains the documented stable entry point; do not assume that all experimental launchers are production ready. Install BLAST+ and HMMER3 separately when those analyses are required, and ensure the tools are accessible on `PATH`. A GUI workstation does not need AF3 weights or GPUs to prepare/export inputs.

**Backend indicators:** After startup, inspect the dependency/backend status panel. An available tool is typically marked in green. A missing optional backend does not prevent all other workflows from running; it disables or limits the corresponding operation.

### 2.2 Installation checks

```bash
python --version
python -c 'import numpy, matplotlib; print("Core Python packages OK")'
python -c 'import PyQt6; print("PyQt6 OK")'
blastp -version      # optional
hmmsearch -h         # optional
```

For graphics errors on a headless cluster, run the desktop GUI on a supported graphical session or workstation. CLI and batch components must be tested separately before claiming a fully headless workflow.

## 3. Tutorial dataset: *Xanthomonas citri* strain 306

The preliminary guide uses the chromosome reference **AE008923.1**. Obtain the nucleotide FASTA from the NCBI Nucleotide record; the source guide demonstrates the **Send to → File → FASTA → Create File** workflow. Record the accession and download date in the project metadata. Subsequent counts depend on genome version, prediction mode, software version and parameters.

### 3.1 Open the genome

1. Select **File → Open FASTA**.
2. Choose the downloaded nucleotide FASTA.
3. Inspect the **Genome** information panel for sequence length, GC content, number of predicted ORFs and annotation state.
4. Save the project using **File → Save Project** to preserve inputs and intermediate results. Verify that the project contains the genome data and project metadata before beginning large analyses.

The preliminary guide illustrates an initial ORF count of zero before translation. Avoid interpreting this as zero coding sequences in the biological genome.

### 3.2 Predict ORFs

Use **Translate genome** and choose **Pyrodigal**, **Hybrid** or **Automatic**. The preliminary tutorial uses Pyrodigal. The example run reported **4,303 predicted ORFs**, which is a tutorial observation and not a universal expected count. Record prediction method and thresholds in your analysis notes.

After completion, confirm that the ORF table is populated and the genome map renders. If table/map synchronization appears delayed, allow the display to finish and avoid repeated rapid zoom operations. The original guide describes reopening the saved project as a recovery action.

### 3.3 Explore the map and ORF table

- Select an ORF in the table to locate it on the map; select a map feature to inspect the related record.
- Use zoom and navigation controls carefully on large genomes.
- Right-click an ORF for contextual actions, including sequence extraction, annotation and interaction-candidate preparation. Exact action labels depend on the installed build.
- Export relevant protein/DNA sequences and retain the original ORF identifiers to preserve traceability.

## 4. Find proteins of interest

### 4.1 BLAST similarity search

Paste a protein query into the **BLAST Query** panel. Configure E-value and, where appropriate, scoring parameters. Perform the search against predicted ORFs; inspect hits and return to the genome map to view their local context. Similarity suggests homology or related function and does not by itself establish an interaction.

### 4.2 HMMER and domain profiles

Load one or more `.hmm` profiles (for example, curated profiles from InterPro/Pfam resources or profiles constructed from alignments). Run the search against predicted proteins, inspect sequence/domain E-values, and link the results back to the ORF table and map. The preliminary guide illustrates **Search all ORFs**, **Annotate HMM**, and separate sequence/domain reporting thresholds. Confirm the exact names in the installed interface.

### 4.3 Genomic neighbourhood

Inspect the ORFs surrounding a candidate. Adjacency, shared orientation, and putative operon-like context can improve candidate prioritization. These observations must be reported separately from evidence for a physical protein–protein interface.

## 5. Build AlphaFold 3 prediction jobs

In the AlphaFold tab or the ORF context menu, add candidate proteins to the selection. The preliminary guide documents these selection strategies:

| Mode | Use |
| --- | --- |
| Add to AF3 list | Build a queue without immediate prediction |
| Hit versus neighbours | Examine genomic proximity around a candidate |
| Selected versus selected | Select explicit pairs or a defined protein set |
| All versus all | Enumerate pairwise combinations in the selected set |
| Homodimer | Test a candidate's self-association |
| HMM-hit combinations | Prioritize sequences matching selected profiles |
| Custom complex | Define a multichain assembly and stoichiometry |

**Worked example:** The original manual selected **16 HMM-annotated T4SS candidates**, yielding **120 unordered non-self protein pairs** (`16 × 15 / 2`). Those 120 pairs are prediction jobs, not 120 confirmed physical interactions. Add homodimers separately where scientifically justified.

### 5.1 Export format matters

The **AlphaFold Server** JSON export and the **local AlphaFold 3** JSON dialect are related but are not necessarily interchangeable. Use the exporter appropriate to the intended backend and validate the JSON before submission. The repository has both server-JSON generation and a local AF3 JSON writer; the latter initially handles **protein targets only**, using the `alphafold3` dialect and model seed configuration.

The upstream AlphaFold Server may change account quotas and accepted batch sizes. **Do not rely on the preliminary manual's old numeric quotas**; consult the current official submission interface before uploading.

## 6. Run predictions: choose an execution target

ppigFinder is **not dependent on the DaVinci cluster**. The desktop application prepares predictions independently of where AF3 is installed. The actual inference requires a separately installed and licensed/authorized backend, appropriate model/database access, and suitable computational resources.

### 6.1 Option A — Web AlphaFold Server

Export server-compatible JSON, submit it manually through the server interface, follow the server's current policy, and download the outputs. Import the output folder into ppigFinder for interpretation. No SSH or Slurm connection is required.

### 6.2 Option B — Generic Linux HPC with Slurm

1. Prepare predictions locally and export an AF3 batch.
2. Select a target HPC profile; for unknown clusters start with `generic_slurm`.
3. Ask the cluster administrator for SSH address, user, optional jump host, work/scratch paths, partition, GPU request syntax, walltime policy, runtime environment and AF3 invocation.
4. Review generated job manifests and `.slurm` files, especially input paths, resources and environment initialization.
5. Transfer the prepared **entire job directory** to a cluster-accessible location. A generated script often contains absolute paths: generating locally and copying it without path adaptation will not work.
6. Connect using your institution's permitted SSH mechanism and run `sbatch <script>` **only after the above checks**. Save the returned job ID.
7. Monitor with `squeue -j <jobid>` and, where configured, `sacct -j <jobid>`; inspect stderr/stdout and final output files.
8. Transfer completed results back to the workstation and import them into ppigFinder.

> **Current implementation boundary:** The repository generates local AF3 inputs, manifests and Slurm arrays. Its modular `hpc_service.py`, `ssh_client.py` and `job_monitor.py` were found empty on 18 November 2026. Therefore, remote upload, submission and monitoring should not be described as a verified automatic end-to-end GUI feature. Legacy GUI code may contain older functionality; test it rather than assuming equivalence.

### 6.3 Option C — DaVinci (optional example)

DaVinci is a named institutional profile, not a requirement. The repository maps resource hints to partitions such as `max50` and `max90` and to site-specific GPU/GRES values. These are **illustrative defaults** and need an administrator check before submission. A documented dry-run validated compatibility of one local-AF3 JSON with the DaVinci wrapper, but **did not execute real inference**.

The generated Slurm job calls `af3 --executor local` **inside an already allocated Slurm job** so that the AF3 wrapper does not submit a nested job. The `af3` executable must be available on the target cluster; when necessary, initialize it through the cluster's module/conda/container configuration.

See [HPC portability and integration](../developer/hpc_portability.md) and the existing [DaVinci validation note](../developer/structure_prediction/davinci_af3_validation.md).

## 7. Interpret and export structural results

Use the application's result-import workflow to read completed AF3 output directories. The repository contains an `AlphaFoldResultsService` with sequential/parallel parsing and optional caching. Report the input pair IDs and chain mapping alongside result paths, model seeds, ranking and confidence metrics.

Typical reported quantities include **ipTM**, **pTM**, interchain PAE, contact percentages and interaction classifications. Do not claim that a high score proves biological binding; assess interfaces, chain-specific confidence, sequence evidence and experimental support. Keep track of failed or incomplete jobs separately from successfully parsed predictions.

Export TSV/CSV result tables, FASTA files, project snapshots and HTML reports as needed. Keep the exact ppigFinder release/commit and external backend versions in the publication record.

## 8. Troubleshooting

| Problem | Check |
| --- | --- |
| GUI does not start | Python/Qt installation, display/X11, application tracebacks |
| BLAST/HMM unavailable | Independent binary installation, PATH, permissions |
| ORF map slow or unsynchronized | Wait for rendering, save/reopen project, reduce rapid zoom |
| Incorrect AF3 JSON | Export dialect, sequence validity, chain IDs, model seeds |
| `af3: command not found` | Remote module/conda/container and `PATH` |
| Slurm rejects job | Partition, account, QoS, GRES, walltime, node resources |
| Job cannot find inputs | Local vs remote absolute paths; transfer and layout |
| AF3 out-of-memory | Adjust backend resources, shorten/separate complexes, retry |
| Import finds no results | Exact prediction directory, completed backend stage, file integrity |

## 9. Reproducibility checklist

Before communicating results, record: genome accession and version; prediction method and thresholds; BLAST database/version; HMM profiles/version; candidate-selection rationale; proteins/chains and stoichiometry; JSON dialect; model seeds; cluster and scheduler details; AF3 version/environment; result timestamps; missing/failed predictions; software commit and source license.

## 10. Source-to-guide map

This user guide preserves the tutorial sequence from **Gabriel U. Oka's 27-page preliminary manual**:

- Pages 1–7: purpose, requirements, Windows/Spyder launch and backend indicators → Sections 1–2.
- Pages 8–15: *X. citri*, FASTA, project saving, ORFs and map navigation → Section 3.
- Pages 16–19: BLAST and HMM annotation → Section 4.
- Pages 20–25: AF3 pair/complex builders and 16→120 candidate-screen example → Section 5.
- Page 26: server export → Section 6.1, with numeric quota claims deliberately omitted pending official verification.
- Page 27: incomplete DaVinci heading → Sections 6.2–6.3, explicitly generalized to other clusters.

Original screenshots should be version-labeled when added to the repository, and any new screenshot must be captured from the actual released GUI. Do not represent old screenshots as verified v2.9.14 layouts.
