#!/usr/bin/env python3
"""
ppigFinder — Protein-Protein Interaction Genomic Finder
========================================================
Version  : 2.00 — v2.0
Released : 2026
License  : MIT (see LICENSE section below)

DESCRIPTION
-----------
ppigFinder is a standalone desktop application for the discovery of novel
protein-protein interactions (PPIs) in bacterial genomes.  Starting from a
raw nucleotide sequence (FASTA, GenBank or SnapGene), the tool identifies
all open reading frames (ORFs), annotates them through homology search
(BLASTp) and Hidden Markov Model (HMM) profile scanning, visualises them
on an interactive genomic map, and generates AlphaFold 3 (AF3) batch jobs
for structural interaction prediction.  Jobs can be submitted to any
remote HPC cluster (SLURM, PBS, LSF) directly from the GUI via SSH/SFTP.

MAIN FEATURES
-------------
  • 6-frame ORF prediction (start/stop codon search or Pyrodigal/Prodigal)
  • BLASTp homology search against all predicted ORFs (local BLAST+ or
    built-in Smith-Waterman / k-mer filter)
  • HMM profile search using HMMER3 or built-in PSSM scanner
    (compatible with Pfam, TIGRFAM and custom .hmm profiles)
  • Interactive zoomable genomic map with colour-coded ORF arrows
  • Genomic neighbourhood analysis (configurable window, FASTA export)
  • AlphaFold 3 job builder: AF3 JSON / ColabFold FASTA export
  • HPC server integration: SSH/SFTP upload and submission to SLURM,
    PBS/Torque or LSF schedulers; job monitoring and result download
  • AF3 results analysis (Interaction Results tab):
      – PAE heatmap (ChimeraX colour scheme) and pLDDT plots
      – Global and focal inter-chain metrics: ipTM, ptm, ranking_score,
        PAE_inter (global mean), PAE_min ★ (chain_pair_pae_min off-diag),
        cp_ipTM ★ (chain_pair_iptm off-diag), Contact% ★ (PAE < 5 Å)
      – Two-tier interaction motif detector: connected-component
        segmentation of the off-diagonal PAE quadrant validated by
        contact_probs and B↔A reciprocity (novel v1.16 algorithm)
      – Confidence classification: HIGH (PAE_min < 4 Å & cp_ipTM ≥ 0.50),
        MED (PAE_min 4–8 Å), LOW (PAE_min ≥ 8 Å)
  • Genomic PPI Map tab: linear genome representation with Bezier arcs
    connecting predicted interaction pairs, colour-coded by PAE_min;
    zoom/pan, click-to-select, SVG/TSV export
  • Project save / load (JSON workspace), multi-language UI (EN/PT/ES/
    FR/ZH/JA), export to GenBank, SnapGene .dna, PDF, TSV

DEPENDENCIES
------------
  Required (install via pip):
    Python       >= 3.8
    PyQt6        >= 6.4   (or PyQt5 >= 5.15 as automatic fallback)
                          https://www.riverbankcomputing.com/software/pyqt/
    matplotlib   >= 3.5   https://matplotlib.org
                          Hunter, J.D. (2007). Matplotlib: A 2D graphics
                          environment. Computing in Science & Engineering,
                          9(3), 90-95. DOI:10.1109/MCSE.2007.55
    numpy        >= 1.21  https://numpy.org
                          Harris, C.R. et al. (2020). Array programming with
                          NumPy. Nature, 585, 357-362.
                          DOI:10.1038/s41586-020-2649-2

  Optional — automatically detected at runtime:
    pyrodigal    >= 2.0   pip install pyrodigal
                          (ML-based prokaryotic gene caller; fallback to
                          built-in codon-based ORF finder if absent)
                          [1] Larralde, M. (2022). Pyrodigal: Python bindings
                          and interface to Prodigal. J. Open Source Softw.,
                          7(72), 4296. DOI:10.21105/joss.04296
                          [2] Hyatt, D. et al. (2010). Prodigal: prokaryotic
                          gene recognition and translation initiation site
                          identification. BMC Bioinformatics, 11, 119.
                          DOI:10.1186/1471-2105-11-119

    NCBI BLAST+  >= 2.12  https://ftp.ncbi.nlm.nih.gov/blast/executables/
                          (local installation; also detectable via WSL on
                          Windows; fallback to built-in k-mer/SW aligner)
                          [3] Camacho, C. et al. (2009). BLAST+: architecture
                          and applications. BMC Bioinformatics, 10, 421.
                          DOI:10.1186/1471-2105-10-421

    HMMER3       >= 3.3   http://hmmer.org
                          conda install -c bioconda hmmer
                          (local or via WSL; fallback to built-in PSSM scanner)
                          [4] Eddy, S.R. (2011). Accelerated Profile HMM
                          Searches. PLoS Comput. Biol., 7(10), e1002195.
                          DOI:10.1371/journal.pcbi.1002195

    paramiko     >= 2.9   pip install paramiko
                          (SSH/SFTP; required for HPC server submission)
                          https://www.paramiko.org

  AlphaFold 3 structural prediction (external service / local install):
                          [5] Abramson, J. et al. (2024). Accurate structure
                          prediction of biomolecular interactions with
                          AlphaFold 3. Nature, 630, 493-500.
                          DOI:10.1038/s41586-024-07487-w

REFERENCES
----------
  [1] Larralde, M. (2022). Pyrodigal: Python bindings and interface to
      Prodigal, an efficient method for gene prediction in prokaryotes.
      Journal of Open Source Software, 7(72), 4296.
      https://doi.org/10.21105/joss.04296

  [2] Hyatt, D., Chen, G.-L., LoCascio, P.F., Land, M.L., Larimer, F.W.,
      & Hauser, L.J. (2010). Prodigal: prokaryotic gene recognition and
      translation initiation site identification.
      BMC Bioinformatics, 11, 119.
      https://doi.org/10.1186/1471-2105-11-119

  [3] Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J.,
      Bealer, K., & Madden, T.L. (2009). BLAST+: architecture and
      applications. BMC Bioinformatics, 10, 421.
      https://doi.org/10.1186/1471-2105-10-421

  [4] Eddy, S.R. (2011). Accelerated Profile HMM Searches.
      PLoS Computational Biology, 7(10), e1002195.
      https://doi.org/10.1371/journal.pcbi.1002195

  [5] Abramson, J., Adler, J., Dunger, J., Evans, R., Green, T.,
      Pritzel, A., ... Jumper, J.M. (2024). Accurate structure prediction
      of biomolecular interactions with AlphaFold 3.
      Nature, 630(8016), 493-500.
      https://doi.org/10.1038/s41586-024-07487-w

  [6] Hunter, J.D. (2007). Matplotlib: A 2D graphics environment.
      Computing in Science & Engineering, 9(3), 90-95.
      https://doi.org/10.1109/MCSE.2007.55

  [7] Harris, C.R., Millman, K.J., van der Walt, S.J. et al. (2020).
      Array programming with NumPy. Nature, 585, 357-362.
      https://doi.org/10.1038/s41586-020-2649-2

  [8] Evans, R., O'Neill, M., Pritzel, A., Antropova, N., Senior, A.,
      Green, T., ... Hassabis, D. (2022). Protein complex prediction with
      AlphaFold-Multimer. bioRxiv 2021.10.04.463034.
      https://doi.org/10.1101/2021.10.04.463034
      [ipTM metric, chain_pair_iptm, chain_pair_pae_min definitions]

  [9] Jumper, J., Evans, R., Pritzel, A., Green, T., Figurnov, M.,
      Ronneberger, O., ... Hassabis, D. (2021). Highly accurate protein
      structure prediction with AlphaFold. Nature, 596, 583-589.
      https://doi.org/10.1038/s41586-021-03819-2
      [PAE (Predicted Aligned Error) metric definition]

  [10] Humphreys, I.R., Pei, J., Baek, M., Krishnakumar, A., Anishchenko, I.,
       Ovchinnikov, S., ... Baker, D. (2021). Computed structures of core
       eukaryotic protein complexes. Science, 374, eabm4805.
       https://doi.org/10.1126/science.abm4805
       [ipTM-based PPI screening at proteome scale; methodology basis]

  [11] Bryant, P., Pozzati, G., & Elofsson, A. (2022). Improved prediction
       of protein-protein interactions using AlphaFold2. Nature Communications,
       13, 1265. https://doi.org/10.1038/s41467-022-28865-w
       [PAE inter-chain analysis for PPI confidence scoring]

  [12] Mirdita, M., Schütze, K., Moriwaki, Y., Heo, L., Ovchinnikov, S.,
       & Steinegger, M. (2022). ColabFold: making protein folding accessible
       to all. Nature Methods, 19, 679-682.
       https://doi.org/10.1038/s41592-022-01488-1
       [ColabFold FASTA format; batch MSA for PPI screening]

  [13] Virtanen, P., Gommers, R., Oliphant, T.E. et al. (2020). SciPy 1.0:
       Fundamental algorithms for scientific computing in Python.
       Nature Methods, 17, 261-272. https://doi.org/10.1038/s41592-019-0686-2
       [scipy.ndimage used in motif detection (connected-component labelling)]

  [14] Chou, T.-F., Bulfer, S.L., Weihl, C.C., Li, K., Leman, L.J., Ghadiri,
       M.R., ... Deshaies, R.J. (2011). Specific inhibition of p97/VCP ATPase
       and kinetic analysis demonstrate interaction between D1 and D2 ATPase
       domains. Journal of Molecular Biology, 406(3), 432-450.
       [Genomic co-localization as PPI evidence in bacterial operons]

  METHODOLOGY NOTE — Focal Interaction Metrics (v1.18):
  The global PAE_inter metric (mean of the entire off-diagonal PAE quadrant)
  is diluted by disordered protein regions and does not reliably identify
  proteins that interact through a single domain. ppigFinder v1.18+
  introduces two focal metrics extracted directly from AF3 summary JSONs:
    • PAE_min (chain_pair_pae_min): minimum PAE in the off-diagonal quadrant.
      A value < 4 Å indicates that AF3 predicts at least one contact point
      with near-atomic confidence, regardless of global disorder.
    • cp_ipTM (chain_pair_iptm): chain-pair ipTM for the interface only,
      removing the contribution of intra-chain folding quality.
  Classification threshold: PAE_min < 4 Å AND cp_ipTM ≥ 0.50 = HIGH
  confidence focal interaction (analogous to the ipTM > 0.75 criterion of
  Evans et al. 2022 but applicable to domain-limited contacts).

INSTALLATION
------------
  pip install PyQt6 matplotlib numpy            # core
  pip install pyrodigal paramiko                # optional recommended

  Then run:
    python ppigfinder.py

SUPPORTED FILE FORMATS
----------------------
  Input  : FASTA (.fasta .fa .fna), GenBank (.gb .gbk), SnapGene (.dna)
  Output : FASTA, GenBank, SnapGene .dna, PDF/EPS, TSV, JSON (AF3),
           FASTA (ColabFold), ppigFinder project (.json)

LICENSE
-------
  MIT License

  Copyright (c) 2026 ppigFinder Contributors

  Permission is hereby granted, free of charge, to any person obtaining a
  copy of this software and associated documentation files (the "Software"),
  to deal in the Software without restriction, including without limitation
  the rights to use, copy, modify, merge, publish, distribute, sublicense,
  and/or sell copies of the Software, and to permit persons to whom the
  Software is furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included
  in all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
  DEALINGS IN THE SOFTWARE.

  ATTRIBUTION REQUIREMENT (MIT with Credit Clause):
  Any publication, software, or derivative work that uses ppigFinder or
  substantial portions of its source code must give appropriate credit to
  the original authors by retaining this notice and/or citing the original
  repository.

CITATION
--------
  If you use ppigFinder in your research, please cite:
    ppigFinder: Protein-Protein Interaction Genomic Finder (2026).
    https://github.com/goka-lab/ppigfinder

CHANGELOG
---------

  v2.0 (2026) — Public release
    • Consolidated all v1.x features into stable public release.
  v1.16 (2026) — Interaction motif detection in the off-diagonal PAE
    • New algorithm: two-tier connected-component segmentation of the
      off-diagonal PAE quadrants, validated by contact_probs and
      reciprocity. Extracts biologically interpretable interaction
      motifs (A[res X-Y] × B[res W-Z]) from each AF3 pair.
        – Tier 1 (core):      PAE < pae_core AND contact_probs ≥ min_contact
        – Tier 2 (extended):  PAE < pae_ext (used to grow the bbox)
        – Morphological closing bridges 1-2 pixel gaps before labelling.
        – Reciprocity: motif must also appear in the B↔A quadrant at
          ≥ 50 % pixel overlap (rejects one-sided ghost motifs).
        – Score 0-100 combining PAE confidence (35 %), size (20 %),
          density (20 %), pLDDT (15 %), reciprocity (10 %).
    • Per-motif metrics surfaced: residue ranges on both chains, size,
      mean / min PAE, density, mean / max contact_probs, pLDDT A / B,
      reciprocity overlap, combined score.
    • scipy.ndimage is used when available for labelling + morphology;
      a pure-numpy BFS + dilate/erode fallback is used otherwise.
    • UI: second toolbar row in the AlphaFold Analysis tab with
      spinboxes for core PAE, ext PAE, min contact_probs, min size,
      'reciprocal' checkbox, Rerun button, Motifs TSV export.
    • UI: motif table appended below PAE/pLDDT plots on row select,
      with colour-coded score / pLDDT cells. Clicking a motif row
      flashes a red highlight rectangle on both off-diagonal quadrants.
    • UI: PAE heatmap overlay now draws a coloured numbered rectangle
      (green / amber / red by score) around every detected motif in the
      A↔B quadrant and a dashed mirror in B↔A.
    • Validated against the uploaded ORF2601 × ORF50 example:
        – synthetic motif at A[30-44] × B[50-69] → detected exactly
          there with mean PAE 3.2 Å, 100 % reciprocity.
        – real ORF2601 × ORF50 → 0 motifs (correct: chain_pair_pae_min
          is 8.70 Å, no contiguous ≥ 5×5 region below 8 Å).
    • Defaults calibrated from real data: pae_core=8, pae_ext=15,
      min_contact=0.05, min_size=5.

  v1.15 (2026) — AlphaFold Analysis tab: bug fixes + AF3-server support
    — Bugs fixed —
    • Critical: selecting a row after sorting the table loaded the
      WRONG job into the PAE/pLDDT plots (view row index was used as a
      data-list index). Index is now stored in UserRole on column 0
      and retrieved via item.data(UserRole).
    • Critical: numeric columns (ipTM, ptm, pLDDT, ranking, PAE_inter)
      were sorted lexicographically ("100.0" < "2.5"). Fixed with a
      _NumericItem subclass with numeric __lt__.
    • Bug: pLDDT extraction missed the AF3 server's per-atom layout.
      AF3 server output stores 'atom_plddts' (per-atom, 9132 values
      for this pair) with NO 'token_plddts'. Parser now aggregates
      atoms → residues using the atom→residue map parsed from the
      model .cif file (whose B_iso column equals the atom pLDDT).
      Verified against the uploaded example: 1170 tokens, 0 missing.
    • Bug: contact markers on the PAE heatmap were drawn as a single
      vertical stripe at the centre of each off-diagonal block.
      Markers now appear at the actual contact residue positions in
      both AB and BA quadrants.
    • Bug: changing the contact threshold did not refresh the "Best
      contact pair" cell in the table. Fixed.
    • Bug: parse errors were silently swallowed to stdout. Failed
      jobs are now collected and surfaced in a single dialog.

    — AF3 server output support —
    • Content-based JSON classification replaces the old filename
      globs. Each .json in a candidate folder is probed for its top-
      level keys and tagged as summary / confidences / data / unknown.
    • Huge input-data JSONs (sequences + MSA, 80+ MB) are skipped by
      size before being opened.
    • seed-*_sample-* subfolders are recognised as per-sample details
      of their parent job (the best model is already promoted to the
      parent) and are no longer scanned as independent jobs.
    • Rich summary fields now extracted: chain_iptm, chain_ptm,
      chain_pair_iptm, chain_pair_pae_min, fraction_disordered,
      has_clash, contact_probs.
    • ranking_scores.csv is parsed and attached to the result; the
      number of diffusion samples is shown in the row tooltip and in
      the TSV export.

    — UX —
    • Filter spinboxes: min ipTM, max PAE_inter.
    • Tooltip on row name cell shows clash, disordered %, pair PAE
      min, per-chain ipTM, and sample count.
    • Plot header now shows disordered %, clash status, and per-chain
      ipTM/pTM in a second line. Per-pair line shows pair ipTM, pair
      PAE min, and pair PAE mean alongside the contact region.
    • TSV export expanded with pae_min_best_pair, pair_iptm_best_pair,
      chain_iptm, chain_ptm, fraction_disordered, has_clash and
      n_diffusion_samples.
    • Double-click a row to open the job folder in the OS file
      manager. High-confidence rows (ipTM > 0.75 and PAEinter < 8 Å)
      are drawn in bold. Progress dialog for large scans. Top-scoring
      row auto-selected after load.

    — AF3 tab Mode dropdown —
    • Removed 'HMM Hits vs Each Other', 'Hit vs All Selected ORFs' and
      'Homodimer (Hit vs Itself)' at user request.
    • Implemented 'Neighbors Interactome' (was silently generating 0
      jobs). Now scans the ENTIRE genome with a sliding window of ±N
      neighbors and deduplicates by (min,max) — produces the canonical
      pair set { (i,j) : 1 ≤ j−i ≤ N }, whose size is N·n − N(N+1)/2.
      Confirms before generating > 5 000 jobs. Validated against the
      user-provided example: for genome={1..6}, N=2 yields exactly
      (1,2),(1,3),(2,3),(2,4),(3,4),(3,5),(4,5),(4,6),(5,6).
    • Implemented 'Trimers (Hit + 2 Neighbors)' (was also silently
      empty). Emits 3-chain AF3 jobs with the Hit and every unordered
      pair of distinct neighbors in the window.

  v1.12 (2026) — File menu cleanup + extended zoom to 1,000,000×
    • Removed from File menu: Open Multi-FASTA, Open SnapGene (.dna),
      Open GenBank (.gb/.gbk), Export as SnapGene (.dna), Export as GenBank
    • Genome map zoom ceiling raised from 200× (20 000%) to 10 000×
      (1 000 000%) to support navigation of large bacterial genomes and
      metagenome-assembled genomes (MAGs > 5 Mbp)
    • Manual, Tutorial, and About dialog fully updated to reflect v1.12
      feature set (FASTA-only input, extended zoom, Download JSONs)
    • ORF-selection-centering logic and map pan/anchor behaviour unchanged

  v1.11 (2026) — AF3 large-batch OOM fix + Download JSONs button
    • "Load session as batch" auto-partitions into chunks of 50 when
      >50 jobs (configurable via _AF3_PARTITION_SIZE), preventing OOM
    • Each partition is staged as a separate JSON and submitted
      sequentially (one job at a time in the queue) — not in parallel
    • Command preview now lists every partition command with job IDs
    • "Download JSONs" button in Submit tab: saves all staged partition
      JSONs to a local folder chosen by the user
    • _dv_on_upload_done handles multi-partition results list
    • Table rows updated per-partition with individual SLURM IDs

  v1.10 (2026) — ORF table copy & export + code cleanup
    • Multi-row selection (Ctrl/Shift+click) in ORF table
    • Ctrl+C copies selected rows as TSV (Excel/Calc ready)
    • Export toolbar: copy selected, copy all, Export menu
    • Export TSV: columns only / full (+DNA+Protein) / annotated only
    • Export FASTA protein/DNA with rich annotation headers
    • Right-click menu expanded with copy + export submenu
    • Removed orphan _show_pyrodigal_dialog (228 lines)
    • Manual and Tutorial fully updated for v1.10
    • QShortcut import fixed for PyQt5

  v1.09 (2026) — Pyrodigal parameters in ORF Analysis Parameters
    • Full parameter dialog: mode, translation table (11/4/25/15),
      min gene size, closed ends, mask N runs
    • Post-prediction start codon filter (ATG/GTG/TTG)
    • Parameters persisted in project save/load

  v1.08 (2026) — [merged into v1.09]

  v1.07 (2026) — EVcouplings tab removed (Python 3.12 incompatible)

  v1.06 (2026) — Definitive project save/load (schema v2)
    • Full ORF annotation state saved and restored
    • hmm_hits_all saved directly; domains re-injected on load
    • AF3 selection table, EVC state, Pyrodigal params persisted
    • AF3 jobs validated; AF3 analysis dir handled gracefully

  v1.05 (2026) — Save/Load critical bug fixes (5 bugs)

  v1.04 (2026) — EVcouplings coevolution tab (test, later removed)

  v1.03 (2026) — DaVinci cluster integration & SLURM defaults
    • SLURM Array anti-OOM export (batches + run_array.sh)
    • max50 partition defaults: 7d walltime, 64G RAM, 64 CPUs
    • Partition tooltips with DaVinci limits

  v1.02 (2026) — Genome-wide interactome scan
    • 'Add All ORFs' button with size/HMM filter dialog
    • Mode 'Interactoma Genômico': selected ORF vs all genome ORFs
    • SLURM array anti-OOM export

  v1.01 (2026) — Initial public release
    • Generic HPC server tab (SSH/SFTP, SLURM/PBS/LSF support)
    • AF3 results analysis tab (PAE heatmap, pLDDT, ipTM scoring)
    • Inter-chain contact detection with configurable PAE threshold
    • Scheduler selector (SLURM / PBS / LSF / None)
    • Environment activation selector (module load / conda / singularity /
      source script / None)
    • Automatic job-name sanitisation (shell-safe characters)
    • Multi-language UI (EN, PT-BR, ES, FR, ZH, JA)
"""

from __future__ import annotations

import sys
import os
import re
import csv
import json
import math
import stat
import random
import base64
import subprocess
import shutil
import tempfile
from pathlib import Path
from collections import defaultdict
from datetime import datetime

# ═══════════════════════════════════════════════════════════════
# PyQt6 / PyQt5 import with fallback
# ═══════════════════════════════════════════════════════════════

try:
    from PyQt6.QtWidgets import (
        QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
        QGridLayout, QSplitter, QTabWidget, QGroupBox, QPushButton,
        QLabel, QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox,
        QCheckBox, QTextEdit, QPlainTextEdit, QTextBrowser, QTableWidget,
        QTableWidgetItem,
        QFileDialog, QMessageBox, QDialog, QDialogButtonBox,
        QMenu, QAction, QToolBar, QScrollArea,
        QFrame, QSizePolicy, QAbstractItemView,
        QInputDialog,
    )
    from PyQt6.QtCore import (
        Qt, QTimer, QThread, pyqtSignal, QObject,
        QPointF, QRectF, QSettings,
    )
    QT_VERSION = 6
    # Qt6 enums
    AlignLeft = Qt.AlignmentFlag.AlignLeft
    AlignRight = Qt.AlignmentFlag.AlignRight
    AlignCenter = Qt.AlignmentFlag.AlignCenter
    AlignTop = Qt.AlignmentFlag.AlignTop
    Horizontal = Qt.Orientation.Horizontal
    Vertical = Qt.Orientation.Vertical
    SelectRows = QAbstractItemView.SelectionBehavior.SelectRows
    SingleSelection = QAbstractItemView.SelectionMode.SingleSelection
    ReadOnly = QLineEdit.EchoMode.Normal  # placeholder
    LeftButton = Qt.MouseButton.LeftButton

except ImportError:
    from PyQt5.QtWidgets import (
        QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
        QGridLayout, QSplitter, QTabWidget, QGroupBox, QPushButton,
        QLabel, QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox,
        QCheckBox, QTextEdit, QPlainTextEdit, QTextBrowser, QTableWidget,
        QTableWidgetItem,
        QFileDialog, QMessageBox, QDialog, QDialogButtonBox,
        QMenu, QAction, QToolBar, QScrollArea,
        QFrame, QSizePolicy, QAbstractItemView,
        QInputDialog,
    )
    from PyQt5.QtCore import (
        Qt, QTimer, QThread, pyqtSignal, QObject,
        QPointF, QSettings,
    )
    from PyQt5.QtGui import (
        QPainter, QPen, QBrush, QColor, QFont, QCursor,
        QPolygonF, QKeySequence,
        QPainterPath, QPainterPathStroker,
    )
    from PyQt5.QtWidgets import QShortcut
    try:
        from PyQt5.QtCore import QRectF
    except ImportError:
        QRectF = None
    QT_VERSION = 5
    AlignLeft = Qt.AlignLeft
    AlignRight = Qt.AlignRight
    AlignCenter = Qt.AlignCenter
    AlignTop = Qt.AlignTop
    Horizontal = Qt.Horizontal
    Vertical = Qt.Vertical
    SelectRows = QAbstractItemView.SelectRows
    SingleSelection = QAbstractItemView.SingleSelection
    LeftButton = Qt.LeftButton

BIOPYTHON_AVAILABLE = False

try:
    import pyrodigal
    PYRODIGAL_AVAILABLE = True
except ImportError:
    PYRODIGAL_AVAILABLE = False

try:
    import paramiko
    PARAMIKO_AVAILABLE = True
except ImportError:
    PARAMIKO_AVAILABLE = False

# Lazy matplotlib loader.
# Importing matplotlib can be slow on HPC/shared filesystems, so it is loaded
# only when plots are actually created.
import importlib.util as _importlib_util

MATPLOTLIB_AVAILABLE = _importlib_util.find_spec("matplotlib") is not None


def _load_matplotlib_objects():
    import matplotlib

    try:
        if QT_VERSION == 6:
            matplotlib.use("Qt6Agg")
            from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as _FigureCanvas
        else:
            matplotlib.use("Qt5Agg")
            from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as _FigureCanvas
    except Exception:
        matplotlib.use("Agg")
        from matplotlib.backends.backend_agg import FigureCanvasAgg as _FigureCanvas

    import matplotlib.pyplot as _plt
    return _plt, _FigureCanvas


class _LazyPyplot:
    def __init__(self):
        self._plt = None

    def _load(self):
        if self._plt is None:
            self._plt, _ = _load_matplotlib_objects()
        return self._plt

    def __getattr__(self, name):
        return getattr(self._load(), name)


class _LazyFigureCanvas:
    def __init__(self):
        self._canvas_cls = None

    def _load(self):
        if self._canvas_cls is None:
            _, self._canvas_cls = _load_matplotlib_objects()
        return self._canvas_cls

    def __call__(self, *args, **kwargs):
        return self._load()(*args, **kwargs)


plt = _LazyPyplot()
FigureCanvas = _LazyFigureCanvas() if MATPLOTLIB_AVAILABLE else None

# Lazy numpy loader.
NUMPY_AVAILABLE = _importlib_util.find_spec("numpy") is not None


class _LazyNumpy:
    def __init__(self):
        self._np = None

    def _load(self):
        if self._np is None:
            import numpy as _np
            self._np = _np
        return self._np

    def __getattr__(self, name):
        return getattr(self._load(), name)


np = _LazyNumpy() if NUMPY_AVAILABLE else None  # type: ignore[assignment]

# scipy.ndimage is used by the motif detection in the AlphaFold
# Analysis tab (connected-component labeling + morphological closing).
# It is OPTIONAL — a pure-numpy BFS fallback is used when it's missing.
# Lazy scipy.ndimage loader.
SCIPY_NDIMAGE_AVAILABLE = _importlib_util.find_spec("scipy.ndimage") is not None


class _LazyScipyNdimage:
    def __init__(self):
        self._mod = None

    def _load(self):
        if self._mod is None:
            import scipy.ndimage as _mod
            self._mod = _mod
        return self._mod

    def __getattr__(self, name):
        return getattr(self._load(), name)


_scipy_ndimage = _LazyScipyNdimage() if SCIPY_NDIMAGE_AVAILABLE else None  # type: ignore[assignment]



# ═══════════════════════════════════════════════════════════════
# EXTERNAL BACKEND DETECTION
# Extracted into ppigfinder.infrastructure.backends
# ═══════════════════════════════════════════════════════════════

try:
    from .infrastructure.backends import BACKENDS, detect_backends
except ImportError:
    from ppigfinder.infrastructure.backends import BACKENDS, detect_backends

try:
    from .bioseq.sequence import (
        gc_content as _bioseq_gc_content,
        reverse_complement as _bioseq_reverse_complement,
        translate_dna as _bioseq_translate_dna,
    )
except ImportError:
    from ppigfinder.bioseq.sequence import (
        gc_content as _bioseq_gc_content,
        reverse_complement as _bioseq_reverse_complement,
        translate_dna as _bioseq_translate_dna,
    )


# ═══════════════════════════════════════════════════════════════
# FILE FORMAT SUPPORT
# Extracted from the legacy monolith into ppigfinder.io
# ═══════════════════════════════════════════════════════════════

try:
    from .io.snapgene import parse_snapgene_dna, write_snapgene_dna
    from .io.genbank import parse_genbank, write_genbank
except ImportError:
    from ppigfinder.io.snapgene import parse_snapgene_dna, write_snapgene_dna
    from ppigfinder.io.genbank import parse_genbank, write_genbank




try:
    from .bioseq.orf_finder import (
        find_orfs as _bioseq_find_orfs,
        find_orfs_pyrodigal as _bioseq_find_orfs_pyrodigal,
        find_orfs_hybrid as _bioseq_find_orfs_hybrid,
    )
except ImportError:
    from ppigfinder.bioseq.orf_finder import (
        find_orfs as _bioseq_find_orfs,
        find_orfs_pyrodigal as _bioseq_find_orfs_pyrodigal,
        find_orfs_hybrid as _bioseq_find_orfs_hybrid,
    )

class AdvancedORFAnalyzer:

    # Class-level error state — always present, no AttributeError on first read
    _last_blast_error: str = ''
    _last_hmm_error:   str = ''

    CODON_TABLE = {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
        'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
        'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
        'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
    }

    # Domains come from HMM search and user reference DB only (no regex false positives)
    DOMAIN_SIGNATURES = {}  # Empty — user builds via HMM profiles

    BLOSUM62 = {
        ('A','A'):4,('A','R'):-1,('A','N'):-2,('A','D'):-2,('A','C'):0,('A','Q'):-1,('A','E'):-1,
        ('A','G'):0,('A','H'):-2,('A','I'):-1,('A','L'):-1,('A','K'):-1,('A','M'):-1,('A','F'):-2,
        ('A','P'):-1,('A','S'):1,('A','T'):0,('A','W'):-3,('A','Y'):-2,('A','V'):0,
        ('R','R'):5,('R','N'):0,('R','D'):-2,('R','C'):-3,('R','Q'):1,('R','E'):0,
        ('R','G'):-2,('R','H'):0,('R','I'):-3,('R','L'):-2,('R','K'):2,('R','M'):-1,('R','F'):-3,
        ('R','P'):-2,('R','S'):-1,('R','T'):-1,('R','W'):-3,('R','Y'):-2,('R','V'):-3,
        ('N','N'):6,('N','D'):1,('N','C'):-3,('N','Q'):0,('N','E'):0,
        ('N','G'):0,('N','H'):1,('N','I'):-3,('N','L'):-3,('N','K'):0,('N','M'):-2,('N','F'):-3,
        ('N','P'):-2,('N','S'):1,('N','T'):0,('N','W'):-4,('N','Y'):-2,('N','V'):-3,
        ('D','D'):6,('D','C'):-3,('D','Q'):0,('D','E'):2,
        ('D','G'):-1,('D','H'):-1,('D','I'):-3,('D','L'):-4,('D','K'):-1,('D','M'):-3,('D','F'):-3,
        ('D','P'):-1,('D','S'):0,('D','T'):-1,('D','W'):-4,('D','Y'):-3,('D','V'):-3,
        ('C','C'):9,('C','Q'):-3,('C','E'):-4,
        ('C','G'):-3,('C','H'):-3,('C','I'):-1,('C','L'):-1,('C','K'):-3,('C','M'):-1,('C','F'):-2,
        ('C','P'):-3,('C','S'):-1,('C','T'):-1,('C','W'):-2,('C','Y'):-2,('C','V'):-1,
        ('Q','Q'):5,('Q','E'):2,('Q','G'):-2,('Q','H'):0,('Q','I'):-3,('Q','L'):-2,('Q','K'):1,
        ('Q','M'):0,('Q','F'):-3,('Q','P'):-1,('Q','S'):0,('Q','T'):-1,('Q','W'):-2,('Q','Y'):-1,('Q','V'):-2,
        ('E','E'):5,('E','G'):-2,('E','H'):0,('E','I'):-3,('E','L'):-3,('E','K'):1,('E','M'):-2,
        ('E','F'):-3,('E','P'):-1,('E','S'):0,('E','T'):-1,('E','W'):-3,('E','Y'):-2,('E','V'):-2,
        ('G','G'):6,('G','H'):-2,('G','I'):-4,('G','L'):-4,('G','K'):-2,('G','M'):-3,('G','F'):-3,
        ('G','P'):-2,('G','S'):0,('G','T'):-2,('G','W'):-2,('G','Y'):-3,('G','V'):-3,
        ('H','H'):8,('H','I'):-3,('H','L'):-3,('H','K'):-1,('H','M'):-2,('H','F'):-1,
        ('H','P'):-2,('H','S'):-1,('H','T'):-2,('H','W'):-2,('H','Y'):2,('H','V'):-3,
        ('I','I'):4,('I','L'):2,('I','K'):-3,('I','M'):1,('I','F'):0,
        ('I','P'):-3,('I','S'):-2,('I','T'):-1,('I','W'):-3,('I','Y'):-1,('I','V'):3,
        ('L','L'):4,('L','K'):-2,('L','M'):2,('L','F'):0,
        ('L','P'):-3,('L','S'):-2,('L','T'):-1,('L','W'):-2,('L','Y'):-1,('L','V'):1,
        ('K','K'):5,('K','M'):-1,('K','F'):-3,('K','P'):-1,('K','S'):0,('K','T'):-1,
        ('K','W'):-3,('K','Y'):-2,('K','V'):-2,
        ('M','M'):5,('M','F'):0,('M','P'):-2,('M','S'):-1,('M','T'):-1,('M','W'):-1,('M','Y'):-1,('M','V'):1,
        ('F','F'):6,('F','P'):-4,('F','S'):-2,('F','T'):-2,('F','W'):1,('F','Y'):3,('F','V'):-1,
        ('P','P'):7,('P','S'):-1,('P','T'):-1,('P','W'):-4,('P','Y'):-3,('P','V'):-2,
        ('S','S'):4,('S','T'):1,('S','W'):-3,('S','Y'):-2,('S','V'):-2,
        ('T','T'):5,('T','W'):-2,('T','Y'):-2,('T','V'):0,
        ('W','W'):11,('W','Y'):2,('W','V'):-3,
        ('Y','Y'):7,('Y','V'):-1,('V','V'):4,
    }

    # Karlin-Altschul parameters for BLOSUM62
    KA_LAMBDA = 0.267
    KA_K = 0.041

    def translate(self, dna_seq):
        return _bioseq_translate_dna(dna_seq)

    def reverse_complement(self, seq):
        return _bioseq_reverse_complement(seq)

    def gc_content(self, seq):
        return _bioseq_gc_content(seq)

    def find_orfs(self, dna_sequence, min_aa=30, start_codons=None):
        return _bioseq_find_orfs(
            dna_sequence,
            min_aa=min_aa,
            start_codons=start_codons,
        )

    def find_orfs_pyrodigal(self, dna_sequence, meta=True, min_aa=30,
                            closed_ends=False, translation_table=11, mask=False):
        return _bioseq_find_orfs_pyrodigal(
            dna_sequence,
            meta=meta,
            min_aa=min_aa,
            closed_ends=closed_ends,
            translation_table=translation_table,
            mask=mask,
        )

    def find_orfs_hybrid(self, dna_sequence,
                         min_aa=30, start_codons=None,
                         pyro_meta=True, pyro_min_aa=30,
                         pyro_closed=False, pyro_translation_table=11,
                         pyro_mask=False, pyro_start_filter=None):
        return _bioseq_find_orfs_hybrid(
            dna_sequence,
            min_aa=min_aa,
            start_codons=start_codons,
            pyro_meta=pyro_meta,
            pyro_min_aa=pyro_min_aa,
            pyro_closed=pyro_closed,
            pyro_translation_table=pyro_translation_table,
            pyro_mask=pyro_mask,
            pyro_start_filter=pyro_start_filter,
        )

    def classify_domains(self, protein_seq):
        found = []
        for name, info in self.DOMAIN_SIGNATURES.items():
            for pattern in info['patterns']:
                for m in re.finditer(pattern, protein_seq, re.IGNORECASE):
                    found.append({'domain': name, 'description': info['desc'],
                                  'system': info['system'], 'role': info['role'],
                                  'start': m.start(), 'end': m.end(), 'match': m.group()})
        return found

    def analyze_neighborhood(self, orfs, target_idx, window_kb=15.0):
        target = orfs[target_idx]
        window_bp = int(window_kb * 1000)
        center = (target['start'] + target['end']) // 2
        ws = max(0, center - window_bp // 2)
        we = center + window_bp // 2
        neighbors = []
        doms_w = defaultdict(int)
        sys_w = defaultdict(int)
        hypo = 0
        for i, orf in enumerate(orfs):
            if i == target_idx:
                continue
            oc = (orf['start'] + orf['end']) // 2
            if ws <= oc <= we:
                d = abs(oc - center)
                neighbors.append({
                    'orf_index': i,
                    'distance_bp': d,
                    'strand': orf['strand'],
                    'length_aa': len(orf['protein'].rstrip('*')),
                    'domains': orf.get('domains', []),
                })
                for dom in orf.get('domains', []):
                    doms_w[dom['domain']] += 1
                    sys_w[dom['system']] += 1
                if not orf.get('domains'):
                    hypo += 1
        cs = 0.0
        if len(neighbors) >= 3:
            cs += 0.2
        if len(sys_w) >= 1:
            cs += 0.3
        if len(doms_w) >= 2:
            cs += 0.2
        return {
            'target_idx': target_idx,
            'window_start': ws,
            'window_end': we,
            'window_kb': window_kb,
            'neighbors': sorted(neighbors, key=lambda x: x['distance_bp']),
            'total_orfs_in_window': len(neighbors),
            'domains_found': dict(doms_w),
            'systems_found': dict(sys_w),
            'hypothetical_proteins': hypo,
            'cluster_score': min(cs, 1.0),
        }

    def calc_evalue(self, score, query_len, db_size):
        """E-value by Karlin-Altschul."""
        if score <= 0: return 999
        try:
            e = self.KA_K * query_len * db_size * math.exp(-self.KA_LAMBDA * score)
            return e
        except OverflowError:
            return 0.0

    def _blosum_score(self, a, b):
        if a == '*' or b == '*': return -4
        if a == 'X' or b == 'X': return -1
        return self.BLOSUM62.get((a, b), self.BLOSUM62.get((b, a), -4))

    # ═══════ MÉTODO 1: K-MER DIAGONAL FILTER (RÁPIDO) ═══════

    def kmer_blast(self, query, subjects, params=None):
        """BLAST rápido: k-mer index + diagonal filter + banded SW."""
        if params is None: params = {}
        threshold = params.get('threshold', 30)
        gap_open = params.get('gap_open', -11)
        gap_extend = params.get('gap_extend', -1)
        word_size = params.get('word_size', 4)
        min_diag_hits = params.get('min_diag_hits', 2)
        evalue_max = params.get('evalue', 0.05)

        query = query.upper().replace('*', '').strip()
        if not query: return []
        q_len = len(query)
        db_size = sum(len(s.replace('*','')) for s in subjects)

        # Fase 1: indexar k-mers do query
        q_index = defaultdict(list)
        for i in range(q_len - word_size + 1):
            kmer = query[i:i+word_size]
            q_index[kmer].append(i)

        hits = []
        band_width = 15

        for idx, subject in enumerate(subjects):
            subj = subject.upper().replace('*', '').strip()
            if not subj or len(subj) < word_size: continue
            s_len = len(subj)

            # Fase 2: contar hits por diagonal
            diag_hits = defaultdict(list)
            for j in range(s_len - word_size + 1):
                kmer = subj[j:j+word_size]
                if kmer in q_index:
                    for i in q_index[kmer]:
                        diag = i - j
                        diag_hits[diag].append((i, j))

            # Keep only diagonals with enough seed hits
            hot_diags = [(d, pts) for d, pts in diag_hits.items()
                         if len(pts) >= min_diag_hits]
            if not hot_diags: continue

            # Phase 3: banded SW on the best diagonals
            best_diag = max(hot_diags, key=lambda x: len(x[1]))
            diag_val, seed_pts = best_diag

            # Determine the band around the best diagonal
            min_i = min(p[0] for p in seed_pts)
            max_i = max(p[0] for p in seed_pts)
            min_j = min(p[1] for p in seed_pts)
            max_j = max(p[1] for p in seed_pts)

            # Expand the alignment window
            q_start = max(0, min_i - band_width)
            q_end = min(q_len, max_i + word_size + band_width)
            s_start = max(0, min_j - band_width)
            s_end = min(s_len, max_j + word_size + band_width)

            q_sub = query[q_start:q_end]
            s_sub = subj[s_start:s_end]
            ql = len(q_sub); sl = len(s_sub)

            if ql < 3 or sl < 3: continue

            # Banded Smith-Waterman with affine gaps
            H = [[0]*(sl+1) for _ in range(ql+1)]
            tb = [[0]*(sl+1) for _ in range(ql+1)]
            max_score = 0; mi = mj = 0

            for i in range(1, ql+1):
                j_center = i - diag_val + s_start - q_start
                j_lo = max(1, j_center - band_width)
                j_hi = min(sl, j_center + band_width)
                for j in range(j_lo, j_hi + 1):
                    diag_s = H[i-1][j-1] + self._blosum_score(q_sub[i-1], s_sub[j-1])
                    # Affine gap: open+extend for new gap, extend for continuation
                    up_open = H[i-1][j] + gap_open + gap_extend
                    left_open = H[i][j-1] + gap_open + gap_extend
                    best = max(0, diag_s, up_open, left_open)
                    H[i][j] = best
                    if best <= 0: tb[i][j] = 0
                    elif best == diag_s: tb[i][j] = 1
                    elif best == up_open: tb[i][j] = 2
                    else: tb[i][j] = 3
                    if best > max_score:
                        max_score = best; mi = i; mj = j

            if max_score <= 0: continue

            # Traceback
            aln_q, aln_s, aln_m = [], [], []
            i, j = mi, mj
            ids = pos = gaps = al = 0
            while i > 0 and j > 0 and H[i][j] > 0:
                t = tb[i][j]
                if t == 1:
                    qa, sa = q_sub[i-1], s_sub[j-1]
                    aln_q.append(qa); aln_s.append(sa)
                    sc = self._blosum_score(qa, sa)
                    if qa == sa: aln_m.append(qa); ids += 1; pos += 1
                    elif sc > 0: aln_m.append('+'); pos += 1
                    else: aln_m.append(' ')
                    al += 1; i -= 1; j -= 1
                elif t == 2:
                    aln_q.append(q_sub[i-1]); aln_s.append('-'); aln_m.append(' ')
                    al += 1; gaps += 1; i -= 1
                elif t == 3:
                    aln_q.append('-'); aln_s.append(s_sub[j-1]); aln_m.append(' ')
                    al += 1; gaps += 1; j -= 1
                else: break

            aln_q.reverse(); aln_s.reverse(); aln_m.reverse()
            id_pct = (ids/al*100) if al > 0 else 0
            pos_pct = (pos/al*100) if al > 0 else 0
            cov = (al/q_len*100) if q_len > 0 else 0
            evalue = self.calc_evalue(max_score, q_len, db_size)

            if id_pct >= threshold and evalue <= evalue_max:
                hits.append({
                    'orf_index': idx, 'identity': round(id_pct, 1),
                    'positives': round(pos_pct, 1), 'score': max_score,
                    'aln_length': al, 'identities_count': ids,
                    'positives_count': pos, 'gaps': gaps,
                    'q_start': q_start + i + 1, 'q_end': q_start + mi,
                    's_start': s_start + j + 1, 's_end': s_start + mj,
                    'coverage': round(cov, 1), 'evalue': evalue,
                    'aln_query': ''.join(aln_q), 'aln_midline': ''.join(aln_m),
                    'aln_subject': ''.join(aln_s),
                })
        return sorted(hits, key=lambda x: x['score'], reverse=True)

    # ═══════ METHOD 2: NCBI BLAST+ (EXTERNAL) ═══════

    def run_ncbi_blast(self, query_protein, orfs, params=None):
        """Run NCBI BLAST+ blastp via subprocess."""
        if not BACKENDS.get('blast+', {}).get('available'):
            return None  # fallback para Python
        if params is None: params = {}
        evalue_thresh = params.get('evalue', 0.05)
        matrix = params.get('matrix', 'BLOSUM62')
        word_size = params.get('word_size', 5)
        max_targets = params.get('max_targets', 100)
        gap_open = abs(params.get('gap_open', -11))
        gap_extend = abs(params.get('gap_extend', -1))
        low_complexity = params.get('low_complexity', True)

        tmpdir = tempfile.mkdtemp(prefix='blast_')
        try:
            # Write query FASTA
            qfile = os.path.join(tmpdir, 'query.fasta')
            with open(qfile, 'w', encoding='utf-8') as f:
                f.write(">query\n")
                for i in range(0, len(query_protein), 80):
                    f.write(query_protein[i:i+80] + "\n")

            # Write subject FASTA (ORF proteins)
            sfile = os.path.join(tmpdir, 'subjects.fasta')
            with open(sfile, 'w', encoding='utf-8') as f:
                for i, orf in enumerate(orfs):
                    prot = orf['protein'].rstrip('*')
                    if prot:
                        f.write(f">ORF{i+1}\n")
                        for j in range(0, len(prot), 80):
                            f.write(prot[j:j+80] + "\n")

            # Run blastp (protein vs ORF proteins)
            outfile = os.path.join(tmpdir, 'results.txt')
            cmd = [
                'blastp',
                '-query', qfile,
                '-subject', sfile,
                '-outfmt', '6 sseqid score pident positive length qstart qend sstart send evalue gaps',
                '-evalue', str(evalue_thresh),
                '-matrix', matrix,
                '-word_size', str(word_size),
                '-gapopen', str(gap_open),
                '-gapextend', str(gap_extend),
                '-max_target_seqs', str(max_targets),
                '-seg', 'yes' if low_complexity else 'no',
                '-out', outfile,
            ]
            r1 = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            if r1.returncode != 0:
                self._last_blast_error = (r1.stderr or r1.stdout or 'blastp failed').strip()

            # Parse tabular output
            hits = []
            if os.path.exists(outfile):
                with open(outfile, encoding='utf-8') as f:
                    for line in f:
                        parts = line.strip().split('\t')
                        if len(parts) >= 11:
                            orf_name = parts[0]
                            orf_idx = int(orf_name.replace('ORF','')) - 1
                            hits.append({
                                'orf_index': orf_idx,
                                'score': int(float(parts[1])),
                                'identity': round(float(parts[2]), 1),
                                'positives': round(float(parts[3]), 1),
                                'aln_length': int(parts[4]),
                                'q_start': int(parts[5]),
                                'q_end': int(parts[6]),
                                's_start': int(parts[7]),
                                's_end': int(parts[8]),
                                'evalue': float(parts[9]),
                                'gaps': int(parts[10]),
                                'identities_count': int(round(float(parts[2])*int(parts[4])/100)),
                                'positives_count': int(round(float(parts[3])*int(parts[4])/100)),
                                'coverage': round(int(parts[4])/len(query_protein)*100, 1),
                                'aln_query': '', 'aln_midline': '', 'aln_subject': '',
                            })

            # Retrieve detailed alignments for top hits
            if hits:
                aln_file = os.path.join(tmpdir, 'results_aln.txt')
                cmd2 = [
                    'blastp', '-query', qfile, '-subject', sfile,
                    '-outfmt', '0', '-evalue', str(evalue_thresh),
                    '-matrix', matrix, '-word_size', str(word_size),
                    '-gapopen', str(gap_open), '-gapextend', str(gap_extend),
                    '-max_target_seqs', str(min(10, max_targets)),
                    '-out', aln_file,
                ]
                r2 = subprocess.run(cmd2, capture_output=True, text=True, timeout=120)
                if r2.returncode != 0 and not self._last_blast_error:
                    self._last_blast_error = (r2.stderr or r2.stdout or 'blastp (aln) failed').strip()
                if os.path.exists(aln_file):
                    self._parse_blast_alignments(aln_file, hits)

            return sorted(hits, key=lambda x: x['score'], reverse=True)

        except subprocess.TimeoutExpired:
            self._last_blast_error = 'BLAST timed out after 120 s.'
            return None
        except FileNotFoundError:
            self._last_blast_error = 'blastp executable not found. Install NCBI BLAST+ or check PATH.'
            return None
        except Exception as e:
            self._last_blast_error = f'{type(e).__name__}: {e}'
            return None
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def _parse_blast_alignments(self, aln_file, hits):
        """Parse alignments from BLAST output format 0."""
        try:
            with open(aln_file) as f:
                content = f.read()
            # blastp -subject writes "> ORF42" (with space); accept both forms
            blocks = re.split(r'>\s*ORF(\d+)', content)
            hit_dict = {h['orf_index']: h for h in hits}
            for i in range(1, len(blocks) - 1, 2):
                try:
                    orf_idx = int(blocks[i]) - 1
                except ValueError:
                    continue
                block = blocks[i + 1]
                q_seqs, s_seqs = [], []
                for qm in re.finditer(r'Query\s+\d+\s+([A-Z\-]+)\s+\d+', block):
                    q_seqs.append(qm.group(1))
                for sm in re.finditer(r'Sbjct\s+\d+\s+([A-Z\-]+)\s+\d+', block):
                    s_seqs.append(sm.group(1))
                if orf_idx in hit_dict and q_seqs and s_seqs:
                    full_q = ''.join(q_seqs)
                    full_s = ''.join(s_seqs)
                    # Extrair midline real do BLAST (linha entre Query e Sbjct)
                    mid_lines = re.findall(
                        r'Query\s+\d+\s+[A-Z\-]+\s+\d+\n([ A-Z\+]*)\nSbjct', block)
                    if mid_lines:
                        raw_mid = ''.join(mid_lines)
                        # Alinhar comprimento com a sequência
                        mid = list(raw_mid[:len(full_q)].ljust(len(full_q)))
                    else:
                        # Fallback: reconstruir via BLOSUM
                        mid = []
                        for qa, sa in zip(full_q, full_s):
                            if qa == sa:
                                mid.append(qa)
                            elif qa != '-' and sa != '-' and self._blosum_score(qa, sa) > 0:
                                mid.append('+')
                            else:
                                mid.append(' ')
                    hit_dict[orf_idx]['aln_query']   = full_q
                    hit_dict[orf_idx]['aln_subject']  = full_s
                    hit_dict[orf_idx]['aln_midline']  = ''.join(mid)
        except Exception:
            pass

    # ═══════ METHOD 3: FULL SMITH-WATERMAN ═══════

    def sw_blast(self, query, subjects, params=None):
        """Full Smith-Waterman alignment (most sensitive, slower)."""
        if params is None: params = {}
        threshold = params.get('threshold', 30)
        gap_open = params.get('gap_open', -11)
        gap_extend = params.get('gap_extend', -1)
        evalue_max = params.get('evalue', 0.05)
        query = query.upper().replace('*', '').strip()
        if not query: return []
        q_len = len(query)
        db_size = sum(len(s.replace('*','')) for s in subjects)
        hits = []
        query_kmers = set()
        for ki in range(q_len - 2):
            query_kmers.add(query[ki:ki+3])

        for idx, subject in enumerate(subjects):
            sc = subject.upper().replace('*', '').strip()
            if not sc: continue
            sl = len(sc)
            # Pre-filtro 3-mer
            found = False
            for si in range(sl - 2):
                if sc[si:si+3] in query_kmers: found = True; break
            if not found: continue

            NEG_INF = float('-inf')
            E_prev = [NEG_INF]*(sl+1)
            H = [[0]*(sl+1) for _ in range(q_len+1)]
            tb = [[0]*(sl+1) for _ in range(q_len+1)]
            ms = mi = mj = 0

            for i in range(1, q_len+1):
                F_val = NEG_INF
                E_curr = [NEG_INF]*(sl+1)
                for j in range(1, sl+1):
                    diag = H[i-1][j-1] + self._blosum_score(query[i-1], sc[j-1])
                    E_curr[j] = max(H[i-1][j]+gap_open+gap_extend, E_prev[j]+gap_extend)
                    f_o = H[i][j-1]+gap_open+gap_extend
                    F_val = max(f_o, F_val+gap_extend)
                    best = max(0, diag, E_curr[j], F_val)
                    H[i][j] = best
                    if best <= 0: tb[i][j] = 0
                    elif best == diag: tb[i][j] = 1
                    elif best == E_curr[j]: tb[i][j] = 2
                    else: tb[i][j] = 3
                    if best > ms: ms = best; mi = i; mj = j
                E_prev = E_curr

            if ms <= 0: continue
            aln_q, aln_s, aln_m = [], [], []
            i, j = mi, mj
            ids = pos = gaps = al = 0
            while i > 0 and j > 0 and H[i][j] > 0:
                t = tb[i][j]
                if t == 1:
                    qa, sa = query[i-1], sc[j-1]
                    aln_q.append(qa); aln_s.append(sa)
                    s = self._blosum_score(qa, sa)
                    if qa == sa: aln_m.append(qa); ids += 1; pos += 1
                    elif s > 0: aln_m.append('+'); pos += 1
                    else: aln_m.append(' ')
                    al += 1; i -= 1; j -= 1
                elif t == 2:
                    aln_q.append(query[i-1]); aln_s.append('-'); aln_m.append(' ')
                    al += 1; gaps += 1; i -= 1
                elif t == 3:
                    aln_q.append('-'); aln_s.append(sc[j-1]); aln_m.append(' ')
                    al += 1; gaps += 1; j -= 1
                else: break
            aln_q.reverse(); aln_s.reverse(); aln_m.reverse()
            id_p = (ids/al*100) if al > 0 else 0
            pos_p = (pos/al*100) if al > 0 else 0
            cov = (al/q_len*100) if q_len > 0 else 0
            ev = self.calc_evalue(ms, q_len, db_size)
            if id_p >= threshold and ev <= evalue_max:
                hits.append({
                    'orf_index': idx, 'identity': round(id_p,1), 'positives': round(pos_p,1),
                    'score': ms, 'aln_length': al, 'identities_count': ids,
                    'positives_count': pos, 'gaps': gaps,
                    'q_start': i+1, 'q_end': mi, 's_start': j+1, 's_end': mj,
                    'coverage': round(cov,1), 'evalue': ev,
                    'aln_query': ''.join(aln_q), 'aln_midline': ''.join(aln_m),
                    'aln_subject': ''.join(aln_s),
                })
        return sorted(hits, key=lambda x: x['score'], reverse=True)

    # ═══════ HMM SEARCH ═══════

    def hmm_scan_orfs(self, hmm_file, orfs, params=None):
        """Search HMM profile against all ORFs. Uses HMMER3 if available."""
        if params is None: params = {}
        if BACKENDS.get('hmmer3', {}).get('available'):
            return self._hmmer3_search(hmm_file, orfs, params)
        return self._pssm_scan(hmm_file, orfs, params)

    def _hmmer3_search(self, hmm_file, orfs, params=None):
        """Chama hmmsearch do HMMER3. Usa --domtblout + -A para coordenadas e alinhamentos."""
        if params is None: params = {}
        hmm_evalue  = params.get('hmm_evalue', 10.0)
        hmm_score   = params.get('hmm_score_thresh', None)
        dom_evalue  = params.get('hmm_dom_evalue', 10.0)
        use_wsl     = BACKENDS.get('hmmer3', {}).get('wsl', False)
        tmpdir      = tempfile.mkdtemp(prefix='hmm_')
        try:
            # Escrever proteínas
            sfile = os.path.join(tmpdir, 'orfs.fasta')
            n_written = 0
            with open(sfile, 'w') as f:
                for i, orf in enumerate(orfs):
                    prot = orf['protein'].rstrip('*')
                    if prot and len(prot) >= 10:
                        f.write(f">ORF{i+1}\n")
                        for j in range(0, len(prot), 80):
                            f.write(prot[j:j+80] + "\n")
                        n_written += 1

            domtbl_file = os.path.join(tmpdir, 'results.domtbl')
            aln_file    = os.path.join(tmpdir, 'results.sto')

            if use_wsl:
                hmm_copy = os.path.join(tmpdir, 'profile.hmm')
                shutil.copy2(hmm_file, hmm_copy)

                def to_wsl_path(p):
                    p = os.path.abspath(p).replace('\\', '/')
                    if len(p) > 1 and p[1] == ':':
                        return f"/mnt/{p[0].lower()}{p[2:]}"
                    return p

                wsl_sfile   = to_wsl_path(sfile)
                wsl_hmm     = to_wsl_path(hmm_copy)
                wsl_domtbl  = to_wsl_path(domtbl_file)
                wsl_aln     = to_wsl_path(aln_file)

                hmm_cmd  = f'hmmsearch --domtblout "{wsl_domtbl}" -A "{wsl_aln}"'
                hmm_cmd += f' -E {hmm_evalue} --domE {dom_evalue}'
                if hmm_score is not None:
                    hmm_cmd += f' -T {hmm_score}'
                hmm_cmd += f' "{wsl_hmm}" "{wsl_sfile}"'
                cmd = ['wsl', 'bash', '-c', hmm_cmd]
            else:
                cmd = ['hmmsearch',
                       '--domtblout', domtbl_file,
                       '-A', aln_file,
                       '-E', str(hmm_evalue),
                       '--domE', str(dom_evalue)]
                if hmm_score is not None:
                    cmd.extend(['-T', str(hmm_score)])
                cmd.extend([hmm_file, sfile])

            result = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
            if result.returncode != 0:
                self._last_hmm_error = (result.stderr or result.stdout or 'hmmsearch failed').strip()

            # ── Parse --domtblout (coordenadas precisas) ──────────────
            hits = []
            if os.path.exists(domtbl_file):
                with open(domtbl_file, encoding='utf-8') as f:
                    for line in f:
                        if line.startswith('#'): continue
                        parts = line.split()
                        if len(parts) < 23: continue
                        orf_name = parts[0]
                        try:
                            orf_idx = int(orf_name.replace('ORF', '')) - 1
                            # domtblout columns (0-based):
                            # 0=target, 2=tlen, 3=query, 5=qlen,
                            # 11=c-evalue, 12=i-evalue, 13=score, 14=bias
                            # 15=hmm_from, 16=hmm_to
                            # 17=ali_from, 18=ali_to
                            # 19=env_from, 20=env_to
                            hmm_from = int(parts[15])
                            hmm_to   = int(parts[16])
                            ali_from = int(parts[17])
                            ali_to   = int(parts[18])
                            score    = float(parts[13])
                            evalue   = float(parts[12])  # i-evalue (independent)
                            tlen     = int(parts[2])
                            qlen     = int(parts[5])
                            hits.append({
                                'orf_index':  orf_idx,
                                'orf_name':   orf_name,
                                'hmm_name':   parts[3],
                                'score':      score,
                                'evalue':     evalue,
                                'bias':       float(parts[14]),
                                'hmm_from':   hmm_from,
                                'hmm_to':     hmm_to,
                                'ali_from':   ali_from,
                                'ali_to':     ali_to,
                                'hmm_len':    qlen,
                                'target_len': tlen,
                                # Região formatada para exibição
                                'match_region': f"HMM:{hmm_from}-{hmm_to}/{qlen}  Prot:{ali_from}-{ali_to}/{tlen}",
                            })
                        except (ValueError, IndexError):
                            continue

            # ── Parse -A Stockholm (alinhamentos) ─────────────────────
            aln_dict = {}
            if os.path.exists(aln_file) and os.path.getsize(aln_file) > 10:
                try:
                    aln_dict = self._parse_stockholm_aln(aln_file)
                except Exception:
                    pass

            # Anexar alinhamento a cada hit
            for h in hits:
                # Stockholm usa nomes como "ORF104/28-251" — normalizar
                orf_name_bare  = h['orf_name']                       # "ORF104"
                orf_name_range = f"{orf_name_bare}/{h.get('ali_from','')}-{h.get('ali_to','')}"

                aln = (aln_dict.get(orf_name_bare)
                    or aln_dict.get(orf_name_range)
                    or next((v for k, v in aln_dict.items()
                             if k.startswith(orf_name_bare + '/')), None)
                    or {})

                if aln:
                    h['aln_hmm']    = aln.get('hmm',    '')
                    h['aln_target'] = aln.get('target', '')
                    h['aln_match']  = aln.get('match',  '')
                else:
                    # Fallback: mostrar a subsequência proteica alinhada
                    # even without a Stockholm alignment file
                    oi = h.get('orf_index', -1)
                    if 0 <= oi < len(orfs):
                        prot      = orfs[oi]['protein'].rstrip('*')
                        ali_from  = h.get('ali_from', 1)
                        ali_to    = h.get('ali_to',   len(prot))
                        try:
                            subseq = prot[int(ali_from) - 1 : int(ali_to)]
                        except (TypeError, ValueError):
                            subseq = prot[:50]
                        h['aln_hmm']    = ''        # HMM consensus não disponível sem -A
                        h['aln_target'] = subseq
                        h['aln_match']  = ''
                    else:
                        h['aln_hmm']    = ''
                        h['aln_target'] = ''
                        h['aln_match']  = ''

            if not hits:
                return [{'error': f"0 hits. (ORFs: {n_written}, return: {result.returncode})"}]
            return sorted(hits, key=lambda x: x['score'], reverse=True)

        except subprocess.TimeoutExpired:
            self._last_hmm_error = 'hmmsearch timed out after 180 s.'
            return [{'error': 'Timeout (180s).'}]
        except FileNotFoundError:
            self._last_hmm_error = 'hmmsearch executable not found. Install HMMER3 or check PATH.'
            return [{'error': self._last_hmm_error}]
        except Exception as e:
            self._last_hmm_error = f'{type(e).__name__}: {e}'
            return [{'error': str(e)}]
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def _parse_stockholm_aln(self, sto_file: str) -> dict:
        """
        Parse Stockholm output from hmmsearch -A.
        Returns dict: orf_name → {'hmm': str, 'target': str, 'match': str}
        """
        result = {}
        with open(sto_file, encoding='utf-8') as f:
            content = f.read()

        # Each alignment block starts with '# STOCKHOLM 1.0' and ends with '//'
        blocks = content.split('//')
        for block in blocks:
            lines = block.strip().splitlines()
            seqs  = {}     # name → sequence fragments
            rf    = ''
            for line in lines:
                if not line.strip() or line.startswith('#=GF'): continue
                if line.startswith('#=GC RF'):
                    rf += line.split()[-1]
                    continue
                if line.startswith('#=GC PP_cons') or line.startswith('#=GC seq_cons'):
                    continue
                if line.startswith('#=GC'):
                    continue
                if line.startswith('#=GR'):
                    # per-residue annotation — skip
                    continue
                if line.startswith('#'):
                    continue
                parts = line.split(None, 1)
                if len(parts) == 2:
                    name = parts[0]
                    seq  = parts[1].strip()
                    seqs.setdefault(name, '')
                    seqs[name] += seq

            # Identify the HMM query (not an ORF)
            hmm_seq    = ''
            target_seqs = {}
            for name, seq in seqs.items():
                if name.startswith('ORF'):
                    target_seqs[name] = seq
                else:
                    hmm_seq  = seq  # name ignored — only sequence is used

            # Build midline by comparing HMM consensus vs target
            for orf_name, tgt_seq in target_seqs.items():
                if not hmm_seq:
                    result[orf_name] = {'hmm': tgt_seq, 'target': tgt_seq, 'match': ''}
                    continue
                mid = []
                for hc, tc in zip(hmm_seq, tgt_seq):
                    if hc == '-' or tc == '-':
                        mid.append(' ')
                    elif hc.upper() == tc.upper():
                        mid.append('|')
                    elif hc.upper() in 'ACDEFGHIKLMNPQRSTVWY' and tc.upper() in 'ACDEFGHIKLMNPQRSTVWY':
                        sc = self._blosum_score(hc.upper(), tc.upper())
                        mid.append('+' if sc > 0 else '.')
                    else:
                        mid.append('.')
                result[orf_name] = {
                    'hmm':    hmm_seq,
                    'target': tgt_seq,
                    'match':  ''.join(mid),
                }
        return result

    def _pssm_scan(self, hmm_file, orfs, params=None):
        """PSSM scan derived from the HMM file (fallback when HMMER3 is unavailable)."""
        if params is None: params = {}
        score_thresh = params.get('hmm_score_thresh', None)
        evalue_thresh = params.get('hmm_evalue', 10.0)
        pssm, aa_order = self._parse_hmm_to_pssm(hmm_file)
        if not pssm: return [{'error': 'Failed to parse the HMM file'}]
        hits = []
        model_len = len(pssm)

        # Calculate null model score (average across all aa per position)
        null_score = 0
        for pos_scores in pssm:
            valid = [v for v in pos_scores.values() if v > -900]
            null_score += sum(valid) / len(valid) if valid else -4
        # Threshold: must beat null by at least 0.5 nats per position
        auto_thresh = null_score + model_len * 0.5
        min_score = score_thresh if score_thresh is not None else auto_thresh

        for i, orf in enumerate(orfs):
            prot = orf['protein'].rstrip('*').upper()
            if len(prot) < model_len: continue
            best_score = float('-inf')
            best_pos = 0
            for start in range(len(prot) - model_len + 1):
                score = 0
                for k in range(model_len):
                    aa = prot[start + k]
                    score += pssm[k].get(aa, -4)
                if score > best_score:
                    best_score = score; best_pos = start
            ev = self.calc_evalue(max(0, best_score - null_score), model_len, len(prot))
            if best_score > min_score and ev <= evalue_thresh:
                hits.append({
                    'orf_index': i, 'orf_name': f'ORF{i+1}',
                    'score': round(best_score, 1), 'position': best_pos,
                    'evalue': ev, 'hmm_name': Path(hmm_file).stem,
                    'match_region': f'{best_pos+1}-{best_pos+model_len}',
                })
        return sorted(hits, key=lambda x: x['score'], reverse=True)

    def _parse_hmm_to_pssm(self, hmm_file):
        """Extract PSSM from HMMER3 Match states. Reads AA order from the file header."""
        try:
            with open(hmm_file, encoding='utf-8') as f:
                content = f.read()
            if 'HMMER3' not in content: return None, None
            lines = content.split('\n')
            in_model = False; pssm = []
            aa_order = list('ACDEFGHIKLMNPQRSTVWY')  # default

            for line in lines:
                # Parse AA order from HMM header line
                if line.strip().startswith('HMM') and not line.strip().startswith('HMMER'):
                    parts = line.split()
                    if len(parts) >= 20:
                        aa_order = [p.strip() for p in parts[1:] if len(p.strip()) == 1 and p.strip().isalpha()]
                        if len(aa_order) < 20:
                            aa_order = list('ACDEFGHIKLMNPQRSTVWY')
                    in_model = True
                    continue
                if not in_model: continue
                if line.startswith('//'): break
                parts = line.split()
                # Match emission lines start with node number
                if len(parts) >= len(aa_order) + 1 and parts[0].isdigit():
                    scores = {}
                    for j, aa in enumerate(aa_order):
                        try:
                            val = parts[j+1]
                            if val == '*':
                                scores[aa] = -999
                            else:
                                # HMMER3: scores are -ln(prob/null), lower=better match
                                # Convert: higher score = better (negate)
                                scores[aa] = -float(val)
                        except (ValueError, IndexError):
                            scores[aa] = -4
                    pssm.append(scores)
            return (pssm, aa_order) if pssm else (None, None)
        except Exception:
            return None, None

    def generate_alphafold_input(self, orfs, indices):
        lines = []
        for idx in indices:
            if idx < len(orfs):
                orf = orfs[idx]; prot = orf['protein'].rstrip('*')
                lines.append(f">ORF{idx+1}|F{orf['frame']}{orf['strand']}|{orf['start']}-{orf['end']}|{len(prot)}aa")
                for i in range(0, len(prot), 80): lines.append(prot[i:i+80])
        return '\n'.join(lines)

    def score_candidate(self, orf, neighborhood):
        s = 0.0
        if orf.get('domains'): s += 0.15 * min(len(orf['domains']), 3)
        if neighborhood: s += 0.2 * neighborhood.get('cluster_score', 0)
        pl = len(orf.get('protein', '').rstrip('*'))
        if 200 <= pl <= 1500: s += 0.1
        elif 100 <= pl <= 200: s += 0.05
        gc = orf.get('gc', 50)
        if 35 <= gc <= 65: s += 0.05
        return min(s, 1.0)


# ═══════════════════════════════════════════════════════════════
# MODULE I: INTERNATIONALISATION (i18n) — v2.00 — English only
# ═══════════════════════════════════════════════════════════════

TRANSLATIONS = {
    # ── English ─────────────────────────────────────────────────
    'en': {
        'app_title':        '🧬 ppigFinder v2.00 — PPI Genomic Finder',
        'menu_file':        '📁 File',
        'menu_params':      '⚙️ Parameters',
        'menu_language':    '🌐 Language',
        'menu_help':        '❓ Help',
        # File menu
        'open_fasta':       '📂 Open FASTA',
        'open_multifasta':  '📂 Open Multi-FASTA',
        'open_snapgene':    '📂 Open SnapGene (.dna)',
        'open_genbank':     '📂 Open GenBank (.gb/.gbk)',
        'load_hmm':         '📎 Load HMM',
        'save_project':     '📦 Save Project',
        'open_project':     '📦 Open Project',
        'save_orfs_fasta':  '💾 Save ORFs (FASTA)',
        'save_cand_json':   '💾 Candidates (JSON)',
        'save_af3':         '💾 AlphaFold3',
        'save_report_tsv':  '📊 Report (TSV)',
        'save_report_html': 'Report (HTML)',
        'export_snapgene':  '🔬 Export as SnapGene (.dna)',
        'export_genbank':   '🔬 Export as GenBank (.gb)',
        'export_map_pdf':   '🖼️ Export Map as PDF',
        'quit':             '❌ Quit',
        # Params menu
        'blast_params':     '🔬 BLAST Parameters...',
        'hmm_params':       '🧬 HMM Parameters...',
        'reset_params':     '🔄 Reset Defaults',
        # Help menu
        'manual':           '📖 Manual',
        'tutorial':         '🎓 Tutorial',
        'install':          '⚙️ Installation Guide',
        'about':            'ℹ️ About',
        # Toolbar
        'btn_open':         '📂 Load a genome file',
        'btn_translate_genome': '🧬 Translate genome',
        'btn_pyrodigal':    '🧪 Pyrodigal',
        'btn_automatic':    '⚙️ Automatic',
        'btn_hybrid':       '🔀 Hybrid',
        'desc_hybrid':      'Pyrodigal (primary) + Automatic ORF scanner (gap-filler) — fills unannotated regions',
        'desc_pyrodigal':   'Gene prediction using Pyrodigal (Prodigal successor) - ML-based prokaryotic gene caller',
        'desc_automatic':   'Simple ORF detection using start/stop codons with size filters (30-5000 aa, frame +/-)',
        'btn_annotate_hmm': '🏷️ Annotate HMM',
        'btn_export_pdf':   '🖼️ Export Map PDF',
        'zoom_label':       'Zoom:',
        # Left panel
        'config_filters':   '⚙️ Config & Filters',
        'min_size_aa':      'Min Size (aa):',
        'start_codons':     'Start Codons:',
        'filters':          'Filters:',
        'search_lbl':       'Search:',
        'frame_lbl':        'Frame:',
        'strand_lbl':       'Strand:',
        'min_filt_aa':      'Min Filter (aa):',
        'apply_btn':        '▶ Apply',
        'genome_lbl':       '🧬 Genome:',
        'blast_lbl':        '🔬 BLAST:',
        'algorithm_lbl':    'Algorithm:',
        'identity_lbl':     'Min Identity (%):',
        'evalue_lbl':       'Max E-value:',
        'run_blast_btn':    '🚀 Run BLASTp',
        # Tabs (notebook right panel)
        'tab_dna':          '🧬 DNA',
        'tab_protein':      '🧪 Protein',
        'tab_domains':      '🏷️ Domains',
        'tab_neighbors':    '🗺️ Neighborhood',
        'tab_blast_query':  '📥 BLAST Query',
        'tab_blast_res':    '📊 BLAST Results',
        'tab_hmm':          '🧬 HMM',
        'tab_af3':          '🔮 AlphaFold',
        # Map
        'map_title':        '🎨 Genome Map (Ctrl+scroll zoom | Shift+drag pan)',
        'locate_orf':       '📍 Locate ORF in Genome',
        'search_btn':       '🔍 Search',
        'orfs_table':       '📋 ORFs',
        # Status
        'ready_status':     '✓ Ready — Load a genome FASTA',
        # Tooltips
        'tip_open':         'Open a genome file — FASTA, GenBank (.gb) or SnapGene (.dna)',
        'tip_analyze':      'Find all ORFs in all 6 reading frames (+/-) across the loaded genome',
        'tip_hmm':          'Apply HMM search results to label ORFs in the table and genome map',
        'tip_pdf':          'Export the map exactly as displayed — zoom, colors and annotations — as PDF/EPS',
        'tip_zoom_in':      'Zoom in on the genome map  (also: Ctrl + scroll wheel)',
        'tip_zoom_out':     'Zoom out on the genome map  (also: Ctrl + scroll wheel)',
        'tip_translate':    'Choose gene prediction method for identifying protein-coding sequences',
        'tip_search':       'Search ORFs by protein sequence or ORF number',
        'tip_frame':        'Filter ORFs by reading frame (+1, +2, +3, -1, -2, -3)',
        'tip_strand':       'Filter ORFs by strand (forward +, reverse -, or both)',
        'tip_min_aa':       'Minimum protein length (amino acids) to include in analysis',
        'tip_apply':        'Apply current search and filter settings to ORF table',
        # Right panel buttons tooltips
        'tip_blast_load':   'Load protein sequence from FASTA file for BLAST search',
        'tip_blast_paste':  'Paste protein sequence directly into the text area',  
        'tip_blast_clear':  'Clear the protein sequence text area',
        'tip_blast_validate': 'Check if the protein sequence format is valid',
        'tip_blast_run':    'Run BLAST search against all ORFs in the genome',
        'tip_blast_copy_hit': 'Copy selected BLAST hit to clipboard',
        'tip_blast_copy_all': 'Copy all BLAST results to clipboard',
        'tip_blast_save':   'Save BLAST results to TSV file',
        'tip_hmm_add':      'Add single HMM profile file for domain annotation',
        'tip_hmm_add_multi': 'Add multiple HMM profile files at once',
        'tip_hmm_search':   'Search all ORFs against loaded HMM profiles',
        'tip_af3_add_sel':  'Add currently selected ORF to AF3 prediction list',
        'tip_af3_add_hmm':  'Add all ORFs with HMM hits to AF3 prediction list',
        'tip_af3_add_all':  'Add ALL ORFs in the genome to AF3 prediction list (genome-wide interactome scan)',
        'tip_chain_copies': 'Number of copies of chain {letter} in the complex',
        'tip_af3_remove':   'Remove selected ORFs from AF3 prediction list',  
        'tip_af3_clear_all': 'Clear all ORFs from AF3 prediction list',
        'tip_af3_generate': 'Generate AF3 jobs based on selected prediction mode',
        'tip_af3_export_cf': 'Export jobs in ColabFold FASTA format',
        'tip_af3_ranking':  'Show ranking of AF3 predictions by confidence scores',
        'tip_af3_clear_jobs': 'Clear all generated AF3 jobs',
        'tip_af3_add_custom': 'Add custom multi-subunit complex job',
        'tip_search_orf':   'Search for specific ORF in the genome map',
        'tip_zoom_minus':   'Zoom out on the genome map',
        'tip_zoom_plus':    'Zoom in on the genome map',
        # Language names
        # ── AF3 tab ──
        'af3_sel_frame':        '📌 Select ORFs for Structure Prediction',
        'af3_add_sel':          '➕ Add Selected ORF',
        'af3_add_hmm':          '📋 Add HMM Hits',
        'af3_add_all':          '🧬 Add All ORFs',
        'af3_remove':           '🗑️ Remove',
        'af3_clear_all':        '🗑️ Clear All',
        'af3_jobs_frame':       '⚡ Generate Jobs',
        'af3_neighbors':        'Neighbors:',
        'af3_mode':             'Mode:',
        'af3_generate':         '⚡ Generate',
        'af3_export_json':      '💾 Export AF3 JSON',
        'af3_export_json_single': '📄 Individual JSONs',
        'af3_export_json_batch':  '📦 Batch JSON',
        'af3_export_slurm_array': '⚡ SLURM Array (anti-OOM)',
        'af3_export_cf':        '🧬 Export ColabFold FASTA',
        'af3_ranking':          '📊 Ranking',
        'af3_clear_jobs':       '🗑️ Clear Jobs',
        'af3_jobs_table':       '📋 Jobs',
        'hmm_add':              '📂 Add HMM Profile',
        'hmm_add_multi':        '📂 Add Multiple Profiles',
        'hmm_remove':           '🗑️ Remove Selected',
        'hmm_edit':             '🎨 Edit Color / Function',
        'hmm_search_all':       '🔍 Search All ORFs',
        'hmm_profiles_frame':   '📋 Loaded HMM Profiles',
        'blast_query_title':    '🔬 BLASTp — Paste your protein query sequence',
        'blast_load_fasta':     '📂 Load FASTA',
        'blast_paste':          '📋 Paste',
        'blast_clear':          '🗑️ Clear',
        'blast_validate':       '✅ Validate',
        'blast_run':            '🚀 Run BLASTp',
        'blast_formats':        'ℹ️ Accepted formats',
        'blast_formats_hint':   '• Raw sequence: MKTLLLFVLALF...   • FASTA: >header\nMKTLLL...',
        'blast_query_frame':    'Query Sequence (Protein)',
        'blast_copy_hit':       '📋 Copy Hit',
        'blast_copy_all':       '📋 Copy All',
        'blast_save':           '💾 Save Results',
        'copy_btn':             '📋 Copy',
        'copy_all_btn':         '📋 Copy All',
        'window_kb':            'Window (kb):',
        'analyze_btn':          '🔬 Analyze',
        'save_btn':             '✓ Save',
        'ok_btn':               '✓ OK',
        'reset_btn':            'Reset Defaults',
        'apply_btn2':           '✓ Apply',
        'color_hex':            'Hex:',
        'hmm_edit_name':        'Name:',
        'hmm_edit_func':        'Function:',
        'hmm_edit_color':       'Color:',
        'orf_obs':              'Observation:',
        'orf_putfunc':          'Putative Function:',
        'orf_gene':             'Gene Name:',
        'orf_notes':            'Notes:',
        'blast_dlg_title':      '🔬 BLAST Parameters',
        'blast_gen_params':     'General Parameters',
        'blast_max_targets':    'Max target sequences:',
        'blast_evalue':         'Expect threshold (E-value):',
        'blast_word_size':      'Word size:',
        'blast_scoring':        'Scoring Parameters',
        'blast_matrix':         'Matrix:',
        'blast_gap_open':       'Gap open cost:',
        'blast_gap_ext':        'Gap extend cost:',
        'blast_filters':        'Filters and Masking',
        'blast_low_complex':    'Filter low-complexity regions',
        'blast_id_filter':      'Identity Filter',
        'blast_min_id':         'Min identity (%):',
        'hmm_dlg_title':        '🧬 HMM Search Parameters',
        'hmm_thresh':           'Reporting Thresholds',
        'hmm_seq_evalue':       'Sequence E-value:',
        'hmm_dom_evalue':       'Domain E-value:',
        'hmm_score_thresh':     'Score threshold:',
        'hmm_backend':          'Backend',
        'hmm_params_passed':    'Parameters will be passed to hmmsearch',
        'hmm_not_found':        '❌ HMMER3 not found',
        'hmm_pssm_fallback':    'Using built-in PSSM scan (lower sensitivity)',
        'hmm_install_hint':     'Install: conda install -c bioconda hmmer',
        # Gene prediction
        'tip_pyrodigal':        'Predict genes using Pyrodigal (Prodigal algorithm) — more accurate gene calling',
        'pyrodigal_mode':       'Pyrodigal Mode:',
        'pyrodigal_meta':       'Metagenomic (no training)',
        'pyrodigal_single':     'Single genome (train on sequence)',
        'pyrodigal_not_avail':  'Pyrodigal not installed. Install: pip install pyrodigal',
        # BLAST Command Preview
        'blast_cmd_preview':    '⌨️ BLAST Command Preview (editable)',
        'blast_cmd_hint':       'Edit the command below before running, or use the parameter menus above.',
        'blast_program':        'Program:',
        'blast_run_custom':     '🚀 Run Custom Command',
        'blast_refresh_cmd':    '🔄 Refresh from Parameters',
    },

}

# Active language
_CURRENT_LANG = ['en']

def t(key: str) -> str:
    """Translate key to current language, fallback to English."""
    lang = _CURRENT_LANG[0]
    return (TRANSLATIONS.get(lang, {}).get(key)
         or TRANSLATIONS['en'].get(key)
         or key)


# ═══════════════════════════════════════════════════════════════
# MODULE II: HELP CONTENT
# ═══════════════════════════════════════════════════════════════

HELP_CONTENT = {
    'manual': {
        'en': """\
ppigFinder v2.00
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
User Manual

OVERVIEW
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
ppigFinder is a standalone desktop application for the discovery of novel
protein-protein interactions (PPIs) in bacterial and archaeal genomes.

The pipeline starts from a raw FASTA nucleotide sequence and proceeds
through five stages:

  1. Gene prediction — identify all protein-coding ORFs using one of
     three complementary prediction engines (see Section 2 below).
  2. Functional annotation — HMM profile scanning (HMMER3 / Pfam /
     TIGRFAM / custom .hmm) and BLASTp homology search.
  3. Visualisation — interactive zoomable genomic map (0.5× – 1,000,000×)
     with colour-coded ORF arrows and genomic neighbourhood analysis.
  4. Interaction prediction — AlphaFold 3 (AF3) batch job builder:
     generates AF3 JSON or ColabFold FASTA files for every desired pair;
     jobs can be submitted directly to an HPC cluster via SSH/SFTP.
  5. Result analysis — PAE heatmaps (ChimeraX colour scheme), per-residue
     pLDDT plots, ipTM / pTM scoring, and inter-chain contact detection
     imported back into the application for ranking and visual inspection.

WHAT'S NEW IN v1.12
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
• Hybrid prediction mode (Pyrodigal + Automatic gap-filler) — see §2.
• File menu streamlined: Multi-FASTA, SnapGene (.dna) and GenBank
  (.gb/.gbk) import/export options removed; FASTA is the sole input.
• Genome map zoom extended to 1,000,000× for large chromosomes / MAGs.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SECTION 1 — LOAD A GENOME
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
File → Open FASTA  —or—  toolbar "Load a genome file"
Accepted format: FASTA (.fasta .fa .fna), single sequence or first
record of a multi-record file.

After loading, the right-panel Genome tab shows sequence name, length,
GC%, and a summary of detected features. The genome map resets to 100%
zoom with the full chromosome visible.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SECTION 2 — GENE PREDICTION MODES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Toolbar → "Translate genome" opens a dropdown with three prediction
engines. The colum "Source" in the ORF table records which engine
produced each entry.

┌─────────────────────────────────────────────────────────────────┐
│ 🧪 PYRODIGAL                                   Source: pyrodigal│
├─────────────────────────────────────────────────────────────────┤
│ Uses the Prodigal / Pyrodigal dynamic-programming algorithm to  │
│ predict protein-coding genes. Evaluates GC-content, ribosome-  │
│ binding site (RBS) motifs, coding potential, and start-codon    │
│ context simultaneously, making it far more selective than a     │
│ simple codon scan.                                              │
│                                                                 │
│ Parameters (Parameters → ORF Analysis Parameters):             │
│   Mode        Metagenomic — uses pre-trained models; ideal for  │
│               short contigs, plasmids, and MAGs.               │
│               Single genome — trains on the input sequence;    │
│               best for closed chromosomes > 100 kb.            │
│   Table       Translation table for start/stop codon logic:    │
│               11 = standard bacteria & archaea (default)       │
│                4 = Mycoplasma / Mollicutes / SR1               │
│               25 = SR1 / Gracilibacteria                       │
│               15 = yeast mitochondria                          │
│   Min size    Minimum predicted protein length (aa). Pyrodigal │
│               internally enforces this as min_gene = aa × 3.  │
│   Closed ends Allow genes that begin or end at a sequence edge.│
│   Mask N runs Soft-mask regions of ambiguous nucleotides.      │
│   Post-filter Optionally restrict results to genes whose first │
│               codon is ATG, GTG, or TTG (any combination),     │
│               applied after Pyrodigal finishes.                │
│                                                                 │
│ When to use: standard bacterial/archaeal genomes; any case     │
│ where annotation accuracy is the top priority.                 │
│ Requires: pip install pyrodigal                                │
└─────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────┐
│ 🔀 HYBRID (Pyrodigal + gap-filler)             Source: mixed    │
├─────────────────────────────────────────────────────────────────┤
│ Addresses a known limitation of Pyrodigal: it occasionally     │
│ misses short genes, frameshifted ORFs, genes with unusual RBS  │
│ signals, or ORFs in low-complexity / GC-extreme regions. In    │
│ those unannotated intervals the 6-frame scanner is applied as  │
│ a gap-filler, recovering candidates that Pyrodigal skipped.    │
│                                                                 │
│ Internal pipeline:                                             │
│   Step 1 — Run Pyrodigal on the full sequence with the         │
│            parameters set in ORF Analysis Parameters.          │
│   Step 2 — Build a merged coverage map: every genomic position │
│            (0-based) spanned by at least one Pyrodigal ORF     │
│            (either strand) is marked as covered.               │
│   Step 3 — Identify contiguous uncovered gaps — runs of        │
│            positions not touched by any Pyrodigal prediction.  │
│   Step 4 — For each gap longer than min_aa × 3 + 3 nt, extract│
│            the subsequence and run the 6-frame start→stop      │
│            codon scanner using the parameters from the ORF     │
│            Detection section (min size and start codons).      │
│            Coordinates are re-mapped from gap-local back to    │
│            global genome space.                                │
│   Step 5 — Merge both ORF sets and sort by genomic start       │
│            position so ORF numbers (ORF1, ORF2 …) always       │
│            increase 5'→3', regardless of source.               │
│                                                                 │
│ Source column values after a Hybrid run:                       │
│   pyrodigal  — predicted by Pyrodigal (primary caller)         │
│   automatic  — found by the 6-frame scanner in a gap           │
│                                                                 │
│ Parameters used:                                               │
│   Pyrodigal section  → controls the primary caller             │
│   ORF Detection section → controls the gap-filler              │
│     (min size in aa, start codons: ATG / GTG / TTG)            │
│                                                                 │
│ When to use: when you suspect important short ORFs or unusual  │
│ genes are being missed by Pyrodigal alone; genomic islands,    │
│ phage insertions, toxin-antitoxin systems, small regulatory    │
│ peptides, or any region of biological interest showing a gap   │
│ in the Pyrodigal annotation.                                   │
│ Requires: pip install pyrodigal (Pyrodigal must be installed)  │
└─────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────┐
│ ⚙️  AUTOMATIC (6-frame scan)                   Source: 6frame   │
├─────────────────────────────────────────────────────────────────┤
│ Exhaustive 6-frame start→stop codon scan across the entire     │
│ sequence on both strands (frames +0, +1, +2, −0, −1, −2).     │
│ Every in-frame interval from a start codon to the next in-     │
│ frame stop codon that meets the minimum size threshold is      │
│ reported. No machine-learning, no training data required.      │
│                                                                 │
│ Parameters (Parameters → ORF Analysis Parameters):             │
│   Min ORF size   Minimum protein length in amino acids.        │
│                  Default: 30 aa (= 90 nt + stop codon).        │
│                  Lower values recover very short peptides but  │
│                  dramatically increase false positives.        │
│   Start codons   ATG (AUG) — canonical methionine start        │
│                  GTG (GUG) — valine start, common in bacteria  │
│                  TTG (UUG) — leucine start, less frequent      │
│                  All three are enabled by default; uncheck to  │
│                  restrict to canonical ATG-only starts.        │
│                                                                 │
│ When to use: quick exploratory scans; sequences from organisms │
│ with non-standard genetic codes; when Pyrodigal is not         │
│ installed; when you want to see every possible reading frame.  │
│ No external dependency needed.                                 │
└─────────────────────────────────────────────────────────────────┘

Recommendation for PPI discovery workflows:
  Use Hybrid mode as the default for bacterial genomes.
  It combines the biological accuracy of Pyrodigal with the
  exhaustiveness of the 6-frame scanner, ensuring no annotatable
  region is left blank.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SECTION 3 — ORF TABLE
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
All predicted genes are listed with the following columns:

  ID           ORF number (ORF1, ORF2 …), assigned by genomic
               start position — always increases 5'→3' regardless
               of prediction mode or strand.
  Frame        Reading frame (0–2 for +strand, 3–5 for −strand).
  Strand       + (sense) or − (antisense).
  Start / End  0-based genomic coordinates (half-open interval).
  Size(aa)     Protein length in amino acids (stop codon excluded).
  GC%          GC content of the coding DNA.
  HMM          Best HMM profile match name and domain region.
  Score        HMM bit-score or Pyrodigal confidence score (0–100).
  Source       Prediction engine: pyrodigal / automatic / 6frame.
               In Hybrid mode both values appear in the same table.
  Obs          Free-text observation / annotation field.
  AF3          ✅ if at least one AF3 prediction result is available.
  Partner      Best interaction partner from AF3 analysis.
  ipTM         Interface predicted TM-score of the best AF3 result.
  PAE_inter    Mean inter-chain PAE (Å) of the best AF3 result.

Interactions with the table:
  • Click a row       → select ORF; genome map centers on that gene.
  • Ctrl+click        → add/remove rows to/from selection.
  • Shift+click       → range selection.
  • Ctrl+C            → copy selected rows as TSV (Excel-ready).
  • Right-click       → context menu: annotate, color, copy sequences,
                        export FASTA, add to AlphaFold list.
  • Export table btn  → TSV (columns only) / TSV full (+DNA+Protein) /
                        TSV annotated only / FASTA protein / FASTA DNA.

Filtering bar:
  Search field  — filter by ORF ID or any protein substring.
  Frame / Strand / Min aa / Source  — additional combo filters.
  Apply button  — refresh the table with the current filter set.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SECTION 4 — HMM ANNOTATION
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
HMM tab (right panel):
  1. "Add HMM Profile" — load a single .hmm file.
     "Add Multiple Profiles" — load an entire folder of .hmm files.
  2. Assign a display colour and a short function label to each profile.
  3. "Search All ORFs" — runs HMMER3 hmmscan (or built-in PSSM scanner
     if HMMER3 is absent) against every predicted ORF protein sequence.
  4. Toolbar "Annotate HMM" — transfers hit information to the ORF table
     (HMM and Score columns) and recolours matching arrows on the map.

Compatible databases: Pfam-A, TIGRFAM, custom domain libraries.
Parameters (Parameters → HMM Parameters): E-value cutoff, bit-score
threshold, domain overlap handling.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SECTION 5 — BLAST HOMOLOGY SEARCH
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
BLAST Query tab:
  Paste or load a protein sequence (FASTA or raw) and click "Run BLAST".
  ppigFinder runs BLASTp against the local ORF database (all predicted
  proteins from the current genome). Hits appear in the BLAST Results tab
  with alignment details. Clicking any ORF link centers the genome map.

Parameters (Parameters → BLAST Parameters):
  Algorithm, scoring matrix (BLOSUM62 default), E-value cutoff,
  word size, gap penalties.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SECTION 6 — GENOMIC NEIGHBOURHOOD
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Neighborhood tab → set the window size (kb, each side) → Analyze.
Displays all ORFs within the window around the selected gene in a
scrollable FASTA list. Useful for operon context analysis and for
manually identifying candidate interaction partners encoded in the
same genomic region.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SECTION 7 — ALPHAFOLD 3 JOB BUILDER
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
AlphaFold tab — four sub-steps:

A) Build the ORF selection list:
   "Add Selected ORF"  — add only the currently highlighted ORF.
   "Add HMM Hits"      — add all ORFs with at least one HMM hit.
   "Add All ORFs"      — opens a filter dialog (size range, HMM
                         annotation required); auto-switches to
                         Interactoma Genômico mode.

B) Choose an interaction mode:
   Pares (hit vs vizinho)  — pair each selected ORF with its N
                              nearest genomic neighbours.
   Pares + Homodímeros     — as above, plus each ORF with itself.
   Trímeros                — three-way combinations.
   All vs All              — all pairwise combinations among selected.
   Hits HMM entre si       — pair ORFs that share an HMM profile.
   Hit vs all selected     — one query ORF against all others.
   Homodímero              — single ORF folded as a homodimer.
   Interactoma Genômico    — every selected ORF vs every genome ORF
                             (full interaction screen).

C) Generate jobs:
   Click "Generate" — the jobs table populates with pair names,
   total residue count, and initial status "pending".

D) Export:
   Export AF3 JSON     — individual JSONs (one per pair) or a
                         single Batch JSON file.
   Export ColabFold    — multi-sequence FASTA for ColabFold batch.
   SLURM Array         — for large interactomes: auto-splits into
                         batches of 50 pairs to prevent GPU OOM;
                         generates run_array.sh for sbatch submission.
   Download JSONs      — save all staged partition JSON files to a
                         local folder (useful for manual HPC transfer).

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SECTION 8 — HPC CLUSTER SUBMISSION
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"Submit AF3 via Server" tab — four sub-tabs:

Connect       Fill SSH credentials: hostname, username, port (22),
              password or key, remote base directory, AF3 run command,
              and module-load string. Click "Test connection".

Submit jobs   "Upload only" — transfer JSON files via SFTP.
              "Upload + Submit all" — transfer and sbatch in one step.
              Advanced: SLURM partition, model seeds, extra flags.
              Large batches are partitioned and submitted sequentially
              (one job per partition) to avoid queue limits.

Monitor       Live job status poll from SLURM / PBS / LSF.
              Click "Refresh" to update the status column.

Results       Download completed prediction folders via SFTP.
              "Auto-import" parses ipTM and pLDDT from AF3 output
              JSONs and updates the AlphaFold Analysis ranking.

DaVinci (ICB/USP) quick reference:
  --partition=basic   max 72 h,   16 CPU,  100 GB RAM,  0 GPU
  --partition=max50   max  8 d,   64 CPU,  500 GB RAM,  1 GPU  ← AF3
  --partition=max90   max 15 d,  110 CPU,    1 TB RAM,  4 GPU
  AF3 module: module load alphafold3
  Run command: af3_run

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SECTION 9 — ALPHAFOLD 3 RESULT ANALYSIS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
AlphaFold Analysis tab → "Load AF3 results folder".

ppigFinder recursively scans for AF3 output JSON files and parses:
  ipTM         — interface predicted TM-score (key interaction metric)
  pTM          — predicted TM-score for the full complex
  mean_pLDDT   — mean per-residue local distance difference test score
  PAE_inter    — mean predicted aligned error across inter-chain pairs
  Contact pair — residue pair with lowest inter-chain PAE value

All predictions are listed in a sortable table. Clicking a row:
  • Renders the PAE heatmap (ChimeraX colour scheme: blue = confident,
    yellow/red = uncertain)
  • Renders the per-residue pLDDT bar plot
  • Centers the genome map on the query ORF
  • Lists inter-chain contacts above the configurable PAE threshold (Å)

Confidence guide:
  ipTM >= 0.75        high-confidence interaction — likely physical
  ipTM 0.50 – 0.75   moderate confidence — worth experimental follow-up
  ipTM < 0.40         low confidence — unlikely direct interaction
  PAE_inter < 10 Å   confident inter-chain contact geometry
  mean_pLDDT > 70    well-structured complex model

Export: "Export plots PDF" saves all PAE heatmaps and pLDDT plots
as a multi-page PDF for reporting.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SECTION 10 — GENOME MAP NAVIGATION
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
The genome map renders ORFs as directional arrows on a horizontal
backbone. Arrow colour reflects HMM annotation (if Annotate HMM was
run) or strand (+ / −).

Zoom range: 0.5× (full chromosome overview) to 1,000,000×
(single-nucleotide resolution for fine inspection).

  Ctrl + scroll wheel    smooth zoom anchored at cursor position
  Toolbar − / +          step zoom ×0.8 (out) / ×1.2 (in)
  Shift + drag           horizontal pan
  Search box             jump to ORF number or protein substring
  Click ORF arrow        select ORF; table scrolls to row; right
                         panel shows DNA / Protein / Domains tabs

When an ORF row is selected in the table, the map always re-centers
on that gene at the current zoom level — this behaviour is preserved
in all prediction modes including Hybrid.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SECTION 11 — PROJECT SAVE / LOAD
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
File → Save Project  saves a single JSON workspace (schema v2)
containing the complete session state:
  • Full genome sequence
  • All ORFs with every annotation field: gene name, putative function,
    observation, notes, custom colour, source (pyrodigal / automatic)
  • HMM profiles loaded and all hit results
  • AlphaFold selection list and all generated jobs
  • AF3 analysis results: PAE matrices, pLDDT arrays, ipTM scores
  • BLAST results and query history
  • Pyrodigal parameters and ORF scanner settings
  • HPC server connection settings
  • UI state: active filters, zoom level, start codon checkboxes

File → Open Project   restores the entire session.
File → Save Project As  makes a fully independent snapshot copy.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SECTION 12 — EXPORT OPTIONS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
ORF table:
  Ctrl+C              copy selected rows as TSV
  Export table menu:
    TSV (columns)     ID, coordinates, size, GC%, HMM, score, source…
    TSV full          above + full DNA and protein sequences
    TSV annotated     only rows with a non-empty Obs / gene name field
  Export FASTA:
    Protein FASTA     rich header: >ORFn|frame|start-end|Naa|source
    DNA FASTA         nucleotide coding sequence with same header

Genome map:
  "Export Map PDF" toolbar button → saves as PNG, PDF, or SVG

AF3 analysis:
  "Export plots PDF" → multi-page PDF of all PAE + pLDDT figures

Project and reports:
  File → Save Project         full JSON workspace
  File → Report (TSV)         legacy ranked ORF table export

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SECTION 13 — KEY SHORTCUTS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Ctrl+C              copy selected ORF table rows as TSV
Ctrl+scroll         zoom genome map
Shift+drag          pan genome map
Right-click (ORF)   annotate / colour / copy / export / add to AF3
Right-click (AF3)   delete selected job
Delete/Backspace    delete selected AF3 jobs (table must be focused)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SECTION 14 — DEPENDENCIES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Required (core GUI and plotting):
  PyQt6 >= 6.4       pip install PyQt6     (PyQt5 >= 5.15 fallback)
  matplotlib >= 3.5  pip install matplotlib
  numpy >= 1.21      pip install numpy

Recommended (gene prediction):
  pyrodigal >= 2.0   pip install pyrodigal
  Required for Pyrodigal and Hybrid modes.

Optional (annotation and submission):
  BLAST+ >= 2.12     https://ftp.ncbi.nlm.nih.gov/blast/executables/
                     Also detectable via WSL on Windows.
  HMMER3 >= 3.3      conda install -c bioconda hmmer
                     Or local install; also via WSL.
  paramiko >= 2.9    pip install paramiko
                     Required for HPC SSH/SFTP submission.

Backend status is shown in real time in the Genome tab (right panel)
and in the toolbar badge (✅ / ❌ per backend).
""",
    },

    'tutorial': {
        'en': """\
ppigFinder v2.00 — Step-by-Step Tutorial
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
A practical guide for novel bacterial PPI discovery using the full
ppigFinder pipeline: gene prediction → HMM / BLAST annotation →
genomic neighbourhood → AlphaFold 3 submission → result analysis.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
STEP 1 — Load a genome
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
1. Click "Load a genome file" in the toolbar  —or—  File → Open FASTA.
2. Select a FASTA file (.fasta .fa .fna).
3. The right-panel Genome tab immediately shows: sequence name, total
   length (bp), GC%, and number of ambiguous nucleotides. The genome
   map renders the backbone at 100% zoom (full sequence visible).

  Tip: If your sequence is in GenBank or SnapGene format, export it to
  FASTA first using any sequence editor (Benchling, Geneious, SnapGene
  Viewer, Biopython: SeqIO.write, or command-line: seqret).

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
STEP 2 — Configure prediction parameters
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Before running any prediction, set the parameters once via:
  Parameters → ORF Analysis Parameters

This dialog has two independent sections that feed two different engines:

  Section A — 6-frame ORF scanner (used by Automatic and Hybrid gap-filler)
  ────────────────────────────────────────────────────────────────────────
  Min ORF size (aa)  The scanner reports only ORFs whose translated
                     protein is at least this many amino acids long
                     (stop codon not counted). Default: 30 aa.
                     Raise to 50–100 to reduce noise; lower to 15–20
                     to capture very small peptides (more false positives).
  Start codons       Check ATG, GTG, and/or TTG. ATG is canonical;
                     GTG and TTG are common alternative starts in
                     bacteria (e.g. many E. coli genes start with GTG).
                     Uncheck GTG/TTG for strict ATG-only scanning.

  Section B — Pyrodigal (used by Pyrodigal and Hybrid primary caller)
  ────────────────────────────────────────────────────────────────────
  Mode        Metagenomic — recommended for most use cases. Uses
              built-in pre-trained models; works on any contig size.
              Single genome — re-trains on the input sequence; requires
              a reasonably complete chromosome (> 100 kb recommended)
              for reliable model training.
  Table       Choose translation table matching your organism:
              11 = bacteria/archaea standard (most common)
               4 = Mycoplasma, Mollicutes, some phages
              25 = SR1 / Gracilibacteria (opal = Trp)
              15 = yeast mitochondria (TAG = Gln)
  Min gene    Minimum protein length for Pyrodigal predictions.
              Independent of the 6-frame scanner min size — in Hybrid
              mode you can set a shorter gap-filler threshold to
              recover small ORFs that Pyrodigal deliberately skips.
  Closed ends If checked, allows Pyrodigal to predict genes that
              run off the edges of the input sequence (partial genes
              at contig boundaries).
  Mask N runs Soft-masks runs of N/n nucleotides before prediction
              to avoid spurious gene calls in assembly gaps.
  Post-filter Optional: after Pyrodigal finishes, keep only genes
              whose first codon matches the checked subset (ATG /
              GTG / TTG / All). "All" (default) makes no restriction.

Click OK to save. Parameters persist for the entire session and are
stored in the project file.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
STEP 3 — Predict genes (choose a mode)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Click "Translate genome" in the toolbar to open the mode dropdown.

───────────────────────────────────────
 🧪 PYRODIGAL
───────────────────────────────────────
Best for: complete or nearly complete bacterial/archaeal chromosomes
where annotation accuracy is the priority.

Pyrodigal applies dynamic programming to score every possible gene on
both strands simultaneously, using GC-content, RBS motif probability,
coding potential, and start-codon context. This produces a clean,
biologically motivated gene set with very low false-positive rates.

All ORFs are tagged Source = pyrodigal.

After running: the ORF table fills and the genome map shows coloured
arrows. Status bar reports: "Pyrodigal: N genes | mode=… table=…".

───────────────────────────────────────
 🔀 HYBRID  ← recommended default
───────────────────────────────────────
Best for: any bacterial genome where completeness matters — especially
when searching for novel interaction partners that may include small or
atypical genes that Pyrodigal tends to miss.

Why Hybrid? Pyrodigal is highly accurate but conservative: it skips
regions with weak RBS signals, very short ORFs (< 60 aa), genes in
repetitive or GC-extreme contexts, and sometimes frameshifted or
overlapping genes. These "dark zones" may encode relevant proteins
(toxin-antitoxin components, small regulatory peptides, phage-related
proteins, signal peptides, etc.).

How it works:
  1. Pyrodigal runs first on the full sequence → primary ORF set.
  2. A merged coverage map is built from Pyrodigal intervals.
     Every nucleotide position covered by at least one Pyrodigal ORF
     (on either strand) is marked. Intervals are merged so that
     overlapping predictions count as a single covered block.
  3. Uncovered gaps (contiguous regions with zero Pyrodigal coverage)
     are identified. Only gaps longer than min_aa × 3 + 3 nt are
     processed (smaller gaps cannot contain a valid ORF anyway).
  4. The 6-frame scanner is applied to each gap subsequence using
     the ORF Detection parameters (min size, start codons).
     Coordinates are translated back to global genome space.
  5. Both ORF sets are merged and sorted by genomic start position,
     so ORF numbering (ORF1, ORF2 …) always follows the 5'→3' order
     of the chromosome regardless of prediction source.

Result: a unified ORF table where Pyrodigal genes carry
Source = pyrodigal and gap-filled genes carry Source = automatic.
You can filter the table by Source to inspect each subset separately.

Status bar after Hybrid run:
  "✓ Hybrid: N ORFs total — X pyrodigal + Y gap-fill (automatic)"

───────────────────────────────────────
 ⚙️ AUTOMATIC
───────────────────────────────────────
Best for: quick exploratory scans, organisms with non-standard genetic
codes, validation runs, or when Pyrodigal is not installed.

Exhaustive 6-frame scan: every interval from a start codon to the next
in-frame stop codon (all 6 reading frames, both strands) that meets the
minimum size threshold is reported. No biological context is used —
purely syntactic. Expect a higher number of ORFs than Pyrodigal, many
of which will be spurious. Use the HMM and BLAST tools to triage.

All ORFs are tagged Source = 6frame.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
STEP 4 — Explore the ORF table and genome map
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
After prediction, the ORF table on the left lists every gene. The
genome map draws directional arrows on a horizontal chromosome backbone.

Navigation:
  • Click any row → genome map centers on that arrow; right panel
    shows the DNA sequence (DNA tab), translated protein (Protein tab),
    and predicted domains (Domains tab).
  • Ctrl+scroll or toolbar − / + → zoom the map (0.5× – 1,000,000×).
  • Shift+drag → pan left/right.
  • Type in the search box → filter by ORF ID or protein substring.

Annotation via right-click on a row:
  • Annotate  — opens a dialog to set gene name, putative function,
                free-text observation, and notes. These fields are
                saved in the project file and exported in TSV/FASTA.
  • Color     — assigns a custom colour to the arrow on the genome map.
  • Copy      — copies the protein sequence, DNA sequence, or full
                FASTA header + sequence to the clipboard.
  • Add to AlphaFold  — queues the ORF for interaction prediction.

Multi-selection: Ctrl+click or Shift+click to select multiple rows.
Ctrl+C copies all selected rows as a tab-separated table (Excel-ready).

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
STEP 5 — HMM annotation
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
HMM scanning matches predicted ORFs against domain profile databases
(Pfam, TIGRFAM, custom), revealing functional families and guiding
which ORFs are worth submitting to AlphaFold.

Workflow:
  1. Open the HMM tab (right panel).
  2. Click "Add HMM Profile" (one .hmm file) or "Add Multiple Profiles"
     (entire folder — useful for full Pfam-A or TIGRFAM libraries).
  3. Assign a colour and short label (e.g. "ToxIN", "ABC transporter")
     to each profile for visual identification on the map.
  4. Click "Search All ORFs" — runs HMMER3 hmmscan against all proteins.
     If HMMER3 is absent, the built-in PSSM scanner is used as fallback.
  5. When the search finishes, click "Annotate HMM" in the toolbar.
     Matching ORFs are recoloured on the map and the HMM and Score
     columns in the ORF table are filled.

  Tip: run HMM annotation after a Hybrid prediction — both pyrodigal
  and automatic ORFs are scanned equally. Hits in gap-filled ORFs are
  particularly informative as evidence of genuine missed genes.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
STEP 6 — BLAST homology search
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
BLAST Query tab:
  1. Paste or load a protein sequence (FASTA or raw amino-acid string)
     into the query box.
  2. Click "Run BLAST". ppigFinder runs BLASTp against all predicted
     ORF protein sequences from the current genome.
  3. Results appear in the BLAST Results tab with e-value, bit-score,
     identity%, and aligned positions. Clicking an ORF link selects
     that row in the table and centers the genome map on it.

  Use case: paste a known interaction partner (e.g. a two-component
  system kinase) to find candidate receiver domains or HAMP linkers
  in the same genome.

Parameters → BLAST Parameters: algorithm (blastp/tblastn), scoring
matrix, E-value threshold, word size, gap open/extend penalties.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
STEP 7 — Build AlphaFold 3 jobs
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
AlphaFold tab — build the list of protein pairs to predict.

A) Populate the ORF selection list:
   "Add Selected ORF"  — adds only the currently highlighted ORF.
   "Add HMM Hits"      — adds all ORFs with at least one HMM match.
                         Ideal after HMM annotation to focus on
                         functionally annotated candidates.
   "Add All ORFs"      — opens a filter dialog (min/max size, require
                         HMM hit). Automatically switches mode to
                         Interactoma Genômico for genome-wide screen.

B) Choose the interaction generation mode:
   Pares (hit vs vizinho) — each ORF paired with its N nearest
                             genomic neighbours (configurable N).
                             Best for operon-level PPI hypotheses.
   Pares + Homodímeros    — as above, plus each ORF paired with
                             itself (homodimer prediction).
   Trímeros               — all three-way combinations.
   All vs All             — all pairwise combinations among selected
                             ORFs. Scales as N²/2 — use with caution
                             for large selection sets.
   Hits HMM entre si      — pair ORFs that share the same HMM profile
                             (e.g. all ABC-ATPase domains with each
                             other).
   Hit vs all selected    — one query ORF against every other
                             selected ORF.
   Homodímero             — a single ORF folded as a homodimer.
   Interactoma Genômico   — every selected ORF vs every ORF in the
                             genome. Full interaction screen — can
                             produce thousands of jobs for large
                             genomes; use SLURM Array export.

C) Generate and review:
   Click "Generate". The jobs table fills with pair names, total
   residue count (sum of both chains), and initial status "pending".
   Review the list — delete unwanted pairs via right-click or
   Delete / Backspace key (table must be focused).

D) Export jobs:
   Export AF3 JSON       — individual JSON files (one per pair) or a
                           single Batch JSON (all pairs in one file).
   Export ColabFold      — multi-sequence FASTA for ColabFold batch.
   SLURM Array (anti-OOM)— auto-splits large batches into partitions
                           of 50 pairs; generates run_array.sh for
                           single sbatch command; prevents GPU OOM.
   Download JSONs        — save all staged partition JSON files to a
                           local folder of your choice, bypassing the
                           SSH upload (manual HPC transfer workflow).

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
STEP 8 — Submit to an HPC cluster
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"Submit AF3 via Server" tab — four sub-tabs:

Connect
  Fill in: hostname (e.g. davinci.icb.usp.br), username, port (22),
  password, remote base directory (will hold all job subdirectories),
  AF3 run command (e.g. af3_run), module-load string
  (e.g. alphafold3). Click "Test connection" to verify SSH/SFTP access.

Submit jobs
  "Upload only"         — transfers JSON files to the remote directory
                          via SFTP without submitting to the scheduler.
  "Upload + Submit all" — transfers and calls sbatch / qsub / bsub in
                          a single click.
  Advanced options: --partition, model seeds, extra SLURM flags.
  Large batches are submitted as sequential partitions (one job per
  partition in the queue) to respect scheduler limits.

Monitor
  Click "Refresh" to poll the scheduler (squeue / qstat / bjobs)
  and update the status column for each submitted job.

Results
  Browse remote directories, select completed job folders, and
  click "Download". ppigFinder downloads via SFTP and optionally
  auto-imports ipTM and pLDDT values into the Ranking tab.

DaVinci (ICB/USP) cluster reference:
  --partition=basic   max  72 h,  16 CPU,  100 GB,  0 GPU
  --partition=max50   max   8 d,  64 CPU,  500 GB,  1 GPU  ← use for AF3
  --partition=max90   max  15 d, 110 CPU,    1 TB,  4 GPU
  Module:  module load alphafold3
  Command: af3_run

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
STEP 9 — Analyse AlphaFold 3 predictions
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
AlphaFold Analysis tab → "Load AF3 results folder".

ppigFinder recursively scans the selected folder for AF3 output JSON
files and automatically parses:
  ipTM          interface predicted TM-score — the primary metric
                for assessing whether two proteins physically interact.
  pTM           predicted TM-score for the entire complex model.
  mean_pLDDT    average per-residue pLDDT across both chains —
                measures overall model confidence.
  PAE_inter     mean predicted aligned error between residues of
                different chains — measures confidence in the
                relative positions of the two proteins.
  Contact pair  the single inter-chain residue pair with the
                lowest PAE value (most geometrically confident
                contact in the interface).

Interpreting results:
  ipTM >= 0.75       high-confidence interaction — structural model
                     suggests direct physical binding.
  ipTM 0.50 – 0.75  moderate confidence — candidate worth experimental
                     validation (pull-down, Y2H, co-IP).
  ipTM < 0.40        low confidence — direct interaction unlikely in
                     this structural context.
  PAE_inter < 10 Å   the inter-chain geometry is well-defined.
  mean_pLDDT > 70    the complex model is overall well-structured.

Clicking a result row:
  • The PAE heatmap renders (ChimeraX colour scheme: blue = low PAE =
    high confidence; yellow / red = high PAE = uncertain geometry).
    The diagonal blocks represent intra-chain PAE; off-diagonal blocks
    represent inter-chain PAE — the key interaction signal.
  • The per-residue pLDDT bar plot renders; low-pLDDT stretches
    indicate disordered or poorly modelled regions.
  • The genome map centers on the query ORF of that prediction.
  • Inter-chain contacts above the threshold (Å, adjustable) are listed.

Export: "Export plots PDF" → multi-page PDF of all PAE and pLDDT
figures for the selected predictions.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
STEP 10 — Navigate the genome map
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
The genome map is central to the ppigFinder workflow — it provides
spatial context for every ORF and supports rapid visual inspection.

Zoom range: 0.5× (full chromosome) → 1,000,000× (nucleotide level).
This range is sufficient for any bacterial or archaeal chromosome,
including large genomes > 5 Mbp and closed plasmids.

  Ctrl + scroll wheel   smooth zoom centred on cursor position.
  Toolbar − / +         step zoom ×0.8 / ×1.2 per click.
  Shift + drag          pan the chromosome left or right.
  Search box            type an ORF number (e.g. "ORF42") or any
                        substring of a protein sequence to jump to
                        matching ORFs. Click Search or press Enter.
  Click an ORF arrow    selects that ORF in the table; right panel
                        switches to show sequences and domains.

Selecting a row in the ORF table always re-centers the map on that
gene at the current zoom level. This is true for all prediction modes
(Pyrodigal, Hybrid, Automatic) — the centering logic is source-agnostic.

  Tip (Hybrid mode): after running Hybrid, zoom to a region of interest
  and look for arrows coloured differently from the Pyrodigal set —
  these are the gap-filled (automatic) ORFs. Check the Source column
  to confirm, then run HMM / BLAST on them to assess biological relevance.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
STEP 11 — Save the project
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
File → Save Project   — saves a single .json workspace.
File → Open Project   — restores the complete session.
File → Save Project As — independent snapshot copy (all data inlined).

The project file preserves:
  ✓ Genome sequence and name
  ✓ All ORFs with every annotation field (gene name, function,
    observation, notes, custom colour, source tag)
  ✓ HMM profiles and all search results
  ✓ AlphaFold selection list and generated job definitions
  ✓ AF3 analysis results: PAE matrices, pLDDT arrays, ipTM / pTM scores
  ✓ BLAST query history and result hits
  ✓ Pyrodigal parameters and ORF scanner settings (both sections)
  ✓ HPC server credentials and connection settings
  ✓ UI state: current zoom level, active Source/Frame/Strand filters,
    start codon checkboxes, min aa spinner

  Tip: save the project immediately after a Hybrid prediction so you
  can restore the full ORF set (including gap-filled ORFs) without
  re-running the prediction.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
STEP 12 — Export data
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
ORF table exports:
  Ctrl+C                  copy selected rows as TSV (paste to Excel).
  Export table → TSV      all visible columns for every ORF (or
    (columns only)        filtered subset).
  Export table →          above + full DNA and protein sequences
    TSV full              appended as extra columns.
  Export table →          only rows where Obs / gene name is filled;
    TSV annotated         useful for reporting confirmed candidates.
  Export FASTA (protein)  >ORFn|Fframe|start–end|Naa|source — one
                          record per ORF; stop codon excluded.
  Export FASTA (DNA)      same header, nucleotide coding sequence.

Genome map:
  "Export Map PDF" toolbar → saves the current map view as PNG, PDF,
  or SVG; resolution and DPI configurable.

AF3 plots:
  "Export plots PDF" → multi-page PDF of selected PAE heatmaps and
  pLDDT bar plots.

Reports:
  File → Report (TSV)     legacy ranked ORF export for spreadsheets.
  File → Save Project     full JSON workspace (all data).
""",
    },

    'install': {
        'en': """\
ppigFinder v2.00 — Installation Guide
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
How to install ppigFinder and all its dependencies on a regular
personal computer — Windows, macOS, and Linux.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
DEPENDENCY OVERVIEW
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
The table below lists every component, whether it is required or
optional, and what breaks without it.

  Package / Tool    Ver      Type      Without it
  ──────────────────────────────────────────────────────────────
  Python            >= 3.8   REQUIRED  app cannot run
  PyQt6             >= 6.4   REQUIRED  no graphical interface
    (PyQt5 >= 5.15 is accepted as automatic fallback)
  matplotlib        >= 3.5   REQUIRED  PAE/pLDDT plots disabled
  numpy             >= 1.21  REQUIRED  PAE/pLDDT plots disabled
  pyrodigal         >= 2.0   OPTIONAL  Pyrodigal & Hybrid modes
                                       unavailable; Automatic mode
                                       still works
  paramiko          >= 2.9   OPTIONAL  HPC SSH/SFTP submission
                                       unavailable; all other
                                       features still work
  NCBI BLAST+       >= 2.12  OPTIONAL  uses built-in k-mer/SW
                                       aligner as fallback
  HMMER3            >= 3.3   OPTIONAL  uses built-in PSSM scanner
                                       as fallback

  On Windows, BLAST+ and HMMER3 are most easily installed through
  WSL (Windows Subsystem for Linux) — see Section 3 below.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SECTION 1 — WINDOWS (without Spyder / Anaconda)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
This route uses the standard Python installer from python.org and
installs packages through pip in a virtual environment.

── 1.1  Install Python ────────────────────────────────────────────
1. Open your browser and go to  https://www.python.org/downloads/
2. Download the latest Python 3.x installer (e.g. Python 3.12.x).
3. Run the installer.
   IMPORTANT: on the first screen, check the box
   "Add Python to PATH" before clicking Install Now.
4. When installation finishes, open the Start menu, search for
   "Command Prompt" or "PowerShell" and open it.
5. Verify the installation:
     python --version
   Expected output: Python 3.12.x  (or similar 3.8+)

── 1.2  Create a virtual environment (recommended) ────────────────
A virtual environment keeps ppigFinder's packages isolated from
other Python projects on your machine.

In the Command Prompt / PowerShell:
  cd %USERPROFILE%\\Desktop
  python -m venv ppigfinder_env
  ppigfinder_env\\Scripts\\activate

Your prompt will now show (ppigfinder_env) at the beginning.
You must activate this environment every time before running the app.

── 1.3  Install core Python packages ──────────────────────────────
With the environment active:
  pip install --upgrade pip
  pip install PyQt6 matplotlib numpy

  Optional but strongly recommended:
  pip install pyrodigal paramiko

  Verify:
  python -c "import PyQt6; import matplotlib; import numpy; print('OK')"
  python -c "import pyrodigal; print('pyrodigal', pyrodigal.__version__)"

── 1.4  Run ppigFinder ────────────────────────────────────────────
Place the file ppigfinderv1_12.py anywhere convenient (e.g. Desktop).
With the virtual environment active:
  python C:\\Users\\YourName\\Desktop\\ppigfinderv1_12.py

To avoid typing this every time, create a one-line batch file
(run_ppigfinder.bat) on your Desktop containing:
  @echo off
  call %USERPROFILE%\\Desktop\\ppigfinder_env\\Scripts\\activate
  python %USERPROFILE%\\Desktop\\ppigfinderv1_12.py

Double-clicking that .bat file will open the app.

── 1.5  Install BLAST+ on Windows ─────────────────────────────────
Option A — Native Windows installer (simplest):
  1. Go to: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
  2. Download the .exe installer for Windows
     (e.g. ncbi-blast-2.16.0+-win64.exe).
  3. Run the installer; accept defaults.
  4. IMPORTANT: the installer should add BLAST to your PATH.
     To verify, open a NEW Command Prompt and type:
       blastp -version
     Expected: blastp: 2.16.0+

Option B — Via WSL (see Section 3): recommended if you also need HMMER3.

── 1.6  Install HMMER3 on Windows (requires WSL) ──────────────────
HMMER3 does not have a native Windows build. Install WSL first
(Section 3), then inside WSL run:
  sudo apt-get update
  sudo apt-get install -y hmmer
  hmmscan -h | head -2     # verify

ppigFinder automatically detects HMMER3 running through WSL.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SECTION 2 — WINDOWS (with Spyder / Anaconda)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Anaconda ships Python, Spyder IDE, and conda package manager together.
Use this route if you already use Anaconda for data science work.

── 2.1  Install Anaconda ──────────────────────────────────────────
If Anaconda is not yet installed:
  1. Go to: https://www.anaconda.com/download
  2. Download the Windows installer and run it.
  3. During installation, select "Add Anaconda to PATH" (or use
     Anaconda Prompt for all commands below).

── 2.2  Create a dedicated conda environment ──────────────────────
Open "Anaconda Prompt" from the Start menu:
  conda create -n ppigfinder python=3.11 -y
  conda activate ppigfinder

Your prompt should now show (ppigfinder).
You must run "conda activate ppigfinder" each time before using the app.

── 2.3  Install Python packages ───────────────────────────────────
  conda install -c conda-forge pyqt matplotlib numpy -y
  pip install pyrodigal paramiko

  Verify:
  python -c "from PyQt6.QtWidgets import QApplication; print('PyQt6 OK')"
  python -c "import pyrodigal; print('pyrodigal', pyrodigal.__version__)"

  If PyQt6 is not found via conda-forge, fall back to pip:
  pip install PyQt6

── 2.4  Install BLAST+ and HMMER3 via conda ───────────────────────
  conda install -c bioconda blast hmmer -y
  blastp -version    # verify
  hmmscan -h | head -2

  Note: bioconda packages work natively on Windows only in certain
  conda configurations. If the install fails, use WSL (Section 3)
  or the native BLAST+ installer (Section 1.5).

── 2.5  Run ppigFinder from Spyder ────────────────────────────────
1. Open Spyder from the Anaconda Navigator (make sure the ppigfinder
   environment is selected in the top-right environment dropdown).
   — or — from Anaconda Prompt with (ppigfinder) active:
   conda install spyder -y && spyder

2. In Spyder, open ppigfinderv1_12.py (File → Open).
3. Press F5 (Run file) or click the green ▶ button.
4. The ppigFinder window will open alongside Spyder.

   Tip: In Spyder Preferences → Run → "Run in external terminal"
   will open ppigFinder as a separate window rather than inside
   Spyder's IPython console, which avoids event-loop conflicts.

── 2.6  Run ppigFinder from Anaconda Prompt ───────────────────────
  conda activate ppigfinder
  python C:\\path\\to\\ppigfinderv1_12.py

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SECTION 3 — WINDOWS SUBSYSTEM FOR LINUX (WSL)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
WSL runs a real Linux environment inside Windows. ppigFinder detects
BLAST+ and HMMER3 installed in WSL automatically via the backend
detection code ("via WSL" label in the Backends panel).

── 3.1  Enable WSL ────────────────────────────────────────────────
Windows 10 (version 2004+) and Windows 11:
1. Open PowerShell as Administrator (right-click → Run as administrator).
2. Run:
     wsl --install
3. Restart your computer when prompted.
4. After restart, Ubuntu will finish installing and ask you to create
   a Linux username and password (this is independent of your Windows
   account).

To verify WSL is working:
  wsl --list --verbose     (should show Ubuntu running)

── 3.2  Install BLAST+ inside WSL ────────────────────────────────
Open the Ubuntu app from the Start menu (or type "wsl" in Powershell):
  sudo apt-get update
  sudo apt-get install -y ncbi-blast+
  blastp -version    # should print blastp: 2.x.x+

── 3.3  Install HMMER3 inside WSL ────────────────────────────────
  sudo apt-get install -y hmmer
  hmmscan -h | head -2    # should print HMMER 3.x

── 3.4  Verify ppigFinder detects WSL tools ───────────────────────
Launch ppigFinder on Windows (Sections 1 or 2). In the right panel,
look at the Backends section:
  BLAST+     ✅  via WSL
  HMMER3     ✅  via WSL

If you see ❌, open a Command Prompt and test:
  wsl bash -c "blastp -version"
  wsl bash -c "hmmscan -h"
These commands must return output for ppigFinder to detect them.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SECTION 4 — macOS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
── 4.1  Install Homebrew (package manager) ────────────────────────
Homebrew is the most convenient way to install system tools on macOS.
Open Terminal (Applications → Utilities → Terminal) and run:
  /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

Follow the on-screen instructions. On Apple Silicon (M1/M2/M3) Macs
the installer may ask you to add Homebrew to your PATH:
  echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> ~/.zprofile
  eval "$(/opt/homebrew/bin/brew shellenv)"

Verify:
  brew --version

── 4.2  Install Python ────────────────────────────────────────────
macOS ships a system Python (often outdated). Install a fresh one:
  brew install python@3.12

Then verify:
  python3 --version    # should be 3.12.x

── 4.3  Create a virtual environment ──────────────────────────────
  cd ~/Desktop
  python3 -m venv ppigfinder_env
  source ppigfinder_env/bin/activate

Your prompt will show (ppigfinder_env). Activate it every session:
  source ~/Desktop/ppigfinder_env/bin/activate

── 4.4  Install Python packages ───────────────────────────────────
  pip install --upgrade pip
  pip install PyQt6 matplotlib numpy
  pip install pyrodigal paramiko    # optional, strongly recommended

  Verify:
  python -c "import PyQt6; import matplotlib; print('GUI packages OK')"
  python -c "import pyrodigal; print(pyrodigal.__version__)"

── 4.5  Install BLAST+ ────────────────────────────────────────────
  brew install blast
  blastp -version    # verify

── 4.6  Install HMMER3 ────────────────────────────────────────────
  brew install hmmer
  hmmscan -h | head -2    # verify

── 4.7  Run ppigFinder ────────────────────────────────────────────
  source ~/Desktop/ppigfinder_env/bin/activate
  python ~/Desktop/ppigfinderv1_12.py

To create a clickable launcher, save the following as
run_ppigfinder.command on your Desktop:
  #!/bin/bash
  source ~/Desktop/ppigfinder_env/bin/activate
  python ~/Desktop/ppigfinderv1_12.py

Then make it executable:
  chmod +x ~/Desktop/run_ppigfinder.command

Double-clicking run_ppigfinder.command in Finder will open the app.

── 4.8  macOS with Anaconda / Miniconda ───────────────────────────
If you prefer conda (similar to the Windows Anaconda route):
  conda create -n ppigfinder python=3.11 -y
  conda activate ppigfinder
  conda install -c conda-forge pyqt matplotlib numpy -y
  conda install -c bioconda blast hmmer pyrodigal -y
  pip install paramiko
  python /path/to/ppigfinderv1_12.py

── 4.9  macOS notes ───────────────────────────────────────────────
• On Apple Silicon (M1/M2/M3), all packages listed above have native
  ARM builds — no Rosetta emulation needed.
• If PyQt6 gives a "This app is not optimised" warning on first launch,
  go to System Preferences → Privacy & Security and allow it.
• Gatekeeper may block the first run of BLAST+ executables. Fix with:
    xattr -dr com.apple.quarantine $(brew --prefix)/bin/blastp

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SECTION 5 — LINUX (Ubuntu / Debian)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
── 5.1  Install system dependencies ───────────────────────────────
  sudo apt-get update
  sudo apt-get install -y python3 python3-pip python3-venv
  sudo apt-get install -y ncbi-blast+ hmmer
  sudo apt-get install -y libgl1 libglib2.0-0    # Qt runtime libs

── 5.2  Create a virtual environment and install packages ──────────
  python3 -m venv ~/ppigfinder_env
  source ~/ppigfinder_env/bin/activate
  pip install --upgrade pip
  pip install PyQt6 matplotlib numpy pyrodigal paramiko

── 5.3  Run ppigFinder ────────────────────────────────────────────
  source ~/ppigfinder_env/bin/activate
  python ~/ppigfinderv1_12.py

── 5.4  Linux with conda ──────────────────────────────────────────
  conda create -n ppigfinder python=3.11 -y
  conda activate ppigfinder
  conda install -c conda-forge pyqt matplotlib numpy -y
  conda install -c bioconda blast hmmer pyrodigal -y
  pip install paramiko
  python /path/to/ppigfinderv1_12.py

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SECTION 6 — VERIFY ALL BACKENDS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
After launching ppigFinder, look at the Backends panel in the
right-side Genome tab. Every installed component shows ✅:

  BLAST+     ✅  blastp: 2.x.x+   (or "via WSL" on Windows)
  HMMER3     ✅  via WSL / 3.x     (or direct path on Mac/Linux)
  Pyrodigal  ✅  3.x.x

If any backend shows ❌, the feature degrades gracefully:
  BLAST+   ❌ → built-in k-mer aligner is used (less sensitive)
  HMMER3   ❌ → built-in PSSM scanner is used (less sensitive)
  Pyrodigal❌ → Pyrodigal and Hybrid modes show an error message;
                use Automatic mode instead

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SECTION 7 — QUICK INSTALL REFERENCE
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Minimum working installation (all platforms, pip):
  pip install PyQt6 matplotlib numpy

Recommended full installation:
  pip install PyQt6 matplotlib numpy pyrodigal paramiko

Full installation via conda:
  conda install -c conda-forge pyqt matplotlib numpy
  conda install -c bioconda blast hmmer pyrodigal
  pip install paramiko

External tools (outside pip/conda):
  Windows  BLAST+  → NCBI .exe installer or WSL (Section 1.5 / 3.2)
  Windows  HMMER3  → WSL only (Section 3.3)
  macOS    BLAST+  → brew install blast
  macOS    HMMER3  → brew install hmmer
  Linux    BLAST+  → apt-get install ncbi-blast+
  Linux    HMMER3  → apt-get install hmmer

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SECTION 8 — TROUBLESHOOTING
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Problem: "No module named PyQt6" when running the app
  Fix: pip install PyQt6
       If that fails, try: pip install PyQt5
       (ppigFinder auto-detects either version)

Problem: The app window opens but plots are blank / grey
  Fix: pip install matplotlib numpy --upgrade

Problem: Pyrodigal install fails on Windows with a compiler error
  Fix: pip install --upgrade pip setuptools wheel
       pip install pyrodigal
       If still failing, install Visual C++ Build Tools from:
       https://visualstudio.microsoft.com/visual-cpp-build-tools/

Problem: BLAST+ not detected even though it is installed
  Check 1: open a terminal and type "blastp -version". If it works
            there but not in ppigFinder, BLAST+ is not on your PATH.
  Windows fix: add the BLAST bin folder to System Environment Variables
               → PATH (e.g. C:\\Program Files\\NCBI\\blast-2.16.0+\\bin)
  macOS/Linux fix: add to ~/.zshrc or ~/.bashrc:
               export PATH="/usr/local/ncbi/blast/bin:$PATH"
               then: source ~/.zshrc

Problem: HMMER3 not detected on Windows
  Fix: install WSL (Section 3) and then:
       wsl bash -c "sudo apt-get install -y hmmer"

Problem: "qt.qpa.plugin: Could not load the Qt platform plugin"
  Linux fix: sudo apt-get install -y libgl1 libxcb-xinerama0
  Conda fix: conda install -c conda-forge libstdcxx-ng

Problem: app opens but all text appears as boxes / garbled
  Fix: install a Unicode font:
       Linux: sudo apt-get install fonts-noto
       macOS: already included (San Francisco / Helvetica Neue)
       Windows: already included (Segoe UI)

Problem: "Permission denied" when running the .command file on macOS
  Fix: chmod +x ~/Desktop/run_ppigfinder.command
""",
    },
}




# ═══════════════════════════════════════════════════════════════
# MODULE 2: GRAPHICAL INTERFACE v9.0
# ═══════════════════════════════════════════════════════════════



# ═══════════════════════════════════════════════════════════════
# WORKER THREAD — thread-safe async operations via Qt signals
# ═══════════════════════════════════════════════════════════════

class AnalysisWorker(QThread):
    """Generic worker thread that emits result or error."""
    finished = pyqtSignal(object)   # result payload
    error    = pyqtSignal(str)      # error message
    progress = pyqtSignal(str)      # status updates

    def __init__(self, fn, *args, **kwargs):
        super().__init__()
        self.fn = fn
        self.args = args
        self.kwargs = kwargs

    def run(self):
        try:
            result = self.fn(*self.args, **self.kwargs)
            self.finished.emit(result)
        except Exception as e:
            self.error.emit(str(e))


# ═══════════════════════════════════════════════════════════════
# GENOME MAP WIDGET — custom QWidget replacing tk.Canvas
# ═══════════════════════════════════════════════════════════════

class GenomeMapWidget(QWidget):
    """Custom painted widget that shows ORFs as directional arrows on a backbone."""
    orf_clicked = pyqtSignal(int)  # emits ORF index
    zoom_changed = pyqtSignal(float)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.dna_length = 0
        self.dna_sequence = ""  # for DNA display at high zoom
        self.orfs = []
        self.hmm_profiles = []
        self.zoom_level = 1.0
        self.pan_offset = 0
        self._dragging = False
        self._drag_start = 0
        self.highlight_idx = -1
        self.setMinimumHeight(220)
        self.setMouseTracking(True)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self._orf_rects = []  # [(QRectF, orf_index), ...]

    def set_data(self, dna_length, orfs, hmm_profiles=None, dna_sequence=""):
        self.dna_length = dna_length
        self.dna_sequence = dna_sequence
        self.orfs = orfs
        self.hmm_profiles = hmm_profiles or []
        self._orf_rects = []
        self.update()

    def set_zoom(self, level, anchor_x: float | None = None):
        """Set zoom level.  If *anchor_x* (pixel) is given, the genomic
        position under that pixel is kept in place after the zoom."""
        bw = max(1, self.width() - 80)   # drawable width (margin=40 each side)
        old_gw = max(1, int(bw * self.zoom_level))
        new_level = max(0.5, min(10000.0, level))

        if anchor_x is not None:
            # Genomic fraction under the cursor
            frac = (self.pan_offset + anchor_x - 40) / old_gw
            new_gw = int(bw * new_level)
            # After zoom, pan so that the same fraction is still under cursor
            self.pan_offset = int(frac * new_gw - (anchor_x - 40))
        self.zoom_level = new_level
        self._clamp_pan()
        self.zoom_changed.emit(self.zoom_level)
        self.update()

    def _clamp_pan(self):
        """Keep pan_offset inside [0, max_overflow].  At zoom ≤ 1 the whole
        genome fits in the widget, so pan is forced to 0."""
        bw = max(1, self.width() - 80)
        gw  = int(bw * self.zoom_level)
        max_pan = max(0, gw - bw)
        self.pan_offset = max(0, min(self.pan_offset, max_pan))

    def wheelEvent(self, event):
        try:
            delta = event.angleDelta().y()
        except AttributeError:
            delta = event.delta()
        ctrl = (Qt.KeyboardModifier.ControlModifier if QT_VERSION == 6
                else Qt.ControlModifier)
        if event.modifiers() & ctrl:
            pos = event.position() if hasattr(event, 'position') else event.pos()
            factor = 1.15 if delta > 0 else 0.87
            self.set_zoom(self.zoom_level * factor, anchor_x=pos.x())
        else:
            self.pan_offset -= delta // 4
            self._clamp_pan()
            self.update()

    def mousePressEvent(self, event):
        if event.button() == LeftButton:
            # Check if clicked on an ORF
            pos = event.position() if hasattr(event, 'position') else event.pos()
            for rect, idx in self._orf_rects:
                if rect.contains(pos):
                    self.orf_clicked.emit(idx)
                    self.highlight_idx = idx
                    self.update()
                    return
            self._dragging = True
            self._drag_start = pos.x()

    def mouseReleaseEvent(self, event):
        self._dragging = False

    def mouseMoveEvent(self, event):
        if self._dragging:
            pos = event.position() if hasattr(event, 'position') else event.pos()
            dx = pos.x() - self._drag_start
            self.pan_offset -= int(dx)
            self._clamp_pan()
            self._drag_start = pos.x()
            self.update()

    def resizeEvent(self, event):
        """Re-clamp pan after resize so we don't get stranded."""
        super().resizeEvent(event)
        self._clamp_pan()
        self.update()

    def paintEvent(self, event):
        if not self.dna_length or not self.orfs:
            return
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)

        w = self.width()
        h = self.height()
        mg = 40
        y_center = h // 2
        bw = w - 2 * mg
        sl = self.dna_length
        gw = int(bw * self.zoom_level)

        # Build HMM color map
        orf_hmm_colors = {}
        for profile in self.hmm_profiles:
            for hit in profile.get('hits', []):
                oi = hit.get('orf_index', -1)
                if oi >= 0:
                    orf_hmm_colors[oi] = profile['color']

        # Backbone
        painter.setPen(QPen(QColor(0, 0, 0), 2))
        x_start = mg - self.pan_offset
        x_end = x_start + gw
        painter.drawLine(int(max(0, x_start)), y_center, int(min(w, x_end)), y_center)

        # Scale markers
        n_marks = max(6, int(self.zoom_level * 6))
        painter.setPen(QPen(QColor(80, 80, 80), 1))
        font_small = QFont('Arial', 7)
        painter.setFont(font_small)
        for i in range(n_marks + 1):
            x = x_start + (i / n_marks) * gw
            if -50 < x < w + 50:
                painter.drawLine(int(x), y_center - 5, int(x), y_center + 5)
                bp_pos = int(i / n_marks * sl)
                painter.drawText(int(x) - 20, y_center + 18, 40, 14,
                                 AlignCenter, f"{bp_pos:,}")

        # DNA sequence display at high zoom
        bp_per_px = sl / gw if gw > 0 else 999
        if bp_per_px < 0.5 and self.dna_sequence:
            # Show individual nucleotides
            font_dna = QFont('Courier New', max(6, min(12, int(1.0 / bp_per_px * 0.8))))
            painter.setFont(font_dna)
            # QFontMetrics and char_w removed — not currently used

            # Calculate visible range
            vis_bp_start = max(0, int((self.pan_offset - mg) / gw * sl))
            vis_bp_end = min(sl, int((self.pan_offset + w - mg) / gw * sl))

            dna_colors = {'A': QColor('#2E7D32'), 'T': QColor('#C62828'),
                          'G': QColor('#1565C0'), 'C': QColor('#FF8F00'),
                          'N': QColor('#757575')}

            for bp in range(vis_bp_start, min(vis_bp_end, len(self.dna_sequence))):
                x = x_start + (bp / sl) * gw
                if 0 <= x <= w:
                    base = self.dna_sequence[bp] if bp < len(self.dna_sequence) else 'N'
                    painter.setPen(QPen(dna_colors.get(base, QColor('#333'))))
                    painter.drawText(int(x), y_center - 8, base.upper())

        # Lane assignment
        LANE_SPACING = 24
        BASE_OFFSET = 22
        MAX_LANES = 4
        GAP_PX = 4

        orf_px = []
        for orf in self.orfs:
            ox1 = x_start + (orf['start'] / sl) * gw
            ox2 = x_start + (orf['end'] / sl) * gw
            orf_px.append((ox1, ox2))

        sorted_all = sorted(range(len(self.orfs)), key=lambda i: self.orfs[i]['start'])
        plus_lane_ends = []
        minus_lane_ends = []
        orf_lanes = {}

        for i in sorted_all:
            orf = self.orfs[i]
            ox1, ox2 = orf_px[i]
            lane_ends = plus_lane_ends if orf['strand'] == '+' else minus_lane_ends
            assigned = False
            for li, last_end in enumerate(lane_ends):
                if ox1 >= last_end + GAP_PX:
                    lane_ends[li] = ox2
                    orf_lanes[i] = li
                    assigned = True
                    break
            if not assigned:
                new_lane = len(lane_ends)
                if new_lane < MAX_LANES:
                    lane_ends.append(ox2)
                    orf_lanes[i] = new_lane
                else:
                    orf_lanes[i] = MAX_LANES - 1

        # Draw ORFs
        self._orf_rects = []
        font_orf = QFont('Arial', 6, QFont.Weight.Bold)
        painter.setFont(font_orf)

        for i in sorted_all:
            x1, x2 = orf_px[i]
            if x2 < -10 or x1 > w + 10:
                continue
            orf = self.orfs[i]

            # Color selection
            if orf.get('custom_color'):
                c = QColor(orf['custom_color'])
            elif i in orf_hmm_colors:
                c = QColor(orf_hmm_colors[i])
            elif orf.get('domains'):
                c = QColor('#607D8B')
                dom_name = orf['domains'][0].get('domain', '')
                for p in self.hmm_profiles:
                    if p['name'] == dom_name:
                        c = QColor(p['color']); break
            else:
                c = QColor('#7986CB') if orf.get('source') == 'pyrodigal' else QColor('#B2BEC3')

            lane = orf_lanes.get(i, 0)
            lane_offset = BASE_OFFSET + lane * LANE_SPACING
            yp = y_center + lane_offset if orf['strand'] == '+' else y_center - lane_offset
            arrow_h = 9
            arrow_tip = max(4, min(12, (x2 - x1) * 0.2))

            if x1 < x2 and (x2 - x1) > 2:
                poly = QPolygonF()
                if orf['strand'] == '+':
                    body_end = x2 - arrow_tip
                    if body_end <= x1: body_end = (x1 + x2) / 2
                    for px, py in [(x1, yp - arrow_h), (body_end, yp - arrow_h),
                                    (x2, yp), (body_end, yp + arrow_h), (x1, yp + arrow_h)]:
                        poly.append(QPointF(px, py))
                else:
                    body_start = x1 + arrow_tip
                    if body_start >= x2: body_start = (x1 + x2) / 2
                    for px, py in [(x2, yp - arrow_h), (body_start, yp - arrow_h),
                                    (x1, yp), (body_start, yp + arrow_h), (x2, yp + arrow_h)]:
                        poly.append(QPointF(px, py))

                is_hmm = i in orf_hmm_colors
                is_highlight = (i == self.highlight_idx)
                outline = QColor('red') if is_highlight else (QColor('black') if is_hmm else QColor('gray'))
                outline_w = 3 if is_highlight else (2 if is_hmm else 1)
                painter.setPen(QPen(outline, outline_w))
                painter.setBrush(QBrush(c))
                painter.drawPolygon(poly)

                # Store rect for click detection
                try:
                    from PyQt6.QtCore import QRectF
                except ImportError:
                    from PyQt5.QtCore import QRectF
                self._orf_rects.append((QRectF(min(x1,x2), yp - arrow_h,
                                                abs(x2-x1), arrow_h*2), i))

                # Label
                mid_x = (x1 + x2) / 2
                if (x2 - x1) > 25:
                    painter.setPen(QPen(QColor('white')))
                    painter.drawText(int(mid_x) - 15, int(yp) - 6, 30, 12,
                                     AlignCenter, f"O{i+1}")

        painter.end()


# ═══════════════════════════════════════════════════════════════
# MAIN MODULE: ppigFinderApp (PyQt6/5)
# ═══════════════════════════════════════════════════════════════

class _PpiArcMapWidget(QWidget):
    """Custom widget that draws ORFs as arrows on a genome backbone and
    overlays cubic Bezier arcs for each predicted protein-protein interaction.

    Arc colour encodes confidence (PAE_min or ipTM).
    Arc height encodes genomic distance (or score, depending on mode).
    Supports scroll-zoom, drag-pan, click on arc/ORF.
    """

    arc_clicked  = pyqtSignal(int)   # arc index
    orf_clicked  = pyqtSignal(int)   # ORF index
    arc_hovered  = pyqtSignal(int)   # arc index (-1 = none)

    # Colour thresholds for PAE_min-based colouring
    _PAE_HIGH = 4.0
    _PAE_MED  = 8.0

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMouseTracking(True)
        self.setMinimumHeight(280)
        # Use module-level QSizePolicy (already imported at top of file)
        try:
            self.setSizePolicy(QSizePolicy.Policy.Expanding,
                               QSizePolicy.Policy.Expanding)
        except AttributeError:
            self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        self._dna_length   = 0
        self._orfs         = []
        self._hmm_profiles = []
        self._arcs         = []   # list of dicts from _ppi_arc_map_refresh

        self._zoom  = 1.0
        self._pan   = 0.0          # fraction offset (0..1-1/zoom)
        self._drag  = False
        self._drag_x = 0

        self._hovered_arc = -1
        self._selected_arc = -1

        # Pre-computed geometry (rebuilt in paintEvent when needed)
        self._orf_rects  = []   # (QRectF, orf_index)
        self._arc_paths  = []   # QPainterPath per arc
        self._dirty = True

    def set_data(self, dna_length, orfs, hmm_profiles, arcs):
        self._dna_length   = dna_length
        self._orfs         = orfs
        self._hmm_profiles = hmm_profiles
        self._arcs         = arcs
        self._dirty        = True
        self.update()

    # ── coordinate helpers ────────────────────────────────────────────────

    def _genomic_to_x(self, pos, w):
        """Convert a genomic position to a pixel X coordinate."""
        if not self._dna_length:
            return 0
        frac = pos / self._dna_length
        visible_start = self._pan
        visible_width = 1.0 / self._zoom
        rel = (frac - visible_start) / visible_width
        margin = 40
        return margin + rel * (w - 2 * margin)

    def _x_to_genomic(self, x, w):
        margin = 40
        if w <= 2 * margin or not self._dna_length:
            return 0
        rel = (x - margin) / (w - 2 * margin)
        visible_start = self._pan
        visible_width = 1.0 / self._zoom
        return int((visible_start + rel * visible_width) * self._dna_length)

    # ── colour helpers ────────────────────────────────────────────────────

    @staticmethod
    def _arc_colour(arc):
        """Return a QColor for the arc based on its score."""
        pae = arc.get('pae_min')
        if pae is None:
            return QColor('#888780')
        if pae < _PpiArcMapWidget._PAE_HIGH:
            return QColor('#1D9E75')   # green — HIGH
        if pae < _PpiArcMapWidget._PAE_MED:
            return QColor('#BA7517')   # amber — MED
        return QColor('#E24B4A')       # red   — LOW

    @staticmethod
    def _orf_colour(orf_idx, orfs, hmm_profiles):
        """Return fill and border colours for an ORF based on HMM hits."""
        if orf_idx < 0 or orf_idx >= len(orfs):
            return QColor('#D3D1C7'), QColor('#888780')
        # Check if this ORF has an HMM hit
        for prof in hmm_profiles:
            for hit in prof.get('hits', []):
                if hit.get('orf_index') == orf_idx:
                    c = prof.get('color', '#1D9E75')
                    return QColor(c).lighter(160), QColor(c)
        return QColor('#D3D1C7'), QColor('#888780')

    # ── paint ─────────────────────────────────────────────────────────────

    def paintEvent(self, event):
        painter = QPainter(self)
        QP = type(painter)   # alias for readability below
        painter.setRenderHint(QP.RenderHint.Antialiasing
                               if hasattr(QP, 'RenderHint') else QP.Antialiasing)

        w = self.width()
        h = self.height()
        backbone_y = h - 80   # Y position of genome backbone

        # Background
        painter.fillRect(0, 0, w, h, self.palette().window())

        if not self._dna_length:
            painter.drawText(w // 2 - 120, h // 2,
                             "Load a genome and AF3 results, then click Refresh")
            return

        # ── Genome backbone ──────────────────────────────────────────────
        pen = QPen(QColor('#B4B2A9'), 1.5)
        painter.setPen(pen)
        painter.drawLine(40, backbone_y, w - 10, backbone_y)

        # Tick marks and position labels
        n_ticks = 6
        for i in range(n_ticks + 1):
            frac = i / n_ticks
            vis_start = self._pan
            vis_width = 1.0 / self._zoom
            genomic_pos = int((vis_start + frac * vis_width) * self._dna_length)
            px = self._genomic_to_x(genomic_pos, w)
            painter.setPen(QPen(QColor('#888780'), 1))
            painter.drawLine(int(px), backbone_y - 5, int(px), backbone_y + 5)
            label = f"{genomic_pos/1e6:.3f}M"
            painter.setFont(QFont('Arial', 7))
            painter.drawText(int(px) - 18, backbone_y + 16, label)

        # ── ORF arrows ───────────────────────────────────────────────────
        self._orf_rects = []
        orf_h = 16
        arrow_tip = 8

        for i, orf in enumerate(self._orfs):
            x1 = self._genomic_to_x(orf['start'], w)
            x2 = self._genomic_to_x(orf['end'],   w)
            if x2 - x1 < 2:
                continue
            if x2 < 40 or x1 > w - 10:
                continue

            fwd = orf.get('strand', '+') == '+'
            y_top = backbone_y - orf_h - 2 if fwd else backbone_y + 3

            fill, border = _PpiArcMapWidget._orf_colour(i, self._orfs, self._hmm_profiles)

            # Draw arrow polygon
            tip_w = min(arrow_tip, x2 - x1)
            if fwd:
                pts = [QPointF(x1, y_top),
                       QPointF(x2 - tip_w, y_top),
                       QPointF(x2, y_top + orf_h / 2),
                       QPointF(x2 - tip_w, y_top + orf_h),
                       QPointF(x1, y_top + orf_h)]
            else:
                pts = [QPointF(x2, y_top),
                       QPointF(x1 + tip_w, y_top),
                       QPointF(x1, y_top + orf_h / 2),
                       QPointF(x1 + tip_w, y_top + orf_h),
                       QPointF(x2, y_top + orf_h)]

            painter.setPen(QPen(border, 1.0))
            painter.setBrush(QBrush(fill))
            painter.drawPolygon(QPolygonF(pts))

            # Label (only if wide enough)
            if x2 - x1 > 30:
                painter.setPen(QPen(border.darker(200), 1))
                painter.setFont(QFont('Arial', 7))
                cy = y_top + orf_h / 2 + 3
                painter.drawText(int(x1 + 2), int(cy),
                                 int(x2 - x1 - 4), 12,
                                 0x0004 | 0x0080,   # AlignHCenter | AlignVCenter
                                 f"ORF{i+1}")

            self._orf_rects.append((QRectF(x1, y_top, x2-x1, orf_h), i))

        # ── Interaction arcs ──────────────────────────────────────────────
        self._arc_paths = []
        max_arc_h = backbone_y - 20   # maximum arc height from backbone

        # Sort so HIGH arcs are drawn on top
        sorted_arcs = sorted(enumerate(self._arcs),
                             key=lambda x: (x[1].get('pae_min') or 99), reverse=True)

        for arc_idx, arc in sorted_arcs:
            cx_a = (self._genomic_to_x(arc['start_a'], w) +
                    self._genomic_to_x(arc['end_a'],   w)) / 2
            cx_b = (self._genomic_to_x(arc['start_b'], w) +
                    self._genomic_to_x(arc['end_b'],   w)) / 2

            if max(cx_a, cx_b) < 30 or min(cx_a, cx_b) > w - 5:
                self._arc_paths.append(None)
                continue

            # Arc height
            hmode = arc.get('height_mode', '')
            if 'distance' in hmode:
                dist = abs(arc['start_b'] - arc['start_a'])
                frac = min(dist / max(self._dna_length, 1), 1.0)
                arc_h = int(30 + frac * (max_arc_h - 30))
            elif 'score' in hmode:
                pae = arc.get('pae_min') or 15
                arc_h = int(30 + (pae / 20) * (max_arc_h - 30))
            else:
                arc_h = max_arc_h // 2

            arc_h = max(25, min(arc_h, max_arc_h))

            path = QPainterPath()
            path.moveTo(cx_a, backbone_y - 2)
            ctrl_x = (cx_a + cx_b) / 2
            ctrl_y = backbone_y - arc_h
            path.cubicTo(cx_a, ctrl_y, cx_b, ctrl_y, cx_b, backbone_y - 2)
            self._arc_paths.append(path)

            colour = _PpiArcMapWidget._arc_colour(arc)
            selected = (arc_idx == self._selected_arc)
            hovered  = (arc_idx == self._hovered_arc)

            pae = arc.get('pae_min') or 99
            is_high = pae < self._PAE_HIGH
            is_med  = pae < self._PAE_MED

            lw = 3.0 if selected else (2.5 if hovered else
                 2.2 if is_high else 1.5)

            if is_high:
                pen = QPen(colour, lw)
            elif is_med:
                pen = QPen(colour, lw)
                try:
                    pen.setStyle(Qt.PenStyle.DashLine)
                except AttributeError:
                    pen.setStyle(Qt.DashLine)
            else:
                pen = QPen(colour, lw)
                try:
                    pen.setStyle(Qt.PenStyle.DotLine)
                except AttributeError:
                    pen.setStyle(Qt.DotLine)

            if selected:
                pen.setColor(QColor('#7F77DD'))
                pen.setWidth(3)

            painter.setPen(pen)
            try:
                painter.setBrush(Qt.BrushStyle.NoBrush)
            except AttributeError:
                painter.setBrush(Qt.NoBrush)
            painter.drawPath(path)

            # Score label near apex
            lbl_x = int(ctrl_x) - 12
            lbl_y = int(ctrl_y) - 2
            pae_v = arc.get('pae_min')
            if pae_v is not None:
                painter.setPen(QPen(colour.darker(130), 1))
                painter.setFont(QFont('Arial', 7))
                painter.drawText(lbl_x, lbl_y, f"{pae_v:.1f}Å")

        painter.end()

    # ── mouse events ──────────────────────────────────────────────────────

    def _arc_hit(self, pos):
        """Return index of arc whose path is within 5px of pos, or -1."""
        pt = QPointF(pos.x(), pos.y()) if hasattr(pos, 'x') else QPointF(pos)
        for i, path in enumerate(self._arc_paths):
            if path is None:
                continue
            stroker = QPainterPathStroker()
            stroker.setWidth(10)
            wide = stroker.createStroke(path)
            if wide.contains(pt):
                return i
        return -1

    def _orf_hit(self, pos):
        """Return ORF index under pos, or -1."""
        pt = QPointF(pos.x(), pos.y())
        for rect, orf_idx in self._orf_rects:
            if rect.contains(pt):
                return orf_idx
        return -1

    def mouseMoveEvent(self, event):
        pos = event.pos()
        if self._drag:
            dx = pos.x() - self._drag_x
            self._drag_x = pos.x()
            vis_width = 1.0 / self._zoom
            delta = -dx / max(self.width() - 80, 1) * vis_width
            self._pan = max(0.0, min(self._pan + delta,
                                     1.0 - vis_width))
            self.update()
            return

        arc_idx = self._arc_hit(pos)
        if arc_idx != self._hovered_arc:
            self._hovered_arc = arc_idx
            self.arc_hovered.emit(arc_idx)
            self.update()

    def mousePressEvent(self, event):
        try:
            LMB = Qt.MouseButton.LeftButton
        except AttributeError:
            LMB = Qt.LeftButton
        if event.button() == LMB:
            arc_idx = self._arc_hit(event.pos())
            if arc_idx >= 0:
                self._selected_arc = arc_idx
                self.arc_clicked.emit(arc_idx)
                self.update()
                return
            orf_idx = self._orf_hit(event.pos())
            if orf_idx >= 0:
                self.orf_clicked.emit(orf_idx)
                return
            # Start pan
            self._drag   = True
            self._drag_x = event.pos().x()
            try:
                self.setCursor(Qt.CursorShape.ClosedHandCursor)
            except AttributeError:
                self.setCursor(Qt.ClosedHandCursor)

    def mouseReleaseEvent(self, event):
        self._drag = False
        try:
            self.setCursor(Qt.CursorShape.ArrowCursor)
        except AttributeError:
            self.setCursor(Qt.ArrowCursor)

    def wheelEvent(self, event):
        try:
            delta = event.angleDelta().y()
        except Exception:
            delta = event.delta()
        factor = 1.15 if delta > 0 else 1 / 1.15
        new_zoom = max(1.0, min(self._zoom * factor, 200.0))
        # Keep genome position under cursor fixed
        anchor_x = event.pos().x()
        w = self.width()
        margin = 40
        if w > 2 * margin:
            frac_under = self._pan + ((anchor_x - margin) / (w - 2 * margin)) / self._zoom
            vis_width_new = 1.0 / new_zoom
            self._pan = max(0.0, min(frac_under - (anchor_x - margin) / (w - 2 * margin) * vis_width_new,
                                     1.0 - vis_width_new))
        self._zoom = new_zoom
        self.update()



class _OrfNumericItem(QTableWidgetItem):
    """QTableWidgetItem whose sort key is numeric, not lexicographic.

    For ORF IDs ('ORF1' … 'ORF4303') it extracts the trailing integer.
    For positional / size / score columns it parses the float value.
    Falls back to string comparison for any non-numeric content.
    """
    _cache: dict = {}   # text → numeric key (avoids repeated re.search)

    def __init__(self, text: str):
        super().__init__(text)
        self._num = _OrfNumericItem._numeric_key(text)

    @staticmethod
    def _numeric_key(text: str) -> float:
        import re as _re
        # Fast path for 'ORF<N>' pattern
        _m = _re.match(r'ORF(\d+)$', text, _re.IGNORECASE)
        if _m:
            return float(_m.group(1))
        # Strip units (Å, %) and commas then try float
        _clean = text.replace(',', '').replace(' Å', '').replace('Å', '').replace('%', '').strip()
        try:
            return float(_clean)
        except ValueError:
            pass
        # Trailing number in any other string
        _m2 = _re.search(r'(\d+(?:\.\d+)?)$', text)
        if _m2:
            return float(_m2.group(1))
        return float('inf')   # non-numeric values sort last

    def __lt__(self, other: 'QTableWidgetItem') -> bool:
        if isinstance(other, _OrfNumericItem):
            if self._num != other._num:
                return self._num < other._num
        return super().__lt__(other)


# ═══════════════════════════════════════════════════════════════
# RESPONSIVE WINDOW MANAGEMENT
# ═══════════════════════════════════════════════════════════════
# Provides relative-sized windows, persistent geometry between
# sessions (QSettings), quick-size presets, custom-size dialog,
# and child-window registration for floating panels/dialogs.
# ═══════════════════════════════════════════════════════════════

# Pre-defined relative-size presets (fraction of available screen)
WIN_PRESETS_RELATIVE = {
    "Small (60%)":   (0.60, 0.65),
    "Medium (75%)":  (0.75, 0.80),
    "Large (85%)":   (0.85, 0.90),    # default
    "X-Large (95%)": (0.95, 0.95),
    "Maximize":      None,
}


class WindowManager:
    """
    Per-app window manager with responsive sizing.

    Handles:
      • Relative sizing (% of available screen)
      • Persistent geometry across sessions (QSettings)
      • Sub-window registration with relative sizing
      • Custom-size dialog (Ctrl+Shift+W)
    """

    def __init__(self, main_window: QMainWindow,
                 app_name: str = "ppigFinder",
                 organisation: str = "ppigFinder"):
        self.main = main_window
        self.app_name = app_name
        self.org = organisation
        self.settings = QSettings(organisation, app_name)
        self._registered_subwindows = {}

    def _available_screen_rect(self):
        """Return available rect of the screen the main window lives on."""
        try:
            screen = self.main.screen() if hasattr(self.main, "screen") else None
        except Exception:
            screen = None
        if screen is None:
            screen = QApplication.primaryScreen()
        if screen is None:
            from PyQt6.QtCore import QRect  # type: ignore
            return QRect(0, 0, 1280, 720)
        return screen.availableGeometry()

    def apply_default_size(self, width_pct=0.85, height_pct=0.90,
                           min_w=1100, min_h=700, center=True):
        """Apply a relative size (default ~85% × 90% of screen)."""
        scr = self._available_screen_rect()
        w = max(int(scr.width() * width_pct), min_w)
        h = max(int(scr.height() * height_pct), min_h)
        w = min(w, scr.width())
        h = min(h, scr.height())
        self.main.setMinimumSize(min_w, min_h)
        self.main.resize(w, h)
        if center:
            self.center_on_screen()

    def apply_relative(self, width_pct, height_pct, center=True):
        """Resize the main window to a fraction of the screen."""
        scr = self._available_screen_rect()
        w = max(int(scr.width() * width_pct),
                self.main.minimumWidth() or 600)
        h = max(int(scr.height() * height_pct),
                self.main.minimumHeight() or 400)
        w = min(w, scr.width())
        h = min(h, scr.height())
        if self.main.isMaximized() or self.main.isFullScreen():
            self.main.showNormal()
        self.main.resize(w, h)
        if center:
            self.center_on_screen()

    def apply_fixed(self, w, h, center=True):
        """Resize the main window to absolute pixel dimensions."""
        if self.main.isMaximized() or self.main.isFullScreen():
            self.main.showNormal()
        self.main.resize(int(w), int(h))
        if center:
            self.center_on_screen()

    def apply_maximized(self):
        self.main.showMaximized()

    def apply_fullscreen(self):
        self.main.showFullScreen()

    def center_on_screen(self):
        """Center the main window on its current screen."""
        scr = self._available_screen_rect()
        geo = self.main.frameGeometry()
        geo.moveCenter(scr.center())
        self.main.move(geo.topLeft())

    def save_geometry(self):
        """Write window state + geometry to QSettings."""
        try:
            self.settings.setValue("main/geometry", self.main.saveGeometry())
            self.settings.setValue("main/state",    self.main.saveState())
            self.settings.setValue("main/maximized", self.main.isMaximized())
        except Exception:
            pass

    def restore_geometry(self):
        """Restore window state + geometry from QSettings."""
        try:
            geo = self.settings.value("main/geometry")
            sta = self.settings.value("main/state")
            if geo:
                self.main.restoreGeometry(geo)
            if sta:
                self.main.restoreState(sta)
            return bool(geo)
        except Exception:
            return False

    def register_subwindow(self, key, widget, rel_w=0.6, rel_h=0.7,
                           min_w=600, min_h=400):
        """Configure a child dialog/window with relative sizing."""
        self._registered_subwindows[key] = widget
        scr = self._available_screen_rect()
        w = max(int(scr.width() * rel_w), min_w)
        h = max(int(scr.height() * rel_h), min_h)
        widget.setMinimumSize(min_w, min_h)
        widget.resize(w, h)
        if self.main is not None:
            try:
                main_geo = self.main.frameGeometry()
                geo = widget.frameGeometry()
                geo.moveCenter(main_geo.center())
                widget.move(geo.topLeft())
            except Exception:
                pass

    def open_custom_size_dialog(self):
        """Show a CustomSizeDialog and apply user's choice."""
        dlg = CustomSizeDialog(self.main, self)
        accepted = (QDialog.DialogCode.Accepted
                    if QT_VERSION == 6 else QDialog.Accepted)
        if dlg.exec() == accepted:
            mode, w, h = dlg.get_result()
            if mode == "relative":
                self.apply_relative(w, h)
            elif mode == "fixed":
                self.apply_fixed(int(w), int(h))
            elif mode == "maximized":
                self.apply_maximized()
            elif mode == "fullscreen":
                self.apply_fullscreen()


class CustomSizeDialog(QDialog):
    """Lets the user dial-in any window size (% of screen or pixels)."""

    def __init__(self, parent, mgr):
        super().__init__(parent)
        self.mgr = mgr
        self.setWindowTitle("🪟 Custom Window Size")
        self.setModal(True)
        self.setMinimumWidth(420)

        scr = mgr._available_screen_rect()
        self.scr_w = scr.width()
        self.scr_h = scr.height()

        lay = QVBoxLayout(self)

        lay.addWidget(QLabel(
            f"<b>Screen available:</b> {self.scr_w} × {self.scr_h} px"))

        # Preset row
        pre_g = QGroupBox("Relative preset")
        pre_l = QVBoxLayout(pre_g)
        self.preset_combo = QComboBox()
        self.preset_combo.addItem("— Custom —")
        for name in WIN_PRESETS_RELATIVE:
            self.preset_combo.addItem(name)
        self.preset_combo.currentIndexChanged.connect(self._on_preset_changed)
        pre_l.addWidget(self.preset_combo)
        lay.addWidget(pre_g)

        # Relative percentage spinboxes
        rel_g = QGroupBox("Relative (% of screen)")
        rel_l = QVBoxLayout(rel_g)
        rh = QHBoxLayout()
        rh.addWidget(QLabel("Width:"))
        self.w_pct = QSpinBox(); self.w_pct.setRange(20, 100); self.w_pct.setValue(85)
        self.w_pct.setSuffix(" %")
        self.w_pct.valueChanged.connect(self._update_pixel_preview)
        rh.addWidget(self.w_pct)
        rh.addWidget(QLabel("Height:"))
        self.h_pct = QSpinBox(); self.h_pct.setRange(20, 100); self.h_pct.setValue(90)
        self.h_pct.setSuffix(" %")
        self.h_pct.valueChanged.connect(self._update_pixel_preview)
        rh.addWidget(self.h_pct)
        rel_l.addLayout(rh)
        self.preview_lbl = QLabel("")
        self.preview_lbl.setStyleSheet("color:#666; font-size:11px;")
        rel_l.addWidget(self.preview_lbl)
        lay.addWidget(rel_g)

        # Fixed pixel option
        fix_g = QGroupBox("Or fixed pixels")
        fix_g.setCheckable(True)
        fix_g.setChecked(False)
        fl = QHBoxLayout(fix_g)
        self.fix_w = QSpinBox(); self.fix_w.setRange(600, 7680); self.fix_w.setValue(1550)
        self.fix_h = QSpinBox(); self.fix_h.setRange(400, 4320); self.fix_h.setValue(980)
        fl.addWidget(QLabel("W:")); fl.addWidget(self.fix_w)
        fl.addWidget(QLabel("H:")); fl.addWidget(self.fix_h)
        self.fix_g = fix_g
        lay.addWidget(fix_g)

        # Special modes
        spc_g = QGroupBox("Special modes")
        spc_l = QHBoxLayout(spc_g)
        self.btn_max = QPushButton("Maximize")
        self.btn_max.clicked.connect(lambda: self._apply_special("maximized"))
        self.btn_full = QPushButton("Fullscreen")
        self.btn_full.clicked.connect(lambda: self._apply_special("fullscreen"))
        spc_l.addWidget(self.btn_max); spc_l.addWidget(self.btn_full)
        lay.addWidget(spc_g)

        # Buttons
        if QT_VERSION == 6:
            bb = QDialogButtonBox(
                QDialogButtonBox.StandardButton.Ok |
                QDialogButtonBox.StandardButton.Cancel)
        else:
            bb = QDialogButtonBox(
                QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        lay.addWidget(bb)

        self._special_mode = None
        self._update_pixel_preview()

    def _on_preset_changed(self, idx):
        if idx <= 0:
            return
        name = self.preset_combo.currentText()
        val = WIN_PRESETS_RELATIVE.get(name)
        if val is None:
            self._apply_special("maximized")
            self.accept()
            return
        wp, hp = val
        self.w_pct.setValue(int(wp * 100))
        self.h_pct.setValue(int(hp * 100))

    def _update_pixel_preview(self):
        wp = self.w_pct.value() / 100.0
        hp = self.h_pct.value() / 100.0
        self.preview_lbl.setText(
            f"≈ {int(self.scr_w * wp)} × {int(self.scr_h * hp)} pixels")

    def _apply_special(self, mode):
        self._special_mode = mode
        self.accept()

    def get_result(self):
        """Return (mode, width, height)."""
        if self._special_mode:
            return self._special_mode, 0, 0
        if self.fix_g.isChecked():
            return "fixed", self.fix_w.value(), self.fix_h.value()
        return ("relative",
                self.w_pct.value() / 100.0,
                self.h_pct.value() / 100.0)


# ═══════════════════════════════════════════════════════════════
# DETACHABLE TAB WIDGET
# ═══════════════════════════════════════════════════════════════
# QTabWidget extension that lets the user "tear off" any tab into
# a free-floating window with its own resizable geometry. Closing
# the floating window re-attaches the tab back into its original
# position. Works with any kind of QWidget content.
#
# Usage: replace any `QTabWidget()` with `DetachableTabWidget()`.
# Per-tab UX:
#   • A "↗" button is added to the corner of every tab
#   • Right-click on a tab also offers "Detach into window"
#   • Once detached, closing the window re-docks the tab
# ═══════════════════════════════════════════════════════════════


class _DetachedTabWindow(QDialog):
    """Floating window that hosts a detached tab.

    On close, signals the parent DetachableTabWidget to re-dock
    the tab at its original index.
    """
    closedSignal = pyqtSignal(object)  # emits self

    def __init__(self, content_widget, tab_label, tab_icon=None,
                 original_index=0, parent=None):
        # Use the top-level main window as parent (not the tab widget)
        # so that window management treats this as a peer dialog and
        # doesn't crop / clip the floating content.
        try:
            top_level = parent.window() if parent is not None else None
        except Exception:
            top_level = parent
        super().__init__(top_level)
        if QT_VERSION == 6:
            self.setWindowFlags(Qt.WindowType.Window)
        else:
            self.setWindowFlags(Qt.Window)
        self.setWindowTitle(f"ppigFinder — {tab_label}")
        if tab_icon is not None:
            try:
                self.setWindowIcon(tab_icon)
            except Exception:
                pass

        self.content = content_widget
        self.tab_label = tab_label
        self.tab_icon = tab_icon
        self.original_index = original_index
        self._reattach_on_close = True   # set to False to dispose

        # Default size: 70% of screen
        try:
            scr = QApplication.primaryScreen().availableGeometry()
            self.resize(int(scr.width() * 0.70), int(scr.height() * 0.80))
            geo = self.frameGeometry()
            geo.moveCenter(scr.center())
            self.move(geo.topLeft())
        except Exception:
            self.resize(1100, 800)

        # Layout: header bar with "Re-dock" button + content
        main_lay = QVBoxLayout(self)
        main_lay.setContentsMargins(4, 4, 4, 4)
        main_lay.setSpacing(4)

        bar = QHBoxLayout()
        info = QLabel(
            f"<b>📌 Detached:</b> <i>{tab_label}</i> "
            "&nbsp;<span style='color:#666; font-size:10px;'>"
            "(close this window to re-dock)</span>")
        info.setStyleSheet("padding:4px;")
        bar.addWidget(info, 1)
        self.btn_redock = QPushButton("↙ Re-dock")
        self.btn_redock.setFixedHeight(26)
        self.btn_redock.setToolTip(
            "Close this window and re-attach this tab to the main app.")
        self.btn_redock.clicked.connect(self.close)
        bar.addWidget(self.btn_redock)
        main_lay.addLayout(bar)

        # Re-parent the content widget INTO this dialog.
        # IMPORTANT: removeTab() on QTabWidget hides the page widget
        # internally — we must explicitly setVisible(True) and call
        # show() AFTER re-parenting, otherwise the floating window
        # appears empty (just shows the header bar).
        content_widget.setParent(self)
        main_lay.addWidget(content_widget, 1)
        content_widget.setVisible(True)
        content_widget.show()
        # Recursively show all child widgets in case some were hidden
        # by their parent QTabWidget.  This guards against deeply
        # nested QStackedWidget pages remaining hidden.
        try:
            for child in content_widget.findChildren(QWidget):
                if not child.isVisible():
                    # Don't force-show widgets that were intentionally
                    # hidden by the user (e.g. collapsed group boxes).
                    # Only force-show direct children of common containers
                    # which Qt may have hidden during removeTab.
                    pass  # rely on parent show() to propagate
        except Exception:
            pass

    def closeEvent(self, event):
        # Notify parent to re-dock unless explicitly disabled
        if self._reattach_on_close:
            try:
                self.closedSignal.emit(self)
            except Exception:
                pass
        super().closeEvent(event)


class DetachableTabWidget(QTabWidget):
    """QTabWidget where every tab can be torn off into a floating window.

    The floating window has its own resizable geometry; closing it
    re-attaches the tab at its original position.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self._detached_windows = {}  # original_index → _DetachedTabWindow
        # Accept right-click to show context menu
        self.tabBar().setContextMenuPolicy(
            Qt.ContextMenuPolicy.CustomContextMenu if QT_VERSION == 6
            else Qt.CustomContextMenu)
        self.tabBar().customContextMenuRequested.connect(
            self._on_tab_context_menu)
        # Accept double-click on tab to detach
        try:
            self.tabBar().tabBarDoubleClicked.connect(self._on_tab_double_clicked)
        except Exception:
            pass

    def _on_tab_double_clicked(self, idx):
        """Double-click on tab → detach."""
        if idx >= 0:
            self.detach_tab(idx)

    def _on_tab_context_menu(self, pos):
        """Right-click on tab → show detach menu."""
        idx = self.tabBar().tabAt(pos)
        if idx < 0:
            return
        menu = QMenu(self)
        act_detach = menu.addAction("↗  Detach into window")
        act_detach.setToolTip("Open this tab in a separate, resizable window")
        menu.addSeparator()
        act_detach_all = menu.addAction("↗  Detach all tabs")
        act_redock_all = menu.addAction("↙  Re-dock all detached tabs")
        if QT_VERSION == 6:
            chosen = menu.exec(self.tabBar().mapToGlobal(pos))
        else:
            chosen = menu.exec_(self.tabBar().mapToGlobal(pos))
        if chosen == act_detach:
            self.detach_tab(idx)
        elif chosen == act_detach_all:
            # Iterate from highest index → lowest because detach removes
            for i in range(self.count() - 1, -1, -1):
                self.detach_tab(i)
        elif chosen == act_redock_all:
            for win in list(self._detached_windows.values()):
                try:
                    win.close()
                except Exception:
                    pass

    def detach_tab(self, idx):
        """Tear off tab `idx` into a floating window."""
        if idx < 0 or idx >= self.count():
            return
        widget = self.widget(idx)
        if widget is None:
            return
        label = self.tabText(idx)
        icon = self.tabIcon(idx)

        # Compute a stable identifier so we can re-dock at the right slot
        # even when other tabs have been detached/redocked in the meantime.
        # We track by widget identity and the original label.
        original_index = idx

        # Remove the tab.  This does NOT delete the widget, but Qt
        # internally hides it (visible=False) so that it doesn't pop
        # up on screen between removeTab and re-parenting.  The
        # _DetachedTabWindow constructor explicitly re-shows it.
        self.removeTab(idx)
        # Clear hidden state proactively
        try:
            widget.setVisible(True)
        except Exception:
            pass

        # Wrap in a floating window.  We pass `self` so the dialog can
        # find the top-level main window for proper window management.
        win = _DetachedTabWindow(widget, label, icon,
                                 original_index=original_index, parent=self)
        win.closedSignal.connect(self._on_detached_closed)
        self._detached_windows[id(widget)] = win
        win.show()
        try:
            win.raise_()
            win.activateWindow()
        except Exception:
            pass

    def _on_detached_closed(self, win):
        """Called when a detached window is closed → re-dock."""
        try:
            wid = id(win.content)
        except Exception:
            wid = None
        # Restore to the closest original index possible
        target_idx = win.original_index
        if target_idx > self.count():
            target_idx = self.count()
        # Re-parent content back to this tab widget.
        # insertTab() handles visibility, but we explicitly setVisible
        # to make sure the widget shows when the tab becomes current.
        try:
            win.content.setParent(self)
            new_idx = self.insertTab(target_idx, win.content,
                                      win.tab_icon, win.tab_label)
            win.content.setVisible(True)
            self.setCurrentIndex(new_idx)
        except Exception as e:
            print(f"[DetachableTabWidget] re-dock error: {e}")
        # Drop the reference
        if wid is not None and wid in self._detached_windows:
            del self._detached_windows[wid]


# ═══════════════════════════════════════════════════════════════
# FLEXIBLE AF3 SERVER SUBMISSION — Customisable command builder
# ═══════════════════════════════════════════════════════════════
# Lets the user pick from named profiles or write custom command
# templates with placeholders (e.g. {json_path}, {job_name},
# {parent_dir}, {output_dir}, plus any user-defined parameter).
# Profiles are persisted to ~/.ppigfinder/af3_server_profiles.json
# ═══════════════════════════════════════════════════════════════

# Default profiles (shipped with the app, cannot be deleted)
AF3_DEFAULT_PROFILES = [
    # Single minimal "blank" starter.  The user is expected to write
    # their own command (or save several named ones) for whatever
    # server / scheduler / runner they have access to.  All previous
    # server-specific built-ins (SLURM, sbatch, PBS, LSF, ...) were
    # removed — they leak assumptions about the user's infrastructure.
    {
        "id":          "blank_starter",
        "name":        "Blank  (write your own command)",
        "description": "Empty starter.  Type your AF3 / scheduler "
                       "command directly in the terminal below.  Click "
                       "'+ New' to save it as a named profile, or "
                       "'Duplicate' to fork from this blank.",
        "is_builtin":  True,
        "template": (
            "# Type the command(s) to run AF3 on YOUR server.\n"
            "# Always-available placeholders:\n"
            "#   {json_path}  {job_name}  {parent_dir}  {output_dir}\n"
            "#   {prefix}     {date}      {timestamp}\n"
            "#\n"
            "# Example (edit freely):\n"
            "cd {parent_dir} && my_af3_cmd --json_path {json_path} --job-name {job_name}"
        ),
        "params": {},
    },
]


class AF3CommandBuilder:
    """Pure-logic AF3 command template resolver."""

    REQUIRED_CONTEXT = {"json_path", "job_name", "parent_dir"}

    @staticmethod
    def detect_placeholders(template):
        """Return list of unique {placeholder} names in the template."""
        pattern = re.compile(r'(?<!\{)\{([a-zA-Z_][a-zA-Z0-9_]*)\}(?!\})')
        seen = []
        for m in pattern.finditer(template):
            n = m.group(1)
            if n not in seen:
                seen.append(n)
        return seen

    @staticmethod
    def auto_context(prefix="af3_batch", parent_dir="~/af3_predictions",
                     job_name=None, json_path=None):
        """Build a minimal context with always-available keys."""
        now = datetime.now()
        ts = now.strftime('%Y%m%d_%H%M%S')
        date = now.strftime('%Y-%m-%d')
        job_name = job_name or f"{prefix}_{ts}"
        json_path = json_path or f"{prefix}_all_jobs.json"
        job_name = re.sub(r'[()\[\]{}|;&!\s]+', '_', job_name).strip('_')
        return {
            "prefix":     prefix,
            "parent_dir": parent_dir.rstrip('/'),
            "json_path":  json_path,
            "job_name":   job_name,
            "output_dir": f"{parent_dir.rstrip('/')}/{job_name}/output",
            "date":       date,
            "timestamp":  ts,
        }

    @staticmethod
    def resolve(template, context, computed=None):
        """Resolve placeholders.  Missing ones are kept as-is."""
        ctx = dict(context)
        if computed:
            for key, src in computed.items():
                try:
                    fn = eval(src, {"__builtins__": {}})
                    ctx[key] = fn(ctx)
                except Exception as e:
                    ctx[key] = f"<{key}: error {e}>"
        def _sub(m):
            n = m.group(1)
            return str(ctx.get(n, m.group(0)))
        return re.sub(r'(?<!\{)\{([a-zA-Z_][a-zA-Z0-9_]*)\}(?!\})',
                      _sub, template)


class AF3ProfileManager:
    """Loads / saves AF3 server profiles to ~/.ppigfinder/af3_server_profiles.json"""

    def __init__(self, config_path=None):
        if config_path is None:
            cfg_dir = Path.home() / ".ppigfinder"
            try:
                cfg_dir.mkdir(exist_ok=True)
            except Exception:
                pass
            config_path = cfg_dir / "af3_server_profiles.json"
        self.config_path = config_path
        self.profiles = []
        self.load()

    def load(self):
        import copy as _copy
        self.profiles = [_copy.deepcopy(p) for p in AF3_DEFAULT_PROFILES]
        if self.config_path.exists():
            try:
                with open(self.config_path, 'r', encoding='utf-8') as f:
                    user_data = json.load(f)
                if isinstance(user_data, list):
                    builtin_ids = {p['id'] for p in self.profiles}
                    for up in user_data:
                        if (isinstance(up, dict) and up.get('id')
                                and up['id'] not in builtin_ids):
                            up['is_builtin'] = False
                            self.profiles.append(up)
            except Exception as e:
                print(f"[AF3ProfileManager] load error: {e}")

    def save(self):
        user_profiles = [p for p in self.profiles if not p.get('is_builtin')]
        try:
            with open(self.config_path, 'w', encoding='utf-8') as f:
                json.dump(user_profiles, f, indent=2)
        except Exception as e:
            print(f"[AF3ProfileManager] save error: {e}")

    def get(self, profile_id):
        for p in self.profiles:
            if p['id'] == profile_id:
                return p
        return None

    def names(self):
        out = []
        for p in self.profiles:
            tag = "  [built-in]" if p.get('is_builtin') else "  [user]"
            out.append((p['id'], p['name'] + tag))
        return out

    def add_or_update(self, profile):
        if profile.get('is_builtin'):
            return False
        for i, p in enumerate(self.profiles):
            if p['id'] == profile['id']:
                if p.get('is_builtin'):
                    return False
                self.profiles[i] = profile
                self.save()
                return True
        self.profiles.append(profile)
        self.save()
        return True

    def delete(self, profile_id):
        for i, p in enumerate(self.profiles):
            if p['id'] == profile_id and not p.get('is_builtin'):
                self.profiles.pop(i)
                self.save()
                return True
        return False

    def duplicate(self, profile_id, new_name):
        import copy as _copy
        src = self.get(profile_id)
        if src is None:
            return None
        new = _copy.deepcopy(src)
        new['id'] = (re.sub(r'[^a-z0-9_]', '_', new_name.lower())[:40]
                     + f"_{int(datetime.now().timestamp())}")
        new['name'] = new_name
        new['is_builtin'] = False
        self.profiles.append(new)
        self.save()
        return new


class FlexibleAF3SubmitWidget(QGroupBox):
    """
    Replaces the legacy "AF3 Advanced Options" group with a fully
    customisable, profile-based command builder.

    Emits commandChanged() whenever the resolved preview changes.
    """
    commandChanged = pyqtSignal()

    def __init__(self, parent=None, profile_manager=None):
        super().__init__("⌨  AF3 Server Submission — Terminal & Profiles",
                         parent)
        self.setCheckable(True)
        self.setChecked(True)
        self.setToolTip(
            "Type the command to run AlphaFold3 (or any scheduler) on YOUR "
            "server in the terminal below.\n"
            "Save it as a named profile to reuse later — profiles are "
            "stored in ~/.ppigfinder/af3_server_profiles.json\n\n"
            "Always-available placeholders:\n"
            "  {json_path}, {job_name}, {parent_dir}, {output_dir},\n"
            "  {prefix}, {date}, {timestamp}")
        self.mgr = profile_manager or AF3ProfileManager()
        self._param_widgets = {}
        self._current_profile_id = None
        self._building = False
        self._popout_window = None     # set by _on_popout()
        self._build_ui()
        if self.mgr.profiles:
            self._load_profile(self.mgr.profiles[0]['id'])

    def _build_ui(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(8, 14, 8, 8)
        outer.setSpacing(6)
        # Ensure the whole widget has enough vertical space so that
        # the inner Parameters / Template / Preview panels are not
        # squeezed by the parent layout.
        self.setMinimumHeight(420)

        # Profile row
        prof_row = QHBoxLayout()
        prof_row.addWidget(QLabel("<b>Profile:</b>"))
        self.profile_combo = QComboBox()
        self.profile_combo.setMinimumWidth(280)
        self.profile_combo.setToolTip(
            "Pre-built profiles for common AF3 / scheduler setups.\n"
            "Select one and edit the parameters below — or create your own.")
        self.profile_combo.currentIndexChanged.connect(self._on_profile_changed)
        prof_row.addWidget(self.profile_combo, 1)

        self.btn_new = QPushButton("➕ New")
        self.btn_new.setToolTip("Create a new profile from scratch.")
        self.btn_new.clicked.connect(self._on_new_profile)
        self.btn_dup = QPushButton("📑 Duplicate")
        self.btn_dup.setToolTip("Copy current profile under a new name.")
        self.btn_dup.clicked.connect(self._on_duplicate_profile)
        self.btn_save = QPushButton("💾 Save")
        self.btn_save.setToolTip("Save edits to current profile (user only).")
        self.btn_save.clicked.connect(self._on_save_profile)
        self.btn_del = QPushButton("🗑 Delete")
        self.btn_del.setToolTip("Delete current profile (user only).")
        self.btn_del.clicked.connect(self._on_delete_profile)
        self.btn_imp = QPushButton("⤴ Import")
        self.btn_imp.setToolTip("Import profiles from a JSON file.")
        self.btn_imp.clicked.connect(self._on_import_profile)
        self.btn_exp = QPushButton("⤵ Export")
        self.btn_exp.setToolTip("Export current profile to a JSON file.")
        self.btn_exp.clicked.connect(self._on_export_profile)
        # Pop-out button — opens this widget in a separate window with
        # plenty of room (avoids the cramped Parameters panel).
        self.btn_popout = QPushButton("↗ Pop-out")
        self.btn_popout.setToolTip(
            "Open this builder in a separate, resizable window "
            "with more room for the Parameters panel and template editor.")
        self.btn_popout.clicked.connect(self._on_popout)
        for b in (self.btn_new, self.btn_dup, self.btn_save,
                  self.btn_del, self.btn_imp, self.btn_exp,
                  self.btn_popout):
            b.setFixedHeight(26)
            prof_row.addWidget(b)
        outer.addLayout(prof_row)

        self.desc_lbl = QLabel("")
        self.desc_lbl.setWordWrap(True)
        self.desc_lbl.setStyleSheet(
            "color:#444; font-size:11px; padding:3px 0;")
        outer.addWidget(self.desc_lbl)

        # Vertical splitter: params (top) | template (bottom)
        if QT_VERSION == 6:
            split = QSplitter(Qt.Orientation.Vertical)
        else:
            split = QSplitter(Qt.Vertical)

        # Params panel
        self.params_box = QGroupBox(
            "Parameters  (edit values to update command)")
        self.params_layout = QGridLayout(self.params_box)
        self.params_layout.setSpacing(4)
        self.params_layout.setContentsMargins(8, 8, 8, 8)
        params_scroll = QScrollArea()
        params_scroll.setWidgetResizable(True)
        params_scroll.setWidget(self.params_box)
        params_scroll.setMinimumHeight(160)
        split.addWidget(params_scroll)

        # Template editor — styled as an interactive terminal
        # (dark background, monospace, prompt-like cursor).  Holds the
        # raw command(s) the user wants to run on their server.
        tmpl_box = QGroupBox("⌨  Terminal — command template")
        tmpl_l = QVBoxLayout(tmpl_box)
        tmpl_l.setContentsMargins(6, 6, 6, 6)
        cheat = QLabel(
            "<b>Always-available placeholders:</b> "
            "<code>{json_path}</code> · <code>{job_name}</code> · "
            "<code>{parent_dir}</code> · <code>{output_dir}</code> · "
            "<code>{prefix}</code> · <code>{date}</code> · "
            "<code>{timestamp}</code>"
            "<br>Add any other <code>{placeholder}</code> and it shows "
            "up as an editable parameter above.")
        cheat.setStyleSheet("color:#888; font-size:10px; padding:2px;")
        cheat.setWordWrap(True)
        tmpl_l.addWidget(cheat)
        self.tmpl_edit = QPlainTextEdit()
        self.tmpl_edit.setFont(QFont('Courier New', 10))
        self.tmpl_edit.setStyleSheet(
            "QPlainTextEdit {"
            "  background:#0c0c0c;"      # near-black terminal bg
            "  color:#e6e6e6;"           # light grey text
            "  selection-background-color:#264f78;"
            "  border:1px solid #333;"
            "  padding:6px;"
            "}")
        self.tmpl_edit.setPlaceholderText(
            "$ Type your AF3 / scheduler command here…\n"
            "  Lines starting with # are treated as comments and stripped\n"
            "  before submission. Use {placeholders} to inject runtime values.")
        self.tmpl_edit.textChanged.connect(self._on_template_changed)
        self.tmpl_edit.setMinimumHeight(160)
        tmpl_l.addWidget(self.tmpl_edit)
        split.addWidget(tmpl_box)
        split.setStretchFactor(0, 2)
        split.setStretchFactor(1, 4)   # give the terminal more space
        outer.addWidget(split, 1)

        # Resolved preview
        prev_box = QGroupBox("Resolved command preview")
        prev_l = QVBoxLayout(prev_box)
        prev_l.setContentsMargins(6, 6, 6, 6)
        self.preview = QPlainTextEdit()
        self.preview.setReadOnly(True)
        self.preview.setFont(QFont('Courier New', 9))
        self.preview.setStyleSheet("background:#1e1e1e; color:#7ec699;")
        self.preview.setMinimumHeight(100)
        prev_l.addWidget(self.preview)
        outer.addWidget(prev_box)

        self._refresh_profile_combo()

    def _refresh_profile_combo(self):
        self._building = True
        self.profile_combo.clear()
        for pid, label in self.mgr.names():
            self.profile_combo.addItem(label, userData=pid)
        self._building = False

    def _on_profile_changed(self, idx):
        if self._building or idx < 0:
            return
        pid = self.profile_combo.itemData(idx)
        if pid:
            self._load_profile(pid)

    def _load_profile(self, profile_id):
        prof = self.mgr.get(profile_id)
        if prof is None:
            return
        self._current_profile_id = profile_id
        self._building = True
        for i in range(self.profile_combo.count()):
            if self.profile_combo.itemData(i) == profile_id:
                self.profile_combo.setCurrentIndex(i)
                break
        self._building = False
        self.desc_lbl.setText(
            f"<i>{prof.get('description', '')}</i>"
            + (" <b>(built-in — clone to edit)</b>"
               if prof.get('is_builtin') else ""))
        self._building = True
        self.tmpl_edit.setPlainText(prof.get('template', ''))
        self._building = False
        merged_params = dict(prof.get('params', {}))
        detected = AF3CommandBuilder.detect_placeholders(prof.get('template', ''))
        auto_keys = ("json_path", "job_name", "parent_dir", "output_dir",
                     "prefix", "date", "timestamp")
        for n in detected:
            if n in auto_keys:
                continue
            if n in (prof.get('computed') or {}):
                continue
            merged_params.setdefault(n, "")
        self._build_param_widgets(merged_params, prof.get('is_builtin', False))
        self._refresh_preview()

    def _build_param_widgets(self, params, read_only):
        while self.params_layout.count():
            item = self.params_layout.takeAt(0)
            wdg = item.widget()
            if wdg is not None:
                wdg.deleteLater()
        self._param_widgets.clear()
        if not params:
            self.params_layout.addWidget(
                QLabel("<i>(this profile has no editable parameters)</i>"),
                0, 0)
            return
        for r, (name, default) in enumerate(params.items()):
            lbl = QLabel(f"<code>{{{name}}}</code>:")
            lbl.setMinimumWidth(160)
            edt = QLineEdit(str(default))
            edt.setReadOnly(read_only)
            if read_only:
                edt.setStyleSheet("background:#f4f4f4;")
                edt.setToolTip(
                    "Built-in profile values are read-only.\n"
                    "Click 'Duplicate' to create an editable copy.")
            edt.textChanged.connect(self._refresh_preview)
            self.params_layout.addWidget(lbl, r, 0)
            self.params_layout.addWidget(edt, r, 1)
            self._param_widgets[name] = edt
        self.params_layout.setColumnStretch(1, 1)

    def _on_template_changed(self):
        if self._building:
            return
        prof = (self.mgr.get(self._current_profile_id)
                if self._current_profile_id else None)
        is_builtin = prof.get('is_builtin', False) if prof else False
        if is_builtin:
            self._building = True
            self.tmpl_edit.setPlainText(prof['template'])
            self._building = False
            return
        new_template = self.tmpl_edit.toPlainText()
        detected = AF3CommandBuilder.detect_placeholders(new_template)
        current_values = {n: w.text() for n, w in self._param_widgets.items()}
        auto_keys = ("json_path", "job_name", "parent_dir", "output_dir",
                     "prefix", "date", "timestamp")
        new_params = {}
        for n in detected:
            if n in auto_keys:
                continue
            new_params[n] = current_values.get(n, "")
        for n, v in current_values.items():
            if n not in new_params:
                new_params[n] = v
        self._build_param_widgets(new_params, read_only=False)
        self._refresh_preview()

    def _refresh_preview(self, *_):
        prof = (self.mgr.get(self._current_profile_id)
                if self._current_profile_id else None)
        template = self.tmpl_edit.toPlainText()
        ctx = AF3CommandBuilder.auto_context(
            prefix="example_batch",
            parent_dir="~/af3_predictions/" + datetime.now().strftime('%Y-%m-%d'),
            json_path="example_batch_001.json")
        for n, w in self._param_widgets.items():
            ctx[n] = w.text()
        computed = (prof.get('computed') if prof else None) or {}
        cmd = AF3CommandBuilder.resolve(template, ctx, computed)
        self.preview.setPlainText(cmd)
        self.commandChanged.emit()

    def _on_new_profile(self):
        name, ok = QInputDialog.getText(self, "New profile", "Profile name:")
        if not ok or not name.strip():
            return
        new = {
            "id": (re.sub(r'[^a-z0-9_]', '_', name.lower())[:40]
                   + f"_{int(datetime.now().timestamp())}"),
            "name":        name.strip(),
            "description": "User-defined profile.",
            "is_builtin":  False,
            "template":    ("cd {parent_dir} && "
                            "{my_command} --json_path {json_path} "
                            "--job-name {job_name}"),
            "params":      {"my_command": "af3_run"},
        }
        self.mgr.add_or_update(new)
        self._refresh_profile_combo()
        self._load_profile(new['id'])

    def _on_duplicate_profile(self):
        if not self._current_profile_id:
            return
        src = self.mgr.get(self._current_profile_id)
        suggested = (src['name'] + " (copy)") if src else "New profile"
        name, ok = QInputDialog.getText(
            self, "Duplicate profile", "New name:", text=suggested)
        if not ok or not name.strip():
            return
        new = self.mgr.duplicate(self._current_profile_id, name.strip())
        if new:
            self._refresh_profile_combo()
            self._load_profile(new['id'])

    def _on_save_profile(self):
        if not self._current_profile_id:
            return
        prof = self.mgr.get(self._current_profile_id)
        if prof is None or prof.get('is_builtin'):
            QMessageBox.information(
                self, "Save profile",
                "Built-in profiles cannot be modified.\n"
                "Click 'Duplicate' to create an editable copy.")
            return
        prof['template'] = self.tmpl_edit.toPlainText()
        prof['params'] = {n: w.text() for n, w in self._param_widgets.items()}
        ok = self.mgr.add_or_update(prof)
        QMessageBox.information(
            self, "Save profile",
            "✓ Profile saved." if ok else "✗ Failed to save profile.")

    def _on_delete_profile(self):
        if not self._current_profile_id:
            return
        prof = self.mgr.get(self._current_profile_id)
        if prof is None or prof.get('is_builtin'):
            QMessageBox.information(
                self, "Delete profile",
                "Built-in profiles cannot be deleted.")
            return
        if QT_VERSION == 6:
            yes_btn = QMessageBox.StandardButton.Yes
        else:
            yes_btn = QMessageBox.Yes
        if QMessageBox.question(
                self, "Delete profile",
                f"Delete profile '{prof['name']}'?") != yes_btn:
            return
        if self.mgr.delete(self._current_profile_id):
            self._refresh_profile_combo()
            if self.mgr.profiles:
                self._load_profile(self.mgr.profiles[0]['id'])

    def _on_import_profile(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Import AF3 server profile", "", "JSON (*.json)")
        if not path:
            return
        try:
            with open(path, 'r', encoding='utf-8') as f:
                data = json.load(f)
            if isinstance(data, dict):
                data = [data]
            cnt = 0
            for p in data:
                if isinstance(p, dict) and p.get('id'):
                    p['is_builtin'] = False
                    if self.mgr.add_or_update(p):
                        cnt += 1
            self._refresh_profile_combo()
            QMessageBox.information(
                self, "Import", f"✓ Imported {cnt} profile(s).")
        except Exception as e:
            QMessageBox.critical(self, "Import failed", str(e))

    def _on_export_profile(self):
        if not self._current_profile_id:
            return
        prof = self.mgr.get(self._current_profile_id)
        if prof is None:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export profile",
            f"{prof['id']}.json", "JSON (*.json)")
        if not path:
            return
        try:
            with open(path, 'w', encoding='utf-8') as f:
                json.dump(prof, f, indent=2)
            QMessageBox.information(self, "Export", f"✓ Saved to {path}")
        except Exception as e:
            QMessageBox.critical(self, "Export failed", str(e))

    def _on_popout(self):
        """Pop this widget out into a separate, resizable window.

        Useful when the host tab is too cramped to comfortably edit
        the Parameters panel or the command template.  Closing the
        pop-out window re-docks the widget back into its host layout.
        """
        # Already detached?
        if getattr(self, '_popout_window', None) is not None:
            try:
                self._popout_window.raise_()
                self._popout_window.activateWindow()
            except Exception:
                pass
            return

        host = self.parent()
        # Figure out the host layout & index so we can dock back later
        host_layout = None
        host_index  = -1
        if host is not None:
            host_layout = host.layout()
            if host_layout is not None:
                for i in range(host_layout.count()):
                    if host_layout.itemAt(i).widget() is self:
                        host_index = i
                        break

        # Build a floating dialog
        if QT_VERSION == 6:
            flags = Qt.WindowType.Window
        else:
            flags = Qt.Window
        dlg = QDialog(host)
        dlg.setWindowFlags(flags)
        dlg.setWindowTitle("⌨  AF3 Server Submission — Terminal & Profiles")
        try:
            scr = QApplication.primaryScreen().availableGeometry()
            dlg.resize(int(scr.width() * 0.70), int(scr.height() * 0.85))
            geo = dlg.frameGeometry()
            geo.moveCenter(scr.center())
            dlg.move(geo.topLeft())
        except Exception:
            dlg.resize(1100, 800)

        dlg_lay = QVBoxLayout(dlg)
        dlg_lay.setContentsMargins(6, 6, 6, 6)
        dlg_lay.setSpacing(4)

        info = QLabel(
            "<b>📌 Detached AF3 Builder</b> &nbsp; "
            "<span style='color:#666; font-size:10px;'>"
            "(close this window to re-dock the builder back into the "
            "Submit AF3 tab)</span>")
        info.setStyleSheet("padding:4px;")
        dlg_lay.addWidget(info)

        # Re-parent self into the dialog.  Same caveat as
        # _DetachedTabWindow: when a widget is removed from one layout
        # and re-parented to another, Qt may keep visibility=False —
        # we have to explicitly setVisible(True) and show().
        self.setParent(dlg)
        dlg_lay.addWidget(self, 1)
        self.setVisible(True)
        self.show()

        # Hide our own pop-out button while we ARE the popout
        self.btn_popout.setEnabled(False)
        self.btn_popout.setText("↙ (already detached)")

        # Re-dock on close
        def _on_dlg_close(event):
            try:
                if host_layout is not None:
                    self.setParent(host)
                    if host_index >= 0:
                        host_layout.insertWidget(host_index, self)
                    else:
                        host_layout.addWidget(self)
                    self.setVisible(True)
                    self.show()
            except Exception as e:
                print(f"[AF3 popout] re-dock error: {e}")
            self._popout_window = None
            self.btn_popout.setEnabled(True)
            self.btn_popout.setText("↗ Pop-out")
            event.accept()

        dlg.closeEvent = _on_dlg_close
        self._popout_window = dlg
        dlg.show()

    # ── Public API used by ppigFinderApp ────────────────────────
    def build_command(self, context):
        """Resolve current template against runtime context."""
        prof = (self.mgr.get(self._current_profile_id)
                if self._current_profile_id else None)
        template = self.tmpl_edit.toPlainText()
        ctx = dict(context)
        for n, w in self._param_widgets.items():
            ctx.setdefault(n, w.text())
        computed = (prof.get('computed') if prof else None) or {}
        return AF3CommandBuilder.resolve(template, ctx, computed)

    def current_profile_summary(self):
        prof = (self.mgr.get(self._current_profile_id)
                if self._current_profile_id else None)
        if not prof:
            return ""
        return f"# AF3 profile: {prof['name']}"

    def get_state(self):
        """Snapshot for project save."""
        return {
            "profile_id": self._current_profile_id,
            "template":   self.tmpl_edit.toPlainText(),
            "params":     {n: w.text()
                           for n, w in self._param_widgets.items()},
        }

    def set_state(self, state):
        """Restore from project save."""
        if not state:
            return
        pid = state.get("profile_id")
        if pid and self.mgr.get(pid):
            self._load_profile(pid)
            tmpl = state.get("template")
            if tmpl:
                self._building = True
                self.tmpl_edit.setPlainText(tmpl)
                self._building = False
            for n, v in (state.get("params") or {}).items():
                if n in self._param_widgets:
                    self._param_widgets[n].setText(v)
            self._refresh_preview()


# ═══════════════════════════════════════════════════════════════
# END of helper classes.  Main app class follows.
# ═══════════════════════════════════════════════════════════════



try:
    from .io.fasta import (
        choose_longest_record as _io_choose_longest_record,
        read_fasta as _io_read_fasta,
        write_orf_protein_fasta as _io_write_orf_protein_fasta,
    )
except ImportError:
    from ppigfinder.io.fasta import (
        choose_longest_record as _io_choose_longest_record,
        read_fasta as _io_read_fasta,
        write_orf_protein_fasta as _io_write_orf_protein_fasta,
    )


try:
    from .ui.file_opening import open_genome_file_into_window as _ui_open_genome_file_into_window
except ImportError:
    from ppigfinder.ui.file_opening import open_genome_file_into_window as _ui_open_genome_file_into_window


try:
    from .io.html_report import write_basic_report as _io_write_basic_report
except ImportError:
    from ppigfinder.io.html_report import write_basic_report as _io_write_basic_report


try:
    from .ui.icon_provider import make_icon as _ui_make_icon
    from .ui.text_fallback import clean_ui_text as _ui_clean_text
except ImportError:
    from ppigfinder.ui.icon_provider import make_icon as _ui_make_icon
    from ppigfinder.ui.text_fallback import clean_ui_text as _ui_clean_text


try:
    from .services.blast_service import (
        BlastSearchParams as _BlastSearchParams,
        BlastSearchService as _BlastSearchService,
    )
except ImportError:
    from ppigfinder.services.blast_service import (
        BlastSearchParams as _BlastSearchParams,
        BlastSearchService as _BlastSearchService,
    )

class ppigFinderApp(QMainWindow):

    def __init__(self):
        super().__init__()
        self.setWindowTitle('🧬 ppigFinder v2.00 — Protein-Protein Interaction Genomic Finder')
        # ── Responsive window sizing (replaces fixed 1550×980) ──────
        self.win_mgr = WindowManager(self, app_name="ppigFinder",
                                      organisation="ppigFinder")
        # Default: 85% of screen width × 90% of screen height,
        # with 1100×700 minimum.  User-configured size is restored
        # at the end of __init__ from QSettings (see win_mgr.restore_geometry()).
        self.win_mgr.apply_default_size(width_pct=0.85, height_pct=0.90,
                                        min_w=1100, min_h=700)

        self.analyzer = AdvancedORFAnalyzer()
        self.dna_sequence = ""
        self.genome_name = ""
        self.current_fasta_path = ""   # source genome file for project save
        self.orfs = []
        self.filtered_orfs = []
        self.selected_orf = None
        self.selected_orf_idx = -1
        self.zoom_level = 1.0
        self.candidates = []
        self._workers = []  # keep refs to prevent GC

        # HMM
        self.hmm_profiles = []
        self.hmm_hits_all = []
        self.hmm_default_colors = ['#E53935','#1E88E5','#43A047','#FB8C00',
                                    '#8E24AA','#00ACC1','#D81B60','#5E35B1']

        # AlphaFold3
        self.af3_jobs = []
        # Pyrodigal parameters (persisted across runs)
        self._pyro_params = {
            'meta':              True,
            'translation_table': 11,
            'min_aa':            30,
            'closed':            False,
            'mask':              False,
        }
        self.af3_n_neighbors = 5
        self.af3_max_residues = 5000

        # SnapGene
        self.snapgene_features = []
        self.snapgene_primers = []

        # HPC server
        self._ssh_client = None
        self._sftp_client = None
        self._hpc_jobs = []   # [{name, slurm_id, remote_dir, remote_json, status, local_output}]
        self._dv_pending_jobs = []   # jobs staged for upload in the Submit tab
        self._hpc_poll_timer = QTimer()
        self._hpc_poll_timer.timeout.connect(self._hpc_poll_queue)
        self._hpc_workers = []
        self._dv_module_loaded = False   # True after successful module load

        # AF3 Analysis (v2)
        self._af3_analysis_results = []  # list of parsed job dicts
        self._af3_analysis_dir    = ''   # last loaded results folder

        # BLAST params (stored as plain values, no tkinter vars)
        self.blast_matrix = "BLOSUM62"
        self.blast_gap_open = 11
        self.blast_gap_ext = 1
        self.blast_word_size = 5
        self.blast_max_targets = 100
        self.blast_low_complexity = True
        self.blast_program = "blastp"
        self.blast_threshold = 30
        self.blast_evalue = 0.05

        # HMM params
        self.hmm_evalue = 10.0
        self.hmm_dom_evalue = 10.0
        self.hmm_score_thresh = None

        # Config
        self.min_length = 30
        self.start_codons = {'ATG', 'GTG', 'TTG'}
        self.algo_choice = "Auto (best available)"

        self.colors = ["#FF6B6B","#4ECDC4","#45B7D1","#FFA07A",
                        "#98D8C8","#A29BFE","#FD79A8","#00B894"]

        self._setup_ui()

        # ── Add Window menu (size presets, custom-size dialog, F11 fullscreen)
        self._build_window_menu()

        # ── Restore last user-chosen geometry (size + position) ────
        # Call AFTER _setup_ui so the menubar/toolbar are already there
        # and saveState() snapshots the right docks.
        try:
            self.win_mgr.restore_geometry()
        except Exception:
            pass

    # ═══════════════════════════════════════════════════════════
    # UI SETUP
    # ═══════════════════════════════════════════════════════════

    def _setup_ui(self):
        self._create_menus()
        self._create_toolbar()
        self._create_central()
        self._create_statusbar()
        self._normalize_toolbar_actions()

    # ─── MENUS ─────────────────────────────────────────────────

    def _create_menus(self):
        mb = self.menuBar()

        # File
        fm = mb.addMenu(t('menu_file'))
        act = fm.addAction(t('open_fasta'), self.load_fasta)
        act.setShortcut('Ctrl+O')
        fm.addAction(t('load_hmm'), self.load_hmm)
        fm.addSeparator()
        act = fm.addAction(t('save_project'), self.save_project)
        act.setShortcut('Ctrl+S')
        act = fm.addAction('💾 Save Project As (full copy)...', self.save_project_as)
        act.setShortcut('Ctrl+Shift+S')
        act = fm.addAction(t('open_project'), self.load_project)
        act.setShortcut('Ctrl+Shift+O')
        fm.addSeparator()
        fm.addAction(t('save_orfs_fasta'), self.save_fasta)
        fm.addAction(t('save_report_tsv'), self.save_report_tsv)
        fm.addAction(t('save_report_html'), self.export_html_report)
        fm.addSeparator()
        act = fm.addAction(t('quit'), self.close)
        act.setShortcut('Ctrl+Q')

        # Parameters
        pm = mb.addMenu(t('menu_params'))
        pm.addAction("🧬 ORF Analysis Parameters...", self._show_orf_params)
        pm.addAction(t('blast_params'), self._show_blast_params)
        pm.addAction(t('hmm_params'), self._show_hmm_params)

        # Help
        hm = mb.addMenu(t('menu_help'))
        hm.addAction(t('manual'), self._show_manual)
        hm.addAction(t('tutorial'), self._show_tutorial)
        hm.addAction(t('install'), self._show_install)
        hm.addSeparator()
        hm.addAction("📈 Interaction Results — analysis guide",
                     self._show_help_interaction_results)
        hm.addAction("🧬 Genomic PPI Map — guide",
                     self._show_help_ppi_map)
        hm.addSeparator()
        hm.addAction("📚 References & methodology",
                     self._show_help_references)
        hm.addSeparator()
        hm.addAction(t('about'), self._show_about)

    # ─── WINDOW MENU (responsive sizing) ──────────────────────

    def _build_window_menu(self):
        """Add a Window menu with size presets, Custom Size dialog,
        Fullscreen toggle, and re-center / reset actions.

        Idempotent: safe to call multiple times — re-uses existing menu.
        """
        mb = self.menuBar()
        # Look for an existing Window menu first
        win_menu = None
        for act in mb.actions():
            if act.text() in ('&Window', 'Window'):
                win_menu = act.menu()
                break
        if win_menu is None:
            win_menu = mb.addMenu("&Window")
        else:
            win_menu.clear()

        # Quick presets
        for name, val in WIN_PRESETS_RELATIVE.items():
            act = QAction(name, self)
            if val is None:
                act.triggered.connect(self.win_mgr.apply_maximized)
            else:
                wp, hp = val
                act.triggered.connect(
                    lambda _checked=False, w=wp, h=hp:
                        self.win_mgr.apply_relative(w, h)
                )
            win_menu.addAction(act)

        win_menu.addSeparator()

        act_custom = QAction("&Custom size...", self)
        act_custom.setShortcut("Ctrl+Shift+W")
        act_custom.setToolTip(
            "Open a dialog to dial-in any window size, "
            "either as a percentage of the screen or in pixels.")
        act_custom.triggered.connect(self.win_mgr.open_custom_size_dialog)
        win_menu.addAction(act_custom)

        act_fs = QAction("&Toggle Fullscreen", self)
        act_fs.setShortcut("F11")
        def _toggle_fs():
            if self.isFullScreen():
                self.showNormal()
            else:
                self.showFullScreen()
        act_fs.triggered.connect(_toggle_fs)
        win_menu.addAction(act_fs)

        win_menu.addSeparator()
        act_center = QAction("Re-center on screen", self)
        act_center.triggered.connect(self.win_mgr.center_on_screen)
        win_menu.addAction(act_center)

        act_reset = QAction("Reset to &default size  (85% × 90%)", self)
        act_reset.triggered.connect(
            lambda: self.win_mgr.apply_default_size(0.85, 0.90, 1100, 700))
        win_menu.addAction(act_reset)

    def closeEvent(self, event):
        """Persist user's window size/position for next session."""
        try:
            if hasattr(self, 'win_mgr'):
                self.win_mgr.save_geometry()
        except Exception:
            pass
        try:
            super().closeEvent(event)
        except Exception:
            event.accept()

    # ─── TOOLBAR ───────────────────────────────────────────────

    def _create_toolbar(self):
        tb = QToolBar("Main Toolbar")
        tb.setObjectName("main_toolbar")
        tb.setMovable(False)
        self.addToolBar(tb)

        self._btn_open = QPushButton(t('btn_open'))
        self._btn_open.clicked.connect(self.load_fasta)
        self._btn_open.setToolTip(t('tip_open'))
        tb.addWidget(self._btn_open)

        # Translate genome dropdown button
        self._btn_translate = QPushButton(t('btn_translate_genome'))
        translate_menu = QMenu(self._btn_translate)
        
        # Pyrodigal option — runs directly with current _pyro_params
        pyrodigal_action = translate_menu.addAction(
            f"{t('btn_pyrodigal')} — {t('desc_pyrodigal')}")
        pyrodigal_action.triggered.connect(self.analyze_orfs_pyrodigal)

        # Hybrid option — Pyrodigal primary + 6-frame gap-filler
        hybrid_action = translate_menu.addAction(
            f"{t('btn_hybrid')} — {t('desc_hybrid')}")
        hybrid_action.triggered.connect(self.analyze_orfs_hybrid)

        # Automatic option with description
        automatic_action = translate_menu.addAction(
            f"{t('btn_automatic')} — {t('desc_automatic')}")
        automatic_action.triggered.connect(self.analyze_orfs)

        self._btn_translate.setMenu(translate_menu)
        self._btn_translate.setToolTip(
            "Choose gene prediction method:\n"
            f"• {t('btn_pyrodigal')}: {t('desc_pyrodigal')}\n"
            f"• {t('btn_hybrid')}: {t('desc_hybrid')}\n"
            f"• {t('btn_automatic')}: {t('desc_automatic')}")
        tb.addWidget(self._btn_translate)

        self._btn_hmm = QPushButton(t('btn_annotate_hmm'))
        self._btn_hmm.clicked.connect(self.classify_all_domains)
        self._btn_hmm.setToolTip(t('tip_hmm'))
        tb.addWidget(self._btn_hmm)

        btn_pdf = QPushButton("Export Map PDF")
        btn_pdf.clicked.connect(self.export_map_pdf)
        btn_pdf.setToolTip("Export the genome map as PDF/PNG image")
        tb.addWidget(btn_pdf)

        tb.addSeparator()

        tb.addWidget(QLabel(f"  {t('zoom_label')} "))
        btn_zm = QPushButton("−")
        btn_zm.setFixedWidth(30)
        btn_zm.clicked.connect(lambda: self._set_zoom(self.zoom_level * 0.8))
        btn_zm.setToolTip(t('tip_zoom_minus'))
        tb.addWidget(btn_zm)
        btn_zp = QPushButton("+")
        btn_zp.setFixedWidth(30)
        btn_zp.clicked.connect(lambda: self._set_zoom(self.zoom_level * 1.2))
        btn_zp.setToolTip(t('tip_zoom_plus'))
        tb.addWidget(btn_zp)
        self._zoom_label = QLabel("100%")
        tb.addWidget(self._zoom_label)

        tb.addSeparator()

        # Backend status
        blast_st = "✅" if BACKENDS.get('blast+',{}).get('available') else "❌"
        hmmer_st = "✅" if BACKENDS.get('hmmer3',{}).get('available') else "❌"
        pyrod_st = "✅" if BACKENDS.get('pyrodigal',{}).get('available') else "❌"
        lbl = QLabel(f"  BLAST+{blast_st}  HMMER3{hmmer_st}  Pyrodigal{pyrod_st}")
        lbl.setStyleSheet("font-size: 10px; color: #666;")
        tb.addWidget(lbl)

    # ─── STATUS BAR ────────────────────────────────────────────

    def _create_statusbar(self):
        self._status = self.statusBar()
        self._status.showMessage(t('ready_status'))

    # ─── CENTRAL LAYOUT ───────────────────────────────────────


    def _normalize_toolbar_actions(self):
        icon_map = {
            "open": "open",
            "load": "open",
            "genome": "orf",
            "orf": "orf",
            "pyrodigal": "orf",
            "automatic": "settings",
            "hybrid": "settings",
            "hmm": "hmm",
            "blast": "blast",
            "af3": "af3",
            "export": "export",
            "pdf": "export",
            "html": "export",
            "hpc": "hpc",
        }

        try:
            actions = self.findChildren(QAction)
        except Exception:
            return

        for action in actions:
            try:
                original = action.text()
                cleaned = _ui_clean_text(original)
                if cleaned and cleaned != original:
                    action.setText(cleaned)

                lower = cleaned.lower()
                for key, icon_name in icon_map.items():
                    if key in lower:
                        action.setIcon(_ui_make_icon(icon_name))
                        break
            except Exception:
                pass


    def _create_central(self):
        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QHBoxLayout(central)
        main_layout.setContentsMargins(3, 3, 3, 3)

        splitter = QSplitter(Horizontal)
        main_layout.addWidget(splitter)

        # LEFT: Map + ORF table (no config panel — moved to menus)
        center = self._create_center_panel()
        splitter.addWidget(center)

        # RIGHT: Tabs
        right = self._create_right_panel()
        splitter.addWidget(right)

        splitter.setSizes([700, 600])

    # ─── HIDDEN INIT: variables previously in left panel ────────

    def _init_filter_vars(self):
        """Initialize filter/config variables (no longer in a visible left panel)."""
        self._min_length_spin = QSpinBox(); self._min_length_spin.setRange(10, 500); self._min_length_spin.setValue(30)
        self._cb_atg = QCheckBox("ATG"); self._cb_atg.setChecked(True)
        self._cb_gtg = QCheckBox("GTG"); self._cb_gtg.setChecked(True)
        self._cb_ttg = QCheckBox("TTG"); self._cb_ttg.setChecked(True)
        self._info_text = None  # replaced by Genome tab
        # BLAST config
        self._algo_combo = QComboBox()
        algo_choices = ["Auto (best available)"]
        if BACKENDS.get('blast+',{}).get('available'):
            algo_choices.append("NCBI BLAST+ (external)")
        algo_choices.extend(["K-mer Filter (fast)", "Smith-Waterman (sensitive)"])
        self._algo_combo.addItems(algo_choices)
        self._identity_spin = QSpinBox(); self._identity_spin.setRange(0, 100); self._identity_spin.setValue(30)
        self._evalue_edit = QLineEdit("0.05")

    # ─── CENTER PANEL: Map + Filter bar + ORF Table ───────────

    def _create_center_panel(self):
        self._init_filter_vars()
        w = QWidget()
        layout = QVBoxLayout(w)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(2)

        # Search bar
        search_frame = QGroupBox(t('locate_orf'))
        sf_layout = QHBoxLayout(search_frame)
        sf_layout.setContentsMargins(4, 2, 4, 2)
        self._map_search_edit = QLineEdit()
        self._map_search_edit.setFont(QFont('Courier', 9))
        self._map_search_edit.setPlaceholderText("ORF123 / protein sequence / DNA")
        self._map_search_edit.returnPressed.connect(self.search_orf_in_map)
        sf_layout.addWidget(self._map_search_edit)
        btn_search = QPushButton(t('search_btn'))
        btn_search.clicked.connect(self.search_orf_in_map)
        btn_search.setToolTip(t('tip_search_orf'))
        sf_layout.addWidget(btn_search)
        layout.addWidget(search_frame)

        # Genome Map + Hits Legend (side by side)
        map_group = QGroupBox(t('map_title'))
        map_outer = QHBoxLayout(map_group)
        map_outer.setContentsMargins(2, 2, 2, 2)
        map_outer.setSpacing(2)

        self._genome_map = GenomeMapWidget()
        self._genome_map.orf_clicked.connect(self._on_map_orf_click)
        self._genome_map.zoom_changed.connect(self._on_map_zoom_changed)
        map_outer.addWidget(self._genome_map, stretch=1)

        # Hits legend list (right side of map)
        self._hits_legend = QTableWidget()
        self._hits_legend.setColumnCount(2)
        self._hits_legend.setHorizontalHeaderLabels(['Color', 'Hit'])
        self._hits_legend.setColumnWidth(0, 30)
        self._hits_legend.horizontalHeader().setStretchLastSection(True)
        self._hits_legend.setMaximumWidth(200)
        self._hits_legend.setMinimumWidth(140)
        self._hits_legend.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers
                                           if QT_VERSION == 6 else QAbstractItemView.NoEditTriggers)
        self._hits_legend.setSelectionBehavior(SelectRows)
        self._hits_legend.verticalHeader().setVisible(False)
        self._hits_legend.setAlternatingRowColors(True)
        self._hits_legend.setStyleSheet("font-size: 10px;")
        self._hits_legend.selectionModel().selectionChanged.connect(self._on_legend_click)
        map_outer.addWidget(self._hits_legend)
        layout.addWidget(map_group)

        # ── Compact filter bar (replaces old left panel) ──
        filter_bar = QFrame()
        filter_bar.setFrameShape(QFrame.Shape.StyledPanel if QT_VERSION == 6 else QFrame.StyledPanel)
        fb = QHBoxLayout(filter_bar)
        fb.setContentsMargins(4, 2, 4, 2)
        fb.setSpacing(6)

        fb.addWidget(QLabel("🔎"))
        self._search_edit = QLineEdit()
        self._search_edit.setPlaceholderText("Search ORF / protein...")
        self._search_edit.setMaximumWidth(160)
        self._search_edit.textChanged.connect(self.filter_orfs)
        self._search_edit.setToolTip(t('tip_search'))
        fb.addWidget(self._search_edit)

        fb.addWidget(QLabel("Frame:"))
        self._frame_combo = QComboBox()
        self._frame_combo.addItems(["All","0","1","2","3","4","5"])
        self._frame_combo.setMaximumWidth(60)
        self._frame_combo.currentTextChanged.connect(self.filter_orfs)
        self._frame_combo.setToolTip(t('tip_frame'))
        fb.addWidget(self._frame_combo)

        fb.addWidget(QLabel("Strand:"))
        self._strand_combo = QComboBox()
        self._strand_combo.addItems(["All","+","-"])
        self._strand_combo.setMaximumWidth(55)
        self._strand_combo.currentTextChanged.connect(self.filter_orfs)
        self._strand_combo.setToolTip(t('tip_strand'))
        fb.addWidget(self._strand_combo)

        fb.addWidget(QLabel("Min aa:"))
        self._size_filter_spin = QSpinBox()
        self._size_filter_spin.setRange(0, 5000)
        self._size_filter_spin.setMaximumWidth(65)
        self._size_filter_spin.setToolTip(t('tip_min_aa'))
        fb.addWidget(self._size_filter_spin)

        fb.addWidget(QLabel("Source:"))
        self._source_combo = QComboBox()
        self._source_combo.addItems(["All", "6frame", "pyrodigal"])
        self._source_combo.setMaximumWidth(90)
        self._source_combo.currentTextChanged.connect(self.filter_orfs)
        fb.addWidget(self._source_combo)

        self._btn_apply = QPushButton("▶ Apply")
        self._btn_apply.clicked.connect(self.filter_orfs)
        self._btn_apply.setToolTip(t('tip_apply'))
        fb.addWidget(self._btn_apply)
        fb.addStretch()

        self._orf_count_label = QLabel("(0 of 0)")
        self._orf_count_label.setStyleSheet("font-weight: bold;")
        fb.addWidget(self._orf_count_label)

        layout.addWidget(filter_bar)

        # ORF Table
        self._orf_table = QTableWidget()
        cols = ['ID','Frame','Strand','Start','End','Size(aa)','GC%',
                'HMM','Score','Source','Obs',
                'AF3','Partner','ipTM','PAE_inter',
                'PAE_min ★','cp_ipTM ★',
                'Contact_region','User_note']
        self._orf_table.setColumnCount(len(cols))
        self._orf_table.setHorizontalHeaderLabels(cols)
        self._orf_table.setSelectionBehavior(SelectRows)
        self._orf_table.setSelectionMode(
            QAbstractItemView.SelectionMode.ExtendedSelection if QT_VERSION == 6
            else QAbstractItemView.ExtendedSelection)
        self._orf_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers
                                         if QT_VERSION == 6 else QAbstractItemView.NoEditTriggers)
        self._orf_table.horizontalHeader().setStretchLastSection(True)
        self._orf_table.setAlternatingRowColors(True)
        self._orf_table.selectionModel().selectionChanged.connect(self._on_orf_table_select)
        # Ctrl+C copies selected rows
        copy_sc = QShortcut(
            QKeySequence.StandardKey.Copy if QT_VERSION == 6 else QKeySequence.Copy,
            self._orf_table)
        copy_sc.activated.connect(self._orf_table_copy_selection)
        self._orf_table.setContextMenuPolicy(
            Qt.ContextMenuPolicy.CustomContextMenu if QT_VERSION == 6 else Qt.CustomContextMenu)
        self._orf_table.customContextMenuRequested.connect(self._on_orf_right_click)
        self._orf_table.setSortingEnabled(True)

        for i, cw in enumerate([55,42,42,70,70,52,42,180,45,60,90,
                                  30,90,50,70,200,120]):
            self._orf_table.setColumnWidth(i, cw)

        # User_note column is editable
        self._orf_table.setEditTriggers(
            QAbstractItemView.EditTrigger.DoubleClicked
            if QT_VERSION == 6 else QAbstractItemView.DoubleClicked)

        # Export toolbar
        export_bar = QHBoxLayout()
        export_bar.setSpacing(4)
        btn_copy_sel = QPushButton("Copy selected rows")
        btn_copy_sel.setToolTip(
            "Copy selected rows as tab-separated text (paste into Excel/Calc).\n"
            "Select rows with mouse (Ctrl/Shift for multi-select), then click.")
        btn_copy_sel.clicked.connect(self._orf_table_copy_selection)
        export_bar.addWidget(btn_copy_sel)

        btn_copy_all = QPushButton("Copy all rows")
        btn_copy_all.setToolTip("Copy ALL visible rows (current filter) as TSV.")
        btn_copy_all.clicked.connect(lambda: self._orf_table_copy_rows(all_rows=True))
        export_bar.addWidget(btn_copy_all)

        btn_exp = QPushButton("Export table...")
        exp_menu = QMenu(btn_exp)
        exp_menu.addAction(
            "TSV — table columns only",
            lambda: self._export_orf_table(fmt='tsv', include_seqs=False))
        exp_menu.addAction(
            "TSV — full (+ DNA + Protein sequences)",
            lambda: self._export_orf_table(fmt='tsv', include_seqs=True))
        exp_menu.addSeparator()
        exp_menu.addAction(
            "FASTA — protein sequences",
            lambda: self._export_orf_fasta(aa=True))
        exp_menu.addAction(
            "FASTA — DNA sequences",
            lambda: self._export_orf_fasta(aa=False))
        exp_menu.addSeparator()
        exp_menu.addAction(
            "TSV — annotated only (observation/function/gene)",
            lambda: self._export_orf_table(fmt='tsv', annotated_only=True, include_seqs=True))
        btn_exp.setMenu(exp_menu)
        btn_exp.setToolTip("Export ORF table in various formats")
        export_bar.addWidget(btn_exp)
        export_bar.addStretch()
        self._orf_export_info = QLabel("")
        self._orf_export_info.setStyleSheet("color: #555; font-size: 11px;")
        export_bar.addWidget(self._orf_export_info)
        layout.addLayout(export_bar)

        layout.addWidget(self._orf_table, stretch=1)

        return w

    # ─── RIGHT PANEL: Tabs ────────────────────────────────────

    def _create_right_panel(self):
        # DetachableTabWidget lets users tear off any tab into a
        # free-floating, resizable window (right-click on a tab,
        # or double-click it).
        self._tabs = DetachableTabWidget()
        self._tabs.setMinimumWidth(450)
        self._tabs.setToolTip(
            "Tip: right-click any tab (or double-click it) to open it in "
            "a separate window — useful when a tab needs more room.")

        # Tab 0: Genome (NEW — moved from left panel)
        self._create_genome_tab()
        # Tab 1: BLAST Query
        self._create_blast_query_tab()
        # Tab 2: BLAST Results
        self._create_blast_results_tab()
        # Tab 3: DNA
        self._create_dna_tab()
        # Tab 4: Protein
        self._create_protein_tab()
        # Tab 5: Domains
        self._create_domains_tab()
        # Tab 6: Neighborhood
        self._create_neighborhood_tab()
        # Tab 7: HMM
        self._create_hmm_tab()
        # Tab 8: AlphaFold
        self._create_af3_tab()
        # Tab 9: Server
        self._create_hpc_server_tab()
        # Tab 10: AF3 Analysis (v2) — after Server
        self._create_af3_analysis_tab()
        # Tab 11: PPI Genomic Arc Map
        self._create_ppi_arc_map_tab()

        return self._tabs

    def _create_genome_tab(self):
        """Genome info tab — replaces old left-panel genome section."""
        w = QWidget()
        layout = QVBoxLayout(w)
        layout.addWidget(QLabel("<b>📊 Genome Information</b>"))
        self._info_text = QTextEdit()
        self._info_text.setReadOnly(True)
        self._info_text.setFont(QFont('Courier', 9))
        self._info_text.setMaximumHeight(200)
        layout.addWidget(self._info_text)

        # Backends status
        blast_ok = "✅" if BACKENDS.get('blast+',{}).get('available') else "❌"
        hmmer_ok = "✅" if BACKENDS.get('hmmer3',{}).get('available') else "❌"
        pyrod_ok = "✅" if BACKENDS.get('pyrodigal',{}).get('available') else "❌"
        status_txt = QTextEdit()
        status_txt.setReadOnly(True)
        status_txt.setFont(QFont('Courier', 9))
        status_txt.setMaximumHeight(100)
        status_txt.setPlainText(
            f"Backends:\n"
            f"  BLAST+     {blast_ok}  {BACKENDS.get('blast+',{}).get('version','')}\n"
            f"  HMMER3     {hmmer_ok}  {BACKENDS.get('hmmer3',{}).get('version','')}\n"
            f"  Pyrodigal  {pyrod_ok}  {BACKENDS.get('pyrodigal',{}).get('version','')}")
        layout.addWidget(QLabel("<b>🔧 Backends</b>"))
        layout.addWidget(status_txt)
        layout.addStretch()
        self._tabs.addTab(w, "📊 Genome")


    # ═══════ TAB CREATION ═══════

    def _create_blast_query_tab(self):
        w = QWidget()
        layout = QVBoxLayout(w)
        layout.addWidget(QLabel(t('blast_query_title')))

        # Buttons row
        btn_row = QHBoxLayout()
        for text, slot in [(t('blast_load_fasta'), self.load_blast_query_fasta),
                           (t('blast_paste'), self._paste_clipboard),
                           (t('blast_clear'), self._clear_query),
                           (t('blast_validate'), self._validate_query)]:
            b = QPushButton(text); b.clicked.connect(slot); btn_row.addWidget(b)
            # Add tooltips
            if 'load_fasta' in slot.__name__:
                b.setToolTip(t('tip_blast_load'))
            elif 'paste' in slot.__name__:
                b.setToolTip(t('tip_blast_paste'))
            elif 'clear' in slot.__name__:
                b.setToolTip(t('tip_blast_clear'))
            elif 'validate' in slot.__name__:
                b.setToolTip(t('tip_blast_validate'))

        btn_row.addWidget(QLabel(f"  {t('blast_program')}"))
        self._blast_prog_combo = QComboBox()
        self._blast_prog_combo.addItems(["blastp","tblastn","blastn","blastx"])
        self._blast_prog_combo.currentTextChanged.connect(
            lambda: self._refresh_blast_cmd_preview())
        btn_row.addWidget(self._blast_prog_combo)

        btn_run = QPushButton(t('blast_run'))
        btn_run.clicked.connect(self.run_blast)
        btn_run.setStyleSheet("font-weight: bold;")
        btn_run.setToolTip(t('tip_blast_run'))
        btn_row.addWidget(btn_run)
        btn_row.addStretch()
        layout.addLayout(btn_row)

        # Query text
        qf = QGroupBox(t('blast_query_frame'))
        qf_layout = QVBoxLayout(qf)
        self._blast_query_text = QTextEdit()
        self._blast_query_text.setFont(QFont('Courier New', 10))
        self._blast_query_text.setPlaceholderText("Paste your protein sequence here...")
        self._blast_query_text.setMaximumHeight(200)
        qf_layout.addWidget(self._blast_query_text)
        layout.addWidget(qf)

        # Command Preview
        cf = QGroupBox(t('blast_cmd_preview'))
        cf_layout = QVBoxLayout(cf)
        cf_layout.addWidget(QLabel(t('blast_cmd_hint')))
        self._blast_cmd_text = QTextEdit()
        self._blast_cmd_text.setFont(QFont('Courier New', 9))
        self._blast_cmd_text.setMaximumHeight(100)
        self._blast_cmd_text.setStyleSheet(
            "background-color: #1E1E1E; color: #D4D4D4; "
            "selection-background-color: #264F78;")
        cf_layout.addWidget(self._blast_cmd_text)
        cmd_btns = QHBoxLayout()
        b_ref = QPushButton(t('blast_refresh_cmd'))
        b_ref.clicked.connect(self._refresh_blast_cmd_preview)
        cmd_btns.addWidget(b_ref)
        b_cust = QPushButton(t('blast_run_custom'))
        b_cust.clicked.connect(self._run_custom_blast_cmd)
        cmd_btns.addWidget(b_cust)
        cmd_btns.addStretch()
        cf_layout.addLayout(cmd_btns)
        layout.addWidget(cf)

        # Status
        self._query_status = QLabel("Waiting for sequence...")
        layout.addWidget(self._query_status)
        layout.addStretch()

        self._tabs.addTab(w, t('tab_blast_query'))
        self._refresh_blast_cmd_preview()

    def _create_blast_results_tab(self):
        w = QWidget()
        layout = QVBoxLayout(w)
        btn_row = QHBoxLayout()
        for text, slot in [(t('blast_copy_hit'), self._copy_blast_hit),
                           (t('blast_copy_all'), self._copy_blast_all),
                           (t('blast_save'), self._save_blast_results)]:
            b = QPushButton(text); b.clicked.connect(slot); btn_row.addWidget(b)
            # Add tooltips
            if 'copy_hit' in slot.__name__:
                b.setToolTip(t('tip_blast_copy_hit'))
            elif 'copy_all' in slot.__name__:
                b.setToolTip(t('tip_blast_copy_all'))
            elif 'save' in slot.__name__:
                b.setToolTip(t('tip_blast_save'))
        btn_row.addStretch()
        layout.addLayout(btn_row)

        self._blast_results_text = QTextBrowser()
        self._blast_results_text.setFont(QFont('Courier New', 9))
        self._blast_results_text.setReadOnly(True)
        self._blast_results_text.setOpenExternalLinks(False)
        self._blast_results_text.setOpenLinks(False)
        self._blast_results_text.anchorClicked.connect(self._on_text_orf_click)
        layout.addWidget(self._blast_results_text)
        self._tabs.addTab(w, t('tab_blast_res'))

    def _create_dna_tab(self):
        w = QWidget()
        layout = QVBoxLayout(w)
        btn = QPushButton("📋 Copy DNA")
        btn.clicked.connect(lambda: self._copy_to_clipboard(self._dna_text.toPlainText()))
        layout.addWidget(btn)
        self._dna_text = QTextEdit()
        self._dna_text.setFont(QFont('Courier New', 10))
        self._dna_text.setReadOnly(True)
        layout.addWidget(self._dna_text)
        self._tabs.addTab(w, t('tab_dna'))

    def _create_protein_tab(self):
        w = QWidget()
        layout = QVBoxLayout(w)
        btn = QPushButton("📋 Copy Protein")
        btn.clicked.connect(lambda: self._copy_to_clipboard(self._protein_text.toPlainText()))
        layout.addWidget(btn)
        self._protein_text = QTextEdit()
        self._protein_text.setFont(QFont('Courier New', 10))
        self._protein_text.setReadOnly(True)
        layout.addWidget(self._protein_text)
        self._tabs.addTab(w, t('tab_protein'))

    def _create_domains_tab(self):
        w = QWidget()
        layout = QVBoxLayout(w)
        self._domains_text = QTextEdit()
        self._domains_text.setFont(QFont('Courier New', 9))
        self._domains_text.setReadOnly(True)
        layout.addWidget(self._domains_text)
        self._tabs.addTab(w, t('tab_domains'))

    def _create_neighborhood_tab(self):
        w = QWidget()
        layout = QVBoxLayout(w)
        ctrl = QHBoxLayout()
        ctrl.addWidget(QLabel(t('window_kb')))
        self._window_kb_spin = QSpinBox()
        self._window_kb_spin.setRange(1, 100)
        self._window_kb_spin.setValue(15)
        ctrl.addWidget(self._window_kb_spin)
        btn = QPushButton(t('analyze_btn'))
        btn.clicked.connect(self._analyze_neighborhood)
        ctrl.addWidget(btn)
        ctrl.addStretch()
        layout.addLayout(ctrl)
        self._neighborhood_text = QTextEdit()
        self._neighborhood_text.setFont(QFont('Courier New', 9))
        self._neighborhood_text.setReadOnly(True)
        layout.addWidget(self._neighborhood_text)
        self._tabs.addTab(w, t('tab_neighbors'))

    def _create_hmm_tab(self):
        w = QWidget()
        layout = QVBoxLayout(w)
        layout.setContentsMargins(2, 2, 2, 2)
        layout.setSpacing(2)

        # Buttons row
        btn_row = QHBoxLayout()
        btn_row.setSpacing(3)
        for text, slot in [("Add HMM Profile", self._add_hmm_profile),
                           ("Add Multiple Profiles", self._add_hmm_multi),
                           ("🗑️ Remove Selected", self._remove_hmm_profile),
                           ("🎨 Edit Color / Function", self._edit_hmm_profile),
                           ("🔍 Search All ORFs", self._hmm_search_all)]:
            b = QPushButton(text); b.clicked.connect(slot); btn_row.addWidget(b)
            # Add tooltips
            if text == "Add HMM Profile":
                b.setToolTip(t('tip_hmm_add'))
            elif text == "Add Multiple Profiles":
                b.setToolTip(t('tip_hmm_add_multi'))
            elif "Search All ORFs" in text:
                b.setToolTip(t('tip_hmm_search'))
        btn_row.addStretch()
        self._hmm_status_label = QLabel("No HMM profiles loaded")
        self._hmm_status_label.setStyleSheet("color: #666; font-style: italic;")
        btn_row.addWidget(self._hmm_status_label)
        layout.addLayout(btn_row)

        # Vertical splitter: profiles table (top) + results (bottom)
        splitter = QSplitter(Vertical)

        # ── Top: Loaded HMM Profiles table ──
        profiles_group = QGroupBox("📋 Loaded HMM Profiles")
        pg_layout = QVBoxLayout(profiles_group)
        pg_layout.setContentsMargins(4, 4, 4, 4)

        self._hmm_profile_table = QTableWidget()
        self._hmm_profile_table.setColumnCount(5)
        self._hmm_profile_table.setHorizontalHeaderLabels(['Name', 'File', 'Color', 'Function', 'Hits'])
        self._hmm_profile_table.horizontalHeader().setStretchLastSection(True)
        self._hmm_profile_table.setSelectionBehavior(SelectRows)
        self._hmm_profile_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers
                                                  if QT_VERSION == 6 else QAbstractItemView.NoEditTriggers)
        self._hmm_profile_table.setColumnWidth(0, 150)
        self._hmm_profile_table.setColumnWidth(1, 180)
        self._hmm_profile_table.setColumnWidth(2, 70)
        self._hmm_profile_table.setColumnWidth(3, 150)
        self._hmm_profile_table.setColumnWidth(4, 50)
        pg_layout.addWidget(self._hmm_profile_table)
        splitter.addWidget(profiles_group)

        # ── Bottom: HMM search results ──
        results_widget = QWidget()
        rl = QVBoxLayout(results_widget)
        rl.setContentsMargins(0, 0, 0, 0)
        self._hmm_text = QTextBrowser()
        self._hmm_text.setFont(QFont('Courier New', 9))
        self._hmm_text.setReadOnly(True)
        self._hmm_text.setOpenExternalLinks(False)
        self._hmm_text.setOpenLinks(False)
        self._hmm_text.anchorClicked.connect(self._on_text_orf_click)
        rl.addWidget(self._hmm_text)
        splitter.addWidget(results_widget)

        splitter.setSizes([180, 400])
        layout.addWidget(splitter)
        self._tabs.addTab(w, t('tab_hmm'))

    def _update_hmm_profile_table(self):
        """Refresh the HMM profiles table from self.hmm_profiles."""
        self._hmm_profile_table.setRowCount(0)
        for p in self.hmm_profiles:
            row = self._hmm_profile_table.rowCount()
            self._hmm_profile_table.insertRow(row)
            n_hits = len(p.get('hits', []))
            items = [p['name'], Path(p['file']).name, p['color'],
                     p.get('function', ''), str(n_hits)]
            for col, val in enumerate(items):
                item = QTableWidgetItem(val)
                if col == 2:  # Color column — show colored background
                    item.setBackground(QColor(p['color']))
                    item.setForeground(QColor('white'))
                self._hmm_profile_table.setItem(row, col, item)
        self._hmm_status_label.setText(
            f"{len(self.hmm_profiles)} HMM profile(s) loaded"
            if self.hmm_profiles else "No HMM profiles loaded")

    def _remove_hmm_profile(self):
        """Remove selected HMM profile."""
        rows = sorted(set(idx.row() for idx in self._hmm_profile_table.selectedIndexes()), reverse=True)
        for r in rows:
            if r < len(self.hmm_profiles):
                self.hmm_profiles.pop(r)
        self._update_hmm_profile_table()

    def _edit_hmm_profile(self):
        """Edit color and function of selected HMM profile."""
        rows = set(idx.row() for idx in self._hmm_profile_table.selectedIndexes())
        if not rows:
            QMessageBox.information(self, "HMM", "Select a profile in the table first.")
            return
        row = min(rows)
        if row >= len(self.hmm_profiles):
            return
        profile = self.hmm_profiles[row]

        dlg = QDialog(self)
        dlg.setWindowTitle(f"🎨 Edit: {profile['name']}")
        dlg.setFixedSize(380, 220)
        layout = QVBoxLayout(dlg)

        layout.addWidget(QLabel(f"<b>{profile['name']}</b> ({Path(profile['file']).name})"))

        form = QGridLayout()
        form.addWidget(QLabel("Function:"), 0, 0)
        func_edit = QLineEdit(profile.get('function', ''))
        form.addWidget(func_edit, 0, 1)

        form.addWidget(QLabel("Color:"), 1, 0)
        color_lay = QHBoxLayout()
        color_edit = QLineEdit(profile['color'])
        color_edit.setMaximumWidth(80)
        color_lay.addWidget(color_edit)
        # Color palette
        for c in ['#E53935','#1E88E5','#43A047','#FB8C00','#8E24AA',
                   '#00ACC1','#D81B60','#FFD600','#795548','#607D8B']:
            btn = QPushButton()
            btn.setFixedSize(22, 22)
            btn.setStyleSheet(f"background-color: {c}; border: 1px solid #888;")
            btn.clicked.connect(lambda checked, col=c: color_edit.setText(col))
            color_lay.addWidget(btn)
        form.addLayout(color_lay, 1, 1)
        layout.addLayout(form)

        btns = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
                                 if QT_VERSION == 6 else QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(dlg.accept); btns.rejected.connect(dlg.reject)
        layout.addWidget(btns)

        if dlg.exec() if QT_VERSION == 6 else dlg.exec_():
            profile['function'] = func_edit.text()
            profile['color'] = color_edit.text()
            self._update_hmm_profile_table()
            self._update_map()
            self._status.showMessage(f"✓ {profile['name']} updated: {color_edit.text()}")



    def _create_af3_tab(self):
        w = QWidget()
        self._af3_tab_widget = w   # stored so _switch_to_af3_tab() can find it
        layout = QVBoxLayout(w)
        layout.setContentsMargins(2, 2, 2, 2)
        layout.setSpacing(2)

        # Use a vertical splitter so user can drag to resize sections
        splitter = QSplitter(Vertical)

        # ═══ SECTION 1: ORF Selection (compact) ═══
        sel_widget = QWidget()
        sel_layout = QVBoxLayout(sel_widget)
        sel_layout.setContentsMargins(4, 4, 4, 2)
        sel_layout.setSpacing(2)

        sb = QHBoxLayout()
        sb.setSpacing(4)
        for text, slot in [(t('af3_add_sel'), self._af3_add_selected),
                           (t('af3_add_hmm'), self._af3_add_hmm_hits),
                           (t('af3_add_all'), self._af3_add_all_orfs),
                           (t('af3_remove'), self._af3_remove_orf),
                           (t('af3_clear_all'), self._af3_clear_all)]:
            b = QPushButton(text); b.clicked.connect(slot); sb.addWidget(b)
            # Add tooltips
            if 'add_selected' in slot.__name__:
                b.setToolTip(t('tip_af3_add_sel'))
            elif 'add_hmm' in slot.__name__:
                b.setToolTip(t('tip_af3_add_hmm'))
            elif 'add_all' in slot.__name__:
                b.setToolTip(t('tip_af3_add_all'))
                b.setStyleSheet("QPushButton { background-color: #2e7d32; color: white; font-weight: bold; "
                                "border-radius: 4px; padding: 3px 8px; }"
                                "QPushButton:hover { background-color: #388e3c; }")
            elif 'remove' in slot.__name__:
                b.setToolTip(t('tip_af3_remove'))
            elif 'clear_all' in slot.__name__:
                b.setToolTip(t('tip_af3_clear_all'))

        # ── "Predict Selected Pair" button ──────────────────────────────────
        # Enabled only when ≥ 2 rows are selected via Ctrl+click in the table.
        self._af3_predict_pair_btn = QPushButton("⚡ Predict Selected Pair")
        self._af3_predict_pair_btn.setEnabled(False)
        self._af3_predict_pair_btn.setToolTip(
            "Ctrl+click two or more ORFs in the list below, then click this button\n"
            "to immediately create an AF3 pairwise prediction job for those ORFs.\n"
            "This bypasses the Mode dropdown — pairs are created exactly as shown.")
        self._af3_predict_pair_btn.setStyleSheet(
            "QPushButton { background-color: #1565c0; color: white; font-weight: bold; "
            "border-radius: 4px; padding: 3px 8px; }"
            "QPushButton:hover { background-color: #1976d2; }"
            "QPushButton:disabled { background-color: #90a4ae; color: #eceff1; }")
        self._af3_predict_pair_btn.clicked.connect(self._af3_predict_selected_pair)
        sb.addWidget(self._af3_predict_pair_btn)

        sb.addStretch()
        self._af3_sel_count = QLabel("0 ORFs selected")
        self._af3_sel_count.setStyleSheet("font-weight: bold;")
        sb.addWidget(self._af3_sel_count)
        sel_layout.addLayout(sb)

        self._af3_sel_table = QTableWidget()
        self._af3_sel_table.setColumnCount(5)
        self._af3_sel_table.setHorizontalHeaderLabels(['ORF', 'Position', 'Size(aa)', 'HMM', 'Note'])
        self._af3_sel_table.horizontalHeader().setStretchLastSection(True)
        self._af3_sel_table.setSelectionBehavior(SelectRows)
        # ExtendedSelection: Ctrl+click = sparse multi-select; Shift+click = range
        self._af3_sel_table.setSelectionMode(
            QAbstractItemView.SelectionMode.ExtendedSelection if QT_VERSION == 6
            else QAbstractItemView.ExtendedSelection)
        self._af3_sel_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers
                                             if QT_VERSION == 6 else QAbstractItemView.NoEditTriggers)
        # Left-click → center genome map + select in main ORF table
        self._af3_sel_table.cellClicked.connect(self._af3_sel_table_click)
        # Right-click context menu
        self._af3_sel_table.setContextMenuPolicy(
            Qt.ContextMenuPolicy.CustomContextMenu if QT_VERSION == 6 else Qt.CustomContextMenu)
        self._af3_sel_table.customContextMenuRequested.connect(self._af3_sel_table_right_click)
        # Enable/disable Predict Pair button when selection changes
        self._af3_sel_table.itemSelectionChanged.connect(self._af3_update_predict_pair_btn)
        sel_layout.addWidget(self._af3_sel_table)
        splitter.addWidget(sel_widget)

        # ═══ SECTION 2: Job Generation + Jobs Table ═══
        job_widget = QWidget()
        jl = QVBoxLayout(job_widget)
        jl.setContentsMargins(4, 2, 4, 2)
        jl.setSpacing(2)

        # ══ Job-generation control area — 2 rows instead of one crowded line ══
        ctrl_frame = QFrame()
        ctrl_frame.setStyleSheet(
            "QFrame { background: #f5f5f5; border: 0.5px solid #ddd; "
            "border-radius: 4px; }")
        ctrl_vbox = QVBoxLayout(ctrl_frame)
        ctrl_vbox.setContentsMargins(6, 4, 6, 4)
        ctrl_vbox.setSpacing(4)

        # ── ROW 1: Neighbors  |  Mode combo  |  Homodimer checkbox ───────────
        row1 = QHBoxLayout()
        row1.setSpacing(6)

        row1.addWidget(QLabel("Neighbors:"))
        self._af3_nb_spin = QSpinBox()
        self._af3_nb_spin.setRange(1, 15)
        self._af3_nb_spin.setValue(5)
        self._af3_nb_spin.setMaximumWidth(50)
        self._af3_nb_spin.setToolTip(
            "Number of genomic neighbors to include in Pairs / Neighbors Interactome modes.")
        row1.addWidget(self._af3_nb_spin)

        row1.addWidget(QLabel("Mode:"))
        self._af3_mode_combo = QComboBox()
        # ── Mode list ────────────────────────────────────────────────────────
        # v2.0: removed 'Pairs + Homodimers' and 'Trimers' (confusing combined
        # modes). Homodimer is now an independent checkbox below the dropdown.
        # Restored 'All vs All (Selected ORFs)', 'HMM Hits vs Each Other',
        # 'Hit vs All Selected'. Added 'Selected vs Selected (Ctrl+click)'.
        _AF3_MODES = [
            "Pairs (Hit vs Neighbor)",
            "Selected vs Selected (Ctrl+click)",
            "All vs All (Selected ORFs)",
            "HMM Hits vs Each Other",
            "Hit vs All Selected",
            "Neighbors Interactome",
            "Genomic Interactome (Selected ORFs vs All ORFs)",
        ]
        self._af3_mode_combo.addItems(_AF3_MODES)
        self._af3_mode_combo.setMinimumWidth(230)
        self._af3_mode_combo.setSizePolicy(
            QSizePolicy.Policy.Expanding if QT_VERSION == 6 else QSizePolicy.Expanding,
            QSizePolicy.Policy.Fixed    if QT_VERSION == 6 else QSizePolicy.Fixed)
        row1.addWidget(self._af3_mode_combo, stretch=1)

        # ── Homodimer checkbox on the same row, clearly separated ─────────────
        _sep = QFrame()
        _sep.setFrameShape(QFrame.Shape.VLine if QT_VERSION == 6 else QFrame.VLine)
        _sep.setStyleSheet("color: #bbb;")
        row1.addWidget(_sep)

        self._af3_homodimer_cb = QCheckBox("＋ Homodimer for each ORF")
        self._af3_homodimer_cb.setToolTip(
            "When checked, a self-vs-self (homodimer) prediction job is added\n"
            "for every Hit ORF, regardless of the mode chosen above.")
        row1.addWidget(self._af3_homodimer_cb)
        ctrl_vbox.addLayout(row1)

        # ── Per-mode descriptions shown in the yellow info label ─────────────
        self._AF3_MODE_DESCS = {
            "Pairs (Hit vs Neighbor)":
                "🔵  <b>Pairs (Hit vs Neighbor)</b><br>"
                "Each selected Hit ORF is paired with each of its <i>N</i> nearest genomic neighbors.<br>"
                "Classic co-localization screen — best for operon-like gene clusters where physically<br>"
                "adjacent proteins are likely to interact.",

            "Selected vs Selected (Ctrl+click)":
                "🟦  <b>Selected vs Selected</b><br>"
                "Ctrl+click two or more ORFs in the selection list above, then click <b>Generate</b>.<br>"
                "Every pairwise combination among the highlighted rows is created — no neighbor window needed.<br>"
                "Ideal for testing specific hypotheses, e.g. ORF2897 vs ORF2596.",

            "All vs All (Selected ORFs)":
                "🟣  <b>All vs All (Selected ORFs)</b><br>"
                "Every pairwise combination among <i>all</i> ORFs in the selection list is generated.<br>"
                "Scales as N×(N-1)/2 — use with caution for large lists. Best after HMM filtering.",

            "HMM Hits vs Each Other":
                "🟤  <b>HMM Hits vs Each Other</b><br>"
                "Pairs every ORF with an HMM hit against every other ORF with an HMM hit in the selection list.<br>"
                "Useful to screen for interactions within a functional family (e.g. all T4SS components).",

            "Hit vs All Selected":
                "🟡  <b>Hit vs All Selected</b><br>"
                "The topmost ORF in the selection list is used as the query and tested against<br>"
                "all other ORFs in the list. Use when you have one anchor protein of known function.",

            "Neighbors Interactome":
                "🟢  <b>Neighbors Interactome</b><br>"
                "Sliding-window pairwise screen over the <b>entire genome</b>: every ORF <i>i</i> is paired<br>"
                "with every ORF within <i>±N</i> genomic positions of it. Symmetric pairs (i↔j) are<br>"
                "deduplicated automatically — the result is the unique set { (i,j) : 1 ≤ j−i ≤ N }.<br>"
                "For a genome of <i>n</i> ORFs and window <i>N</i>, that's <b>N·n − N(N+1)/2</b> jobs.<br>"
                "Does not use the 'Selected ORFs' list — always scans the full genome.",

            "Genomic Interactome (Selected ORFs vs All ORFs)":
                "🔴  <b>Genomic Interactome (Selected ORFs vs All ORFs)</b><br>"
                "Every <i>selected</i> ORF is tested against <b>every ORF in the full genome table</b>.<br>"
                "Symmetric duplicates (A↔B) are collapsed automatically.<br>"
                "Maximum coverage — high job count. Use size / HMM filters to manage scale.<br>"
                "Jobs &gt; 5 000 will prompt for confirmation before generating.",
        }

        # (mode combo and homodimer checkbox are in row1 above)

        # ── Yellow description label (updates on mode change) ────────────────
        self._af3_mode_desc_label = QLabel()
        self._af3_mode_desc_label.setWordWrap(True)
        self._af3_mode_desc_label.setTextFormat(
            Qt.TextFormat.RichText if QT_VERSION == 6 else Qt.RichText)
        self._af3_mode_desc_label.setStyleSheet(
            "background:#fffde7; color:#333; border:1px solid #f9a825;"
            "border-radius:4px; padding:5px 8px; font-size:11px;")
        self._af3_mode_desc_label.setMinimumHeight(62)

        def _update_mode_desc(idx):
            txt = self._AF3_MODE_DESCS.get(
                self._af3_mode_combo.currentText(), "")
            self._af3_mode_desc_label.setText(txt)

        self._af3_mode_combo.currentIndexChanged.connect(_update_mode_desc)
        _update_mode_desc(0)   # initialise with first item

        # ── Live job-count preview label ──────────────────────────────────────
        self._af3_job_preview_lbl = QLabel("— jobs estimated")
        self._af3_job_preview_lbl.setStyleSheet(
            "color: #1565c0; font-size: 11px; font-style: italic;")
        # (added to jl after ctrl_frame — see below)

        def _update_job_preview():
            n_sel = self._af3_sel_table.rowCount()
            n_nb  = self._af3_nb_spin.value()
            mode  = self._af3_mode_combo.currentText()
            n_genome = len(self.orfs) if hasattr(self, 'orfs') and self.orfs else 0
            if mode.startswith("Pairs (Hit vs Neighbor)"):
                est = n_sel * 2 * n_nb
            elif mode.startswith("Selected vs Selected"):
                rows_sel = len(set(idx.row() for idx in self._af3_sel_table.selectedIndexes()))
                est = rows_sel * (rows_sel - 1) // 2
            elif mode.startswith("All vs All"):
                est = n_sel * (n_sel - 1) // 2
            elif mode.startswith("HMM Hits"):
                n_hmm = sum(1 for r in range(n_sel) if self._af3_sel_table.item(r, 3) and
                            self._af3_sel_table.item(r, 3).text() not in ('-', ''))
                est = n_hmm * (n_hmm - 1) // 2
            elif mode.startswith("Hit vs All"):
                est = max(0, n_sel - 1)
            elif mode.startswith("Neighbors Interactome"):
                est = max(0, n_nb * n_genome - n_nb * (n_nb + 1) // 2)
            elif mode.startswith("Genomic Interactome"):
                est = n_sel * (n_genome - 1)
            else:
                est = 0
            if self._af3_homodimer_cb.isChecked() and not mode.startswith(("Neighbors", "Genomic")):
                est += n_sel
            color = "#c62828" if est > 5000 else "#1565c0" if est > 500 else "#2e7d32"
            self._af3_job_preview_lbl.setStyleSheet(
                f"color: {color}; font-size: 11px; font-style: italic;")
            self._af3_job_preview_lbl.setText(
                f"≈ {est:,} jobs will be generated")

        self._af3_mode_combo.currentIndexChanged.connect(_update_job_preview)
        self._af3_nb_spin.valueChanged.connect(_update_job_preview)
        self._af3_homodimer_cb.stateChanged.connect(_update_job_preview)
        self._af3_sel_table.itemSelectionChanged.connect(_update_job_preview)
        _update_job_preview()   # populate label immediately on tab creation

        # ── ROW 2: action buttons ─────────────────────────────────────────
        row2 = QHBoxLayout()
        row2.setSpacing(4)
        for text, slot in [(t('af3_generate'), self._af3_generate_jobs),
                           (t('af3_export_cf'), self._af3_export_colabfold),
                           ("Ranking", self._af3_show_ranking),
                           ("Clear Jobs", self._af3_clear_jobs)]:
            b = QPushButton(text); b.clicked.connect(slot)
            if 'generate' in slot.__name__:
                b.setToolTip(t('tip_af3_generate'))
                b.setStyleSheet(
                    "QPushButton{background:#2e7d32;color:white;font-weight:bold;"
                    "border-radius:4px;padding:3px 10px;}"
                    "QPushButton:hover{background:#388e3c;}")
            elif 'export_colabfold' in slot.__name__:
                b.setToolTip(t('tip_af3_export_cf'))
            elif 'ranking' in slot.__name__:
                b.setToolTip(t('tip_af3_ranking'))
            elif 'clear_jobs' in slot.__name__:
                b.setToolTip(t('tip_af3_clear_jobs'))
                b.setStyleSheet(
                    "QPushButton{color:#c62828;border:1px solid #c62828;"
                    "border-radius:4px;padding:3px 8px;}"
                    "QPushButton:hover{background:#ffebee;}")
            row2.addWidget(b)
        # AF3 JSON export menu button (individual vs batch)
        export_btn = QPushButton(t('af3_export_json'))
        export_menu = QMenu(export_btn)
        export_menu.addAction(t('af3_export_json_single'), self._af3_export_json)
        export_menu.addAction(t('af3_export_json_batch'), self._af3_export_json_batch)
        export_menu.addSeparator()
        _slurm_act = export_menu.addAction(t('af3_export_slurm_array'), self._af3_export_slurm_array)
        _slurm_act.setToolTip('Export JSONs in numbered batches + one SLURM array script. '
                              'Submit with: sbatch run_array.sh — one command, no OOM.')
        export_btn.setMenu(export_menu)
        row2.addWidget(export_btn)
        row2.addStretch()
        ctrl_vbox.addLayout(row2)
        jl.addWidget(ctrl_frame)
        jl.addWidget(self._af3_job_preview_lbl)
        jl.addWidget(self._af3_mode_desc_label)

        self._af3_jobs_table = QTableWidget()
        self._af3_jobs_table.setColumnCount(8)
        self._af3_jobs_table.setHorizontalHeaderLabels(
            ['Job', 'Hit', 'Partner', 'Residues', 'Status', 'ipTM', 'PAEinter', 'Confidence'])
        self._af3_jobs_table.horizontalHeader().setStretchLastSection(True)
        self._af3_jobs_table.setSelectionBehavior(SelectRows)
        self._af3_jobs_table.setSelectionMode(
            QAbstractItemView.SelectionMode.ExtendedSelection if QT_VERSION == 6
            else QAbstractItemView.ExtendedSelection)
        self._af3_jobs_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers
                                              if QT_VERSION == 6 else QAbstractItemView.NoEditTriggers)
        self._af3_jobs_table.setContextMenuPolicy(
            Qt.ContextMenuPolicy.CustomContextMenu if QT_VERSION == 6 else Qt.CustomContextMenu)
        self._af3_jobs_table.customContextMenuRequested.connect(self._af3_jobs_right_click)
        # Left-click on a job row → center genome map on Hit ORF
        self._af3_jobs_table.cellClicked.connect(self._af3_jobs_table_click)
        # NOTE: _af3_jobs_table is added to jl AFTER the custom section below
        # ── Custom complex builder — below mode desc, above jobs table ──────

        # Thin separator line before custom section
        _custom_sep = QFrame()
        _custom_sep.setFrameShape(QFrame.Shape.HLine if QT_VERSION == 6 else QFrame.HLine)
        _custom_sep.setStyleSheet("color: #ddd;")
        jl.addWidget(_custom_sep)

        # Custom header row
        custom_header = QHBoxLayout()
        custom_header.setSpacing(4)
        custom_header.addWidget(QLabel("⚡ Custom:"))
        custom_header.addWidget(QLabel("Subunits:"))
        self._custom_n_subunits = QSpinBox()
        self._custom_n_subunits.setRange(1, 11)
        self._custom_n_subunits.setValue(2)
        self._custom_n_subunits.setMaximumWidth(52)
        self._custom_n_subunits.setToolTip(
            "Number of subunits (chains) in the complex.\n"
            "Each subunit is assigned a chain letter (A, B, C … K).\n"
            "Maximum: 11 subunits.")
        self._custom_n_subunits.valueChanged.connect(self._af3_rebuild_custom_rows)
        custom_header.addWidget(self._custom_n_subunits)

        btn_add_custom = QPushButton("➕ Add")
        btn_add_custom.clicked.connect(self._af3_add_custom_job)
        btn_add_custom.setToolTip(t('tip_af3_add_custom'))
        custom_header.addWidget(btn_add_custom)
        custom_header.addStretch()
        jl.addLayout(custom_header)

        # Yellow explanatory label
        _custom_info = QLabel(
            "⚡ <b>Custom complex builder</b> — Define any multi-chain assembly manually.<br>"
            "Set the number of <b>Subunits</b> (A, B, C…) and type each chain's ORF name "
            "(e.g. <tt>ORF42</tt>). Set <b>n=</b> for stoichiometry "
            "(e.g. n=2 means two copies of that chain — homodimer, trimer, etc.). "
            "Click <b>➕ Add</b> to append the job to the list above.<br>"
            "This overrides all Mode settings and lets you specify any complex exactly."
        )
        _custom_info.setWordWrap(True)
        _custom_info.setTextFormat(
            Qt.TextFormat.RichText if QT_VERSION == 6 else Qt.RichText)
        _custom_info.setStyleSheet(
            "background:#fffde7; color:#555; border:1px solid #f9a825;"
            "border-radius:4px; padding:5px 8px; font-size:11px;")
        jl.addWidget(_custom_info)

        # Scroll area for dynamic per-subunit rows (A: ORF n=1, B: ORF n=1, …)
        self._custom_scroll = QScrollArea()
        self._custom_scroll.setWidgetResizable(True)
        self._custom_scroll.setFrameShape(QFrame.Shape.NoFrame
                                          if QT_VERSION == 6 else QFrame.NoFrame)
        self._custom_scroll.setMaximumHeight(110)
        self._custom_scroll.setMinimumHeight(55)
        self._custom_rows_widget = QWidget()
        self._custom_rows_layout = QVBoxLayout(self._custom_rows_widget)
        self._custom_rows_layout.setContentsMargins(0, 0, 0, 0)
        self._custom_rows_layout.setSpacing(2)
        self._custom_scroll.setWidget(self._custom_rows_widget)
        jl.addWidget(self._custom_scroll)

        # Internal list rebuilt by _af3_rebuild_custom_rows
        self._custom_subunit_rows: list = []
        self._af3_rebuild_custom_rows(2)

        # Thin separator before jobs table
        _jobs_sep = QFrame()
        _jobs_sep.setFrameShape(QFrame.Shape.HLine if QT_VERSION == 6 else QFrame.HLine)
        _jobs_sep.setStyleSheet("color: #ddd;")
        jl.addWidget(_jobs_sep)

        # Jobs table — shown below custom builder
        jl.addWidget(self._af3_jobs_table)

        splitter.addWidget(job_widget)

        # 2-pane splitter: selection list (top) | everything else (bottom)
        splitter.setSizes([180, 600])
        layout.addWidget(splitter)
        self._tabs.addTab(w, t('tab_af3'))

    # ═══════════════════════════════════════════════════════════
    # CORE ANALYSIS METHODS
    # ═══════════════════════════════════════════════════════════

    def _run_worker(self, fn, on_done, on_error=None, *args, **kwargs):
        """Launch a worker thread for fn, connecting result/error signals."""
        worker = AnalysisWorker(fn, *args, **kwargs)
        worker.finished.connect(on_done)
        worker.error.connect(on_error or (lambda e: QMessageBox.critical(self, "Error", e)))
        self._workers.append(worker)
        worker.finished.connect(lambda: self._workers.remove(worker) if worker in self._workers else None)
        worker.start()

    def analyze_orfs(self):
        if not self.dna_sequence:
            QMessageBox.warning(self, "Warning", "Load a FASTA file first!")
            return
        self._status.showMessage("⏳ Analyzing ORFs...")
        sc = set()
        if self._cb_atg.isChecked(): sc.add('ATG')
        if self._cb_gtg.isChecked(): sc.add('GTG')
        if self._cb_ttg.isChecked(): sc.add('TTG')
        min_aa = self._min_length_spin.value()

        def work():
            return self.analyzer.find_orfs(self.dna_sequence, min_aa=min_aa, start_codons=sc)

        def done(orfs):
            self.orfs = orfs
            self.filtered_orfs = orfs.copy()
            self._update_orfs_list()
            self._update_info()
            self._update_map()
            self._status.showMessage(f"✓ {len(orfs)} ORFs (6-frame)")

        self._run_worker(work, done)


    def analyze_orfs_pyrodigal(self):
        """Run Pyrodigal using current _pyro_params (set via Parameters menu)."""
        if not self.dna_sequence:
            QMessageBox.warning(self, "Warning", "Load a FASTA file first!")
            return
        if not PYRODIGAL_AVAILABLE:
            QMessageBox.critical(self, "Pyrodigal", t('pyrodigal_not_avail'))
            return
        p = self._pyro_params
        mode_str = 'meta' if p.get('meta', True) else 'single'
        self._status.showMessage(
            f"Running Pyrodigal — mode={mode_str}, "
            f"table={p.get('translation_table',11)}, "
            f"min={p.get('min_aa',30)}aa...")

        def work():
            orfs = self.analyzer.find_orfs_pyrodigal(
                self.dna_sequence,
                meta=p.get('meta', True),
                min_aa=p.get('min_aa', 30),
                closed_ends=p.get('closed', False),
                translation_table=p.get('translation_table', 11),
                mask=p.get('mask', False),
            )
            # Post-prediction start codon filter
            sf = p.get('start_filter', {'all': True})
            if not sf.get('all', True):
                allowed = set()
                if sf.get('ATG'): allowed.add('ATG')
                if sf.get('GTG'): allowed.add('GTG')
                if sf.get('TTG'): allowed.add('TTG')
                if allowed:
                    orfs = [
                        o for o in orfs
                        if o.get('dna', '')[:3].upper() in allowed
                    ]
            return orfs

        def done(orfs):
            self.orfs = orfs
            self.filtered_orfs = orfs.copy()
            self._update_orfs_list()
            self._update_info()
            self._update_map()
            sf = p.get('start_filter', {'all': True})
            filter_str = ('all starts' if sf.get('all', True)
                else '+'.join(k for k in ('ATG','GTG','TTG') if sf.get(k)))
            self._status.showMessage(
                f"Pyrodigal: {len(orfs)} genes | "
                f"mode={mode_str}, table={p.get('translation_table',11)}, "
                f"min={p.get('min_aa',30)}aa, starts={filter_str}")

        self._run_worker(work, done)

    def analyze_orfs_hybrid(self):
        """Hybrid mode: Pyrodigal as primary caller, 6-frame scanner fills gaps."""
        if not self.dna_sequence:
            QMessageBox.warning(self, "Warning", "Load a FASTA file first!")
            return
        if not PYRODIGAL_AVAILABLE:
            QMessageBox.critical(self, "Hybrid mode — Pyrodigal required",
                t('pyrodigal_not_avail') + "\n\nHybrid mode requires Pyrodigal "
                "as the primary caller.\nUse 'Automatic' mode instead, or "
                "install Pyrodigal:\n  pip install pyrodigal")
            return

        # Collect Pyrodigal params
        p = self._pyro_params
        mode_str = 'meta' if p.get('meta', True) else 'single'

        # Collect 6-frame scanner params
        sc = set()
        if self._cb_atg.isChecked(): sc.add('ATG')
        if self._cb_gtg.isChecked(): sc.add('GTG')
        if self._cb_ttg.isChecked(): sc.add('TTG')
        min_aa = self._min_length_spin.value()

        self._status.showMessage(
            f"⏳ Hybrid mode — Pyrodigal (mode={mode_str}, "
            f"table={p.get('translation_table',11)}, min={p.get('min_aa',30)}aa) "
            f"+ gap-fill scanner (min={min_aa}aa, starts={','.join(sorted(sc)) or 'ATG'})…")

        def work():
            return self.analyzer.find_orfs_hybrid(
                self.dna_sequence,
                # 6-frame scanner params
                min_aa=min_aa,
                start_codons=sc,
                # Pyrodigal params
                pyro_meta=p.get('meta', True),
                pyro_min_aa=p.get('min_aa', 30),
                pyro_closed=p.get('closed', False),
                pyro_translation_table=p.get('translation_table', 11),
                pyro_mask=p.get('mask', False),
                pyro_start_filter=p.get('start_filter', {'all': True}),
            )

        def done(orfs):
            self.orfs = orfs
            self.filtered_orfs = orfs.copy()
            self._update_orfs_list()
            self._update_info()
            self._update_map()
            n_pyro = sum(1 for o in orfs if o.get('source') == 'pyrodigal')
            n_auto = sum(1 for o in orfs if o.get('source') == 'automatic')
            self._status.showMessage(
                f"✓ Hybrid: {len(orfs)} ORFs total — "
                f"{n_pyro} pyrodigal + {n_auto} gap-fill (automatic)")

        self._run_worker(work, done)

    def classify_all_domains(self):
        if not self.orfs:
            QMessageBox.warning(self, "Warning", "Run ORF analysis first!")
            return
        if not self.hmm_hits_all:
            QMessageBox.information(self, "Annotate Domains",
                "No HMM results available.\n\n"
                "1. Go to the HMM tab\n2. Add HMM profiles\n"
                "3. Click 'Search All ORFs'\n4. Then click 'Annotate HMM' again")
            return
        for orf in self.orfs:
            orf['domains'] = []
        td = 0
        for hit in self.hmm_hits_all:
            oi = hit.get('orf_index', -1)
            if 0 <= oi < len(self.orfs):
                profile_name = hit.get('profile_name', hit.get('hmm_name', '?'))
                profile_func = hit.get('profile_function', '')
                ali_from = hit.get('ali_from', 0); ali_to = hit.get('ali_to', 0)
                domain_entry = {
                    'domain': profile_name,
                    'description': profile_func or f"HMM: {profile_name}",
                    'system': profile_func or 'HMM hit', 'role': 'HMM',
                    'start': ali_from, 'end': ali_to,
                    'ali_region': f"{ali_from}-{ali_to}" if ali_from and ali_to else '',
                    'evalue': hit.get('evalue', 999), 'score': hit.get('score', 0),
                }
                existing = [d['domain'] for d in self.orfs[oi]['domains']]
                if profile_name not in existing:
                    self.orfs[oi]['domains'].append(domain_entry); td += 1
        self.filtered_orfs = self.orfs.copy()
        self._update_orfs_list()
        self._update_map()
        self._status.showMessage(f"✓ {td} HMM domains annotated")

    # ═══════════════════════════════════════════════════════════
    # BLAST
    # ═══════════════════════════════════════════════════════════

    def _build_blast_cmd_string(self):
        prog = self._blast_prog_combo.currentText() if hasattr(self, '_blast_prog_combo') else 'blastp'
        evalue = self._evalue_edit.text() if hasattr(self, '_evalue_edit') else '0.05'
        parts = [prog, '-query <query.fasta>', '-subject <orfs_db.fasta>',
                 '-outfmt "6 qseqid score pident positive length qstart qend sstart send evalue gaps"',
                 f'-evalue {evalue}', f'-matrix {self.blast_matrix}',
                 f'-word_size {self.blast_word_size}',
                 f'-gapopen {self.blast_gap_open}', f'-gapextend {self.blast_gap_ext}',
                 f'-max_target_seqs {self.blast_max_targets}',
                 f'-seg {"yes" if self.blast_low_complexity else "no"}']
        cmd = ' \\\n  '.join(parts)
        thresh = self._identity_spin.value() if hasattr(self, '_identity_spin') else 30
        cmd += f'\n\n# Post-filter: min identity ≥ {thresh}%'
        algo = self._algo_combo.currentText() if hasattr(self, '_algo_combo') else 'Auto'
        cmd += f'\n# Algorithm: {algo}'
        if not BACKENDS.get('blast+',{}).get('available'):
            cmd += '\n\n# ⚠️  NCBI BLAST+ not detected — will use built-in Python algorithm'
        return cmd

    def _refresh_blast_cmd_preview(self):
        if hasattr(self, '_blast_cmd_text'):
            self._blast_cmd_text.setPlainText(self._build_blast_cmd_string())

    def _run_custom_blast_cmd(self):
        self.run_blast()

    def run_blast(self):
        raw = self._blast_query_text.toPlainText().strip()
        if not raw:
            QMessageBox.warning(self, "BLAST", "Paste a protein sequence first!")
            return

        qp = self._parse_fasta_query(raw)
        if not qp or len(qp) < 5:
            QMessageBox.warning(self, "BLAST", "Invalid or too short sequence!")
            return

        if not self.orfs:
            QMessageBox.warning(self, "BLAST", "Run ORF analysis first!")
            return

        algo = self._algo_combo.currentText()

        try:
            evalue = float(self._evalue_edit.text() or "0.05")
        except ValueError:
            evalue = 0.05

        params = _BlastSearchParams(
            threshold=self._identity_spin.value(),
            gap_open=-self.blast_gap_open,
            gap_extend=-self.blast_gap_ext,
            evalue=evalue,
            word_size=self.blast_word_size,
            max_targets=self.blast_max_targets,
            matrix=self.blast_matrix,
            low_complexity=self.blast_low_complexity,
        )

        self._status.showMessage(
            f"BLAST: {len(qp)} aa vs {len(self.orfs)} ORFs..."
        )

        service = _BlastSearchService(self.analyzer)

        def work():
            return service.search(
                query_sequence=qp,
                orfs=self.orfs,
                algorithm=algo,
                params=params,
            )

        def done(result):
            self._show_blast_results(
                result.hits,
                result.query_sequence,
                result.algorithm_used,
            )

            if result.backend_error:
                self._status.showMessage(
                    f"{result.algorithm_used}: {len(result.hits)} hits "
                    f"(fallback used — {result.backend_error[:80]})"
                )
            else:
                self._status.showMessage(
                    f"{result.algorithm_used}: {len(result.hits)} hits"
                )

        self._run_worker(work, done)

    def _show_blast_results(self, hits, query_protein, algo_used):
        """Display BLAST results with colored HTML alignments."""
        self._blast_results_text.clear()
        if not hits:
            self._blast_results_text.setHtml("<h3>No hits found.</h3>")
            return

        html = []
        html.append("<pre style='font-family: Courier New, monospace; font-size: 9pt;'>")
        html.append(f"<b style='color:#1A237E;'>BLASTp — {algo_used}</b>")
        html.append(f"\n{'='*72}")
        html.append(f"\n<span style='color:#006064;'>Query: {len(query_protein)} aa  |  "
                     f"Database: {len(self.orfs)} ORFs  |  Hits: {len(hits)}</span>")
        html.append(f"\n{'='*72}\n")

        # Summary table
        html.append("\n<b style='color:#4A148C;'>SUMMARY TABLE</b>")
        html.append(f"\n<span style='color:#BDBDBD;'>{'-'*80}</span>")
        html.append(f"\n<b>{'#':<4} {'ORF':<9} {'Score':<7} {'E-value':<10} "
                     f"{'Ident%':<8} {'Pos%':<7} {'AlnLen':<7} {'Cov%':<7} Domains</b>")
        html.append(f"\n<span style='color:#BDBDBD;'>{'-'*80}</span>")

        for i, hit in enumerate(hits[:30]):
            oi = hit['orf_index']
            orf = self.orfs[oi] if oi < len(self.orfs) else {}
            doms = ', '.join(d['domain'] for d in orf.get('domains', [])) or '-'
            ev = hit.get('evalue', 999)
            ev_str = f"{ev:.1e}" if ev < 0.01 else f"{ev:.3f}"
            html.append(f"\n{i+1:<4} <a href='orf:{oi}' style='color:#0D47A1;font-weight:bold;'>ORF{oi+1}</a>   "
                         f"{hit['score']:<7} {ev_str:<10} {hit['identity']:>5.1f}%  "
                         f"{hit['positives']:>5.1f}% {hit['aln_length']:<7} "
                         f"{hit.get('coverage',0):>5.0f}%  {doms}")

        # Detailed alignments
        html.append(f"\n\n{'='*72}")
        html.append("\n<b style='color:#1A237E;'>DETAILED ALIGNMENTS (top 10)</b>")
        html.append(f"\n{'='*72}")
        html.append("\nLegend: "
                     "<span style='background:#2E7D32;color:white;'> Identity </span>  "
                     "<span style='background:#1565C0;color:white;'> Positive </span>  "
                     "<span style='background:#FFCDD2;color:#B71C1C;'> Gap </span>  "
                     "Mismatch\n")

        for i, hit in enumerate(hits[:10]):
            oi = hit['orf_index']
            orf = self.orfs[oi] if oi < len(self.orfs) else {}
            pl = len(orf.get('protein', '').rstrip('*'))
            doms = ', '.join(d['domain'] for d in orf.get('domains', [])) or 'hypothetical protein'
            ev = hit.get('evalue', 999)
            ev_str = f"{ev:.1e}" if ev < 0.01 else f"{ev:.3f}"

            html.append(f"\n<span style='color:#BDBDBD;'>{'─'*72}</span>")
            html.append(f"\n<b style='color:#4A148C;'>Hit #{i+1}: "
                         f"<a href='orf:{oi}' style='color:#0D47A1;'>ORF{oi+1}</a></b>")
            html.append(f"\n<span style='color:#006064;'>"
                         f"  Score: {hit['score']}  |  E-value: {ev_str}  |  "
                         f"Identities: {hit.get('identities_count',0)}/{hit['aln_length']} ({hit['identity']:.0f}%)  |  "
                         f"Positives: {hit.get('positives_count',0)}/{hit['aln_length']} ({hit['positives']:.0f}%)")
            html.append(f"\n  Subject: ORF{oi+1} ({pl} aa) | Frame: F{orf.get('frame','?')}{orf.get('strand','?')} | "
                         f"Pos: {orf.get('start',0):,}-{orf.get('end',0):,} | {doms}</span>")

            # Colored alignment blocks
            aq = hit.get('aln_query', '')
            am = hit.get('aln_midline', '')
            asb = hit.get('aln_subject', '')
            if not aq:
                continue

            bs = 50
            qp = hit.get('q_start', 1)
            sp = hit.get('s_start', 1)
            for b in range(0, len(aq), bs):
                qb = aq[b:b+bs]
                mb = am[b:b+bs] if am else ' ' * len(qb)
                sb = asb[b:b+bs]
                qe = qp + len(qb.replace('-', '')) - 1
                se = sp + len(sb.replace('-', '')) - 1

                # Query line with colors
                html.append(f"\n<span style='color:#1B5E20;font-weight:bold;'>  Query  {qp:>5}  </span>")
                html.append(self._color_aln_html(qb, mb))
                html.append(f"<span style='color:#1B5E20;font-weight:bold;'>  {qe}</span>")

                # Midline
                html.append("\n               ")
                html.append(self._color_mid_html(mb))

                # Subject line with colors
                html.append(f"\n<span style='color:#B71C1C;font-weight:bold;'>  Sbjct  {sp:>5}  </span>")
                html.append(self._color_aln_html(sb, mb))
                html.append(f"<span style='color:#B71C1C;font-weight:bold;'>  {se}</span>\n")

                qp = qe + 1
                sp = se + 1

        html.append("</pre>")
        self._blast_results_text.setHtml(''.join(html))
        self._tabs.setCurrentIndex(2)

    def _color_aln_html(self, seq, mid):
        """Color an alignment sequence based on midline match info."""
        parts = []
        for aa, m in zip(seq, mid):
            if aa == '-':
                parts.append(f"<span style='background:#FFCDD2;color:#B71C1C;'>{aa}</span>")
            elif m != ' ' and m != '-':
                if m == '+':
                    parts.append(f"<span style='background:#BBDEFB;color:#0D47A1;'>{aa}</span>")
                else:
                    parts.append(f"<span style='background:#C8E6C9;color:#1B5E20;'>{aa}</span>")
            else:
                parts.append(f"<span style='color:#666;'>{aa}</span>")
        return ''.join(parts)

    def _color_mid_html(self, mid):
        """Color the midline of an alignment."""
        parts = []
        for ch in mid:
            if ch == ' ':
                parts.append(' ')
            elif ch == '+':
                parts.append(f"<span style='color:#1565C0;'>{ch}</span>")
            else:
                parts.append(f"<span style='color:#2E7D32;font-weight:bold;'>{ch}</span>")
        return ''.join(parts)

    def _parse_fasta_query(self, raw):
        lines = raw.strip().split('\n')
        seq_lines = [l.strip() for l in lines if l.strip() and not l.strip().startswith('>')]
        seq = ''.join(seq_lines).upper()
        seq = re.sub(r'[^A-Z]', '', seq)
        # Heuristic: if >50% ATCG, it's probably DNA, not protein
        atcg = sum(1 for c in seq if c in 'ATCGN')
        if len(seq) > 0 and atcg / len(seq) > 0.85:
            return None  # DNA, not protein
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        seq = ''.join(c for c in seq if c in valid_aa or c == 'X')
        return seq if len(seq) >= 5 else None

    # ═══════════════════════════════════════════════════════════
    # FILE I/O
    # ═══════════════════════════════════════════════════════════

    def load_fasta(self):
        f, _ = QFileDialog.getOpenFileName(
            self,
            "Open sequence file",
            "",
            "Sequence files (*.fasta *.fa *.fna *.faa *.gb *.gbk *.genbank *.dna);;All (*)",
        )
        if not f:
            return

        _ui_open_genome_file_into_window(self, f)

    def load_multi_fasta(self):
        self.load_fasta()

    def load_snapgene(self):
        f, _ = QFileDialog.getOpenFileName(self, "Open SnapGene",
            "", "SnapGene (*.dna);;All (*)")
        if not f: return
        result = parse_snapgene_dna(f)
        self.dna_sequence = result['sequence']
        self.genome_name = result['name']
        self.snapgene_features = result.get('features', [])
        self.snapgene_primers = result.get('primers', [])
        self._on_sequence_loaded(f)

    def load_genbank(self):
        f, _ = QFileDialog.getOpenFileName(self, "Open GenBank",
            "", "GenBank (*.gb *.gbk *.genbank);;All (*)")
        if not f: return
        result = parse_genbank(f)
        self.dna_sequence = result['sequence']
        self.genome_name = result['name']
        self.snapgene_features = result.get('features', [])
        self._on_sequence_loaded(f)

    def _on_sequence_loaded(self, filepath):
        self.current_fasta_path = filepath   # remember source for project save
        self._update_info()
        self._genome_map.set_data(len(self.dna_sequence), [], self.hmm_profiles, self.dna_sequence)
        self._status.showMessage(
            f"✓ Loaded: {self.genome_name} ({len(self.dna_sequence):,} bp)")

    def load_hmm(self):
        f, _ = QFileDialog.getOpenFileName(self, "Load HMM",
            "", "HMM (*.hmm);;All (*)")
        if f: self._add_hmm_file(f)

    def _add_hmm_file(self, filepath):
        name = Path(filepath).stem
        color = self.hmm_default_colors[len(self.hmm_profiles) % len(self.hmm_default_colors)]
        self.hmm_profiles.append({
            'file': filepath, 'name': name, 'color': color,
            'function': '', 'hits': []})
        self._update_hmm_profile_table()

    def save_fasta(self):
        if not self.orfs:
            QMessageBox.information(
                self,
                "Save FASTA",
                "No ORFs to export. Run ORF analysis first.",
            )
            return

        f, _ = QFileDialog.getSaveFileName(
            self,
            "Save FASTA",
            "",
            "FASTA (*.fasta)",
        )
        if not f:
            return

        try:
            _io_write_orf_protein_fasta(f, self.orfs)
        except OSError as exc:
            QMessageBox.critical(self, "Save FASTA", f"Could not write file:\n{exc}")
            return

        self._status.showMessage(f"✓ Saved {len(self.orfs)} ORFs")

    def export_html_report(self):
        f, _ = QFileDialog.getSaveFileName(
            self,
            "Export HTML report",
            "",
            "HTML report (*.html)",
        )
        if not f:
            return

        try:
            _io_write_basic_report(
                f,
                title="ppigFinder Report",
                genome_name=getattr(self, "genome_name", ""),
                genome_length=len(getattr(self, "dna_sequence", "") or ""),
                orfs=getattr(self, "orfs", []) or [],
                interaction_results=getattr(self, "af3_results", []) or [],
            )
        except Exception as exc:
            QMessageBox.critical(
                self,
                "Export HTML report",
                f"Could not export HTML report:\n{exc}",
            )
            return

        try:
            self._status.showMessage(f"✓ HTML report exported: {f}")
        except Exception:
            pass


    # ═══════════════════════════════════════════════════════════
    # ORF TABLE — COPY & EXPORT
    # ═══════════════════════════════════════════════════════════

    def _orf_table_headers(self):
        """Return current column header labels."""
        t = self._orf_table
        return [t.horizontalHeaderItem(c).text()
                for c in range(t.columnCount())]

    def _orf_table_rows_tsv(self, row_indices):
        """Return header + given rows as tab-separated string."""
        t = self._orf_table
        headers = self._orf_table_headers()
        lines = ['	'.join(headers)]
        for r in row_indices:
            cells = []
            for c in range(t.columnCount()):
                item = t.item(r, c)
                cells.append(item.text() if item else '')
            lines.append('\t'.join(cells))
        return '\n'.join(lines)

    def _orf_table_copy_selection(self):
        """Copy selected rows as TSV to clipboard."""
        rows = sorted(set(
            idx.row() for idx in self._orf_table.selectedIndexes()))
        if not rows:
            self._status.showMessage("No rows selected — use Ctrl/Shift+click to select rows")
            return
        text = self._orf_table_rows_tsv(rows)
        QApplication.clipboard().setText(text)
        self._status.showMessage(
            f"Copied {len(rows)} row(s) to clipboard (tab-separated — paste into Excel)")

    def _orf_table_copy_rows(self, all_rows=False):
        """Copy all visible rows as TSV to clipboard."""
        n = self._orf_table.rowCount()
        if n == 0:
            self._status.showMessage("No ORFs to copy"); return
        text = self._orf_table_rows_tsv(range(n))
        QApplication.clipboard().setText(text)
        self._status.showMessage(
            f"Copied {n} rows to clipboard (tab-separated — paste into Excel)")

    def _export_orf_table(self, fmt='tsv', include_seqs=False, annotated_only=False):
        """Export ORF table to TSV or tab-formatted TXT.

        Parameters
        ----------
        fmt : 'tsv' | 'txt'
        include_seqs : bool — add DNA and Protein columns
        annotated_only : bool — only ORFs with observation/function/gene_name
        """
        if not self.orfs:
            QMessageBox.warning(self, "Export", "No ORFs to export."); return

        suffix = '.tsv' if fmt == 'tsv' else '.txt'
        filter_str = "TSV (*.tsv);;All (*)" if fmt == 'tsv' else "TXT (*.txt);;All (*)"
        genome_safe = re.sub(r'[^\w\-]', '_', self.genome_name or 'orfs')
        seq_tag = '_full' if include_seqs else ''
        ann_tag = '_annotated' if annotated_only else ''
        default = f"{genome_safe}_orfs{ann_tag}{seq_tag}{suffix}"

        path, _ = QFileDialog.getSaveFileName(
            self, "Export ORF table", default, filter_str)
        if not path: return

        # Build ORF list to export
        orfs_to_export = []
        for i, orf in enumerate(self.orfs):
            if annotated_only:
                if not any(orf.get(k) for k in
                           ('observation', 'putative_function', 'gene_name')):
                    continue
            orfs_to_export.append((i, orf))

        # Column definitions
        base_cols = [
            'ORF_ID', 'Frame', 'Strand', 'Start', 'End', 'Size_aa', 'GC_pct',
            'Source', 'Score', 'HMM_domains', 'HMM_evalue',
            'Gene_name', 'Putative_function', 'Observation', 'Notes',
            'Custom_color', 'RBS_motif', 'Partial',
            'AF3_done', 'AF3_partner', 'ipTM', 'PAE_inter', 'Contact_region',
        ]
        seq_cols = ['DNA_sequence', 'Protein_sequence'] if include_seqs else []
        all_cols = base_cols + seq_cols

        try:
            with open(path, 'w', newline='', encoding='utf-8') as fh:
                w = csv.writer(fh, delimiter='	')
                w.writerow(all_cols)

                for i, orf in orfs_to_export:
                    # HMM info
                    domains = ';'.join(d.get('domain','?') for d in orf.get('domains',[]))
                    evalues = ';'.join(str(d.get('evalue','?')) for d in orf.get('domains',[]))

                    # AF3 info from analysis results
                    orf_label = f"ORF{i+1}"
                    af3_done = partner = iptm_s = pae_s = contact_s = '-'
                    for res in getattr(self, '_af3_analysis_results', []):
                        if orf_label in res.get('orf_names', []):
                            af3_done  = 'yes'
                            partner   = res.get('partner_name', '-')
                            iptm_v    = res.get('iptm')
                            iptm_s    = f"{iptm_v:.4f}" if iptm_v is not None else '-'
                            pae_v     = res.get('pae_inter')
                            pae_s     = f"{pae_v:.2f}" if pae_v is not None else '-'
                            contact_s = res.get('contact_region', '-')
                            break

                    row = [
                        orf_label,
                        orf.get('frame', ''),
                        orf.get('strand', ''),
                        orf.get('start', ''),
                        orf.get('end', ''),
                        len(orf.get('protein', '').rstrip('*')),
                        f"{orf.get('gc', 0):.2f}",
                        orf.get('source', '6frame'),
                        f"{orf.get('candidate_score', 0):.4f}",
                        domains or '-',
                        evalues or '-',
                        orf.get('gene_name', ''),
                        orf.get('putative_function', ''),
                        orf.get('observation', ''),
                        orf.get('notes', ''),
                        orf.get('custom_color', ''),
                        orf.get('rbs_motif', ''),
                        'yes' if orf.get('partial') else 'no',
                        af3_done, partner, iptm_s, pae_s, contact_s,
                    ]
                    if include_seqs:
                        row += [
                            orf.get('dna', ''),
                            orf.get('protein', '').rstrip('*'),
                        ]
                    w.writerow(row)

            n = len(orfs_to_export)
            kb = Path(path).stat().st_size / 1024
            msg = (f"Exported {n} ORFs to {Path(path).name} "
                   f"({len(all_cols)} columns, {kb:.0f} KB)")
            self._status.showMessage(f"Export: {msg}")
            if hasattr(self, '_orf_export_info'):
                self._orf_export_info.setText(msg)

        except OSError as e:
            QMessageBox.critical(self, "Export error", str(e))

    def _export_orf_fasta(self, aa=True):
        """Export ORF sequences as FASTA with full annotation headers."""
        if not self.orfs:
            QMessageBox.warning(self, "Export", "No ORFs to export."); return

        genome_safe = re.sub(r'[^\w\-]', '_', self.genome_name or 'orfs')
        mol = 'protein' if aa else 'dna'
        default = f"{genome_safe}_{mol}.fasta"
        path, _ = QFileDialog.getSaveFileName(
            self, f"Export {'protein' if aa else 'DNA'} FASTA",
            default, "FASTA (*.fasta *.fa);;All (*)")
        if not path: return

        n_written = 0
        try:
            with open(path, 'w', encoding='utf-8') as fh:
                for i, orf in enumerate(self.orfs):
                    seq = orf.get('protein', '').rstrip('*') if aa else orf.get('dna', '')
                    if not seq: continue

                    # Rich FASTA header
                    gene   = orf.get('gene_name', '')
                    func   = orf.get('putative_function', '')
                    obs    = orf.get('observation', '')
                    hmms   = ';'.join(d.get('domain','') for d in orf.get('domains',[]))
                    rbs    = orf.get('rbs_motif', '')
                    parts  = [
                        f"ORF{i+1}",
                        f"loc={orf.get('start',0)}-{orf.get('end',0)}{orf.get('strand','')}",
                        f"frame={orf.get('frame','')}",
                        f"len={len(seq)}{'aa' if aa else 'bp'}",
                        f"gc={orf.get('gc',0):.1f}",
                        f"src={orf.get('source','6frame')}",
                    ]
                    if gene:  parts.append(f"gene={gene}")
                    if func:  parts.append(f"function={func}")
                    if obs:   parts.append(f"obs={obs}")
                    if hmms:  parts.append(f"hmm={hmms}")
                    if rbs:   parts.append(f"rbs={rbs}")
                    if orf.get('partial'): parts.append("partial=yes")

                    header = ' '.join(parts)
                    fh.write(f">{header}\n")
                    # Write 60 chars per line
                    for j in range(0, len(seq), 60):
                        fh.write(seq[j:j+60] + '\n')
                    n_written += 1

            kb = Path(path).stat().st_size / 1024
            self._status.showMessage(
                f"FASTA exported: {n_written} sequences -> {Path(path).name} ({kb:.0f} KB)")
        except OSError as e:
            QMessageBox.critical(self, "Export error", str(e))


    def save_report_tsv(self):
        if not self.orfs:
            QMessageBox.information(self, "Save Report",
                "No ORFs to export. Run ORF analysis first.")
            return
        f, _ = QFileDialog.getSaveFileName(self, "Save Report", "", "TSV (*.tsv)")
        if not f: return
        try:
            with open(f, 'w', newline='', encoding='utf-8') as fh:
                w = csv.writer(fh, delimiter='\t')
                w.writerow(['Rank','ORF','Start','End','Frame','Strand','Len_aa','GC',
                            'Score','Source','Domains','RBS','Protein50'])
                for r, orf in enumerate(sorted(self.orfs,
                        key=lambda x: x.get('candidate_score',0), reverse=True), 1):
                    idx = self.orfs.index(orf)+1
                    ds = ';'.join(d['domain'] for d in orf.get('domains',[]))
                    w.writerow([r, f'ORF{idx}', orf['start'], orf['end'], orf['frame'],
                        orf['strand'], len(orf['protein'].rstrip('*')), f"{orf['gc']:.1f}",
                        f"{orf.get('candidate_score',0):.3f}", orf.get('source','6frame'),
                        ds or '-', orf.get('rbs_motif','') or '-', orf['protein'][:50]])
        except OSError as e:
            QMessageBox.critical(self, "Save Report", f"Could not write file:\n{e}")
            return
        self._status.showMessage(f"✓ Report saved: {Path(f).name}")

    # ───────────────────────────────────────────────────────────
    # PROJECT SAVE / LOAD  (directory-based, v34)
    # ───────────────────────────────────────────────────────────
    PROJECT_MANIFEST = "project.json"
    PROJECT_VERSION  = "v2.00"

    # ─────────────────────────────────────────────────────────
    # PROJECT SAVE / LOAD
    # ─────────────────────────────────────────────────────────

    def _build_manifest(self, proj_dir=None):
        """Build project manifest dict.
        proj_dir: Path — if given, HMM/genome file paths are relative to it.
        Heavy AF3 arrays (pae_matrix, contact_probs, plddt_arr, token_res_ids)
        are NEVER serialised — they are always reloaded on demand from job_dir.
        """

        # HMM manifest
        hmm_manifest = []
        for p in self.hmm_profiles:
            src = Path(p.get('file', ''))
            rel = f"hmm/{src.name}" if (proj_dir and src.is_file()) else ''
            hmm_manifest.append({
                'name': p['name'], 'file': rel,
                'color': p['color'], 'function': p.get('function', ''),
                'hits': p.get('hits', []),
            })

        # AF3 result score JSONs list
        results_manifest = []
        if proj_dir:
            for job in self.af3_jobs:
                if job.get('iptm') is not None:
                    jfn = re.sub(r'[^\w\-.]', '_', job['name']) + '.json'
                    results_manifest.append(f"results/{jfn}")

        # AF3 Analysis results
        af3_analysis_ser = []
        # Heavy keys that are NEVER saved in the manifest because they can be
        # reloaded on demand from job_dir.  Excluding them reduces a typical
        # 100-job project from ~200 MB to <1 MB and makes json.dump instant.
        _HEAVY_KEYS = frozenset({
            'pae_matrix',       # N×N float list — by far the biggest item
            'contact_probs',    # N×N float list — same size as pae_matrix
            'plddt_arr',        # per-residue float list
            'token_res_ids',    # per-residue int list
            'ranking_samples',  # per-sample ranking CSV rows (rarely needed)
        })
        for res in getattr(self, '_af3_analysis_results', []):
            entry = {k: v for k, v in res.items() if k not in _HEAVY_KEYS}
            entry['_lightweight'] = True   # signals reload-on-demand on load

            # JSON requires string keys. pair_metrics uses (chain_A, chain_B)
            # tuple keys → convert to "A-B" strings before serialising.
            pm = entry.get('pair_metrics')
            if isinstance(pm, dict):
                entry['pair_metrics'] = {
                    f"{k[0]}-{k[1]}" if isinstance(k, tuple) else str(k): v
                    for k, v in pm.items()
                }

            # best_pair is also a tuple → store as "A-B" string
            bp = entry.get('best_pair')
            if isinstance(bp, tuple):
                entry['best_pair'] = f"{bp[0]}-{bp[1]}" if len(bp) == 2 else str(bp)

            af3_analysis_ser.append(entry)

        # BLAST state
        blast_query = ''
        blast_html  = ''
        try:
            blast_query = self._blast_query_text.toPlainText().strip()
            blast_html  = self._blast_results_text.toHtml()
        except Exception:
            pass

        # HPC server state
        dv_pwd_enc = ''
        try:
            pwd = self._dv_pwd.text()
            if pwd:
                dv_pwd_enc = base64.b64encode(pwd.encode()).decode('ascii')
        except Exception:
            pass

        hpc_state = {
            'host':       getattr(self, '_dv_host',       None) and self._dv_host.text()       or '',
            'user':       getattr(self, '_dv_user',       None) and self._dv_user.text()       or '',
            'port':       getattr(self, '_dv_port',       None) and self._dv_port.value()      or 22,
            'password':   dv_pwd_enc,
            'base_path':  getattr(self, '_dv_base_path',  None) and self._dv_base_path.text()  or '~/af3_predictions',
            'af3cmd':     getattr(self, '_dv_af3cmd',     None) and self._dv_af3cmd.text()     or 'af3_run',
            'module_cmd': getattr(self, '_dv_module_cmd', None) and self._dv_module_cmd.text() or 'alphafold3',
        }

        # UI state
        ui_state = {
            'min_length':    self._min_length_spin.value(),
            'start_codons':  {'ATG': self._cb_atg.isChecked(),
                              'GTG': self._cb_gtg.isChecked(),
                              'TTG': self._cb_ttg.isChecked()},
            'filter_frame':  self._frame_combo.currentText(),
            'filter_strand': self._strand_combo.currentText(),
            'filter_min_aa': self._size_filter_spin.value(),
            'filter_source': self._source_combo.currentText(),
            'filter_search': self._search_edit.text(),
            'zoom_level':    self.zoom_level,
            'blast_algorithm':     self._algo_combo.currentText(),
            'blast_identity':      self._identity_spin.value(),
            'blast_evalue':        self._evalue_edit.text(),
            'blast_program':       self.blast_program,
            'blast_matrix':        self.blast_matrix,
            'blast_evalue_val':    self.blast_evalue,
            'blast_word_size':     self.blast_word_size,
            'blast_gap_open':      self.blast_gap_open,
            'blast_gap_ext':       self.blast_gap_ext,
            'blast_max_targets':   self.blast_max_targets,
            'blast_threshold':     self.blast_threshold,
            'blast_low_complexity':self.blast_low_complexity,
            'hmm_evalue':          self.hmm_evalue,
            'hmm_dom_evalue':      self.hmm_dom_evalue,
            'hmm_score_thresh':    self.hmm_score_thresh,
            'af3_n_neighbors':     self.af3_n_neighbors,
            'af3_max_residues':    self.af3_max_residues,
            'pyrodigal_params':    getattr(self, '_pyro_params', {}),
        }

        genome_rel = ''
        if proj_dir and self.dna_sequence:
            safe = re.sub(r'[^\w\.\-]', '_', self.genome_name or 'genome')
            genome_rel = f"genome/{safe}.fasta"

                # ── AF3 selection table ────────────────────────────────
        af3_sel = []
        try:
            for r in range(self._af3_sel_table.rowCount()):
                if self._af3_sel_table.item(r, 0):
                    af3_sel.append({
                        'orf_name': self._af3_sel_table.item(r, 0).text(),
                        'position': self._af3_sel_table.item(r, 1).text(),
                        'size_aa':  self._af3_sel_table.item(r, 2).text(),
                        'hmm':      self._af3_sel_table.item(r, 3).text(),
                        'note':     self._af3_sel_table.item(r, 4).text(),
                    })
        except Exception:
            pass

        # ── Sanitise af3_jobs before saving ───────────────────
        safe_af3_jobs = []
        for j in self.af3_jobs:
            jc = dict(j)
            jc.setdefault('sequences', [])
            jc.setdefault('iptm', None)
            jc.setdefault('plddt', None)
            jc.setdefault('status', 'unknown')
            safe_af3_jobs.append(jc)

        return {
            # ── Schema version ──────────────────────────────────
            'version':              self.PROJECT_VERSION,
            'schema':               2,   # bumped: definitive save format
            'saved_at':             datetime.now().isoformat(timespec='seconds'),
            # ── Genome ──────────────────────────────────────────
            'genome_name':          self.genome_name,
            'genome_file':          genome_rel,
            'dna_sequence':         self.dna_sequence,
            # ── ORFs (full annotation state) ────────────────────
            # Every annotation field is preserved inside each orf dict:
            # start, end, strand, frame, protein, dna, gc, source,
            # domains, observation, gene_name, putative_function,
            # custom_color, candidate_score, af3_user_note, notes
            'orfs':                 self.orfs,
            # ── HMM ─────────────────────────────────────────────
            'hmm_profiles':         hmm_manifest,
            'hmm_hits_all':         self.hmm_hits_all,
            # ── AlphaFold jobs ──────────────────────────────────
            'af3_jobs':             safe_af3_jobs,
            'af3_selected_orfs':    af3_sel,
            'result_files':         results_manifest,
            # ── AlphaFold analysis ──────────────────────────────
            'af3_analysis_results': af3_analysis_ser,
            'af3_analysis_dir':     getattr(self, '_af3_analysis_dir', ''),
            # ── BLAST ────────────────────────────────────────────
            'blast_query':          blast_query,
            'blast_results_html':   blast_html,
            # ── Other ────────────────────────────────────────────
            'snapgene':             {'features': self.snapgene_features,
                                    'primers':  self.snapgene_primers},
            'ui_state':             ui_state,
            'hpc_server':           hpc_state,
            'hpc_jobs':             getattr(self, '_hpc_jobs', []),
        }

    def save_project(self):
        """Save project as a single self-contained JSON file.
        Includes ORFs, HMM hits, AF3 jobs and AF3 result metadata.
        Heavy arrays (PAE matrices, pLDDT, contact_probs) are NOT saved —
        they are reloaded on demand from job_dir when the user selects a row.
        File → Save Project  (Ctrl+S style — single JSON, no file copying).
        """
        safe = re.sub(r'[^\w\.\-]', '_', self.genome_name or 'project')
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Project as JSON",
            f"{safe}.json",
            "Project JSON (*.json);;All (*)")
        if not path:
            return

        manifest = self._build_manifest(proj_dir=None)
        # Cap blast HTML to avoid bloating the file
        _bh = manifest.get('blast_results_html', '')
        if len(_bh) > 512 * 1024:
            manifest['blast_results_html'] = (
                _bh[:512*1024] + '<!-- TRUNCATED -->')

        # Write in background thread to keep UI responsive
        import threading as _threading, time as _ttime
        _save_err: list = []; _save_done: list = []
        def _bg_write():
            try:
                data = json.dumps(manifest, separators=(',', ':'),
                                  ensure_ascii=False)
                with open(path, 'w', encoding='utf-8') as fh:
                    fh.write(data)
            except Exception as _e:
                _save_err.append(str(_e))
            finally:
                _save_done.append(True)
        _t = _threading.Thread(target=_bg_write, daemon=True)
        _t.start()
        # Show a non-blocking wait cursor
        try:
            QApplication.setOverrideCursor(
                Qt.CursorShape.WaitCursor if QT_VERSION == 6 else Qt.WaitCursor)
        except Exception:
            pass
        while not _save_done:
            QApplication.processEvents()
            _ttime.sleep(0.02)
        _t.join(timeout=60)
        try:
            QApplication.restoreOverrideCursor()
        except Exception:
            pass
        if _save_err:
            QMessageBox.critical(self, "Save Project",
                                 f"Cannot write file:\n{_save_err[0]}")
            return

        n_af3 = len(manifest.get('af3_analysis_results', []))
        kb = Path(path).stat().st_size / 1024
        self._status.showMessage(
            f"✓ Saved: {Path(path).name}  "
            f"({len(self.orfs)} ORFs, {len(self.hmm_profiles)} HMM, "
            f"{n_af3} AF3 results, {kb:.0f} KB)")
        QMessageBox.information(
            self, "Project Saved",
            f"✓ Project saved!\n\n"
            f"  📄 {path}\n\n"
            f"  ORFs:         {len(self.orfs)}\n"
            f"  HMM profiles: {len(self.hmm_profiles)}\n"
            f"  AF3 jobs:     {len(self.af3_jobs)}\n"
            f"  AF3 analysis: {n_af3} result(s) — metadata only\n"
            f"  File size:    {kb:.0f} KB\n\n"
            f"  PAE/pLDDT arrays are reloaded on demand from job_dir.\n"
            f"  To reopen: File → Open Project → select this JSON file.\n\n"
            f"  ⚠ Server password stored as base64 (not encrypted).")

    def save_project_as(self):
        """Save As — project folder with ALL associated files copied.

        Structure:
            <ProjectName>/
                project.json          ← lightweight manifest + relative links
                genome/               ← FASTA
                hmm/                  ← .hmm profile files
                results/              ← AF3 score JSONs
                blast/                ← BLAST query sequence
                af3_predictions/      ← full AF3 output folders (PAE, model.cif…)
                    orf1303_vs_orf1298_up5/
                        *_confidences.json
                        *_summary_confidences.json
                        *_model.cif
                        ranking_scores.csv
                        seed_XXXX_sample-N/
                            confidences.json
                            model.cif
                            summary_confidences.json
                    orf1303_vs_orf1299_up4/
                        …
        """

        base_dir = QFileDialog.getExistingDirectory(
            self, "Choose folder to save project into")
        if not base_dir:
            return

        safe_name = re.sub(r'[^\w\.\-]', '_', self.genome_name or 'project')
        proj_dir  = Path(base_dir) / safe_name
        try:
            for sub in ('genome', 'hmm', 'results', 'blast', 'af3_predictions'):
                (proj_dir / sub).mkdir(parents=True, exist_ok=True)
        except OSError as e:
            QMessageBox.critical(self, "Save Project As",
                                 f"Cannot create directories:\n{e}")
            return

        # Progress dialog
        try:
            from PyQt6.QtWidgets import QProgressDialog
        except ImportError:
            from PyQt5.QtWidgets import QProgressDialog

        prog = QProgressDialog("Saving project...", "Cancel", 0, 100, self)
        prog.setWindowTitle("Save Project As")
        prog.setMinimumDuration(0)
        prog.setValue(0)
        QApplication.processEvents()

        def _step(pct, msg):
            prog.setValue(pct)
            prog.setLabelText(msg)
            QApplication.processEvents()
            return not prog.wasCanceled()

        # 1. Genome FASTA
        if not _step(5, "Saving genome FASTA..."): return
        genome_fname = f"{safe_name}.fasta"
        if self.dna_sequence:
            try:
                if self.current_fasta_path and Path(self.current_fasta_path).is_file():
                    shutil.copy2(self.current_fasta_path,
                                 proj_dir / "genome" / genome_fname)
                else:
                    with open(proj_dir / "genome" / genome_fname, 'w', encoding='utf-8') as fh:
                        fh.write(f">{self.genome_name}\n")
                        for i in range(0, len(self.dna_sequence), 60):
                            fh.write(self.dna_sequence[i:i+60] + "\n")
            except OSError:
                pass

        # 2. HMM profiles
        if not _step(12, "Copying HMM profiles..."): return
        for p in self.hmm_profiles:
            src = Path(p.get('file', ''))
            if src.is_file():
                try:
                    shutil.copy2(src, proj_dir / "hmm" / src.name)
                except OSError:
                    pass

        # 3. AF3 result score JSONs
        if not _step(18, "Saving AF3 score files..."): return
        for job in self.af3_jobs:
            if job.get('iptm') is not None:
                jfn = re.sub(r'[^\w\.\-]', '_', job['name']) + '.json'
                try:
                    with open(proj_dir / "results" / jfn, 'w', encoding='utf-8') as fh:
                        json.dump({
                            "name":       job['name'],
                            "modelSeeds": [],
                            "sequences":  job.get('sequences', []),
                            "dialect":    "alphafoldserver",
                            "version":    2,
                            "_iptm":      job.get('iptm'),
                            "_plddt":     job.get('plddt'),
                        }, fh, indent=2, ensure_ascii=False)
                except OSError:
                    pass

        # 4. BLAST query
        try:
            bq = self._blast_query_text.toPlainText().strip()
            if bq:
                with open(proj_dir / "blast" / "query.fasta", 'w', encoding='utf-8') as fh:
                    if not bq.startswith('>'):
                        fh.write(">blast_query\n")
                    fh.write(bq + "\n")
        except Exception:
            pass

        # 5. AF3 prediction folders — the big one
        af3_results = getattr(self, '_af3_analysis_results', [])
        n_af3       = len(af3_results)
        af3_pred_dir = proj_dir / 'af3_predictions'
        copied_jobs  = []

        for i, res in enumerate(af3_results):
            if prog.wasCanceled():
                break
            pct = 20 + int(65 * (i / max(n_af3, 1)))
            name = res.get('job_name', f'job_{i+1}')
            if not _step(pct, f"Copying AF3 folder {i+1}/{n_af3}: {name}"):
                break
            src_dir = Path(res.get('job_dir', ''))
            if not src_dir.is_dir():
                continue
            dst_dir = af3_pred_dir / src_dir.name
            try:
                shutil.copytree(str(src_dir), str(dst_dir), dirs_exist_ok=True)
                copied_jobs.append(src_dir.name)
            except Exception as e:
                print(f"[Save As] {src_dir.name}: {e}")


        # 6. Write project.json — stripped of redundant sequences
        # Genome is already saved as FASTA (step 1), so strip dna_sequence
        # and orf["dna"] to avoid serialising 8-15 MB on the main thread
        # (that caused the visible 88% hang on Windows).
        if not _step(88, "Writing project manifest..."): return
        manifest = self._build_manifest(proj_dir=proj_dir)
        manifest['project_copy_mode'] = 'full_copy_light_manifest'
        manifest['af3_analysis_dir'] = 'af3_predictions'
        for entry in manifest.get('af3_analysis_results', []):
            jname = entry.get('job_name', '')
            entry['job_dir'] = f"af3_predictions/{jname}"
        # Strip sequences already saved elsewhere
        manifest['dna_sequence'] = ''
        for _orf in manifest.get('orfs', []):
            _orf.pop('dna', None)
        _bh = manifest.get('blast_results_html', '')
        if len(_bh) > 512 * 1024:
            manifest['blast_results_html'] = _bh[:512*1024] + '<!-- TRUNCATED -->'
        # Serialise + write in a background thread so the progress dialog
        # stays alive on Windows (json.dump blocks the event loop).
        import threading as _threading, time as _ttime
        _err: list = []; _done: list = []
        def _do_write():
            try:
                data = json.dumps(manifest, separators=(',', ':'), ensure_ascii=False)
                with open(proj_dir / self.PROJECT_MANIFEST, 'w', encoding='utf-8') as fh:
                    fh.write(data)
            except Exception as _e:
                _err.append(str(_e))
            finally:
                _done.append(True)
        _t = _threading.Thread(target=_do_write, daemon=True)
        _t.start()
        _pct = 88
        while not _done:
            QApplication.processEvents()
            if prog.wasCanceled(): break
            _pct = min(97, _pct + 1)
            prog.setValue(_pct)
            _ttime.sleep(0.05)
        _t.join(timeout=60)
        if _err:
            QMessageBox.critical(self, "Save Project As",
                                 f"Cannot write manifest:\n{_err[0]}")
            prog.close(); return
        _step(100, "Done!")
        prog.close()

        total_kb = sum(f.stat().st_size for f in proj_dir.rglob('*')
                       if f.is_file()) / 1024
        self._status.showMessage(
            f"✓ Saved As: {proj_dir.name}/  "
            f"({len(copied_jobs)} AF3 folders, {total_kb:.0f} KB total)")
        QMessageBox.information(
            self, "Project Saved As",
            f"✓ Full project copy saved!\n\n"
            f"  📁 {proj_dir}\n\n"
            f"  genome/           — FASTA\n"
            f"  hmm/              — {len(self.hmm_profiles)} HMM profile(s)\n"
            f"  results/          — AF3 score JSONs\n"
            f"  blast/            — BLAST query\n"
            f"  af3_predictions/  — {len(copied_jobs)} prediction folder(s)\n"
            f"    (PAE, pLDDT, model.cif, ranking_scores, seed_*/)\n"
            f"  project.json      — lightweight manifest (AF3 arrays reloaded from folders)\n\n"
            f"  Total size: {total_kb:.0f} KB\n\n"
            f"⚠ Server password stored as base64 (not encrypted).")

    def load_project(self):
        """Open a project saved by save_project (JSON) or save_project_as (folder)."""

        dlg = QDialog(self)
        dlg.setWindowTitle("Open Project")
        dlg.setFixedSize(430, 150)
        dlg_lay = QVBoxLayout(dlg)
        dlg_lay.addWidget(QLabel(
            "<b>How is the project stored?</b><br>"
            "Use <b>JSON file</b> for single-file projects (File → Save Project).<br>"
            "Use <b>Project folder</b> for full-copy projects (File → Save Project As)."))
        btn_row = QHBoxLayout()
        btn_file   = QPushButton("📄  Select JSON file")
        btn_folder = QPushButton("📁  Select project folder")
        btn_cancel = QPushButton("Cancel")
        btn_file.setDefault(True)
        for b in (btn_file, btn_folder, btn_cancel):
            b.setFixedHeight(32); btn_row.addWidget(b)
        dlg_lay.addLayout(btn_row)

        _choice = ['']
        btn_file.clicked.connect(  lambda: (_choice.__setitem__(0, 'file'),   dlg.accept()))
        btn_folder.clicked.connect(lambda: (_choice.__setitem__(0, 'folder'), dlg.accept()))
        btn_cancel.clicked.connect(dlg.reject)

        if not (dlg.exec() if QT_VERSION == 6 else dlg.exec_()):
            return
        mode = _choice[0]

        if mode == 'file':
            f, _ = QFileDialog.getOpenFileName(
                self, "Open Project JSON", "", "JSON (*.json);;All (*)")
            if not f: return
            json_path = Path(f)
            proj_dir  = json_path.parent
        else:
            d = QFileDialog.getExistingDirectory(self, "Select project folder")
            if not d: return
            proj_dir  = Path(d)
            json_path = proj_dir / self.PROJECT_MANIFEST
            if not json_path.is_file():
                QMessageBox.warning(self, "Open Project",
                    f"Folder does not contain '{self.PROJECT_MANIFEST}'.\n"
                    f"Select the root project folder.")
                return

        try:
            with open(json_path, encoding='utf-8') as fh:
                data = json.load(fh)
        except (OSError, json.JSONDecodeError) as e:
            QMessageBox.critical(self, "Open Project", f"Cannot read file:\n{e}")
            return

        ver = data.get('version', 'unknown')

        # ══ DEFINITIVE LOAD (schema v2) ═══════════════════════

        # ── 1. Genome ─────────────────────────────────────────
        self.genome_name   = data.get('genome_name', '')
        self.dna_sequence  = data.get('dna_sequence', '')
        rel_genome = data.get('genome_file', '')
        if rel_genome:
            abs_g = proj_dir / rel_genome
            self.current_fasta_path = str(abs_g) if abs_g.is_file() else ''

        # ── Reload genome FASTA when dna_sequence was stripped during save ──
        # save_project_as strips dna_sequence (stored as FASTA in genome/).
        # Reload it transparently here so the genome map and ORF table work.
        if not self.dna_sequence and rel_genome:
            abs_g = proj_dir / rel_genome
            if abs_g.is_file():
                try:
                    with open(abs_g, encoding='utf-8', errors='replace') as _f:
                        _seqs = {}
                        _cur  = None
                        for _line in _f:
                            _line = _line.rstrip()
                            if _line.startswith('>'):
                                _cur = _line[1:].split()[0]
                                _seqs[_cur] = []
                            elif _cur:
                                _seqs[_cur].append(_line)
                    if _seqs:
                        _name, _parts = next(iter(_seqs.items()))
                        self.dna_sequence = ''.join(_parts).upper()
                        if not self.genome_name:
                            self.genome_name = _name
                except (OSError, StopIteration):
                    pass

        # ── 2. ORFs (full annotation state) ───────────────────
        self.orfs = data.get('orfs', [])
        # Ensure every annotation field exists with safe defaults
        _orf_defaults = {
            'domains': [], 'observation': '', 'gene_name': '',
            'putative_function': '', 'custom_color': '', 'notes': '',
            'af3_user_note': '', 'candidate_score': 0.0,
            'gc': 0.0, 'source': '6frame',
        }
        for orf in self.orfs:
            for k, v in _orf_defaults.items():
                orf.setdefault(k, v)
        self.filtered_orfs = self.orfs.copy()

        # ── Reconstruct orf["dna"] when stripped during save ──────────
        # save_project_as strips orf["dna"] (derivable from dna_sequence).
        # Rebuild now so all downstream code that expects orf["dna"] works.
        if self.dna_sequence:
            _dna = self.dna_sequence
            _dn  = len(_dna)
            for _orf in self.orfs:
                if not _orf.get("dna"):
                    _s = max(0, _orf.get("start", 0))
                    _e = min(_dn, _orf.get("end", _dn))
                    _orf["dna"] = _dna[_s:_e]

        # Re-apply manual annotations saved separately (orf_annotations)
        # This handles projects where orfs were saved WITHOUT annotation fields
        for ann in data.get('orf_annotations', []):
            i = ann.get('idx', -1)
            if 0 <= i < len(self.orfs):
                for k in ('observation', 'putative_function', 'gene_name',
                          'custom_color', 'af3_user_note', 'candidate_score'):
                    if ann.get(k):
                        self.orfs[i][k] = ann[k]

        # ── 3. AF3 jobs ────────────────────────────────────────
        self.af3_jobs = data.get('af3_jobs', [])
        for _j in self.af3_jobs:
            _j.setdefault('sequences', [])
            _j.setdefault('iptm', None)
            _j.setdefault('plddt', None)
            _j.setdefault('status', 'unknown')

        # ── Restore AF3 selection table (Bug 1 fix) ──────────
        self._af3_sel_table.setRowCount(0)
        for entry in data.get('af3_selected_orfs', []):
            row = self._af3_sel_table.rowCount()
            self._af3_sel_table.insertRow(row)
            for col, key in enumerate(
                    ['orf_name', 'position', 'size_aa', 'hmm', 'note']):
                self._af3_sel_table.setItem(
                    row, col, QTableWidgetItem(entry.get(key, '')))
        self._af3_sel_count.setText(
            f"{self._af3_sel_table.rowCount()} ORFs selected")

        # ── HMM profiles ──────────────────────────────────────
        self.hmm_profiles = []
        for p in data.get('hmm_profiles', []):
            rel = p.get('file', '')
            abs_path = str(proj_dir / rel) if rel else ''
            # If abs_path doesn't exist (JSON-only save), keep empty
            if abs_path and not Path(abs_path).is_file():
                abs_path = ''
            self.hmm_profiles.append({
                'name':     p.get('name', ''),
                'file':     abs_path,
                'color':    p.get('color', '#888888'),
                'function': p.get('function', ''),
                'hits':     p.get('hits', []),
            })
        if 'hmm_hits_all' in data and data['hmm_hits_all']:
            self.hmm_hits_all = data['hmm_hits_all']
        else:
            self.hmm_hits_all = []
            for p in self.hmm_profiles:
                for h in p.get('hits', []):
                    self.hmm_hits_all.append(dict(h, profile_name=p['name']))

        # ── Re-inject HMM domains into ORFs (Bug 4 fix) ──────
        if self.hmm_hits_all and self.orfs:
            for orf in self.orfs:
                if 'domains' not in orf:
                    orf['domains'] = []
            for hit in self.hmm_hits_all:
                oi = hit.get('orf_index', -1)
                if 0 <= oi < len(self.orfs):
                    pn = hit.get('profile_name', hit.get('hmm_name', '?'))
                    existing_domains = [d['domain'] for d in self.orfs[oi]['domains']]
                    if pn not in existing_domains:
                        self.orfs[oi]['domains'].append({
                            'domain':      pn,
                            'description': hit.get('profile_function', f'HMM: {pn}'),
                            'system':      hit.get('profile_function', 'HMM hit'),
                            'role':        'HMM',
                            'start':       hit.get('ali_from', 0),
                            'end':         hit.get('ali_to', 0),
                            'evalue':      hit.get('evalue', 999),
                            'score':       hit.get('score', 0),
                        })

        # ── AF3 result JSONs ──────────────────────────────────
        for rel_res in data.get('result_files', []):
            res_path = proj_dir / rel_res
            if not res_path.is_file(): continue
            try:
                with open(res_path, encoding='utf-8') as fh:
                    rdata = json.load(fh)
                rname = rdata.get('name', res_path.stem)
                for job in self.af3_jobs:
                    if job['name'] == rname:
                        job['iptm']  = rdata.get('_iptm')
                        job['plddt'] = rdata.get('_plddt')
                        if job['iptm'] is not None:
                            job['status'] = 'done'
                        break
            except Exception:
                pass

        # ── SnapGene ──────────────────────────────────────────
        sg = data.get('snapgene', {})
        self.snapgene_features = sg.get('features', [])
        self.snapgene_primers  = sg.get('primers', [])

        # ── BLAST ─────────────────────────────────────────────
        try:
            bq = data.get('blast_query', '')
            if bq: self._blast_query_text.setPlainText(bq)
            bh = data.get('blast_results_html', '')
            if bh: self._blast_results_text.setHtml(bh)
        except (AttributeError, TypeError) as e:
            print(f"[load_project] BLAST state restore: {e}")

        # ── AF3 Analysis results (v2) ─────────────────────────
        af3_results_raw = data.get('af3_analysis_results', [])
        self._af3_analysis_results = []
        af3_dir_rel = data.get('af3_analysis_dir', '')
        for entry in af3_results_raw:
            # Resolve job_dir path — may be relative to proj_dir
            jdir = entry.get('job_dir', '')
            if jdir and not Path(jdir).is_absolute():
                abs_jdir = proj_dir / jdir
                entry['job_dir'] = str(abs_jdir) if abs_jdir.is_dir() else jdir
            # Restore pair_metrics: "A-B" string keys → (A, B) tuple keys
            pm = entry.get('pair_metrics')
            if isinstance(pm, dict):
                restored = {}
                for k, v in pm.items():
                    if isinstance(k, str) and '-' in k:
                        parts = k.split('-', 1)
                        restored[(parts[0], parts[1])] = v
                    else:
                        restored[k] = v
                entry['pair_metrics'] = restored
            # Restore best_pair: "A-B" string → (A, B) tuple
            bp = entry.get('best_pair')
            if isinstance(bp, str) and '-' in bp:
                parts = bp.split('-', 1)
                entry['best_pair'] = (parts[0], parts[1])
            self._af3_analysis_results.append(entry)
        if af3_dir_rel:
            abs_af3_dir = proj_dir / af3_dir_rel
            self._af3_analysis_dir = str(abs_af3_dir) if abs_af3_dir.is_dir() else af3_dir_rel

        # Full-copy projects may store AF3 analysis in lightweight form to avoid
        # very large project.json files. In that case, rebuild the full analysis
        # records directly from the copied af3_predictions/ folders.
        needs_af3_rescan = False
        if self._af3_analysis_results:
            needs_af3_rescan = any(r.get('_lightweight') or 'pae_matrix' not in r
                                   for r in self._af3_analysis_results)
        if needs_af3_rescan:
            af3_dir = getattr(self, '_af3_analysis_dir', '')
            if af3_dir and Path(af3_dir).is_dir():
                try:
                    self._af3a_scan_folder(af3_dir)
                except Exception:
                    try: self._af3a_populate_table()
                    except AttributeError as e:
                        print(f"[load_project] AF3 table populate: {e}")
            else:
                try: self._af3a_populate_table()
                except AttributeError as e:
                    print(f"[load_project] AF3 table populate: {e}")
                if af3_dir:
                    self._status.showMessage(
                        f'\u26a0 AF3 analysis dir not found: {af3_dir} '
                        '-- scores loaded, PAE/pLDDT unavailable')
        else:
            # Repopulate the AF3 Analysis table
            try:
                self._af3a_populate_table()
            except AttributeError as e:
                print(f"[load_project] AF3 table populate: {e}")

        # ── UI state ──────────────────────────────────────────
        ui = data.get('ui_state', {})
        try:
            if 'min_length' in ui:
                self._min_length_spin.setValue(ui['min_length'])
                self.min_length = ui['min_length']
            codons = ui.get('start_codons', {})
            if codons:
                self._cb_atg.setChecked(codons.get('ATG', True))
                self._cb_gtg.setChecked(codons.get('GTG', True))
                self._cb_ttg.setChecked(codons.get('TTG', True))
                self.start_codons = {c for c, v in codons.items() if v}
            for attr, combo in [('filter_frame', self._frame_combo),
                                 ('filter_strand', self._strand_combo),
                                 ('filter_source', self._source_combo)]:
                if attr in ui:
                    idx = combo.findText(ui[attr])
                    if idx >= 0: combo.setCurrentIndex(idx)
            if 'filter_min_aa' in ui: self._size_filter_spin.setValue(ui['filter_min_aa'])
            if 'filter_search' in ui: self._search_edit.setText(ui['filter_search'])
            if 'blast_algorithm' in ui:
                idx = self._algo_combo.findText(ui['blast_algorithm'])
                if idx >= 0: self._algo_combo.setCurrentIndex(idx)
                self.algo_choice = ui['blast_algorithm']
            if 'blast_identity' in ui:
                self._identity_spin.setValue(ui['blast_identity'])
                self.blast_threshold = ui['blast_identity']
            if 'blast_evalue' in ui: self._evalue_edit.setText(ui['blast_evalue'])
            for k in ('blast_program', 'blast_matrix', 'blast_evalue_val',
                      'blast_word_size', 'blast_gap_open', 'blast_gap_ext',
                      'blast_max_targets', 'blast_threshold', 'blast_low_complexity',
                      'hmm_evalue', 'hmm_dom_evalue', 'hmm_score_thresh',
                      'af3_n_neighbors', 'af3_max_residues'):
                if k in ui: setattr(self, k, ui[k])
            if 'blast_program' in ui:
                try: self._blast_prog_combo.setCurrentText(ui['blast_program'])
                except AttributeError as e:
                    print(f"[load_project] blast_program widget: {e}")
            if 'af3_n_neighbors' in ui:
                try: self._af3_nb_spin.setValue(ui['af3_n_neighbors'])
                except AttributeError as e:
                    print(f"[load_project] af3_n_neighbors widget: {e}")
            if 'pyrodigal_params' in ui and ui['pyrodigal_params']:
                self._pyro_params.update(ui['pyrodigal_params'])
            if 'zoom_level' in ui:
                self.zoom_level = ui['zoom_level']
                self._zoom_label.setText(f"{int(self.zoom_level * 100)}%")
        except (AttributeError, TypeError, ValueError) as e:
            print(f"[load_project] UI state restore: {e}")

        # ── HPC server ───────────────────────────────────────────
        dv = data.get('hpc_server', {})
        try:
            if dv.get('host'):       self._dv_host.setText(dv['host'])
            if dv.get('user'):       self._dv_user.setText(dv['user'])
            if dv.get('port'):       self._dv_port.setValue(int(dv['port']))
            if dv.get('password'):
                try:
                    self._dv_pwd.setText(
                        base64.b64decode(dv['password'].encode()).decode('utf-8'))
                except (ValueError, UnicodeDecodeError) as e:
                    print(f"[load_project] password decode: {e}")
            if dv.get('base_path'): self._dv_base_path.setText(dv['base_path'])
            if dv.get('af3cmd'):    self._dv_af3cmd.setText(dv['af3cmd'])
            if dv.get('module_cmd') is not None:
                self._dv_module_cmd.setText(dv['module_cmd'])
        except AttributeError as e:
            print(f"[load_project] HPC server state: {e}")
        try:
            self._hpc_jobs = data.get('hpc_jobs', [])
            self._dv_refresh_monitor_table()
        except AttributeError as e:
            print(f"[load_project] HPC jobs restore: {e}")

        # ── Refresh all UI ────────────────────────────────────
        self._update_orfs_list()
        self._update_info()
        self._genome_map.set_data(
            len(self.dna_sequence), self.orfs,
            self.hmm_profiles, self.dna_sequence)
        self._genome_map.set_zoom(self.zoom_level)
        self._update_hmm_profile_table()
        self._af3_update_jobs_table()
        self._refresh_blast_cmd_preview()
        # Refresh Genomic PPI Map if AF3 results are loaded
        if getattr(self, '_af3_analysis_results', []):
            try:
                self._ppi_arc_map_refresh()
            except Exception as _e:
                print(f"[load_project] PPI map refresh: {_e}")

        n_annot    = sum(1 for o in self.orfs if o.get('observation')
                         or o.get('putative_function') or o.get('gene_name')
                         or o.get('custom_color'))
        n_hmm_hits = sum(len(p.get('hits', [])) for p in self.hmm_profiles)
        n_af3a     = len(self._af3_analysis_results)
        n_af3_sel     = self._af3_sel_table.rowCount() if hasattr(self, '_af3_sel_table') else 0

        self._status.showMessage(
            f"✓ {self.genome_name}  "
            f"{len(self.orfs)} ORFs | "
            f"{n_hmm_hits} HMM hits | "
            f"{len(self.af3_jobs)} AF3 jobs ({n_af3_sel} sel) | "
            f"{n_af3a} AF3 analyses | "
            f"{n_annot} annotated  [{ver}]")

        saved_at = data.get('saved_at', '')
        lines = [
            'Projeto carregado com sucesso!',
            '',
            f'  Genoma:          {self.genome_name}',
            f'  ORFs:            {len(self.orfs)} ({n_annot} anotadas manualmente)',
            f'  HMM profiles:    {len(self.hmm_profiles)} ({n_hmm_hits} hits)',
            f'  HMM hits_all:    {len(self.hmm_hits_all)}',
            f'  AF3 selecionados:{n_af3_sel} ORFs na fila',
            f'  AF3 jobs:        {len(self.af3_jobs)}',
            f'  AF3 analyses:    {n_af3a} (com PAE/pLDDT)',
            f'  Formato:         {ver}',
        ]
        if saved_at: lines.append(f'  Salvo em:        {saved_at}')
        QMessageBox.information(self, 'Projeto Carregado', '\n'.join(lines))

    def export_map_pdf(self):
        """Export genome map as PDF or PNG."""
        if not self.dna_sequence:
            QMessageBox.warning(self, "Export", "Load a sequence first!")
            return
        f, _ = QFileDialog.getSaveFileName(self, "Export Genome Map",
            f"{self.genome_name or 'genome_map'}.png",
            "PNG Image (*.png);;PDF Document (*.pdf);;SVG (*.svg);;All (*)")
        if not f:
            return
        try:
            if f.lower().endswith('.pdf'):
                try:
                    if QT_VERSION == 6:
                        from PyQt6.QtPrintSupport import QPrinter
                    else:
                        from PyQt5.QtPrintSupport import QPrinter
                    printer = QPrinter(QPrinter.RenderMode.HighResolution if QT_VERSION == 6
                                        else QPrinter.HighResolution)
                    printer.setOutputFormat(QPrinter.OutputFormat.PdfFormat if QT_VERSION == 6
                                             else QPrinter.PdfFormat)
                    printer.setOutputFileName(f)
                    painter = QPainter(printer)
                    self._genome_map.render(painter)
                    painter.end()
                except ImportError:
                    # Fallback to PNG if print support not available
                    f = f.replace('.pdf', '.png')
                    self._genome_map.grab().save(f)
            elif f.lower().endswith('.svg'):
                try:
                    if QT_VERSION == 6:
                        from PyQt6.QtSvg import QSvgGenerator
                    else:
                        from PyQt5.QtSvg import QSvgGenerator
                    gen = QSvgGenerator()
                    gen.setFileName(f)
                    gen.setSize(self._genome_map.size())
                    painter = QPainter(gen)
                    self._genome_map.render(painter)
                    painter.end()
                except ImportError:
                    f = f.replace('.svg', '.png')
                    self._genome_map.grab().save(f)
            else:
                # PNG
                self._genome_map.grab().save(f)

            sz = os.path.getsize(f) / 1024
            self._status.showMessage(f"✓ Map exported: {Path(f).name} ({sz:.0f} KB)")
        except Exception as e:
            QMessageBox.critical(self, "Export Error", str(e))

    def export_genbank(self):
        if not self.dna_sequence: return
        f, _ = QFileDialog.getSaveFileName(self, "Export GenBank", "", "GenBank (*.gb)")
        if f:
            write_genbank(f, self.dna_sequence, self.orfs,
                         self.snapgene_features, self.genome_name)
            self._status.showMessage("✓ GenBank exported")

    def export_snapgene(self):
        if not self.dna_sequence: return
        f, _ = QFileDialog.getSaveFileName(self, "Export SnapGene", "", "SnapGene (*.dna)")
        if f:
            features = []
            for i, orf in enumerate(self.orfs):
                features.append({'name': f"ORF{i+1}", 'type': 'CDS',
                    'start': orf['start'], 'end': orf['end'],
                    'strand': orf['strand'], 'color': '#99ccff'})
            write_snapgene_dna(f, self.dna_sequence, features,
                              name=self.genome_name)
            self._status.showMessage("✓ SnapGene exported")

    # ═══════════════════════════════════════════════════════════
    # UI UPDATE METHODS
    # ═══════════════════════════════════════════════════════════

    def _update_orfs_list(self):
        # Block all signals and disable sorting during the rebuild.
        # Without this, every insertRow fires selectionChanged which calls
        # _on_orf_table_select → _select_and_center_orf, leaving
        # _orf_table_selecting=True and making subsequent centering silently
        # no-op.  Sorting must also be off during insert or Qt re-orders
        # rows mid-way and corrupts the filtered_orfs→row mapping.
        self._orf_table.setSortingEnabled(False)
        self._orf_table.blockSignals(True)
        self._orf_table_selecting = False   # reset any stuck re-entrancy flag
        self._orf_table.setRowCount(0)
        # Build lookup: orf_name → list of AF3 analysis results
        af3_lookup = {}
        for res in getattr(self, '_af3_analysis_results', []):
            for orf_name in res.get('orf_names', []):
                af3_lookup.setdefault(orf_name, []).append(res)

        for i, orf in enumerate(self.filtered_orfs):
            n = self.orfs.index(orf) + 1 if orf in self.orfs else i + 1
            orf_label = f"ORF{n}"
            hmm_parts = []
            for d in orf.get('domains', []):
                domain_name = d['domain']
                region = d.get('ali_region', '')
                profile_display = "?"
                for profile in self.hmm_profiles:
                    if any(h.get('profile_name') == domain_name or
                           h.get('hmm_name') == domain_name
                           for h in profile.get('hits', [])):
                        profile_display = profile['name']
                        break
                label = f"{profile_display}:{domain_name}"
                if region:
                    label += f" [{region}]"
                hmm_parts.append(label)

            hmm = ', '.join(hmm_parts) or '-'

            # AF3 columns — find best result for this ORF
            af3_results = af3_lookup.get(orf_label, [])
            if af3_results:
                best = max(af3_results,
                           key=lambda r: r.get('iptm', 0) or 0)
                af3_done     = '✅'
                partner      = best.get('partner_name', '-')
                iptm_s       = f"{best.get('iptm', 0):.3f}" \
                               if best.get('iptm') is not None else '-'
                pae_inter_s  = f"{best.get('pae_inter', 0):.1f} Å" \
                               if best.get('pae_inter') is not None else '-'
                pae_min_v    = best.get('pae_min_inter')
                pae_min_s    = f"{pae_min_v:.2f} Å" if pae_min_v is not None else '-'
                cp_iptm_v    = best.get('cp_iptm_inter')
                cp_iptm_s    = f"{cp_iptm_v:.2f}" if cp_iptm_v is not None else '-'
                contact_s    = best.get('contact_region', '-')
            else:
                af3_done     = '-'
                partner      = '-'
                iptm_s       = '-'
                pae_inter_s  = '-'
                pae_min_s    = '-'
                cp_iptm_s    = '-'
                contact_s    = '-'

            user_note = orf.get('af3_user_note', '')

            row = self._orf_table.rowCount()
            self._orf_table.insertRow(row)
            items = [
                orf_label, f"F{orf['frame']}", orf['strand'],
                f"{orf['start']:,}", f"{orf['end']:,}",
                str(len(orf['protein'].rstrip('*'))),
                f"{orf['gc']:.1f}", hmm[:50],
                f"{orf.get('candidate_score',0):.2f}",
                orf.get('source', '6frame'),
                orf.get('observation', '')[:30],
                af3_done, partner, iptm_s, pae_inter_s, pae_min_s, cp_iptm_s, contact_s, user_note,
            ]
            for col, val in enumerate(items):
                # Use numeric-sorting item for columns that contain numbers:
                #   0  = ORF id  (ORF1, ORF2 … → sort by trailing integer)
                #   3  = Start bp
                #   4  = End bp
                #   5  = Size(aa)
                #   8  = Score
                #   13 = ipTM
                #   14 = PAE_inter
                #   15 = PAE_min
                #   16 = cp_ipTM
                if col in (0, 3, 4, 5, 8, 13, 14, 15, 16):
                    item = _OrfNumericItem(str(val))
                else:
                    item = QTableWidgetItem(str(val))
                # Color AF3 columns by quality
                if col == 13 and iptm_s not in ('-', ''):  # ipTM
                    try:
                        v = float(iptm_s)
                        if v >= 0.75:
                            item.setBackground(QColor('#C8E6C9'))
                        elif v >= 0.50:
                            item.setBackground(QColor('#FFF9C4'))
                        else:
                            item.setBackground(QColor('#FFCDD2'))
                    except ValueError:
                        pass
                if col == 15 and pae_min_s not in ('-', ''):  # PAE_min ★
                    try:
                        v = float(pae_min_s.replace(' Å','').replace('Å',''))
                        if v < 4.0:
                            item.setBackground(QColor('#C8E6C9'))
                        elif v < 8.0:
                            item.setBackground(QColor('#FFF9C4'))
                        else:
                            item.setBackground(QColor('#FFCDD2'))
                    except ValueError:
                        pass
                if col == 16 and cp_iptm_s not in ('-', ''):  # cp_ipTM ★
                    try:
                        v = float(cp_iptm_s)
                        if v >= 0.65:
                            item.setBackground(QColor('#C8E6C9'))
                        elif v >= 0.50:
                            item.setBackground(QColor('#FFF9C4'))
                        else:
                            item.setBackground(QColor('#FFCDD2'))
                    except ValueError:
                        pass
                # Make User_note column editable
                if col == 18:
                    item.setFlags(item.flags() |
                                  (Qt.ItemFlag.ItemIsEditable if QT_VERSION == 6
                                   else Qt.ItemIsEditable))
                self._orf_table.setItem(row, col, item)

        self._orf_table.blockSignals(False)
        self._orf_table.setSortingEnabled(True)
        self._orf_count_label.setText(
            f"({len(self.filtered_orfs)} of {len(self.orfs)})")

    def _on_orf_note_changed(self, item):
        """Save edited User_note back to the orf dict."""
        if item.column() != 16:
            return
        row = item.row()
        id_item = self._orf_table.item(row, 0)
        if not id_item:
            return
        orf_label = id_item.text()
        try:
            idx = int(orf_label.replace('ORF', '')) - 1
            if 0 <= idx < len(self.orfs):
                self.orfs[idx]['af3_user_note'] = item.text()
        except (ValueError, IndexError):
            pass

    def _update_info(self):
        if not self.dna_sequence:
            self._info_text.setPlainText("No genome loaded"); return
        nd = sum(1 for o in self.orfs if o.get('domains'))
        n_pyro = sum(1 for o in self.orfs if o.get('source') == 'pyrodigal')
        n_6fr = sum(1 for o in self.orfs if o.get('source','6frame') == '6frame')
        src = f"6frame:{n_6fr} pyro:{n_pyro}" if n_pyro else f"6frame:{n_6fr}"
        self._info_text.setPlainText(
            f"{self.genome_name[:25]}\n{len(self.dna_sequence):,} bp\n"
            f"GC: {self.analyzer.gc_content(self.dna_sequence):.1f}%\n"
            f"ORFs: {len(self.orfs)} ({src})\nDomains: {nd}")

    def _update_map(self):
        self._genome_map.set_data(len(self.dna_sequence), self.orfs,
                                   self.hmm_profiles, self.dna_sequence)
        self._update_hits_legend()
        # Execute any pending center request now that data is loaded
        pending = getattr(self, '_pending_center_idx', -1)
        if pending >= 0:
            self._select_and_center_orf(pending)

    def _set_zoom(self, level):
        target = max(0.5, min(10000.0, level))
        center_x = self._genome_map.width() / 2
        self._genome_map.set_zoom(target, anchor_x=center_x)

    def _on_map_zoom_changed(self, level):
        self.zoom_level = float(level)
        self._zoom_label.setText(f"{int(self.zoom_level * 100)}%")

    def filter_orfs(self):
        self.filtered_orfs = self.orfs.copy()
        s = self._search_edit.text().upper()
        if s:
            self.filtered_orfs = [o for o in self.filtered_orfs
                if s in o['protein'] or s in str(self.orfs.index(o)+1)]
        fr = self._frame_combo.currentText()
        if fr != "All":
            self.filtered_orfs = [o for o in self.filtered_orfs if o['frame'] == int(fr)]
        st = self._strand_combo.currentText()
        if st != "All":
            self.filtered_orfs = [o for o in self.filtered_orfs if o['strand'] == st]
        ms = self._size_filter_spin.value()
        if ms > 0:
            self.filtered_orfs = [o for o in self.filtered_orfs
                if len(o['protein'].rstrip('*')) >= ms]
        src = self._source_combo.currentText()
        if src != "All":
            self.filtered_orfs = [o for o in self.filtered_orfs
                if o.get('source', '6frame') == src]
        self._update_orfs_list()

    def _on_text_orf_click(self, url):
        """Handle clicks on ORF links in BLAST/HMM results."""
        url_str = url.toString() if hasattr(url, 'toString') else str(url)
        # Handle both "orf:123" and full URL like "orf:123" or just the path
        m = re.search(r'orf[:/]?(\d+)', url_str)
        if m:
            idx = int(m.group(1))
            if 0 <= idx < len(self.orfs):
                self._select_and_center_orf(idx)

    def _on_legend_click(self):
        """Click on hits legend → select and center that ORF."""
        rows = self._hits_legend.selectionModel().selectedRows()
        if not rows:
            return
        row = rows[0].row()
        item = self._hits_legend.item(row, 1)
        if item:
            text = item.text()
            m = re.match(r'ORF(\d+)', text)
            if m:
                idx = int(m.group(1)) - 1
                if 0 <= idx < len(self.orfs):
                    self._select_and_center_orf(idx)

    def _update_hits_legend(self):
        """Refresh the hits legend panel next to the genome map."""
        self._hits_legend.setRowCount(0)

        # Collect all annotated ORFs (HMM + custom colored)
        entries = []  # (color, label, orf_index)

        # HMM hits
        for profile in self.hmm_profiles:
            color = profile['color']
            for hit in profile.get('hits', []):
                oi = hit.get('orf_index', -1)
                if 0 <= oi < len(self.orfs):
                    label = f"ORF{oi+1} {profile['name']}"
                    entries.append((color, label, oi))

        # Custom colored ORFs (without HMM)
        for i, orf in enumerate(self.orfs):
            if orf.get('custom_color') and not any(e[2] == i for e in entries):
                entries.append((orf['custom_color'], f"ORF{i+1} (custom)", i))

        # Sort by ORF index
        entries.sort(key=lambda x: x[2])

        # Remove duplicates (same ORF)
        seen = set()
        unique = []
        for color, label, oi in entries:
            if oi not in seen:
                unique.append((color, label, oi))
                seen.add(oi)

        for color, label, oi in unique:
            row = self._hits_legend.rowCount()
            self._hits_legend.insertRow(row)
            # Color swatch
            color_item = QTableWidgetItem("■")
            color_item.setForeground(QColor(color))
            color_item.setFont(QFont('Arial', 14))
            self._hits_legend.setItem(row, 0, color_item)
            self._hits_legend.setRowHeight(row, 18)
            # Label
            self._hits_legend.setItem(row, 1, QTableWidgetItem(label))

    def _on_orf_table_select(self):
        # Guard: _select_and_center_orf calls selectRow which re-triggers this.
        if getattr(self, '_orf_table_selecting', False):
            return
        rows = self._orf_table.selectionModel().selectedRows()
        if not rows:
            return
        row = rows[0].row()
        if row >= len(self.filtered_orfs):
            return
        orf = self.filtered_orfs[row]
        orf_idx = self.orfs.index(orf) if orf in self.orfs else -1
        if orf_idx >= 0:
            self._select_and_center_orf(orf_idx)

    def _on_map_orf_click(self, idx):
        if 0 <= idx < len(self.orfs):
            self._select_and_center_orf(idx)

    @staticmethod
    def _parse_orf_idx_from_text(text: str) -> int:
        """Extract a 0-based ORF index from any string that contains 'ORF{N}'.
        Handles emoji prefixes (e.g. '✅ ORF2588'), plain 'ORF2588', etc.
        Returns -1 if no valid index found.
        """
        m = re.search(r'ORF(\d+)', text, re.IGNORECASE)
        if m:
            return int(m.group(1)) - 1
        # Fallback: any digit sequence in the text
        nums = re.findall(r'\d+', text)
        if nums:
            try:
                return int(nums[-1]) - 1
            except ValueError:
                pass
        return -1

    def _show_orf_details(self, orf):
        idx = self.orfs.index(orf) + 1 if orf in self.orfs else 0
        self._dna_text.setPlainText(
            f">ORF{idx} DNA ({len(orf['dna'])} bp)\n{orf['dna']}")
        self._protein_text.setPlainText(
            f">ORF{idx} Protein ({len(orf['protein'].rstrip('*'))} aa)\n{orf['protein']}")
        if orf.get('domains'):
            txt = f"ORF{idx} Domains:\n"
            for d in orf['domains']:
                txt += f"  • {d['domain']} [{d.get('ali_region','')}] "
                txt += f"score={d.get('score',0):.1f} E={d.get('evalue',0):.1e}\n"
            self._domains_text.setPlainText(txt)
        else:
            self._domains_text.setPlainText(f"ORF{idx}: No domains annotated")

    def search_orf_in_map(self):
        query = self._map_search_edit.text().strip()
        if not query or not self.orfs: return
        m = re.match(r'[Oo][Rr][Ff](\d+)', query)
        if m:
            idx = int(m.group(1)) - 1
            if 0 <= idx < len(self.orfs):
                self._select_and_center_orf(idx)
        else:
            # Search by protein sequence
            query_up = query.upper()
            for i, orf in enumerate(self.orfs):
                if query_up in orf['protein']:
                    self._select_and_center_orf(i)
                    return

    def _select_and_center_orf(self, idx):
        """Select an ORF in table, highlight on map, and center the map view on it."""
        if idx < 0 or idx >= len(self.orfs): return
        orf = self.orfs[idx]
        self.selected_orf = orf
        self.selected_orf_idx = idx
        self._show_orf_details(orf)

        # Highlight on map and center
        self._genome_map.highlight_idx = idx
        # Correct centering formula:
        #   gw = full genome width in pixels at current zoom
        #   orf_pixel = absolute pixel position of ORF centre in full genome
        #   pan_offset = scroll so that orf_pixel appears at screen centre
        if self.dna_sequence:
            sl  = len(self.dna_sequence)
            w   = self._genome_map.width()
            if w < 80:
                # Widget not yet visible (hidden tab / minimised window).
                # Store the pending center request and execute it the next
                # time the genome map widget is resized/shown.
                self._pending_center_idx = idx
            else:
                mg  = 40
                bw  = w - 2 * mg
                gw  = int(bw * self.zoom_level)
                orf_pixel = int((orf['start'] + orf['end']) / 2 / sl * gw)
                new_pan = orf_pixel - (w // 2 - mg)
                self._genome_map.pan_offset = max(0, new_pan)
                self._genome_map._clamp_pan()
                self._pending_center_idx = -1  # fulfilled
        self._genome_map.update()

        # Select in table — set flag to block re-entrancy in _on_orf_table_select
        # Use index-based lookup (not identity) so annotations that modify
        # the orf dict in-place or replace it are handled correctly.
        self._orf_table_selecting = True
        try:
            # Prefer looking up by ORF index (position in self.orfs)
            target_label = f"ORF{idx+1}"
            row = -1
            for r in range(self._orf_table.rowCount()):
                cell = self._orf_table.item(r, 0)
                if cell is not None and cell.text().strip() == target_label:
                    row = r
                    break
            if row >= 0:
                self._orf_table.selectRow(row)
                self._orf_table.scrollTo(
                    self._orf_table.model().index(row, 0))
        finally:
            self._orf_table_selecting = False

        self._status.showMessage(
            f"✓ ORF{idx+1} selected ({orf['start']:,}–{orf['end']:,})")

    # ═══════════════════════════════════════════════════════════
    # RIGHT-CLICK CONTEXT MENU (ORF Table)
    # ═══════════════════════════════════════════════════════════

    def _on_orf_right_click(self, pos):
        """Right-click on ORF table.

        Key behaviour change (v2.0):
        - If the click lands on a row that is already part of a multi-selection
          (Ctrl+click), the multi-selection is preserved and the AF3 submenu
          shows group actions for all selected ORFs.
        - If the click lands on an unselected row, it selects that single row
          (existing behaviour).
        """
        row = self._orf_table.rowAt(pos.y())
        if row < 0 or row >= len(self.filtered_orfs):
            return

        # ── Preserve multi-selection if the click is inside it ──────────────
        sel_rows = sorted(set(
            idx.row() for idx in self._orf_table.selectedIndexes()))
        if row not in sel_rows:
            # Click is outside the current selection → select just this row
            self._orf_table.selectRow(row)
            sel_rows = [row]

        # Resolve the single ORF under the cursor (for single-item actions)
        orf      = self.filtered_orfs[row]
        orf_idx  = self.orfs.index(orf) if orf in self.orfs else -1
        n_sel    = len(sel_rows)

        menu = QMenu(self)

        # ════════════════════════════════════════════════════════
        # A) Single-ORF actions (annotation, color, copy)
        # ════════════════════════════════════════════════════════
        if n_sel == 1:
            menu.addAction(f"📝 Annotate ORF{orf_idx+1}…",
                           lambda: self._annotate_orf(orf, orf_idx))
            menu.addAction(f"🎨 Color ORF{orf_idx+1}…",
                           lambda: self._color_orf(orf, orf_idx))
            menu.addSeparator()
            menu.addAction("📋 Copy protein (FASTA)",
                           lambda: self._copy_to_clipboard(
                               f">ORF{orf_idx+1}|{orf['start']}-{orf['end']}\n"
                               f"{orf['protein'].rstrip('*')}"))
            menu.addAction("📋 Copy DNA",
                           lambda: self._copy_to_clipboard(orf['dna']))
            menu.addAction("📋 Copy raw protein sequence",
                           lambda: self._copy_to_clipboard(orf['protein'].rstrip('*')))
            menu.addSeparator()

        # ════════════════════════════════════════════════════════
        # B) AlphaFold 3 prediction actions — adapt to selection size
        # ════════════════════════════════════════════════════════
        menu.addSeparator()

        if n_sel == 1:
            # Single ORF — simple add + quick-predict sub-menu
            af3_sub = menu.addMenu("🔮 AlphaFold 3…")
            af3_sub.addAction(
                f"➕ Add ORF{orf_idx+1} to AF3 list",
                lambda: self._af3_add_orf_by_index(orf_idx))
            af3_sub.addSeparator()
            af3_sub.addAction(
                f"⚡ Add + Predict vs {self._af3_nb_spin.value()} neighbors",
                lambda: self._orf_table_af3_quick_neighbors([orf_idx]))
            af3_sub.addAction(
                "⚡ Add + Predict homodimer",
                lambda: self._orf_table_af3_quick_homodimer(orf_idx))

        else:
            # Multi-ORF — compute budget hints for display
            resolved = self._orf_table_resolve_sel_indices(sel_rows)
            n_orfs   = len(resolved)
            n_pairs  = n_orfs * (n_orfs - 1) // 2
            avg_res  = self._orf_table_avg_residues(resolved)
            # Smart mode suggestion
            suggestion = self._orf_table_af3_suggest_mode(resolved)

            # ── Budget preview header (disabled item used as label) ──────────
            hdr = QAction(
                f"🔮  {n_orfs} ORFs selected  →  {n_pairs} pairwise job(s)  "
                f"│  avg {avg_res} aa/chain", menu)
            hdr.setEnabled(False)
            menu.addAction(hdr)

            if suggestion:
                hint = QAction(f"💡 Suggested: {suggestion}", menu)
                hint.setEnabled(False)
                menu.addAction(hint)

            menu.addSeparator()

            # ── Multi-ORF actions ────────────────────────────────────────────
            menu.addAction(
                f"➕ Add {n_orfs} ORFs to AF3 list",
                lambda _r=resolved: self._orf_table_af3_add_multi(_r, switch_tab=False))

            menu.addAction(
                f"⚡ Add {n_orfs} ORFs + predict all pairs ({n_pairs} jobs)",
                lambda _r=resolved: self._orf_table_af3_predict_allvsall(_r))

            menu.addAction(
                f"⚡ Add {n_orfs} ORFs + predict vs {self._af3_nb_spin.value()} neighbors each",
                lambda _r=resolved: self._orf_table_af3_quick_neighbors(_r))

            # HMM sub-mode: only shown when ≥2 ORFs have HMM hits
            hmm_orfs = [i for i in resolved
                        if any(h.get('orf_index') == i for h in self.hmm_hits_all)]
            if len(hmm_orfs) >= 2:
                menu.addAction(
                    f"🧬 Predict {len(hmm_orfs)} HMM hits vs each other",
                    lambda _h=hmm_orfs: self._orf_table_af3_predict_allvsall(_h,
                        label="hmmhits"))

            menu.addSeparator()
            menu.addAction(
                f"➕ Add {n_orfs} ORFs to AF3 list (no generate)",
                lambda _r=resolved: self._orf_table_af3_add_multi(_r, switch_tab=True))

        # ════════════════════════════════════════════════════════
        # C) Table copy / export (always shown)
        # ════════════════════════════════════════════════════════
        menu.addSeparator()
        menu.addAction(
            f"Copy {n_sel} selected row(s) as TSV",
            self._orf_table_copy_selection)
        menu.addAction(
            "Copy all visible rows as TSV",
            lambda: self._orf_table_copy_rows(all_rows=True))
        menu.addSeparator()
        exp_sub = menu.addMenu("Export table…")
        exp_sub.addAction("TSV — table columns only",
            lambda: self._export_orf_table(fmt='tsv', include_seqs=False))
        exp_sub.addAction("TSV — full (+ DNA + Protein)",
            lambda: self._export_orf_table(fmt='tsv', include_seqs=True))
        exp_sub.addAction("TSV — annotated only",
            lambda: self._export_orf_table(fmt='tsv', annotated_only=True, include_seqs=True))
        exp_sub.addSeparator()
        exp_sub.addAction("FASTA — protein",
            lambda: self._export_orf_fasta(aa=True))
        exp_sub.addAction("FASTA — DNA",
            lambda: self._export_orf_fasta(aa=False))

        menu.exec(self._orf_table.viewport().mapToGlobal(pos))

    # ──────────────────────────────────────────────────────────────────────────
    # ORF-table → AF3  helper methods  (added v2.0)
    # ──────────────────────────────────────────────────────────────────────────

    def _switch_to_af3_tab(self):
        """Switch the right-panel tab widget to the AlphaFold prediction tab."""
        if hasattr(self, '_af3_tab_widget'):
            self._tabs.setCurrentWidget(self._af3_tab_widget)

    def _orf_table_resolve_sel_indices(self, sel_rows: list) -> list:
        """Convert selection rows (in filtered_orfs) to global ORF indices."""
        indices = []
        for r in sel_rows:
            if 0 <= r < len(self.filtered_orfs):
                orf = self.filtered_orfs[r]
                try:
                    idx = self.orfs.index(orf)
                    indices.append(idx)
                except ValueError:
                    pass
        return indices

    def _orf_table_avg_residues(self, orf_indices: list) -> int:
        """Return average protein length (aa) for a list of ORF indices."""
        if not orf_indices:
            return 0
        sizes = [len(self.orfs[i]['protein'].rstrip('*')) for i in orf_indices
                 if 0 <= i < len(self.orfs)]
        return int(sum(sizes) / len(sizes)) if sizes else 0

    def _orf_table_af3_suggest_mode(self, orf_indices: list) -> str:
        """Heuristically suggest the best AF3 mode for a given set of ORFs."""
        if not orf_indices or not self.orfs:
            return ""
        # Check if ORFs are genomic neighbors (sorted by position, max gap = N)
        orfs_by_pos = sorted(range(len(self.orfs)), key=lambda i: self.orfs[i]['start'])
        pos_ranks   = {idx: rank for rank, idx in enumerate(orfs_by_pos)}
        ranks = sorted(pos_ranks.get(i, -1) for i in orf_indices)
        max_gap  = max(ranks[k+1] - ranks[k] for k in range(len(ranks)-1)) if len(ranks) > 1 else 0
        n_nb     = self._af3_nb_spin.value()
        is_nbrs  = max_gap <= n_nb

        # Check if ORFs all have HMM hits
        hmm_count = sum(1 for i in orf_indices
                        if any(h.get('orf_index') == i for h in self.hmm_hits_all))

        if is_nbrs:
            return f"Pairs (Hit vs Neighbor) — ORFs are within {n_nb} genomic positions of each other"
        if hmm_count == len(orf_indices):
            return "HMM Hits vs Each Other — all selected ORFs have domain annotations"
        return "All vs All — mixed ORF set"

    def _orf_table_af3_add_multi(self, orf_indices: list, switch_tab: bool = True):
        """Add a list of ORF indices to the AF3 selection table.
        Silently skips duplicates. Optionally switches to the AF3 tab."""
        existing = set()
        for r in range(self._af3_sel_table.rowCount()):
            item = self._af3_sel_table.item(r, 0)
            if item:
                existing.add(item.text())
        added = 0
        for idx in orf_indices:
            if not (0 <= idx < len(self.orfs)):
                continue
            name = f"ORF{idx+1}"
            if name in existing:
                continue
            orf = self.orfs[idx]
            hmm_names = [p['name'] for p in self.hmm_profiles
                         for h in p.get('hits', [])
                         if h.get('orf_index') == idx]
            row = self._af3_sel_table.rowCount()
            self._af3_sel_table.insertRow(row)
            for col, val in enumerate([
                    name,
                    f"{orf['start']:,}-{orf['end']:,}",
                    str(len(orf['protein'].rstrip('*'))),
                    ', '.join(hmm_names) or '-',
                    '']):
                self._af3_sel_table.setItem(row, col, QTableWidgetItem(val))
            existing.add(name)
            added += 1
        self._af3_sel_count.setText(f"{self._af3_sel_table.rowCount()} ORFs selected")
        self._status.showMessage(f"✓ {added} ORF(s) added to AF3 list")
        if switch_tab:
            self._switch_to_af3_tab()

    def _orf_table_af3_predict_allvsall(self, orf_indices: list, label: str = "selected"):
        """Add ORFs to AF3 list, generate all-vs-all pairwise jobs, switch tab."""
        if not self.orfs:
            QMessageBox.warning(self, "AF3", "Run ORF analysis first!")
            return
        self._orf_table_af3_add_multi(orf_indices, switch_tab=False)
        existing_names = {j['name'] for j in self.af3_jobs}
        added = 0
        seen_pairs: set = set()
        for i_a in range(len(orf_indices)):
            for i_b in range(i_a + 1, len(orf_indices)):
                hi = orf_indices[i_a]; pi = orf_indices[i_b]
                key = (min(hi, pi), max(hi, pi))
                if key in seen_pairs:
                    continue
                seen_pairs.add(key)
                hn = f"ORF{hi+1}"; pn = f"ORF{pi+1}"
                job_name = f"{hn}_vs_{pn}_{label}"
                if job_name in existing_names:
                    continue
                hp = self.orfs[hi]['protein'].rstrip('*')
                pp = self.orfs[pi]['protein'].rstrip('*')
                tr = len(hp) + len(pp)
                self.af3_jobs.append({
                    'name':            job_name,
                    'hit_orf_idx':     hi,
                    'partner_orf_idx': pi,
                    'hit_name':        hn,
                    'partner_name':    pn,
                    'total_residues':  tr,
                    'paeinter':        None,
                    'status':          ('pending' if tr <= self.af3_max_residues
                                        else f'>{self.af3_max_residues}!'),
                    'iptm':   None,
                    'plddt':  None,
                    'sequences': [
                        {'proteinChain': {'sequence': hp, 'count': 1}},
                        {'proteinChain': {'sequence': pp, 'count': 1}}],
                })
                existing_names.add(job_name)
                added += 1
        self._af3_update_jobs_table()
        self._switch_to_af3_tab()
        self._status.showMessage(
            f"✓ {added} pairwise job(s) added from {len(orf_indices)} selected ORFs "
            f"→ AlphaFold tab")

    def _orf_table_af3_quick_neighbors(self, orf_indices: list):
        """Add ORFs to AF3 list, generate Pairs-vs-Neighbor jobs, switch tab."""
        if not self.orfs:
            QMessageBox.warning(self, "AF3", "Run ORF analysis first!")
            return
        self._orf_table_af3_add_multi(orf_indices, switch_tab=False)
        n_nb = self._af3_nb_spin.value()
        orfs_by_pos = sorted(enumerate(self.orfs), key=lambda x: x[1]['start'])
        pos_to_rank = {idx: rank for rank, (idx, _) in enumerate(orfs_by_pos)}
        existing_names = {j['name'] for j in self.af3_jobs}
        added = 0
        for hi in orf_indices:
            if not (0 <= hi < len(self.orfs)):
                continue
            ho = self.orfs[hi]; hp = ho['protein'].rstrip('*'); hn = f"ORF{hi+1}"
            hr = pos_to_rank.get(hi, 0)
            for d in range(-n_nb, n_nb + 1):
                if d == 0:
                    continue
                nr = hr + d
                if not (0 <= nr < len(orfs_by_pos)):
                    continue
                ni, no = orfs_by_pos[nr]
                np_s = no['protein'].rstrip('*')
                nn   = f"ORF{ni+1}"
                tr   = len(hp) + len(np_s)
                job_name = f"{hn}_vs_{nn}_{'up' if d < 0 else 'dn'}{abs(d)}_tbl"
                if job_name in existing_names:
                    continue
                self.af3_jobs.append({
                    'name':            job_name,
                    'hit_orf_idx':     hi,
                    'partner_orf_idx': ni,
                    'hit_name':        hn,
                    'partner_name':    nn,
                    'total_residues':  tr,
                    'paeinter':        None,
                    'status':          ('pending' if tr <= self.af3_max_residues
                                        else f'>{self.af3_max_residues}!'),
                    'iptm':  None,
                    'plddt': None,
                    'sequences': [
                        {'proteinChain': {'sequence': hp,   'count': 1}},
                        {'proteinChain': {'sequence': np_s, 'count': 1}}],
                })
                existing_names.add(job_name)
                added += 1
        self._af3_update_jobs_table()
        self._switch_to_af3_tab()
        self._status.showMessage(
            f"✓ {added} neighbor jobs from {len(orf_indices)} ORF(s), N={n_nb} → AlphaFold tab")

    def _orf_table_af3_quick_homodimer(self, orf_idx: int):
        """Add one ORF to AF3 list and generate its homodimer job, switch tab."""
        if not (0 <= orf_idx < len(self.orfs)):
            return
        self._orf_table_af3_add_multi([orf_idx], switch_tab=False)
        orf  = self.orfs[orf_idx]
        hp   = orf['protein'].rstrip('*')
        hn   = f"ORF{orf_idx+1}"
        name = f"{hn}_homodimer_tbl"
        if any(j['name'] == name for j in self.af3_jobs):
            self._switch_to_af3_tab()
            return
        tr = len(hp) * 2
        self.af3_jobs.append({
            'name':            name,
            'hit_orf_idx':     orf_idx,
            'partner_orf_idx': orf_idx,
            'hit_name':        hn,
            'partner_name':    hn,
            'total_residues':  tr,
            'paeinter':        None,
            'status':          ('pending' if tr <= self.af3_max_residues
                                else f'>{self.af3_max_residues}!'),
            'iptm':  None,
            'plddt': None,
            'sequences': [{'proteinChain': {'sequence': hp, 'count': 2}}],
        })
        self._af3_update_jobs_table()
        self._switch_to_af3_tab()
        self._status.showMessage(f"✓ Homodimer job added for {hn} → AlphaFold tab")

    def _annotate_orf(self, orf, orf_idx):
        """Dialog to annotate an ORF with observation, function, gene name."""
        dlg = QDialog(self)
        dlg.setWindowTitle(f"📝 Annotate ORF{orf_idx+1}")
        dlg.setFixedSize(450, 320)
        layout = QVBoxLayout(dlg)

        layout.addWidget(QLabel(f"<b>ORF{orf_idx+1}</b> ({orf['start']:,}-{orf['end']:,}) | "
                                 f"{len(orf['protein'].rstrip('*'))} aa | GC: {orf['gc']:.1f}%"))

        form = QGridLayout()
        form.addWidget(QLabel("Observation:"), 0, 0)
        obs_edit = QLineEdit(orf.get('observation', ''))
        form.addWidget(obs_edit, 0, 1)

        form.addWidget(QLabel("Putative Function:"), 1, 0)
        func_edit = QLineEdit(orf.get('putative_function', ''))
        form.addWidget(func_edit, 1, 1)

        form.addWidget(QLabel("Gene Name:"), 2, 0)
        gene_edit = QLineEdit(orf.get('gene_name', ''))
        form.addWidget(gene_edit, 2, 1)

        form.addWidget(QLabel("Notes:"), 3, 0)
        notes_edit = QTextEdit()
        notes_edit.setMaximumHeight(80)
        notes_edit.setPlainText(orf.get('notes', ''))
        form.addWidget(notes_edit, 3, 1)
        layout.addLayout(form)

        btns = QDialogButtonBox(QDialogButtonBox.StandardButton.Save | QDialogButtonBox.StandardButton.Cancel
                                 if QT_VERSION == 6 else QDialogButtonBox.Save | QDialogButtonBox.Cancel)
        btns.accepted.connect(dlg.accept)
        btns.rejected.connect(dlg.reject)
        layout.addWidget(btns)

        if dlg.exec() if QT_VERSION == 6 else dlg.exec_():
            orf['observation'] = obs_edit.text()
            orf['putative_function'] = func_edit.text()
            orf['gene_name'] = gene_edit.text()
            orf['notes'] = notes_edit.toPlainText()
            self._update_orfs_list()
            self._status.showMessage(f"✓ ORF{orf_idx+1} annotated")

    def _color_orf(self, orf, orf_idx):
        """Dialog to set a custom color for an ORF."""
        dlg = QDialog(self)
        dlg.setWindowTitle(f"🎨 Color ORF{orf_idx+1}")
        dlg.setFixedSize(350, 180)
        layout = QVBoxLayout(dlg)
        layout.addWidget(QLabel(f"<b>Color for ORF{orf_idx+1}</b>"))

        current = orf.get('custom_color', '#B2BEC3')
        color_var = [current]

        # Color palette
        palette = QHBoxLayout()
        for c in ['#E53935','#1E88E5','#43A047','#FB8C00','#8E24AA',
                   '#00ACC1','#D81B60','#FFD600','#B2BEC3','#795548']:
            btn = QPushButton()
            btn.setFixedSize(26, 26)
            btn.setStyleSheet(f"background-color: {c}; border: 1px solid #888;")
            btn.clicked.connect(lambda checked, col=c: (color_var.__setitem__(0, col),
                                                          hex_edit.setText(col)))
            palette.addWidget(btn)
        layout.addLayout(palette)

        hex_lay = QHBoxLayout()
        hex_lay.addWidget(QLabel("Hex:"))
        hex_edit = QLineEdit(current)
        hex_edit.setMaximumWidth(100)
        hex_lay.addWidget(hex_edit)
        hex_lay.addStretch()
        layout.addLayout(hex_lay)

        btns = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
                                 if QT_VERSION == 6 else QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(dlg.accept)
        btns.rejected.connect(dlg.reject)
        layout.addWidget(btns)

        if dlg.exec() if QT_VERSION == 6 else dlg.exec_():
            orf['custom_color'] = hex_edit.text()
            self._update_map()
            self._status.showMessage(f"✓ ORF{orf_idx+1} colored: {hex_edit.text()}")

    # ═══════════════════════════════════════════════════════════
    # HELPER / STUB METHODS
    # ═══════════════════════════════════════════════════════════

    def _copy_to_clipboard(self, text):
        QApplication.clipboard().setText(text)
        self._status.showMessage(f"✓ Copied ({len(text)} chars)")

    def _paste_clipboard(self):
        text = QApplication.clipboard().text()
        if text: self._blast_query_text.setPlainText(text)

    def _clear_query(self):
        self._blast_query_text.clear()

    def _validate_query(self):
        raw = self._blast_query_text.toPlainText().strip()
        qp = self._parse_fasta_query(raw)
        if qp:
            self._query_status.setText(f"✓ Valid protein: {len(qp)} aa")
            self._query_status.setStyleSheet("color: green;")
        else:
            self._query_status.setText("❌ Invalid sequence (DNA detected or empty)")
            self._query_status.setStyleSheet("color: red;")

    def load_blast_query_fasta(self):
        f, _ = QFileDialog.getOpenFileName(self, "Load Query FASTA",
            "", "FASTA (*.fasta *.fa *.faa);;All (*)")
        if f:
            with open(f, encoding='utf-8') as fh: self._blast_query_text.setPlainText(fh.read())

    def _copy_blast_hit(self):
        if self.selected_orf:
            self._copy_to_clipboard(self.selected_orf['protein'])

    def _copy_blast_all(self):
        self._copy_to_clipboard(self._blast_results_text.toPlainText())

    def _save_blast_results(self):
        f, _ = QFileDialog.getSaveFileName(self, "Save BLAST Results", "", "Text (*.txt)")
        if f:
            with open(f, 'w', encoding='utf-8') as fh: fh.write(self._blast_results_text.toPlainText())

    def _analyze_neighborhood(self):
        if self.selected_orf_idx < 0:
            QMessageBox.warning(self, "Warning", "Select an ORF first!"); return
        idx = self.selected_orf_idx
        orf = self.orfs[idx]
        wk = self._window_kb_spin.value()
        result = self.analyzer.analyze_neighborhood(self.orfs, idx, wk)
        score = self.analyzer.score_candidate(orf, result)
        orf['candidate_score'] = score

        ws = result['window_start']; we = result['window_end']

        txt = f"""GENOMIC NEIGHBORHOOD ANALYSIS
{'='*70}
Target ORF: ORF{idx+1} ({orf['start']:,}-{orf['end']:,}) | Frame: F{orf['frame']}{orf['strand']}
Protein: {orf['protein'][:80]}...
Window: {wk} kb ({ws:,}-{we:,})

CLUSTER:
  ORFs in window: {result['total_orfs_in_window']}
  Hypothetical: {result['hypothetical_proteins']}
  Cluster score: {result['cluster_score']:.2f}
  Candidate score: {score:.2f}

DOMAINS IN NEIGHBORHOOD:
"""
        if result['domains_found']:
            for d, c in result['domains_found'].items(): txt += f"  {d}: {c}x\n"
        else: txt += "  (none)\n"
        txt += "\nSYSTEMS:\n"
        if result['systems_found']:
            for s, c in result['systems_found'].items(): txt += f"  {s}: {c}\n"
        else: txt += "  (none)\n"

        # Neighboring ORFs table
        txt += "\nNEIGHBORING ORFs:\n"
        txt += f"{'ORF':<10} {'Dist(bp)':<10} {'Strand':<7} {'Size(aa)':<9} {'HMM':<25} {'Obs'}\n"
        txt += "-"*75 + "\n"
        for n in result['neighbors'][:30]:
            ni = n['orf_index']
            ds = ', '.join(d['domain'] for d in n['domains']) or '-'
            obs_n = self.orfs[ni].get('observation', '') if ni < len(self.orfs) else ''
            txt += f"ORF{ni+1:<7} {n['distance_bp']:<10,} {n['strand']:<7} {n['length_aa']:<8} {ds:<25} {obs_n}\n"

        # DNA sequence of the genomic window
        txt += f"\n{'='*70}\n"
        txt += f"GENOMIC WINDOW DNA SEQUENCE ({ws:,}-{we:,} = {we-ws:,} bp)\n"
        txt += f"{'='*70}\n"
        if self.dna_sequence:
            dna_window = self.dna_sequence[ws:min(we, len(self.dna_sequence))]
            for i in range(0, len(dna_window), 80):
                txt += dna_window[i:i+80] + "\n"

        # Protein sequences (FASTA) for all ORFs in window
        txt += f"\n{'='*70}\n"
        txt += "ORF PROTEIN SEQUENCES IN WINDOW (FASTA)\n"
        txt += f"{'='*70}\n"

        # Target ORF first
        target_prot = orf['protein'].rstrip('*')
        txt += f">ORF{idx+1}|TARGET|F{orf['frame']}{orf['strand']}|{orf['start']:,}-{orf['end']:,}|{len(target_prot)}aa"
        hmm_t = ', '.join(d['domain'] for d in orf.get('domains', []))
        if hmm_t: txt += f"|{hmm_t}"
        txt += "\n"
        for j in range(0, len(target_prot), 80):
            txt += target_prot[j:j+80] + "\n"

        # Neighbor proteins
        for n in result['neighbors'][:30]:
            ni = n['orf_index']
            if ni >= len(self.orfs): continue
            norf = self.orfs[ni]
            nprot = norf['protein'].rstrip('*')
            ds = ', '.join(d['domain'] for d in n['domains'])
            txt += f">ORF{ni+1}|F{norf['frame']}{norf['strand']}|{norf['start']:,}-{norf['end']:,}|{len(nprot)}aa"
            if ds: txt += f"|{ds}"
            txt += "\n"
            for j in range(0, len(nprot), 80):
                txt += nprot[j:j+80] + "\n"

        self._neighborhood_text.setPlainText(txt)

    def _add_hmm_profile(self):
        f, _ = QFileDialog.getOpenFileName(self, "Add HMM Profile", "", "HMM (*.hmm)")
        if f: self._add_hmm_file(f)

    def _add_hmm_multi(self):
        files, _ = QFileDialog.getOpenFileNames(self, "Add Multiple HMM Profiles",
            "", "HMM files (*.hmm);;All files (*)")
        if files:
            for f in files:
                self._add_hmm_file(f)
            self._status.showMessage(f"✓ {len(self.hmm_profiles)} HMM profile(s) loaded")

    def _add_hmm_from_folder(self):
        """Scan a folder (and subfolders) for all .hmm files."""
        d = QFileDialog.getExistingDirectory(self, "Select folder with HMM profiles")
        if not d:
            return
        found = list(Path(d).rglob('*.hmm'))
        if not found:
            QMessageBox.information(self, "HMM",
                f"No .hmm files found in:\n{d}\n\n"
                f"(searched subfolders too)")
            return
        for f in sorted(found):
            self._add_hmm_file(str(f))
        self._hmm_text.append(f"\n✓ Found {len(found)} .hmm files in {Path(d).name}/")
        self._status.showMessage(f"✓ {len(self.hmm_profiles)} HMM profile(s) loaded")

    def _hmm_search_all(self):
        if not self.orfs:
            QMessageBox.warning(self, "HMM", "Run ORF analysis first!"); return
        if not self.hmm_profiles:
            QMessageBox.warning(self, "HMM", "Add HMM profiles first!"); return
        self._status.showMessage("⏳ Running HMM search...")
        self.hmm_hits_all = []

        def work():
            all_hits = []
            errors = []
            params = {'hmm_evalue': self.hmm_evalue, 'hmm_dom_evalue': self.hmm_dom_evalue}
            for profile in self.hmm_profiles:
                hits = self.analyzer.hmm_scan_orfs(profile['file'], self.orfs, params)
                error_hits = [h for h in hits if 'error' in h]
                valid_hits = [h for h in hits if 'error' not in h]
                if error_hits:
                    errors.append(f"{profile['name']}: {error_hits[0]['error']}")
                for h in valid_hits:
                    h['profile_name'] = profile['name']
                    h['profile_color'] = profile['color']
                    h['profile_function'] = profile.get('function', '')
                profile['hits'] = valid_hits
                all_hits.extend(valid_hits)
            return (all_hits, errors)

        def done(result):
            all_hits, errors = result
            self.hmm_hits_all = all_hits
            self._show_hmm_results(errors)
            self._update_map()
            self._status.showMessage(f"✓ HMM search: {len(all_hits)} hits")

        self._run_worker(work, done)

    def _show_hmm_results(self, errors=None):
        """Display HMM results with colored alignments in HTML."""
        method = "HMMER3" if BACKENDS.get('hmmer3', {}).get('available') else "PSSM (built-in)"
        wsl_tag = " (via WSL)" if BACKENDS.get('hmmer3', {}).get('wsl') else ""
        n_hits = len(self.hmm_hits_all)

        html = []
        html.append("<pre style='font-family: Courier New, monospace; font-size: 9pt;'>")
        html.append(f"<b style='color:#0D47A1;'>HMM SEARCH — {method}{wsl_tag}</b>")
        html.append(f"\n{'='*70}")
        html.append(f"\n<span style='color:#424242;'>Profiles: {len(self.hmm_profiles)}  |  "
                     f"Database: {len(self.orfs)} ORFs  |  Total hits: {n_hits}</span>")
        html.append(f"\n{'='*70}\n")

        if errors:
            for err in errors:
                html.append(f"\n<span style='color:red;'>⚠ {err}</span>")
            html.append("\n")

        for profile in self.hmm_profiles:
            hits = [h for h in profile.get('hits', []) if 'error' not in h]
            color = profile.get('color', '#333')

            html.append(f"\n<b style='color:{color};'>┌─ {profile['name']}  "
                         f"({Path(profile['file']).name})</b>")
            html.append(f"\n<span style='color:#424242;'>│  Color: {color}  |  "
                         f"Function: {profile.get('function', '-')}</span>")
            html.append(f"\n<span style='color:#424242;'>│  Hits: {len(hits)}</span>")

            if hits:
                html.append("\n│")
                html.append(f"\n<b style='color:#1B5E20;'>│  {'#':<3} {'ORF':<9} {'Score':<8} "
                             f"{'E-value':<12} {'HMM region':<18} {'Protein region'}</b>")
                html.append(f"\n<span style='color:#9E9E9E;'>│  {'─'*70}</span>")

                for idx_h, h in enumerate(hits[:30]):
                    ev = h.get('evalue', 999)
                    ev_s = f"{ev:.1e}" if ev < 0.01 else f"{ev:.3f}"
                    oi = h.get('orf_index', -1)
                    hmm_from = h.get('hmm_from', '?')
                    hmm_to = h.get('hmm_to', '?')
                    hmm_len = h.get('hmm_len', '?')
                    ali_from = h.get('ali_from', '?')
                    ali_to = h.get('ali_to', '?')
                    tgt_len = h.get('target_len', '?')
                    hmm_reg = f"{hmm_from}-{hmm_to}/{hmm_len}"
                    prot_reg = f"{ali_from}-{ali_to}/{tgt_len}"

                    html.append(f"\n│  {idx_h+1:<3} "
                                 f"<a href='orf:{oi}' style='color:#0D47A1;font-weight:bold;'>"
                                 f"{h.get('orf_name','?')}</a>    "
                                 f"{h['score']:<8.1f} {ev_s:<12} {hmm_reg:<18} {prot_reg}")

                    # Show alignment if available
                    aln_hmm = h.get('aln_hmm', '')
                    aln_tgt = h.get('aln_target', '')
                    aln_mid = h.get('aln_match', '')

                    if aln_hmm and aln_tgt:
                        html.append("\n│")
                        block_w = 60
                        h_pos = int(hmm_from) if str(hmm_from).isdigit() else 1
                        t_pos = int(ali_from) if str(ali_from).isdigit() else 1

                        for b in range(0, len(aln_hmm), block_w):
                            hb = aln_hmm[b:b+block_w]
                            mb = aln_mid[b:b+block_w] if aln_mid else ' ' * len(hb)
                            tb = aln_tgt[b:b+block_w]
                            h_adv = len(hb.replace('-', '').replace('.', ''))
                            t_adv = len(tb.replace('-', ''))
                            h_end = h_pos + h_adv - 1
                            t_end = t_pos + t_adv - 1

                            # HMM consensus line
                            hmm_colored = []
                            for ch in hb:
                                if ch == '-':
                                    hmm_colored.append(f"<span style='color:#B71C1C;'>{ch}</span>")
                                elif ch.isupper():
                                    hmm_colored.append(f"<span style='background:#E8F5E9;color:#1B5E20;'>{ch}</span>")
                                else:
                                    hmm_colored.append(f"<span style='background:#FFF3E0;color:#E65100;'>{ch}</span>")
                            html.append(f"\n│  <span style='color:#424242;'>HMM    {h_pos:>5}  </span>"
                                         f"{''.join(hmm_colored)}"
                                         f"<span style='color:#424242;'>  {h_end}</span>")

                            # Match line
                            mid_colored = []
                            for ch in mb:
                                if ch == '|':
                                    mid_colored.append(f"<span style='color:#1B5E20;font-weight:bold;'>{ch}</span>")
                                elif ch == '+':
                                    mid_colored.append(f"<span style='color:#E65100;'>{ch}</span>")
                                else:
                                    mid_colored.append(f"<span style='color:#BDBDBD;'>{ch}</span>")
                            html.append(f"\n│         {'':>5}  {''.join(mid_colored)}")

                            # Target protein line
                            tgt_colored = []
                            for ch in tb:
                                if ch == '-':
                                    tgt_colored.append(f"<span style='color:#B71C1C;'>{ch}</span>")
                                else:
                                    tgt_colored.append(f"<span style='color:#424242;'>{ch}</span>")
                            html.append(f"\n│  <span style='color:#424242;'>Prot   {t_pos:>5}  </span>"
                                         f"{''.join(tgt_colored)}"
                                         f"<span style='color:#424242;'>  {t_end}</span>")
                            html.append("\n│")

                            h_pos = h_end + 1
                            t_pos = t_end + 1

            html.append(f"\n<span style='color:#9E9E9E;'>└{'─'*70}</span>\n")

        html.append("</pre>")
        self._hmm_text.setHtml(''.join(html))
        self._update_hmm_profile_table()

    # ═══════ AF3 METHODS ═══════

    def _af3_add_selected(self):
        if not self.selected_orf or self.selected_orf_idx < 0:
            QMessageBox.information(self, "AF3", "Select an ORF from the table first."); return
        self._af3_add_orf_by_index(self.selected_orf_idx)

    def _af3_add_orf_by_index(self, orf_idx):
        if orf_idx < 0 or orf_idx >= len(self.orfs): return
        orf = self.orfs[orf_idx]
        for r in range(self._af3_sel_table.rowCount()):
            if self._af3_sel_table.item(r, 0) and self._af3_sel_table.item(r, 0).text() == f"ORF{orf_idx+1}":
                self._status.showMessage(f"ORF{orf_idx+1} already in AF3 list"); return
        hmm_names = [p['name'] for p in self.hmm_profiles for h in p.get('hits',[]) if h.get('orf_index')==orf_idx]
        row = self._af3_sel_table.rowCount(); self._af3_sel_table.insertRow(row)
        for col, val in enumerate([f"ORF{orf_idx+1}", f"{orf['start']:,}-{orf['end']:,}",
                                    str(len(orf['protein'].rstrip('*'))), ', '.join(hmm_names) or '-', '']):
            self._af3_sel_table.setItem(row, col, QTableWidgetItem(val))
        self._af3_sel_count.setText(f"{self._af3_sel_table.rowCount()} ORFs selected")

    def _af3_add_hmm_hits(self):
        if not self.hmm_hits_all:
            QMessageBox.information(self, "AF3", "No HMM hits. Run HMM search first."); return
        existing = set()
        for r in range(self._af3_sel_table.rowCount()):
            if self._af3_sel_table.item(r,0): existing.add(self._af3_sel_table.item(r,0).text())
        added = 0
        for hit in self.hmm_hits_all:
            oi = hit.get('orf_index',-1)
            if oi < 0 or oi >= len(self.orfs): continue
            name = f"ORF{oi+1}"
            if name in existing: continue
            orf = self.orfs[oi]; row = self._af3_sel_table.rowCount(); self._af3_sel_table.insertRow(row)
            for col, val in enumerate([name, f"{orf['start']:,}-{orf['end']:,}",
                                        str(len(orf['protein'].rstrip('*'))), hit.get('profile_name','?'),
                                        f"score={hit.get('score',0):.1f}"]):
                self._af3_sel_table.setItem(row, col, QTableWidgetItem(val))
            existing.add(name); added += 1
        self._af3_sel_count.setText(f"{self._af3_sel_table.rowCount()} ORFs selected")
        self._status.showMessage(f"✓ {added} HMM hits added")

    def _af3_add_all_orfs(self):
        """Add ALL predicted ORFs in the genome to the AF3 selection list.

        This enables a genome-wide interactome scan: every ORF is a candidate
        for pair-wise structural interaction prediction with AlphaFold 3.
        A size-filter dialog is shown first so users can cap the total number
        of jobs that would be generated downstream.
        """
        if not self.orfs:
            QMessageBox.information(self, "Add All ORFs",
                "No ORFs found. Run ORF analysis (DNA tab) first."); return

        n_total = len(self.orfs)

        # ── Confirmation / filter dialog ──────────────────────────────────
        dlg = QDialog(self)
        dlg.setWindowTitle("🧬 Add All ORFs — Genome-wide Interactome")
        dlg.setMinimumWidth(420)
        dl = QVBoxLayout(dlg)

        info = QLabel(
            f"<b>{n_total} ORFs</b> detected in the genome.<br>"
            "All will be added to the AF3 selection list for genome-wide<br>"
            "interactome scanning.<br><br>"
            f"<i>Tip: Mode <b>Genomic Interactome</b> generates up to "
            f"<b>{n_total*(n_total-1)//2:,}</b> pairs.<br>"
            f"Mode <b>Neighbors Interactome</b> (sliding window N) is "
            f"much lighter: <b>N·n − N(N+1)/2</b> pairs — e.g. "
            f"{5*n_total - 15:,} pairs at N=5.<br>"
            "Apply size filters below to keep the job count manageable.</i>")
        info.setWordWrap(True)
        dl.addWidget(info)

        # Min / max size filters
        fbox = QHBoxLayout()
        fbox.addWidget(QLabel("Min size (aa):"))
        min_spin = QSpinBox(); min_spin.setRange(1, 10000); min_spin.setValue(30)
        fbox.addWidget(min_spin)
        fbox.addSpacing(12)
        fbox.addWidget(QLabel("Max size (aa):"))
        max_spin = QSpinBox(); max_spin.setRange(1, 100000); max_spin.setValue(2000)
        fbox.addWidget(max_spin)
        dl.addLayout(fbox)

        # HMM-only checkbox
        hmm_only_cb = QCheckBox("Only ORFs with HMM / BLAST annotation")
        hmm_only_cb.setToolTip(
            "When checked, only ORFs that have at least one HMM hit or BLAST "
            "annotation are added — dramatically reduces job count while keeping "
            "functionally characterised proteins.")
        dl.addWidget(hmm_only_cb)

        # Live counter label
        counter_lbl = QLabel()
        counter_lbl.setStyleSheet("color: #1565c0; font-weight: bold;")
        dl.addWidget(counter_lbl)

        def _update_count():
            mn = min_spin.value(); mx = max_spin.value()
            hmm_set = set()
            if hmm_only_cb.isChecked():
                for h in self.hmm_hits_all:
                    hmm_set.add(h.get('orf_index', -1))
                for bi in range(self._blast_results_table.rowCount()
                                if hasattr(self, '_blast_results_table') else 0):
                    try:
                        idx = int(self._blast_results_table.item(bi, 0).text().replace('ORF',''))-1
                        hmm_set.add(idx)
                    except Exception:
                        pass
            count = 0
            for i, o in enumerate(self.orfs):
                sz = len(o['protein'].rstrip('*'))
                if sz < mn or sz > mx: continue
                if hmm_only_cb.isChecked() and i not in hmm_set: continue
                count += 1
            pairs = count*(count-1)//2
            counter_lbl.setText(
                f"→ {count} ORFs will be added  |  "
                f"All-vs-All would produce {pairs:,} job pairs")

        min_spin.valueChanged.connect(_update_count)
        max_spin.valueChanged.connect(_update_count)
        hmm_only_cb.stateChanged.connect(_update_count)
        _update_count()

        btns = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
            if QT_VERSION == 6 else
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(dlg.accept)
        btns.rejected.connect(dlg.reject)
        dl.addWidget(btns)

        if dlg.exec() != QDialog.DialogCode.Accepted if QT_VERSION == 6 else dlg.exec() != QDialog.Accepted:
            return

        mn = min_spin.value(); mx = max_spin.value()
        hmm_only = hmm_only_cb.isChecked()

        # Build set of annotated ORF indices
        annotated = set()
        if hmm_only:
            for h in self.hmm_hits_all:
                annotated.add(h.get('orf_index', -1))

        # Collect existing entries to avoid duplicates
        existing = set()
        for r in range(self._af3_sel_table.rowCount()):
            if self._af3_sel_table.item(r, 0):
                existing.add(self._af3_sel_table.item(r, 0).text())

        added = 0
        self._af3_sel_table.setUpdatesEnabled(False)
        try:
            for i, orf in enumerate(self.orfs):
                sz = len(orf['protein'].rstrip('*'))
                if sz < mn or sz > mx: continue
                if hmm_only and i not in annotated: continue
                name = f"ORF{i+1}"
                if name in existing: continue

                # Collect HMM annotations for display
                hmm_names = []
                for hit in self.hmm_hits_all:
                    if hit.get('orf_index') == i:
                        pn = hit.get('profile_name', hit.get('hmm_name', ''))
                        if pn and pn not in hmm_names:
                            hmm_names.append(pn)

                row = self._af3_sel_table.rowCount()
                self._af3_sel_table.insertRow(row)
                for col, val in enumerate([
                        name,
                        f"{orf['start']:,}-{orf['end']:,}",
                        str(sz),
                        ', '.join(hmm_names) or '-',
                        'genome-wide']):
                    self._af3_sel_table.setItem(row, col, QTableWidgetItem(val))
                existing.add(name)
                added += 1
        finally:
            self._af3_sel_table.setUpdatesEnabled(True)

        total_sel = self._af3_sel_table.rowCount()
        self._af3_sel_count.setText(f"{total_sel} ORFs selected")
        self._status.showMessage(
            f"✓ {added} ORFs added (genome-wide) — "
            f"{total_sel} total | Use mode 'Genomic Interactome' to generate all pairs")

        # Auto-switch mode combo to Interactoma Genômico for convenience
        idx = self._af3_mode_combo.findText("Genomic Interactome (Selected ORFs vs All ORFs)")
        if idx >= 0:
            self._af3_mode_combo.setCurrentIndex(idx)

    def _af3_remove_orf(self):
        rows = sorted(set(idx.row() for idx in self._af3_sel_table.selectedIndexes()), reverse=True)
        for r in rows: self._af3_sel_table.removeRow(r)
        self._af3_sel_count.setText(f"{self._af3_sel_table.rowCount()} ORFs selected")

    def _af3_clear_all(self):
        self._af3_sel_table.setRowCount(0)
        self._af3_jobs_table.setRowCount(0)
        self.af3_jobs = []
        self._af3_sel_count.setText("0 ORFs selected")

    def _af3_generate_jobs(self):
        n_sel = self._af3_sel_table.rowCount()
        if n_sel == 0: QMessageBox.warning(self, "AF3", "Select ORFs first!"); return
        if not self.orfs: QMessageBox.warning(self, "AF3", "Run ORF analysis first!"); return
        n_nb = self._af3_nb_spin.value(); mode = self._af3_mode_combo.currentText(); self.af3_jobs = []
        sel_indices = []
        for r in range(n_sel):
            try:
                idx = self._parse_orf_idx_from_text(self._af3_sel_table.item(r, 0).text())
                if 0 <= idx < len(self.orfs): sel_indices.append(idx)
            except (AttributeError, ValueError, TypeError): continue
        orfs_by_pos = sorted(enumerate(self.orfs), key=lambda x: x[1]['start'])
        pos_to_rank = {idx: rank for rank, (idx, _) in enumerate(orfs_by_pos)}

        # ── Genomic Interactome: each SELECTED ORF vs ALL genome ORFs ──────
        if mode.startswith("Genomic Interactome"):
            n_sel_orfs = len(sel_indices)
            n_genome   = len(self.orfs)
            estimated  = n_sel_orfs * (n_genome - 1)
            if estimated > 5000:
                ans = QMessageBox.question(
                    self, "Genomic Interactome — large job set",
                    f"This will generate up to <b>{estimated:,}</b> pairwise AF3 jobs<br>"
                    f"(<b>{n_sel_orfs}</b> selected ORF(s) × <b>{n_genome}</b> genome ORFs).<br><br>"
                    "Symmetric duplicates (A↔B) are automatically collapsed.<br>"
                    "This may produce very large output files. Continue?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
                    if QT_VERSION == 6 else QMessageBox.Yes | QMessageBox.No)
                if ans != (QMessageBox.StandardButton.Yes if QT_VERSION == 6 else QMessageBox.Yes):
                    return
            seen_pairs = set()
            for hi in sel_indices:
                ho = self.orfs[hi]; hp = ho['protein'].rstrip('*'); hn = f"ORF{hi+1}"
                for gj, go in enumerate(self.orfs):
                    if gj == hi: continue
                    pair_key = (min(hi, gj), max(hi, gj))
                    if pair_key in seen_pairs: continue
                    seen_pairs.add(pair_key)
                    gp = go['protein'].rstrip('*'); gn = f"ORF{gj+1}"; tr = len(hp) + len(gp)
                    self.af3_jobs.append({
                        'name': f"{hn}_vs_{gn}_interactome",
                        'hit_orf_idx': hi, 'partner_orf_idx': gj,
                        'hit_name': hn, 'partner_name': gn,
                        'total_residues': tr, 'paeinter': None,
                        'status': 'pending' if tr <= self.af3_max_residues else f'>{self.af3_max_residues}!',
                        'iptm': None, 'plddt': None,
                        'sequences': [
                            {'proteinChain': {'sequence': hp, 'count': 1}},
                            {'proteinChain': {'sequence': gp, 'count': 1}}]})
            self._af3_update_jobs_table()
            self._status.showMessage(
                f"✓ {len(self.af3_jobs)} interactome jobs — "
                f"{n_sel_orfs} selected ORF(s) vs {n_genome} genome ORFs "
                f"({len(seen_pairs)} unique pairs)")
            return

        # ── Neighbors Interactome: genome-wide sliding window ───────────────
        if mode.startswith("Neighbors Interactome"):
            ordered   = [idx for idx, _ in orfs_by_pos]
            n_gen     = len(ordered)
            estimated = sum(1 for ra in range(n_gen)
                            for d in range(1, n_nb+1) if ra+d < n_gen)
            if estimated > 5000:
                ans = QMessageBox.question(
                    self, "Neighbors Interactome — large job set",
                    f"This will generate <b>{estimated:,}</b> pairwise AF3 jobs "
                    f"over the entire genome<br>"
                    f"(<b>{n_gen}</b> ORFs, window N=<b>{n_nb}</b>).<br><br>"
                    "Symmetric duplicates (A↔B) are already removed.<br>Continue?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
                    if QT_VERSION == 6 else QMessageBox.Yes | QMessageBox.No)
                if ans != (QMessageBox.StandardButton.Yes if QT_VERSION == 6 else QMessageBox.Yes):
                    return
            seen_pairs = set()
            for rank_a in range(n_gen):
                ia = ordered[rank_a]; oa = self.orfs[ia]
                sa = oa['protein'].rstrip('*'); na = f"ORF{ia+1}"
                for d in range(1, n_nb+1):
                    rank_b = rank_a + d
                    if rank_b >= n_gen: break
                    ib = ordered[rank_b]; ob = self.orfs[ib]
                    sb = ob['protein'].rstrip('*'); nb = f"ORF{ib+1}"
                    key = (ia, ib) if ia < ib else (ib, ia)
                    if key in seen_pairs: continue
                    seen_pairs.add(key)
                    tr = len(sa) + len(sb)
                    self.af3_jobs.append({
                        'name': f"{na}_vs_{nb}_nbr{d}",
                        'hit_orf_idx': ia, 'partner_orf_idx': ib,
                        'hit_name': na, 'partner_name': nb,
                        'total_residues': tr, 'paeinter': None,
                        'status': ('pending' if tr <= self.af3_max_residues
                                   else f'>{self.af3_max_residues}!'),
                        'iptm': None, 'plddt': None,
                        'sequences': [
                            {'proteinChain': {'sequence': sa, 'count': 1}},
                            {'proteinChain': {'sequence': sb, 'count': 1}}]})
            self._af3_update_jobs_table()
            self._status.showMessage(
                f"✓ {len(self.af3_jobs)} Neighbors Interactome jobs — "
                f"genome-wide sliding window N={n_nb} "
                f"({len(seen_pairs)} unique pairs over {n_gen} ORFs)")
            return

        # ── Selected vs Selected (Ctrl+click highlighted rows) ──────────────
        if mode.startswith("Selected vs Selected"):
            highlighted = sorted(set(
                idx.row() for idx in self._af3_sel_table.selectedIndexes()))
            hl_indices  = []
            for r in highlighted:
                try:
                    idx = self._parse_orf_idx_from_text(self._af3_sel_table.item(r, 0).text())
                    if 0 <= idx < len(self.orfs): hl_indices.append(idx)
                except (AttributeError, ValueError, TypeError): continue
            if len(hl_indices) < 2:
                QMessageBox.information(
                    self, "AF3 — Selected vs Selected",
                    "Ctrl+click at least 2 ORFs in the selection list above first.")
                return
            seen_pairs = set()
            for i_a in range(len(hl_indices)):
                for i_b in range(i_a + 1, len(hl_indices)):
                    hi = hl_indices[i_a]; pi = hl_indices[i_b]
                    key = (min(hi, pi), max(hi, pi))
                    if key in seen_pairs: continue
                    seen_pairs.add(key)
                    hn = f"ORF{hi+1}"; pn = f"ORF{pi+1}"
                    hp = self.orfs[hi]['protein'].rstrip('*')
                    pp = self.orfs[pi]['protein'].rstrip('*')
                    tr = len(hp) + len(pp)
                    self.af3_jobs.append({
                        'name': f"{hn}_vs_{pn}_selected",
                        'hit_orf_idx': hi, 'partner_orf_idx': pi,
                        'hit_name': hn, 'partner_name': pn,
                        'total_residues': tr, 'paeinter': None,
                        'status': ('pending' if tr <= self.af3_max_residues
                                   else f'>{self.af3_max_residues}!'),
                        'iptm': None, 'plddt': None,
                        'sequences': [
                            {'proteinChain': {'sequence': hp, 'count': 1}},
                            {'proteinChain': {'sequence': pp, 'count': 1}}]})
            self._af3_update_jobs_table()
            self._status.showMessage(
                f"✓ {len(self.af3_jobs)} job(s) from {len(hl_indices)} Ctrl+selected ORFs")
            return

        # ── All vs All (all ORFs in selection list) ──────────────────────────
        if mode.startswith("All vs All"):
            n_s = len(sel_indices); est = n_s * (n_s - 1) // 2
            if est > 5000:
                ans = QMessageBox.question(
                    self, "All vs All — large job set",
                    f"This will generate <b>{est:,}</b> pairwise jobs "
                    f"from <b>{n_s}</b> selected ORFs. Continue?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
                    if QT_VERSION == 6 else QMessageBox.Yes | QMessageBox.No)
                if ans != (QMessageBox.StandardButton.Yes if QT_VERSION == 6 else QMessageBox.Yes):
                    return
            seen_pairs = set()
            for i_a in range(len(sel_indices)):
                for i_b in range(i_a + 1, len(sel_indices)):
                    hi = sel_indices[i_a]; pi = sel_indices[i_b]
                    key = (min(hi, pi), max(hi, pi))
                    if key in seen_pairs: continue
                    seen_pairs.add(key)
                    hn = f"ORF{hi+1}"; pn = f"ORF{pi+1}"
                    hp = self.orfs[hi]['protein'].rstrip('*')
                    pp = self.orfs[pi]['protein'].rstrip('*')
                    tr = len(hp) + len(pp)
                    self.af3_jobs.append({
                        'name': f"{hn}_vs_{pn}_allvsall",
                        'hit_orf_idx': hi, 'partner_orf_idx': pi,
                        'hit_name': hn, 'partner_name': pn,
                        'total_residues': tr, 'paeinter': None,
                        'status': ('pending' if tr <= self.af3_max_residues
                                   else f'>{self.af3_max_residues}!'),
                        'iptm': None, 'plddt': None,
                        'sequences': [
                            {'proteinChain': {'sequence': hp, 'count': 1}},
                            {'proteinChain': {'sequence': pp, 'count': 1}}]})
            self._af3_update_jobs_table()
            self._status.showMessage(f"✓ {len(self.af3_jobs)} All vs All jobs generated")
            return

        # ── HMM Hits vs Each Other ───────────────────────────────────────────
        if mode.startswith("HMM Hits vs Each Other"):
            hmm_indices = [i for i in sel_indices
                           if any(h.get('orf_index') == i for h in self.hmm_hits_all)]
            if len(hmm_indices) < 2:
                QMessageBox.information(
                    self, "AF3 — HMM Hits vs Each Other",
                    "Need at least 2 ORFs with HMM hits in the selection list.\n"
                    "Run HMM search first, then use 'Add HMM Hits'.")
                return
            seen_pairs = set()
            for i_a in range(len(hmm_indices)):
                for i_b in range(i_a + 1, len(hmm_indices)):
                    hi = hmm_indices[i_a]; pi = hmm_indices[i_b]
                    key = (min(hi, pi), max(hi, pi))
                    if key in seen_pairs: continue
                    seen_pairs.add(key)
                    hn = f"ORF{hi+1}"; pn = f"ORF{pi+1}"
                    hp = self.orfs[hi]['protein'].rstrip('*')
                    pp = self.orfs[pi]['protein'].rstrip('*')
                    tr = len(hp) + len(pp)
                    self.af3_jobs.append({
                        'name': f"{hn}_vs_{pn}_hmmhits",
                        'hit_orf_idx': hi, 'partner_orf_idx': pi,
                        'hit_name': hn, 'partner_name': pn,
                        'total_residues': tr, 'paeinter': None,
                        'status': ('pending' if tr <= self.af3_max_residues
                                   else f'>{self.af3_max_residues}!'),
                        'iptm': None, 'plddt': None,
                        'sequences': [
                            {'proteinChain': {'sequence': hp, 'count': 1}},
                            {'proteinChain': {'sequence': pp, 'count': 1}}]})
            self._af3_update_jobs_table()
            self._status.showMessage(
                f"✓ {len(self.af3_jobs)} jobs from {len(hmm_indices)} HMM-hit ORFs")
            return

        # ── Hit vs All Selected ──────────────────────────────────────────────
        if mode.startswith("Hit vs All Selected"):
            if len(sel_indices) < 2:
                QMessageBox.information(
                    self, "AF3 — Hit vs All Selected",
                    "Need at least 2 ORFs in the selection list.\n"
                    "The topmost row is used as the query ORF.")
                return
            query_idx = sel_indices[0]
            qn = f"ORF{query_idx+1}"
            qp = self.orfs[query_idx]['protein'].rstrip('*')
            for pi in sel_indices[1:]:
                pn = f"ORF{pi+1}"; pp = self.orfs[pi]['protein'].rstrip('*')
                tr = len(qp) + len(pp)
                self.af3_jobs.append({
                    'name': f"{qn}_vs_{pn}_hitall",
                    'hit_orf_idx': query_idx, 'partner_orf_idx': pi,
                    'hit_name': qn, 'partner_name': pn,
                    'total_residues': tr, 'paeinter': None,
                    'status': ('pending' if tr <= self.af3_max_residues
                               else f'>{self.af3_max_residues}!'),
                    'iptm': None, 'plddt': None,
                    'sequences': [
                        {'proteinChain': {'sequence': qp, 'count': 1}},
                        {'proteinChain': {'sequence': pp, 'count': 1}}]})
            self._af3_update_jobs_table()
            self._status.showMessage(
                f"✓ {len(self.af3_jobs)} jobs: {qn} vs {len(sel_indices)-1} selected ORFs")
            return

        # ── Pairs (Hit vs Neighbor) — default neighbor-walk mode ────────────
        for hi in sel_indices:
            ho = self.orfs[hi]; hr = pos_to_rank.get(hi, 0)
            hp = ho['protein'].rstrip('*'); hn = f"ORF{hi+1}"
            nbs = []
            for d in range(-n_nb, n_nb + 1):
                if d == 0: continue
                nr = hr + d
                if 0 <= nr < len(orfs_by_pos):
                    ni, no = orfs_by_pos[nr]; nbs.append((ni, no, d))
            for ni, no, d in nbs:
                np_s = no['protein'].rstrip('*'); tr = len(hp) + len(np_s)
                self.af3_jobs.append({
                    'name': f"{hn}_vs_ORF{ni+1}_{'up' if d < 0 else 'down'}{abs(d)}",
                    'hit_orf_idx': hi, 'partner_orf_idx': ni,
                    'hit_name': hn, 'partner_name': f"ORF{ni+1}",
                    'total_residues': tr, 'paeinter': None,
                    'status': 'pending' if tr <= self.af3_max_residues else f'>{self.af3_max_residues}!',
                    'iptm': None, 'plddt': None,
                    'sequences': [
                        {'proteinChain': {'sequence': hp, 'count': 1}},
                        {'proteinChain': {'sequence': np_s, 'count': 1}}]})

        # ── Homodimer checkbox — appended after any mode ─────────────────────
        if self._af3_homodimer_cb.isChecked():
            existing_names = {j['name'] for j in self.af3_jobs}
            for hi in sel_indices:
                ho = self.orfs[hi]; hp = ho['protein'].rstrip('*'); hn = f"ORF{hi+1}"
                hname = f"{hn}_homodimer"
                if hname in existing_names: continue
                tr2 = len(hp) * 2
                self.af3_jobs.append({
                    'name': hname,
                    'hit_orf_idx': hi, 'partner_orf_idx': hi,
                    'hit_name': hn, 'partner_name': hn,
                    'total_residues': tr2, 'paeinter': None,
                    'status': 'pending' if tr2 <= self.af3_max_residues else f'>{self.af3_max_residues}!',
                    'iptm': None, 'plddt': None,
                    'sequences': [{'proteinChain': {'sequence': hp, 'count': 2}}]})

        self._af3_update_jobs_table()
        self._status.showMessage(f"✓ {len(self.af3_jobs)} AF3 jobs generated")

    def _af3_update_jobs_table(self):
        self._af3_jobs_table.setRowCount(0)
        for j in self.af3_jobs:
            row = self._af3_jobs_table.rowCount()
            self._af3_jobs_table.insertRow(row)
            iptm_s   = f"{j['iptm']:.3f}" if j.get('iptm') is not None else '-'
            paeinter = j.get('paeinter')
            pae_s    = f"{paeinter:.1f}" if paeinter is not None else '-'
            # Confidence classification (ipTM > 0.75 AND PAEinter < 8 Å = HIGH)
            iptm_val = j.get('iptm')
            if iptm_val is not None and paeinter is not None:
                if iptm_val >= 0.75 and paeinter < 8.0:
                    conf_str = "HIGH ★"
                    conf_color = "#1b5e20"
                    conf_bg    = "#e8f5e9"
                elif iptm_val >= 0.5 or paeinter < 15.0:
                    conf_str = "MED"
                    conf_color = "#e65100"
                    conf_bg    = "#fff3e0"
                else:
                    conf_str = "LOW"
                    conf_color = "#b71c1c"
                    conf_bg    = "#ffebee"
            else:
                conf_str = "-"
                conf_color = None
                conf_bg    = None

            vals = [j['name'], j['hit_name'], j['partner_name'],
                    str(j['total_residues']), j['status'],
                    iptm_s, pae_s, conf_str]
            high_confidence = (iptm_val is not None and paeinter is not None and
                                iptm_val >= 0.75 and paeinter < 8.0)
            for col, val in enumerate(vals):
                item = QTableWidgetItem(val)
                if high_confidence:
                    f = item.font(); f.setBold(True); item.setFont(f)
                if col == 7 and conf_bg:
                    try:
                        from PyQt6.QtGui import QColor as _QColor
                    except ImportError:
                        from PyQt5.QtGui import QColor as _QColor
                    item.setBackground(_QColor(conf_bg))
                    item.setForeground(_QColor(conf_color))
                self._af3_jobs_table.setItem(row, col, item)

    def _af3_export_json(self):
        if not self.af3_jobs: QMessageBox.warning(self,"AF3","Generate jobs first!"); return
        folder = QFileDialog.getExistingDirectory(self, "Select output folder")
        if not folder: return
        errors = []
        for j in self.af3_jobs:
            safe_name = re.sub(r'[^\w.\-]', '_', j['name'])
            data = {"name": j['name'], "modelSeeds": [], "sequences": j['sequences'],
                    "dialect": "alphafoldserver", "version": 2}
            try:
                with open(os.path.join(folder, f"{safe_name}.json"), 'w',
                          encoding='utf-8') as f:
                    json.dump(data, f, indent=2, ensure_ascii=False)
            except OSError as e:
                errors.append(f"{safe_name}.json: {e}")
        if errors:
            QMessageBox.warning(self, "AF3 Export",
                f"{len(errors)} file(s) could not be written:\n" + "\n".join(errors))
        self._status.showMessage(
            f"✓ {len(self.af3_jobs) - len(errors)} AF3 JSONs exported")

    def _af3_export_json_batch(self):
        """Export all AF3 jobs as a single batch JSON file for AlphaFold3."""
        if not self.af3_jobs: 
            QMessageBox.warning(self, "AF3", "Generate jobs first!")
            return
        
        # Ask for single file instead of folder
        default_name = f"{self.genome_name or 'batch'}_af3_jobs.json"
        filepath, _ = QFileDialog.getSaveFileName(
            self, "Export AF3 Batch JSON", default_name, "JSON (*.json);;All (*)")
        if not filepath:
            return

        # Create batch array - each job as individual AF3 format
        batch_jobs = []
        for j in self.af3_jobs:
            job_data = {
                "name": j['name'],
                "modelSeeds": [],
                "sequences": j['sequences'],
                "dialect": "alphafoldserver", 
                "version": 1  # AF3 API uses version 1 for batch
            }
            batch_jobs.append(job_data)
        
        try:
            with open(filepath, 'w', encoding='utf-8') as f:
                json.dump(batch_jobs, f, indent=2, ensure_ascii=False)
            
            file_size = Path(filepath).stat().st_size / 1024  # KB
            self._status.showMessage(
                f"✓ AF3 batch exported: {len(self.af3_jobs)} jobs → "
                f"{Path(filepath).name} ({file_size:.0f} KB)")
        except OSError as e:
            QMessageBox.critical(self, "Export Error", f"Failed to write batch file:\n{e}")

    def _af3_export_colabfold(self):
        if not self.af3_jobs: QMessageBox.warning(self,"AF3","Generate jobs first!"); return
        folder = QFileDialog.getExistingDirectory(self, "Select output folder")
        if not folder: return
        csv_lines = ["name,sequence"]
        errors = []
        for j in self.af3_jobs:
            safe_name = re.sub(r'[^\w.\-]', '_', j['name'])
            chains = [s['proteinChain']['sequence'] for s in j['sequences']]
            joined = ':'.join(chains)
            csv_lines.append(f"{safe_name},{joined}")
            try:
                with open(os.path.join(folder, f"{safe_name}.fasta"), 'w',
                          encoding='utf-8') as f:
                    f.write(f">{j['name']}\n{joined}\n")
            except OSError as e:
                errors.append(f"{safe_name}.fasta: {e}")
        safe_genome = re.sub(r'[^\w.\-]', '_', self.genome_name or 'batch')
        try:
            with open(os.path.join(folder, f"{safe_genome}_batch.csv"), 'w',
                      encoding='utf-8') as f:
                f.write('\n'.join(csv_lines))
        except OSError as e:
            errors.append(f"{safe_genome}_batch.csv: {e}")
        if errors:
            QMessageBox.warning(self, "ColabFold Export",
                f"{len(errors)} file(s) could not be written:\n" + "\n".join(errors))
        self._status.showMessage(f"✓ ColabFold exported: {len(self.af3_jobs)} jobs")

    def _af3_export_slurm_array(self):
        """Export AF3 jobs as numbered JSON batches + a ready-to-submit SLURM array script.

        One sbatch command submits all batches automatically — no manual loop needed.
        Each array task runs its own AF3 process, so RAM is fully released between batches
        (prevents OUT_OF_MEMORY on large interactome scans).

        Output layout:
            <output_dir>/
                batches/
                    batch_001/  job_001.json  job_002.json  …  (batch_size jobs)
                    batch_002/  …
                    …
                run_array.sh    ← sbatch run_array.sh  (one command!)
                submit_all.sh   ← fallback sequential submitter
        """
        if not self.af3_jobs:
            QMessageBox.warning(self, "SLURM Array Export", "Generate AF3 jobs first!"); return

        # ── Config dialog ──────────────────────────────────────────────────
        dlg = QDialog(self)
        dlg.setWindowTitle("⚡ Export SLURM Array — anti-OOM batch splitter")
        dlg.setMinimumWidth(480)
        dl = QVBoxLayout(dlg)

        n_jobs = len(self.af3_jobs)
        info = QLabel(
            f"<b>{n_jobs} AF3 jobs</b> will be split into batches.<br>"
            "Each batch runs as a SLURM array task — its own process,<br>"
            "RAM fully released between batches. <b>One sbatch command.</b>")
        info.setWordWrap(True)
        dl.addWidget(info)

        grid = QGridLayout()
        grid.setSpacing(6)

        grid.addWidget(QLabel("Jobs per batch (batch size):"), 0, 0)
        batch_spin = QSpinBox(); batch_spin.setRange(1, 500); batch_spin.setValue(50)
        batch_spin.setToolTip("50 is safe for most clusters. Reduce for large proteins.")
        grid.addWidget(batch_spin, 0, 1)

        n_batches_lbl = QLabel()
        grid.addWidget(n_batches_lbl, 0, 2)

        grid.addWidget(QLabel("RAM per task (--mem):"), 1, 0)
        mem_edit = QLineEdit("64G")
        grid.addWidget(mem_edit, 1, 1)

        grid.addWidget(QLabel("Time per task (--time):"), 2, 0)
        time_edit = QLineEdit("7-00:00:00")
        time_edit.setToolTip("SLURM format: D-HH:MM:SS\nDaVinci: basic=3d | max50=8d | max90=15d")
        grid.addWidget(time_edit, 2, 1)

        grid.addWidget(QLabel("Partition:"), 3, 0)
        part_edit = QLineEdit("max50")
        part_edit.setToolTip("DaVinci partitions:\n  basic : 72h  | 16 CPUs | 100 GB | 0 GPU\n  max50 : 8d   | 64 CPUs | 500 GB | 1 GPU  (recomendada para AF3)\n  max90 : 15d  | 110 CPUs| 1 TB   | 4 GPUs (jobs muito grandes)")
        grid.addWidget(part_edit, 3, 1)

        grid.addWidget(QLabel("GPUs per task (--gres):"), 4, 0)
        gpu_edit = QLineEdit("gpu:1")
        grid.addWidget(gpu_edit, 4, 1)

        grid.addWidget(QLabel("CPUs per task (--cpus):"), 5, 0)
        cpu_spin = QSpinBox(); cpu_spin.setRange(1, 64); cpu_spin.setValue(16)
        cpu_spin.setToolTip("DaVinci max50: até 64 CPUs | max90: até 110 CPUs")
        grid.addWidget(cpu_spin, 5, 1)

        grid.addWidget(QLabel("AF3 command on cluster:"), 6, 0)
        af3cmd_edit = QLineEdit()
        try: af3cmd_edit.setText(self._dv_af3cmd.text() or "af3_run")
        except Exception: af3cmd_edit.setText("af3_run")
        af3cmd_edit.setToolTip("AF3 command or path on the server")
        grid.addWidget(af3cmd_edit, 6, 1, 1, 2)

        grid.addWidget(QLabel("Base folder on cluster:"), 7, 0)
        remote_edit = QLineEdit()
        try: remote_edit.setText(self._dv_base_path.text() or "~/af3_predictions")
        except Exception: remote_edit.setText("~/af3_predictions")
        grid.addWidget(remote_edit, 7, 1, 1, 2)

        dl.addLayout(grid)

        # Live batch counter
        def _upd():
            bs = batch_spin.value()
            nb = (n_jobs + bs - 1) // bs
            n_batches_lbl.setText(f"→ {nb} arrays")
        batch_spin.valueChanged.connect(_upd); _upd()

        btns = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
            if QT_VERSION == 6 else QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(dlg.accept); btns.rejected.connect(dlg.reject)
        dl.addWidget(btns)

        if (dlg.exec() if QT_VERSION == 6 else dlg.exec_()) != (
                QDialog.DialogCode.Accepted if QT_VERSION == 6 else QDialog.Accepted):
            return

        batch_size  = batch_spin.value()
        mem         = mem_edit.text().strip() or "32G"
        walltime    = time_edit.text().strip() or "7-00:00:00"
        partition   = part_edit.text().strip() or "max50"
        gres        = gpu_edit.text().strip() or "gpu:1"
        ncpus       = cpu_spin.value()
        af3cmd      = af3cmd_edit.text().strip() or "af3_run"
        remote_base = remote_edit.text().strip() or "~/af3_predictions"

        # ── Choose output folder ────────────────────────────────────────────
        out_dir = QFileDialog.getExistingDirectory(self, "Select output folder for SLURM array export")
        if not out_dir: return
        out_path = Path(out_dir)
        batches_dir = out_path / "batches"

        # ── Split jobs into batches ─────────────────────────────────────────
        batches = []
        for i in range(0, len(self.af3_jobs), batch_size):
            batches.append(self.af3_jobs[i:i + batch_size])
        n_batches = len(batches)

        prog_dlg = None
        try:
            from PyQt6.QtWidgets import QProgressDialog
        except ImportError:
            try:
                from PyQt5.QtWidgets import QProgressDialog
            except ImportError:
                pass

        if QProgressDialog:
            prog_dlg = QProgressDialog(
                f"Writing {n_jobs} JSONs into {n_batches} batches...", "Cancel", 0, n_batches, self)
            prog_dlg.setWindowTitle("SLURM Array Export")
            prog_dlg.setMinimumDuration(0); prog_dlg.setValue(0)
            QApplication.processEvents()

        try:
            for bi, batch in enumerate(batches):
                if prog_dlg and prog_dlg.wasCanceled(): return
                bdir = batches_dir / f"batch_{bi+1:03d}"
                bdir.mkdir(parents=True, exist_ok=True)
                for j in batch:
                    jdata = {"name": j['name'], "modelSeeds": [],
                             "sequences": j.get('sequences', []),
                             "dialect": "alphafoldserver", "version": 2}
                    safe_name = re.sub(r'[^\w\-.]', '_', j['name'])
                    with open(bdir / f"{safe_name}.json", 'w', encoding='utf-8') as fh:
                        json.dump(jdata, fh, indent=2, ensure_ascii=False)
                if prog_dlg: prog_dlg.setValue(bi + 1); QApplication.processEvents()
        finally:
            if prog_dlg: prog_dlg.close()

        # ── Generate run_array.sh ───────────────────────────────────────────
        genome_safe = re.sub(r'[^\w\-]', '_', self.genome_name or 'interactome')
        array_script = f"""#!/bin/bash
#SBATCH --job-name={genome_safe}_array
#SBATCH --array=1-{n_batches}
#SBATCH --mem={mem}
#SBATCH --time={walltime}
#SBATCH --partition={partition}
#SBATCH --gres={gres}
#SBATCH --cpus-per-task={ncpus}
#SBATCH --output={remote_base}/logs/array_%A_%a.out
#SBATCH --error={remote_base}/logs/array_%A_%a.err

# ─────────────────────────────────────────────────────────────────────────────
# ppigFinder SLURM Array — {genome_safe}
# {n_jobs} AF3 jobs / {batch_size} per batch = {n_batches} array tasks
# Submit with:  sbatch run_array.sh
# Each task runs independently → RAM fully released between batches (anti-OOM)
# ─────────────────────────────────────────────────────────────────────────────

BATCH_DIR="{remote_base}/batches/batch_$(printf '%03d' $SLURM_ARRAY_TASK_ID)"
OUT_DIR="{remote_base}/results/batch_$(printf '%03d' $SLURM_ARRAY_TASK_ID)"
mkdir -p "$OUT_DIR"
mkdir -p "{remote_base}/logs"

echo "[$(date)] Starting array task $SLURM_ARRAY_TASK_ID / {n_batches}"
echo "  Batch dir : $BATCH_DIR"
echo "  Output dir: $OUT_DIR"
echo "  Node      : $(hostname)"
echo "  GPU       : $CUDA_VISIBLE_DEVICES"

# Run AF3 for every JSON in this batch
for json_file in "$BATCH_DIR"/*.json; do
    job_name=$(basename "$json_file" .json)
    job_out="$OUT_DIR/$job_name"
    mkdir -p "$job_out"
    echo "  → Predicting: $job_name"
    {af3cmd} \
        --json_path="$json_file" \
        --output_dir="$job_out"
done

echo "[$(date)] Task $SLURM_ARRAY_TASK_ID complete."
"""

        array_path = out_path / "run_array.sh"
        with open(array_path, 'w', encoding='utf-8') as fh:
            fh.write(array_script)

        # ── Generate submit_all.sh (sequential fallback) ────────────────────
        seq_script = f"""#!/bin/bash
# Sequential fallback — submits each batch as a separate job (not an array).
# Prefer run_array.sh when your cluster supports job arrays.
for i in $(seq -w 1 {n_batches}); do
    sbatch --job-name={genome_safe}_b${{i}} \
           --mem={mem} --time={walltime} \
           --partition={partition} --gres={gres} \
           --cpus-per-task={ncpus} \
           --wrap="{af3cmd} --json_path={remote_base}/batches/batch_${{i}} --output_dir={remote_base}/results/batch_${{i}}"
done
echo "Submitted {n_batches} jobs."
"""
        with open(out_path / "submit_all.sh", 'w', encoding='utf-8') as fh:
            fh.write(seq_script)

        # ── Summary ─────────────────────────────────────────────────────────
        total_kb = sum(
            f.stat().st_size for f in batches_dir.rglob("*.json")
        ) // 1024

        self._status.showMessage(
            f"✓ SLURM array: {n_jobs} jobs → {n_batches} batches × {batch_size} | "
            f"run_array.sh ready — sbatch run_array.sh to submit all at once")

        _msg = "\n".join([
            "Export completo!",
            "",
            f"  Pasta: {out_path}",
            f"  batches/  ({n_batches} pastas, {total_kb:,} KB de JSONs)",
            "  run_array.sh   <- sbatch run_array.sh",
            "  submit_all.sh  <- fallback sequencial",
            "",
            f"  {n_jobs} jobs | {batch_size} por batch | {n_batches} arrays",
            f"  RAM per task: {mem} | Time: {walltime}",
            "",
            "Como submeter (UM unico comando):",
            "  1. Copie a pasta 'batches/' para o cluster",
            "  2. sbatch run_array.sh",
            "",
            f"O SLURM dispara os {n_batches} tasks automaticamente.",
            "RAM released between each batch — no OOM risk.",
        ])
        QMessageBox.information(self, "SLURM Array Export", _msg)

    def _af3_show_ranking(self):
        """Show ipTM ranking in a floating dialog (no longer uses the removed _af3_text panel)."""
        done = [j for j in self.af3_jobs if j.get('iptm') is not None]
        if not done:
            QMessageBox.information(self, "AF3 Ranking", "No results imported yet.")
            return
        ranked = sorted(done, key=lambda x: x['iptm'], reverse=True)
        txt = f"RANKING — ipTM\n{'='*62}\nipTM > 0.75 + PAEinter < 8 Å = HIGH confidence\n\n"
        txt += f"{'#':<4}  {'ipTM':>6}  {'pLDDT':>6}  {'Hit':<12}  {'Partner':<12}  {'Job name'}\n"
        txt += "-" * 62 + "\n"
        for i, j in enumerate(ranked):
            p  = j.get('plddt', 0) or 0
            pi = j.get('paeinter')
            stars = "★★★" if j['iptm'] >= 0.8 else "★★" if j['iptm'] >= 0.6 else "★" if j['iptm'] >= 0.4 else ""
            pae_str = f"  PAEi={pi:.1f}" if pi is not None else ""
            txt += (f"{i+1:<4}  {j['iptm']:>6.3f}  {p:>6.1f}"
                    f"  {j['hit_name']:<12}  {j['partner_name']:<12}"
                    f"  {j['name']} {stars}{pae_str}\n")

        dlg = QDialog(self)
        dlg.setWindowTitle(f"AF3 Ranking — {len(ranked)} result(s)")
        dlg.setMinimumSize(700, 420)
        dl = QVBoxLayout(dlg)
        te = QTextEdit()
        te.setFont(QFont('Courier New', 9))
        te.setReadOnly(True)
        te.setPlainText(txt)
        dl.addWidget(te)
        btns = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Close if QT_VERSION == 6
            else QDialogButtonBox.Close)
        btns.rejected.connect(dlg.reject)
        dl.addWidget(btns)
        if QT_VERSION == 6:
            dlg.exec()
        else:
            dlg.exec_()

    def _af3_clear_jobs(self):
        self._af3_jobs_table.setRowCount(0)
        self.af3_jobs = []

    # ──────────────────────────────────────────────────────────
    # Dynamic multi-subunit custom job helpers
    # ──────────────────────────────────────────────────────────
    _CHAIN_LETTERS = list("ABCDEFGHIJK")   # A … K  (11 subunidades max)

    def _af3_rebuild_custom_rows(self, n_subunits: int | None = None):
        """Rebuild the per-subunit input rows (ORF field + n= spinbox) for
        the given number of subunits.  Preserves existing field values."""
        if n_subunits is None:
            n_subunits = self._custom_n_subunits.value()

        # Snapshot current values so we can restore them after rebuild
        prev_vals = [(row[0].text(), row[1].value())
                     for row in self._custom_subunit_rows]

        # Remove all existing widgets from the layout
        while self._custom_rows_layout.count():
            item = self._custom_rows_layout.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()
        self._custom_subunit_rows.clear()

        # Build n_subunits rows
        for i in range(n_subunits):
            letter = self._CHAIN_LETTERS[i] if i < len(self._CHAIN_LETTERS) else str(i + 1)
            row_widget = QWidget()
            row_layout = QHBoxLayout(row_widget)
            row_layout.setContentsMargins(0, 0, 0, 0)
            row_layout.setSpacing(3)

            lbl = QLabel(f"  <b>{letter}</b>:")
            lbl.setTextFormat(Qt.TextFormat.RichText if QT_VERSION == 6
                              else Qt.RichText)
            lbl.setMinimumWidth(22)
            row_layout.addWidget(lbl)

            orf_edit = QLineEdit()
            placeholder = f"ORF{i+1}"
            orf_edit.setPlaceholderText(placeholder)
            orf_edit.setMinimumWidth(70)
            orf_edit.setMaximumWidth(100)
            # restore previous value if it existed
            if i < len(prev_vals):
                orf_edit.setText(prev_vals[i][0])
            row_layout.addWidget(orf_edit)

            row_layout.addWidget(QLabel("n="))

            n_spin = QSpinBox()
            n_spin.setRange(1, 20)
            n_spin.setValue(prev_vals[i][1] if i < len(prev_vals) else 1)
            n_spin.setMaximumWidth(50)
            n_spin.setToolTip(t("tip_chain_copies").format(letter=letter))
            row_layout.addWidget(n_spin)

            row_layout.addStretch()
            self._custom_rows_layout.addWidget(row_widget)
            self._custom_subunit_rows.append((orf_edit, n_spin))

        self._custom_rows_layout.addStretch()

    def _af3_add_custom_job(self):
        """Add a manually defined custom job with N subunits and stoichiometry."""

        def _parse_orf_idx(text: str) -> int:
            """Return 0-based ORF index from 'ORF123' text, or -1 on failure."""
            m = re.match(r'ORF(\d+)', text.strip().upper())
            return int(m.group(1)) - 1 if m else -1

        if not self._custom_subunit_rows:
            QMessageBox.warning(self, "Custom Job", "No subunits defined.")
            return

        # ── Validate & collect all subunits ──
        subunits = []    # list of (orf_idx, n_copies, chain_letter)
        errors = []
        for i, (orf_edit, n_spin) in enumerate(self._custom_subunit_rows):
            letter = self._CHAIN_LETTERS[i] if i < len(self._CHAIN_LETTERS) else str(i + 1)
            raw = orf_edit.text().strip()
            if not raw:
                errors.append(f"Chain {letter}: ORF field is empty.")
                continue
            idx = _parse_orf_idx(raw)
            if idx < 0 or idx >= len(self.orfs):
                errors.append(f"Chain {letter}: '{raw}' not found in ORF list.")
                continue
            subunits.append((idx, n_spin.value(), letter))

        if errors:
            QMessageBox.warning(self, "Custom Job", "\n".join(errors))
            return
        if not subunits:
            QMessageBox.warning(self, "Custom Job", "No valid subunits defined.")
            return

        # ── Build AF3 sequences list (merge identical ORFs) ──
        sequences = []
        seen: dict[int, dict] = {}   # orf_idx → existing sequence entry
        total_residues = 0
        for orf_idx, n_copies, letter in subunits:
            prot = self.orfs[orf_idx]['protein'].rstrip('*')
            total_residues += len(prot) * n_copies
            if orf_idx in seen:
                # AlphaFold3 supports multiple entries for the same sequence;
                # add a new entry (allows independent chain IDs & stoichiometry).
                sequences.append({'proteinChain': {'sequence': prot, 'count': n_copies}})
            else:
                entry = {'proteinChain': {'sequence': prot, 'count': n_copies}}
                sequences.append(entry)
                seen[orf_idx] = entry

        # ── Build job name  e.g.  ORF1x2_ORF3x1_ORF5x3_custom ──
        name_parts = [f"ORF{idx+1}x{n}" for idx, n, _ in subunits]
        name = "_".join(name_parts) + "_custom"

        # ── Summary string for status bar ──
        chain_summary = ", ".join(
            f"{lt}=ORF{idx+1}×{n}" for idx, n, lt in subunits)

        # Use first and last ORF as hit/partner for table display
        hit_idx    = subunits[0][0]
        part_idx   = subunits[-1][0]

        self.af3_jobs.append({
            'name':            name,
            'hit_orf_idx':     hit_idx,
            'partner_orf_idx': part_idx,
            'hit_name':        f"ORF{hit_idx+1}",
            'partner_name':    f"ORF{part_idx+1}",
            'total_residues':  total_residues,
            'status':          ('pending'
                                if total_residues <= self.af3_max_residues
                                else f'>{self.af3_max_residues}!'),
            'iptm':            None,
            'plddt':           None,
            'paeinter':        None,
            'sequences':       sequences,
        })
        self._af3_update_jobs_table()
        n_chains = len(subunits)
        self._status.showMessage(
            f"✓ Custom job: {n_chains} subunit(s) [{chain_summary}]  "
            f"({total_residues} residues)")

    def _af3_jobs_right_click(self, pos):
        """Right-click on AF3 jobs table → enhanced context menu."""
        rows = sorted(set(idx.row() for idx in self._af3_jobs_table.selectedIndexes()))
        if not rows:
            return
        menu = QMenu(self)
        n = len(rows)

        # ── Navigation ──────────────────────────────────────────────────────
        if len(rows) == 1 and rows[0] < len(self.af3_jobs):
            j = self.af3_jobs[rows[0]]
            hi = j.get('hit_orf_idx', -1)
            pi = j.get('partner_orf_idx', -1)
            if hi >= 0:
                menu.addAction(
                    f"🗺  Focus genome → {j['hit_name']}",
                    lambda _hi=hi: self._select_and_center_orf(_hi))
            if pi >= 0 and pi != hi:
                menu.addAction(
                    f"🗺  Focus genome → {j['partner_name']}",
                    lambda _pi=pi: self._select_and_center_orf(_pi))
            menu.addSeparator()

        # ── Export selected jobs only ────────────────────────────────────────
        menu.addAction(
            f"⬇  Export {n} selected job{'s' if n > 1 else ''} as JSON…",
            lambda: self._af3_export_selected_jobs_json(rows))
        menu.addSeparator()

        # ── Clipboard ────────────────────────────────────────────────────────
        menu.addAction(
            "📋  Copy job names",
            lambda: self._copy_to_clipboard(
                '\n'.join(self.af3_jobs[r]['name']
                          for r in rows if r < len(self.af3_jobs))))
        menu.addSeparator()

        # ── Delete ───────────────────────────────────────────────────────────
        menu.addAction(
            f"🗑  Delete {n} selected job{'s' if n > 1 else ''}",
            self._af3_delete_selected_jobs)

        menu.exec(self._af3_jobs_table.viewport().mapToGlobal(pos))

    # ──────────────────────────────────────────────────────────────────────────
    # AF3 selection-table helpers (NEW in v2.0)
    # ──────────────────────────────────────────────────────────────────────────

    def _af3_sel_table_click(self, row, col):
        """Left-click a row in the AF3 selection table →
        center genome map on that ORF and select it in the main ORF table."""
        item = self._af3_sel_table.item(row, 0)
        if item is None:
            return
        idx = self._parse_orf_idx_from_text(item.text())
        if idx >= 0:
            self._select_and_center_orf(idx)

    def _af3_jobs_table_click(self, row, col):
        """Left-click a row in the AF3 jobs table →
        center genome map on the Hit ORF of that job."""
        if row < 0 or row >= len(self.af3_jobs):
            return
        j = self.af3_jobs[row]
        hi = j.get('hit_orf_idx', -1)
        if hi >= 0:
            self._select_and_center_orf(hi)

    def _af3_update_predict_pair_btn(self):
        """Enable/disable the 'Predict Selected Pair' button based on how many
        rows are highlighted in the AF3 selection table."""
        n = len(set(idx.row() for idx in self._af3_sel_table.selectedIndexes()))
        enabled = n >= 2
        self._af3_predict_pair_btn.setEnabled(enabled)
        if enabled:
            self._af3_predict_pair_btn.setText(
                f"⚡ Predict Selected Pair ({n} ORFs)")
        else:
            self._af3_predict_pair_btn.setText("⚡ Predict Selected Pair")

    def _af3_predict_selected_pair(self):
        """Create AF3 pairwise jobs for the rows currently highlighted in the
        AF3 selection table via Ctrl+click (≥ 2 rows required).
        All pairwise combinations among the highlighted ORFs are generated and
        appended to the jobs list immediately — no mode dropdown involved."""
        sel_rows = sorted(set(idx.row() for idx in self._af3_sel_table.selectedIndexes()))
        if len(sel_rows) < 2:
            QMessageBox.information(
                self, "Predict Selected Pair",
                "Ctrl+click at least 2 ORFs in the selection list first.")
            return
        if not self.orfs:
            QMessageBox.warning(self, "AF3", "Run ORF analysis first!")
            return

        # Collect ORF indices for the highlighted rows
        pair_indices = []
        for r in sel_rows:
            item = self._af3_sel_table.item(r, 0)
            if item is None:
                continue
            try:
                idx = int(item.text().replace('ORF', '')) - 1
                if 0 <= idx < len(self.orfs):
                    pair_indices.append(idx)
            except (ValueError, TypeError):
                continue

        if len(pair_indices) < 2:
            QMessageBox.warning(self, "AF3", "Could not resolve ORF indices.")
            return

        added = 0
        existing_names = {j['name'] for j in self.af3_jobs}
        seen_pairs: set = set()
        for i_a in range(len(pair_indices)):
            for i_b in range(i_a + 1, len(pair_indices)):
                hi = pair_indices[i_a]
                pi = pair_indices[i_b]
                key = (min(hi, pi), max(hi, pi))
                if key in seen_pairs:
                    continue
                seen_pairs.add(key)
                hn = f"ORF{hi + 1}"; pn = f"ORF{pi + 1}"
                job_name = f"{hn}_vs_{pn}_selected"
                if job_name in existing_names:
                    continue
                hp = self.orfs[hi]['protein'].rstrip('*')
                pp = self.orfs[pi]['protein'].rstrip('*')
                tr = len(hp) + len(pp)
                self.af3_jobs.append({
                    'name':            job_name,
                    'hit_orf_idx':     hi,
                    'partner_orf_idx': pi,
                    'hit_name':        hn,
                    'partner_name':    pn,
                    'total_residues':  tr,
                    'status':          ('pending'
                                        if tr <= self.af3_max_residues
                                        else f'>{self.af3_max_residues}!'),
                    'iptm':      None,
                    'plddt':     None,
                    'paeinter':  None,
                    'sequences': [
                        {'proteinChain': {'sequence': hp, 'count': 1}},
                        {'proteinChain': {'sequence': pp, 'count': 1}}],
                })
                existing_names.add(job_name)
                added += 1

        self._af3_update_jobs_table()
        self._status.showMessage(
            f"✓ {added} pairwise job(s) added from Ctrl+click selection "
            f"({len(pair_indices)} ORFs → {added} pair(s))")

    def _af3_sel_table_right_click(self, pos):
        """Right-click on the AF3 selection table → context menu."""
        row = self._af3_sel_table.rowAt(pos.y())
        if row < 0:
            return
        sel_rows = sorted(set(idx.row() for idx in self._af3_sel_table.selectedIndexes()))
        menu = QMenu(self)

        # ── Navigation (single-row click) ────────────────────────────────────
        item = self._af3_sel_table.item(row, 0)
        if item:
            try:
                orf_idx = int(item.text().replace('ORF', '')) - 1
            except (ValueError, TypeError):
                orf_idx = -1
            orf_name = item.text()
            if orf_idx >= 0:
                menu.addAction(
                    f"🗺  Focus genome on {orf_name}",
                    lambda _i=orf_idx: self._select_and_center_orf(_i))
                menu.addAction(
                    f"🔍  Select {orf_name} in ORF table",
                    lambda _i=orf_idx: self._select_and_center_orf(_i))
                menu.addSeparator()

                # ── Quick pair prediction ────────────────────────────────────
                n_sel = len(sel_rows)
                if n_sel >= 2:
                    menu.addAction(
                        f"⚡  Generate prediction for {n_sel} selected ORFs",
                        self._af3_predict_selected_pair)
                else:
                    menu.addAction(
                        f"⚡  Predict {orf_name} vs neighbors (Pairs mode)",
                        lambda _i=orf_idx: self._af3_predict_single_vs_neighbors(_i))
                menu.addSeparator()

                # ── Sequence copy ────────────────────────────────────────────
                if 0 <= orf_idx < len(self.orfs):
                    orf = self.orfs[orf_idx]
                    menu.addAction(
                        "📋  Copy protein sequence (FASTA)",
                        lambda _o=orf, _n=orf_name: self._copy_to_clipboard(
                            f">{_n}\n{_o['protein'].rstrip('*')}"))
                    menu.addAction(
                        "📋  Copy DNA sequence",
                        lambda _o=orf, _n=orf_name: self._copy_to_clipboard(
                            f">{_n}\n{self.dna_sequence[_o['start']:_o['end']]}"
                            if self.dna_sequence else "(no DNA loaded)"))
                menu.addSeparator()

        # ── Remove ───────────────────────────────────────────────────────────
        n = len(sel_rows) if sel_rows else 1
        menu.addAction(
            f"✖  Remove {n} ORF{'s' if n > 1 else ''} from AF3 list",
            self._af3_remove_orf)

        menu.exec(self._af3_sel_table.viewport().mapToGlobal(pos))

    def _af3_predict_single_vs_neighbors(self, orf_idx: int):
        """Quick-predict one ORF against its genomic neighbors using the
        current N-neighbors setting — without touching the mode combo."""
        if not self.orfs:
            return
        n_nb = self._af3_nb_spin.value()
        orf  = self.orfs[orf_idx]
        hp   = orf['protein'].rstrip('*')
        hn   = f"ORF{orf_idx + 1}"
        orfs_by_pos = sorted(enumerate(self.orfs), key=lambda x: x[1]['start'])
        pos_to_rank = {idx: rank for rank, (idx, _) in enumerate(orfs_by_pos)}
        hr = pos_to_rank.get(orf_idx, 0)
        existing_names = {j['name'] for j in self.af3_jobs}
        added = 0
        for d in range(-n_nb, n_nb + 1):
            if d == 0:
                continue
            nr = hr + d
            if not (0 <= nr < len(orfs_by_pos)):
                continue
            ni, no = orfs_by_pos[nr]
            np_s = no['protein'].rstrip('*')
            nn   = f"ORF{ni + 1}"
            tr   = len(hp) + len(np_s)
            job_name = f"{hn}_vs_{nn}_{'up' if d < 0 else 'down'}{abs(d)}_quick"
            if job_name in existing_names:
                continue
            self.af3_jobs.append({
                'name':            job_name,
                'hit_orf_idx':     orf_idx,
                'partner_orf_idx': ni,
                'hit_name':        hn,
                'partner_name':    nn,
                'total_residues':  tr,
                'status':          ('pending'
                                    if tr <= self.af3_max_residues
                                    else f'>{self.af3_max_residues}!'),
                'iptm':     None,
                'plddt':    None,
                'paeinter': None,
                'sequences': [
                    {'proteinChain': {'sequence': hp,   'count': 1}},
                    {'proteinChain': {'sequence': np_s, 'count': 1}}],
            })
            existing_names.add(job_name)
            added += 1
        self._af3_update_jobs_table()
        self._status.showMessage(
            f"✓ {added} quick-predict jobs added for {hn} vs {n_nb} neighbors")

    def _af3_export_selected_jobs_json(self, rows: list):
        """Export only the selected rows from the jobs table as individual JSONs."""
        jobs_to_export = [self.af3_jobs[r] for r in rows if r < len(self.af3_jobs)]
        if not jobs_to_export:
            return
        folder = QFileDialog.getExistingDirectory(self, "Select output folder for selected jobs")
        if not folder:
            return
        errors = []
        for j in jobs_to_export:
            safe_name = re.sub(r'[^\w.\-]', '_', j['name'])
            data = {"name": j['name'], "modelSeeds": [],
                    "sequences": j['sequences'],
                    "dialect": "alphafoldserver", "version": 2}
            try:
                with open(os.path.join(folder, f"{safe_name}.json"),
                          'w', encoding='utf-8') as f:
                    json.dump(data, f, indent=2, ensure_ascii=False)
            except OSError as e:
                errors.append(f"{safe_name}.json: {e}")
        if errors:
            QMessageBox.warning(self, "AF3 Export",
                f"{len(errors)} file(s) could not be written:\n" + "\n".join(errors))
        self._status.showMessage(
            f"✓ {len(jobs_to_export) - len(errors)} selected job(s) exported as JSON")

    def _af3_delete_selected_jobs(self):
        """Delete selected rows from AF3 jobs table and internal list."""
        rows = sorted(set(idx.row() for idx in self._af3_jobs_table.selectedIndexes()), reverse=True)
        if not rows: return
        for r in rows:
            if r < len(self.af3_jobs):
                self.af3_jobs.pop(r)
        self._af3_update_jobs_table()
        self._status.showMessage(f"✓ {len(rows)} job(s) removed — {len(self.af3_jobs)} remaining")

    def keyPressEvent(self, event):
        """Handle Delete key for AF3 jobs table."""
        key = event.key()
        if key in (Qt.Key.Key_Delete if QT_VERSION == 6 else Qt.Key_Delete,
                   Qt.Key.Key_Backspace if QT_VERSION == 6 else Qt.Key_Backspace):
            if self._af3_jobs_table.hasFocus() and self._af3_jobs_table.selectedIndexes():
                self._af3_delete_selected_jobs()
                return
        super().keyPressEvent(event)

    # ═══════════════════════════════════════════════════════════
    # DIALOGS
    # ═══════════════════════════════════════════════════════════

    def _show_orf_params(self):
        """Dialog for ORF analysis parameters — 6-frame scanner and Pyrodigal."""
        dlg = QDialog(self)
        dlg.setWindowTitle("ORF Analysis Parameters")
        dlg.setMinimumWidth(460)
        layout = QVBoxLayout(dlg)
        layout.setSpacing(8)

        # ── Section 1: 6-frame scanner ──────────────────────────
        gf = QGroupBox("6-frame ORF scanner (Translate genome → Automatic)")
        gf_l = QGridLayout(gf)
        gf_l.setSpacing(6)

        gf_l.addWidget(QLabel("Min ORF size (aa):"), 0, 0)
        min_spin = QSpinBox(); min_spin.setRange(10, 500)
        min_spin.setValue(self._min_length_spin.value())
        min_spin.setToolTip("Minimum protein length for 6-frame scanner.")
        gf_l.addWidget(min_spin, 0, 1)

        gf_l.addWidget(QLabel("Start codons:"), 1, 0)
        cb_atg = QCheckBox("ATG"); cb_atg.setChecked(self._cb_atg.isChecked())
        cb_gtg = QCheckBox("GTG"); cb_gtg.setChecked(self._cb_gtg.isChecked())
        cb_ttg = QCheckBox("TTG"); cb_ttg.setChecked(self._cb_ttg.isChecked())
        for cb in (cb_atg, cb_gtg, cb_ttg):
            cb.setToolTip("Start codons recognised by the 6-frame scanner.")
        codon_lay = QHBoxLayout()
        codon_lay.addWidget(cb_atg); codon_lay.addWidget(cb_gtg); codon_lay.addWidget(cb_ttg)
        codon_lay.addStretch()
        gf_l.addLayout(codon_lay, 1, 1)
        layout.addWidget(gf)

        # ── Section 2: Pyrodigal ────────────────────────────────
        pyro_ok = PYRODIGAL_AVAILABLE
        pyro_ver = pyrodigal.__version__ if pyro_ok else "not installed"
        pf = QGroupBox(f"Pyrodigal {pyro_ver} (Translate genome → Pyrodigal)")
        pf_l = QGridLayout(pf)
        pf_l.setSpacing(6)

        status_lbl = QLabel(
            "Status: " + ("Installed" if pyro_ok else "NOT installed — pip install pyrodigal"))
        status_lbl.setStyleSheet(
            "color: #2e7d32; font-weight:bold;" if pyro_ok else "color: #c62828; font-weight:bold;")
        pf_l.addWidget(status_lbl, 0, 0, 1, 4)

        # Mode
        pf_l.addWidget(QLabel("Mode:"), 1, 0)
        mode_combo = QComboBox()
        mode_combo.addItems([
            "Metagenomic (meta=True)  — recommended",
            "Single genome (meta=False)  — trains on sequence",
        ])
        mode_combo.setCurrentIndex(0 if self._pyro_params.get('meta', True) else 1)
        mode_combo.setToolTip(
            "Metagenomic: pre-trained models, any contig length.\n"
            "Single genome: trains Prodigal on the sequence (needs >= 100 kb).")
        mode_combo.setEnabled(pyro_ok)
        pf_l.addWidget(mode_combo, 1, 1, 1, 3)

        # Translation table
        pf_l.addWidget(QLabel("Translation table:"), 2, 0)
        tt_combo = QComboBox()
        tt_items = [
            ("11 — Bacteria / Archaea (standard)", 11),
            ("4  — Mycoplasma / Spiroplasma  (UGA→Trp)", 4),
            ("25 — SR1 / Gracilibacteria  (UGA→Gly)", 25),
            ("15 — Yeast mitochondria", 15),
        ]
        for lbl, val in tt_items:
            tt_combo.addItem(lbl, val)
        cur_tt = self._pyro_params.get('translation_table', 11)
        for i, (_, v) in enumerate(tt_items):
            if v == cur_tt:
                tt_combo.setCurrentIndex(i); break
        tt_combo.setToolTip(
            "Genetic code for translation.\n"
            "Table 11 is correct for most bacteria.\n"
            "Table 4: Mycoplasma, Spiroplasma, Phytoplasma.\n"
            "Table 25: SR1/Gracilibacteria.")
        tt_combo.setEnabled(pyro_ok)
        pf_l.addWidget(tt_combo, 2, 1, 1, 3)

        # Min gene size (shared with 6-frame but independent)
        pf_l.addWidget(QLabel("Min gene size (aa):"), 3, 0)
        pyro_min = QSpinBox(); pyro_min.setRange(10, 500)
        pyro_min.setValue(self._pyro_params.get('min_aa', 30))
        pyro_min.setSuffix(" aa")
        pyro_min.setToolTip("Minimum protein length for Pyrodigal predictions.")
        pyro_min.setEnabled(pyro_ok)
        pf_l.addWidget(pyro_min, 3, 1)

        # Closed ends
        closed_cb = QCheckBox("Closed ends")
        closed_cb.setChecked(self._pyro_params.get('closed', False))
        closed_cb.setToolTip(
            "Checked: genes must start AND end within the sequence.\n"
            "Unchecked: allows partial genes at contig edges (recommended for drafts).")
        closed_cb.setEnabled(pyro_ok)
        pf_l.addWidget(closed_cb, 3, 2)

        # Mask N runs
        mask_cb = QCheckBox("Mask N runs")
        mask_cb.setChecked(self._pyro_params.get('mask', False))
        mask_cb.setToolTip("Mask regions with runs of N before prediction (draft genomes).")
        mask_cb.setEnabled(pyro_ok)
        pf_l.addWidget(mask_cb, 3, 3)

        # Start codons note for Pyrodigal
        note = QLabel(
            "Start codons: Pyrodigal selects ATG / GTG / TTG automatically\n"
            "based on the genetic code and training data — no manual override.\n"
            "Post-filter: check below to keep only specific starts after prediction.")
        note.setStyleSheet("color: #555; font-size: 11px;")
        note.setWordWrap(True)
        pf_l.addWidget(note, 4, 0, 1, 4)

        # Post-prediction start codon filter
        pf_l.addWidget(QLabel("Post-filter starts:"), 5, 0)
        pyro_cb_atg = QCheckBox("ATG")
        pyro_cb_gtg = QCheckBox("GTG")
        pyro_cb_ttg = QCheckBox("TTG")
        pyro_cb_all = QCheckBox("All (no filter)")
        # Load from params
        pf_opts = self._pyro_params.get('start_filter', {'ATG':True,'GTG':True,'TTG':True,'all':True})
        pyro_cb_atg.setChecked(pf_opts.get('ATG', True))
        pyro_cb_gtg.setChecked(pf_opts.get('GTG', True))
        pyro_cb_ttg.setChecked(pf_opts.get('TTG', True))
        pyro_cb_all.setChecked(pf_opts.get('all', True))
        for cb in (pyro_cb_atg, pyro_cb_gtg, pyro_cb_ttg, pyro_cb_all):
            cb.setEnabled(pyro_ok)
        # "All" disables/enables individual checkboxes
        def _toggle_filter(state):
            all_checked = pyro_cb_all.isChecked()
            for c in (pyro_cb_atg, pyro_cb_gtg, pyro_cb_ttg):
                c.setEnabled(pyro_ok and not all_checked)
        pyro_cb_all.stateChanged.connect(_toggle_filter)
        _toggle_filter(None)
        starts_lay = QHBoxLayout()
        starts_lay.addWidget(pyro_cb_atg); starts_lay.addWidget(pyro_cb_gtg)
        starts_lay.addWidget(pyro_cb_ttg); starts_lay.addWidget(pyro_cb_all)
        starts_lay.addStretch()
        pf_l.addLayout(starts_lay, 5, 1, 1, 3)

        layout.addWidget(pf)

        # ── Info box ────────────────────────────────────────────
        info_box = QLabel()
        info_box.setWordWrap(True)
        info_box.setStyleSheet(
            "background:#f3f8ff;border:1px solid #bbdefb;"
            "border-radius:4px;padding:5px;color:#0d47a1;font-size:11px;")
        layout.addWidget(info_box)

        def _upd_info():
            is_meta = mode_combo.currentIndex() == 0
            tt = tt_combo.currentData() if pyro_ok else 11
            seq_kb = len(self.dna_sequence) / 1000 if self.dna_sequence else 0
            if not pyro_ok:
                info_box.setText("Pyrodigal not installed. Install: pip install pyrodigal")
                return
            warn = ""
            if not is_meta and seq_kb < 100:
                warn = f"Warning: sequence is {seq_kb:.0f} kb — single-genome needs >= 100 kb. "
            info_box.setText(
                f"{warn}Table {tt} | Mode: {'meta' if is_meta else 'single'} | "
                f"Seq: {seq_kb:.0f} kb")
        if pyro_ok:
            mode_combo.currentIndexChanged.connect(_upd_info)
            tt_combo.currentIndexChanged.connect(_upd_info)
        _upd_info()

        # ── Buttons ─────────────────────────────────────────────
        btns = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
            if QT_VERSION == 6 else QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(dlg.accept); btns.rejected.connect(dlg.reject)
        layout.addWidget(btns)

        if dlg.exec() if QT_VERSION == 6 else dlg.exec_():
            # Save 6-frame params
            self._min_length_spin.setValue(min_spin.value())
            self._cb_atg.setChecked(cb_atg.isChecked())
            self._cb_gtg.setChecked(cb_gtg.isChecked())
            self._cb_ttg.setChecked(cb_ttg.isChecked())
            # Save Pyrodigal params
            self._pyro_params.update({
                'meta':              mode_combo.currentIndex() == 0,
                'translation_table': tt_combo.currentData() if pyro_ok else 11,
                'min_aa':            pyro_min.value(),
                'closed':            closed_cb.isChecked(),
                'mask':              mask_cb.isChecked(),
                'start_filter': {
                    'ATG': pyro_cb_atg.isChecked(),
                    'GTG': pyro_cb_gtg.isChecked(),
                    'TTG': pyro_cb_ttg.isChecked(),
                    'all': pyro_cb_all.isChecked(),
                },
            })
            self._status.showMessage(
                f"Parameters saved — 6-frame: min={min_spin.value()}aa | "
                f"Pyrodigal: mode={'meta' if self._pyro_params['meta'] else 'single'}, "
                f"table={self._pyro_params['translation_table']}, "
                f"min={self._pyro_params['min_aa']}aa")

    def _show_blast_params(self):
        dlg = QDialog(self)
        dlg.setWindowTitle(t('blast_dlg_title'))
        dlg.setFixedSize(420, 580)
        layout = QVBoxLayout(dlg)

        # Algorithm
        af = QGroupBox("Algorithm")
        af_l = QGridLayout(af)
        af_l.addWidget(QLabel("Algorithm:"), 0, 0)
        algo_cb = QComboBox()
        algo_choices = ["Auto (best available)"]
        if BACKENDS.get('blast+',{}).get('available'):
            algo_choices.append("NCBI BLAST+ (externo)")
        algo_choices.extend(["K-mer Filter (rápido)", "Smith-Waterman (sensível)"])
        algo_cb.addItems(algo_choices)
        algo_cb.setCurrentText(self._algo_combo.currentText())
        af_l.addWidget(algo_cb, 0, 1)
        af_l.addWidget(QLabel("Min identity (%):"), 1, 0)
        id_spin = QSpinBox(); id_spin.setRange(0,100); id_spin.setValue(self._identity_spin.value())
        af_l.addWidget(id_spin, 1, 1)
        af_l.addWidget(QLabel("Max E-value:"), 2, 0)
        ev_edit = QLineEdit(self._evalue_edit.text())
        af_l.addWidget(ev_edit, 2, 1)
        layout.addWidget(af)

        # Program
        pg = QGroupBox(t('blast_program'))
        pg_l = QGridLayout(pg)
        pg_l.addWidget(QLabel(t('blast_program')), 0, 0)
        prog_cb = QComboBox()
        prog_cb.addItems(["blastp","tblastn","blastn","blastx"])
        prog_cb.setCurrentText(self.blast_program)
        pg_l.addWidget(prog_cb, 0, 1)
        layout.addWidget(pg)

        # General
        gf = QGroupBox(t('blast_gen_params'))
        gf_l = QGridLayout(gf)
        gf_l.addWidget(QLabel("Max targets:"), 0, 0)
        max_t = QSpinBox(); max_t.setRange(10,1000); max_t.setValue(self.blast_max_targets)
        gf_l.addWidget(max_t, 0, 1)
        gf_l.addWidget(QLabel("E-value:"), 1, 0)
        ev = QLineEdit(str(self.blast_evalue))
        gf_l.addWidget(ev, 1, 1)
        gf_l.addWidget(QLabel("Word size:"), 2, 0)
        ws = QComboBox(); ws.addItems(["3","4","5","6"]); ws.setCurrentText(str(self.blast_word_size))
        gf_l.addWidget(ws, 2, 1)
        layout.addWidget(gf)

        # Scoring
        sf = QGroupBox(t('blast_scoring'))
        sf_l = QGridLayout(sf)
        sf_l.addWidget(QLabel("Matrix:"), 0, 0)
        mx = QComboBox(); mx.addItems(["BLOSUM45","BLOSUM62","BLOSUM80","PAM30","PAM70"])
        mx.setCurrentText(self.blast_matrix)
        sf_l.addWidget(mx, 0, 1)
        sf_l.addWidget(QLabel("Gap open:"), 1, 0)
        go = QComboBox(); go.addItems([str(i) for i in range(7,14)])
        go.setCurrentText(str(self.blast_gap_open))
        sf_l.addWidget(go, 1, 1)
        sf_l.addWidget(QLabel("Gap extend:"), 2, 0)
        ge = QComboBox(); ge.addItems(["1","2","3"]); ge.setCurrentText(str(self.blast_gap_ext))
        sf_l.addWidget(ge, 2, 1)
        layout.addWidget(sf)

        # Buttons
        btns = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
                                 if QT_VERSION == 6 else QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(dlg.accept)
        btns.rejected.connect(dlg.reject)
        layout.addWidget(btns)

        if dlg.exec() if QT_VERSION == 6 else dlg.exec_():
            self._algo_combo.setCurrentText(algo_cb.currentText())
            self._identity_spin.setValue(id_spin.value())
            self._evalue_edit.setText(ev_edit.text())
            self.blast_program = prog_cb.currentText()
            self.blast_max_targets = max_t.value()
            self.blast_evalue = float(ev.text() or '0.05')
            self.blast_word_size = int(ws.currentText())
            self.blast_matrix = mx.currentText()
            self.blast_gap_open = int(go.currentText())
            self.blast_gap_ext = int(ge.currentText())
            self._refresh_blast_cmd_preview()

    def _show_hmm_params(self):
        dlg = QDialog(self)
        dlg.setWindowTitle(t('hmm_dlg_title'))
        dlg.setFixedSize(380, 250)
        layout = QVBoxLayout(dlg)
        gf = QGroupBox(t('hmm_thresh'))
        gf_l = QGridLayout(gf)
        gf_l.addWidget(QLabel("Seq E-value:"), 0, 0)
        ev = QLineEdit(str(self.hmm_evalue)); gf_l.addWidget(ev, 0, 1)
        gf_l.addWidget(QLabel("Dom E-value:"), 1, 0)
        dev = QLineEdit(str(self.hmm_dom_evalue)); gf_l.addWidget(dev, 1, 1)
        layout.addWidget(gf)
        btns = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
                                 if QT_VERSION == 6 else QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(dlg.accept); btns.rejected.connect(dlg.reject)
        layout.addWidget(btns)
        if dlg.exec() if QT_VERSION == 6 else dlg.exec_():
            self.hmm_evalue = float(ev.text() or '10.0')
            self.hmm_dom_evalue = float(dev.text() or '10.0')

    def _show_manual(self):
        self._show_help_dlg('manual', 'manual')

    def _show_tutorial(self):
        self._show_help_dlg('tutorial', 'tutorial')

    def _show_install(self):
        self._show_help_dlg('install', 'install', width=800, height=640)

    def _show_help_dlg(self, title_key, content_key, width=720, height=560):
        lang = _CURRENT_LANG[0]
        content = HELP_CONTENT.get(content_key, {}).get(lang) or \
                  HELP_CONTENT.get(content_key, {}).get('en', '(no content)')
        dlg = QDialog(self)
        dlg.setWindowTitle(t(title_key))
        dlg.resize(width, height)
        layout = QVBoxLayout(dlg)
        txt = QTextEdit()
        txt.setFont(QFont('Courier', 9))
        txt.setReadOnly(True)
        txt.setPlainText(content)
        layout.addWidget(txt)
        btns = QDialogButtonBox(QDialogButtonBox.StandardButton.Close if QT_VERSION == 6
                                 else QDialogButtonBox.Close)
        btns.rejected.connect(dlg.close)
        layout.addWidget(btns)
        dlg.exec() if QT_VERSION == 6 else dlg.exec_()

    def _show_help_interaction_results(self):
        """Help dialog: Interaction Results tab — full analysis guide."""
        dlg = QDialog(self)
        dlg.setWindowTitle("📈 Interaction Results — Analysis Guide")
        dlg.resize(820, 680)
        lay = QVBoxLayout(dlg)
        txt = QTextEdit()
        txt.setFont(QFont('Courier', 9))
        txt.setReadOnly(True)
        txt.setPlainText("""\
INTERACTION RESULTS TAB — ppigFinder v2.0
══════════════════════════════════════════════════════════════════

PURPOSE
───────
The Interaction Results tab imports AlphaFold 3 output folders and
displays ranked, annotated results for every predicted protein pair.
It is the primary decision-support tool for identifying real
protein-protein interactions from large AF3 batch runs.

HOW TO USE
──────────
  1. Click "Load AF3 results folder" → select a folder containing one
     or more AF3 job subfolders (each must contain a
     *_summary_confidences_0.json AND a *_full_data_0.json).
  2. Results appear in the table, sorted by ranking_score by default.
  3. Click any row to view the PAE heatmap, pLDDT plot, and motif table
     below.
  4. Adjust filter controls to focus on high-confidence pairs.
  5. Export to TSV for downstream analysis.

TABLE COLUMNS (left → right)
─────────────────────────────
  Job name         Folder / job identifier
  Chains (n)       Number of protein chains (usually 2 for pairwise)
  Chain IDs        ORF names and sizes for each chain
  ipTM             Interface pTM score (0–1). Measures how well AF3
                   predicts the RELATIVE positions of the two chains.
                   > 0.75 = high confidence  |  0.50–0.75 = likely
  ptm              Predicted TM-score (single-chain quality; less
                   informative for interaction assessment)
  mean_pLDDT       Mean per-residue confidence (0–100). Both chains.
  ranking_score    AF3's own combined confidence score
  PAE_inter (Å)    GLOBAL mean of the entire off-diagonal PAE quadrant.
                   ⚠ Diluted by disordered residues — use PAE_min ★
  PAE_min ★ (Å)   Minimum PAE in the off-diagonal quadrant.
                   = chain_pair_pae_min from summary_confidences.json.
                   < 4 Å = at least one contact predicted with near-
                   atomic confidence regardless of global disorder.
                   KEY METRIC for domain-limited interactions.
  cp_ipTM ★        Chain-pair ipTM for the interface only.
                   = chain_pair_iptm[A][B] from summary_confidences.json.
                   Removes intra-chain folding quality contribution.
  Contact% ★       Fraction of residue pairs with PAE < 5 Å.
                   Even small values (5–10%) confirm a focal interface.
  Best contact     Focal residue ranges in each chain where PAE < 5 Å.

CONFIDENCE CLASSIFICATION
──────────────────────────
  HIGH ★  (bold green row)
          PAE_min < 4 Å  AND  cp_ipTM ≥ 0.50
          → Focal domain contact confirmed by AF3 with high precision.
          Recommended for follow-up (mutagenesis, pull-down, etc.).

  MED     (amber row)
          PAE_min 4–8 Å  (any cp_ipTM)
          → Possible interaction; check the motif detector results and
          pLDDT plots. May be a real but flexible interface.

  LOW     (pink row)
          PAE_min ≥ 8 Å
          → No confident contact predicted. Could still be correct if
          the proteins require a cofactor or membrane environment.

FILTER CONTROLS
───────────────
  Contact threshold (Å)   Residues with inter-chain PAE below this
                           value are counted as contacts. Affects the
                           contact-region string and n_contacts display.
                           Does NOT change ipTM or PAE_min.
                           Typical: 5 Å (atomic contact); 10 Å (proximity).

  min ipTM                 Hide rows with global ipTM below threshold.
                           Set to 0.00 to show all results.

  max PAE_inter            Hide rows with global PAE_inter above threshold.
                           ⚠ Use with caution — PAE_inter is diluted by
                           disorder. Prefer filtering by PAE_min ★ column.

MOTIF DETECTION (second toolbar row)
──────────────────────────────────────
  core PAE (Å)    Tier-1 threshold: pixels below this form the motif
                  core (high-confidence contacts). Typical: 5–8 Å.
                  Lower = stricter, fewer but higher-confidence motifs.

  ext PAE (Å)     Tier-2 threshold: used to grow the core bounding box
                  to include peripheral contacts. Must be > core PAE.
                  Typical: 12–18 Å.

  min contact ≥   Minimum AF3 contact_probs value for a pixel to join
                  the core. contact_probs is AF3's internal probability
                  that two residues are in physical contact (independent
                  of PAE). Set to 0.00 to disable (use PAE alone).

  min size        Minimum bounding-box dimension (residues) for a motif
                  to survive. Rejects isolated noise pixels.
                  Recommended: 3–5 for domain contacts; 2 for peptides.

  reciprocal ☑    Require the motif to appear in BOTH PAE quadrants
                  (A→B and B→A) with ≥ 50% pixel overlap. Eliminates
                  one-sided ghost motifs caused by disordered alignment.
                  Strongly recommended: keep checked.

  Rerun           Re-run motif detection on all loaded results with
                  the current parameter values.

REFERENCES
──────────
  [Evans et al. 2022]   ipTM, chain_pair_iptm, chain_pair_pae_min
                        (AlphaFold-Multimer paper)
  [Abramson et al. 2024] AF3 contact_probs, ranking_score
  [Jumper et al. 2021]  PAE (Predicted Aligned Error) metric definition
  [Humphreys et al. 2021] ipTM-based proteome-scale PPI screening
  [Bryant et al. 2022]  PAE inter-chain analysis for PPI confidence
  See Help → References & methodology for full citations.
""")
        lay.addWidget(txt)
        btns = QDialogButtonBox(QDialogButtonBox.StandardButton.Close
                                if QT_VERSION == 6 else QDialogButtonBox.Close)
        btns.rejected.connect(dlg.close)
        lay.addWidget(btns)
        dlg.exec() if QT_VERSION == 6 else dlg.exec_()

    def _show_help_ppi_map(self):
        """Help dialog: Genomic PPI Map tab guide."""
        dlg = QDialog(self)
        dlg.setWindowTitle("🧬 Genomic PPI Map — Guide")
        dlg.resize(780, 560)
        lay = QVBoxLayout(dlg)
        txt = QTextEdit()
        txt.setFont(QFont('Courier', 9))
        txt.setReadOnly(True)
        txt.setPlainText("""\
GENOMIC PPI MAP TAB — ppigFinder v2.0
══════════════════════════════════════════════════════════════════

PURPOSE
───────
Displays all AF3-predicted protein-protein interactions as curved
arcs drawn above a linear genomic map — analogous to published operon
interaction diagrams (e.g. T4SS component interaction figures).

This view makes it immediately apparent whether interacting partners
are genomic neighbours (likely operon members) or distant co-operating
proteins elsewhere in the chromosome.

READING THE MAP
───────────────
  Horizontal axis   Genomic position (bp). Scroll to pan; mouse wheel
                    to zoom (position under cursor is preserved).
  ORF arrows        Directional arrows on the backbone; colour encodes
                    HMM annotation family (green = HMM hit, orange =
                    custom/VirD4, gray = no annotation).
  Arcs above ORFs   Each arc connects a predicted interaction pair.
    ─── green solid    PAE_min < 4 Å  — HIGH confidence
    ─── amber dashed   PAE_min 4–8 Å  — MED confidence
    ··· red dotted     PAE_min > 8 Å  — LOW / no contact
  Arc height        By default, proportional to genomic distance:
                    neighbouring ORFs = low flat arc (like operon diagrams);
                    distant ORFs = tall arc.
  Score label       PAE_min value shown near the arc apex.

INTERACTION
───────────
  Click an arc      → Shows interaction details in the info bar AND
                    selects the corresponding row in the Interaction
                    Results table.
  Click an ORF      → Centers the main genome map on that ORF.
  Drag              → Pan the view.
  Mouse wheel       → Zoom in/out (genomic position preserved).

FILTER CONTROLS
───────────────
  Show              Filter by confidence: All / HIGH only / HIGH+MED /
                    Focal hits (PAE_min < 4 Å AND cp_ipTM ≥ 0.50).
  Color by          Metric used to determine arc colour (PAE_min ★
                    recommended; also ipTM, cp_ipTM, Contact%).
  Arc height        "By genomic distance" (default) reproduces the
                    classic operon interaction diagram style;
                    "By score" puts high-confidence arcs lower/closer;
                    "Fixed" gives equal height to all arcs.

EXPORT
──────
  Export SVG        Saves the entire arc map as a scalable vector
                    graphic (SVG) suitable for publication figures.
  Export TSV        Saves the currently displayed interactions as a
                    tab-separated file with all confidence metrics.

NOTE
────
The Genomic PPI Map reads from the same data as the Interaction Results
table. Click "Refresh" after loading new AF3 results or changing filters
in the Interaction Results tab.
""")
        lay.addWidget(txt)
        btns = QDialogButtonBox(QDialogButtonBox.StandardButton.Close
                                if QT_VERSION == 6 else QDialogButtonBox.Close)
        btns.rejected.connect(dlg.close)
        lay.addWidget(btns)
        dlg.exec() if QT_VERSION == 6 else dlg.exec_()

    def _show_help_references(self):
        """Help dialog: full references and methodology note."""
        dlg = QDialog(self)
        dlg.setWindowTitle("📚 References & Methodology")
        dlg.resize(820, 700)
        lay = QVBoxLayout(dlg)
        txt = QTextEdit()
        txt.setFont(QFont('Courier', 9))
        txt.setReadOnly(True)
        txt.setPlainText("""\
REFERENCES & METHODOLOGY — ppigFinder v2.0
══════════════════════════════════════════════════════════════════

CORE DEPENDENCIES
─────────────────
[1]  Larralde, M. (2022). Pyrodigal: Python bindings and interface to
     Prodigal. Journal of Open Source Software, 7(72), 4296.
     https://doi.org/10.21105/joss.04296

[2]  Hyatt, D. et al. (2010). Prodigal: prokaryotic gene recognition and
     translation initiation site identification.
     BMC Bioinformatics, 11, 119.
     https://doi.org/10.1186/1471-2105-11-119

[3]  Camacho, C. et al. (2009). BLAST+: architecture and applications.
     BMC Bioinformatics, 10, 421.
     https://doi.org/10.1186/1471-2105-10-421

[4]  Eddy, S.R. (2011). Accelerated Profile HMM Searches.
     PLoS Computational Biology, 7(10), e1002195.
     https://doi.org/10.1371/journal.pcbi.1002195

[5]  Abramson, J. et al. (2024). Accurate structure prediction of
     biomolecular interactions with AlphaFold 3.
     Nature, 630(8016), 493-500.
     https://doi.org/10.1038/s41586-024-07487-w
     ► Used for: AF3 JSON generation, contact_probs, ranking_score.

[6]  Hunter, J.D. (2007). Matplotlib: A 2D graphics environment.
     Computing in Science & Engineering, 9(3), 90-95.
     https://doi.org/10.1109/MCSE.2007.55

[7]  Harris, C.R. et al. (2020). Array programming with NumPy.
     Nature, 585, 357-362.
     https://doi.org/10.1038/s41586-020-2649-2

INTERACTION SCORING METHODOLOGY
────────────────────────────────
[8]  Evans, R. et al. (2022). Protein complex prediction with
     AlphaFold-Multimer. bioRxiv 2021.10.04.463034.
     https://doi.org/10.1101/2021.10.04.463034
     ► Defines ipTM, chain_pair_iptm, and chain_pair_pae_min.
       These are extracted directly from AF3 summary_confidences.json.
       ipTM > 0.75 is the canonical high-confidence threshold.

[9]  Jumper, J. et al. (2021). Highly accurate protein structure
     prediction with AlphaFold. Nature, 596, 583-589.
     https://doi.org/10.1038/s41586-021-03819-2
     ► Defines PAE (Predicted Aligned Error): the expected positional
       error of residue j when the structure is aligned on residue i.
       Off-diagonal PAE < 5 Å = confident relative positioning of
       two residues on different chains.

[10] Humphreys, I.R. et al. (2021). Computed structures of core
     eukaryotic protein complexes. Science, 374, eabm4805.
     https://doi.org/10.1126/science.abm4805
     ► Demonstrates proteome-scale ipTM-based PPI screening.
       Methodology basis for ppigFinder's batch AF3 workflow.

[11] Bryant, P., Pozzati, G., & Elofsson, A. (2022). Improved prediction
     of protein-protein interactions using AlphaFold2.
     Nature Communications, 13, 1265.
     https://doi.org/10.1038/s41467-022-28865-w
     ► Systematic analysis of PAE inter-chain metrics for PPI
       confidence scoring. Validates the use of min PAE (rather
       than mean) as a more sensitive detector of focal interfaces.

[12] Mirdita, M. et al. (2022). ColabFold: making protein folding
     accessible to all. Nature Methods, 19, 679-682.
     https://doi.org/10.1038/s41592-022-01488-1
     ► ColabFold FASTA format used by ppigFinder for batch MSA
       generation and local AF2-Multimer runs.

MOTIF DETECTION
───────────────
[13] Virtanen, P. et al. (2020). SciPy 1.0: Fundamental algorithms for
     scientific computing in Python. Nature Methods, 17, 261-272.
     https://doi.org/10.1038/s41592-019-0686-2
     ► scipy.ndimage used for connected-component labelling and
       morphological operations in the motif detector (v1.16+).
       A pure-NumPy BFS fallback is used when SciPy is absent.

NOVEL METHODOLOGY — ppigFinder v1.18: Focal Interaction Metrics
────────────────────────────────────────────────────────────────
The standard PAE_inter metric (mean of the entire off-diagonal PAE
quadrant) is mathematically correct but biologically misleading for
proteins that interact through a single domain embedded in a longer,
disordered sequence. For example:

  ORF2588 (431 aa, XVIPCD domain at N-terminus) × ORF2589 (267 aa):
    PAE_inter_global = 22.0 Å  → reported as "not relevant"
    PAE_min          =  2.1 Å  → actual focal contact confirmed

Only 6.6% of residue pairs are in contact; the global mean is dominated
by 93.4% of pairs involving ORF2588 residues outside the domain.

ppigFinder v1.18 introduces three focal metrics extracted directly
from AF3 summary_confidences.json (no PAE matrix parsing required):

  PAE_min ★  = min(chain_pair_pae_min[A][B], chain_pair_pae_min[B][A])
               Minimum PAE in the off-diagonal quadrant.
               < 4 Å = at least one contact point predicted with
               near-atomic positional certainty.

  cp_ipTM ★  = chain_pair_iptm[A][B]  (off-diagonal)
               Chain-pair ipTM removing intra-chain folding contribution.

  Contact%   = fraction of residue pairs with PAE < 5 Å
               Even small values (5–10%) confirm a real interface.

Classification rule (analogous to ipTM > 0.75 of Evans et al. 2022
but applicable to domain-limited contacts):
  PAE_min < 4 Å  AND  cp_ipTM ≥ 0.50  →  HIGH confidence
  PAE_min 4–8 Å                         →  MED (check motifs)
  PAE_min ≥ 8 Å                         →  LOW

This methodology is implemented in ppigFinder v1.18 and validated
on T4SS toxin-antitoxin pairs (XVIPCD × xac2610, VirB11 × VirB10
homologs) from Xanthomonas axonopodis pv. citri MAFF311018 genome.
""")
        lay.addWidget(txt)
        btns = QDialogButtonBox(QDialogButtonBox.StandardButton.Close
                                if QT_VERSION == 6 else QDialogButtonBox.Close)
        btns.rejected.connect(dlg.close)
        lay.addWidget(btns)
        dlg.exec() if QT_VERSION == 6 else dlg.exec_()

    def _show_about(self):
        blast_ok = "✅" if BACKENDS.get('blast+',{}).get('available') else "❌"
        hmmer_ok = "✅" if BACKENDS.get('hmmer3',{}).get('available') else "❌"
        pyrod_ok = "✅" if BACKENDS.get('pyrodigal',{}).get('available') else "❌"
        mat_ok   = "✅" if MATPLOTLIB_AVAILABLE else "❌"
        n_res    = len(getattr(self, '_af3_analysis_results', []))
        n_jobs   = len(getattr(self, 'af3_jobs', []))
        n_orfs   = len(getattr(self, 'orfs', []))
        QMessageBox.about(self, t('about'),
            f"🧬 ppigFinder — Protein-Protein Interaction Genomic Finder\n"
            f"Version 2.0 — v2.0  |  MIT License\n\n"
            f"Discovery of novel bacterial PPIs via ORF prediction\n"
            f"(Pyrodigal / 6-frame scan), HMM/BLAST annotation,\n"
            f"genomic neighbourhood analysis, and AlphaFold 3\n"
            f"structural interaction prediction.\n\n"
            f"Key features (v2.0):\n"
            f"  • Interaction Results: PAE_min ★, cp_ipTM ★, Contact% ★\n"
            f"    focal metrics + HIGH / MED / LOW confidence classification\n"
            f"  • Column-header explanations on hover (all result tables)\n"
            f"  • Main ORF table shows PAE_min ★ and cp_ipTM ★ per ORF\n"
            f"  • Two-tier interaction motif detector (v1.16+)\n"
            f"  • Genomic PPI Map: arc diagram over genome backbone (v2.0)\n"
            f"  • Help → References: full methodology & 14 citations\n\n"
            f"Current session: {n_orfs:,} ORFs  |  {n_jobs:,} AF3 jobs  |  "
            f"{n_res:,} results loaded\n\n"
            f"https://github.com/goka-lab/ppigfinder\n\n"
            f"Backends:\n"
            f"  BLAST+     {blast_ok}\n"
            f"  HMMER3     {hmmer_ok}\n"
            f"  Pyrodigal  {pyrod_ok}\n"
            f"  Matplotlib {mat_ok}\n\n"
            f"PyQt{QT_VERSION} | Python {sys.version.split()[0]}\n\n"
            f"Cite: ppigFinder (2026) github.com/goka-lab/ppigfinder\n"
            f"Refs: Evans 2022; Abramson 2024; Jumper 2021; Bryant 2022")



    # ═══════════════════════════════════════════════════════════
    # AF3 ANALYSIS TAB  (v2)
    # ═══════════════════════════════════════════════════════════
    #
    # v1.15 fixes (AlphaFold Analysis tab):
    #   • Bug: row selection after column sort loaded the wrong job's
    #     data into the plots. Fixed by storing the result-list index
    #     as UserRole data on column 0 of each row.
    #   • Bug: numeric columns (ipTM, ptm, pLDDT, ranking, PAE_inter)
    #     were sorted lexicographically ("100.0" < "2.5"). Fixed via
    #     a _NumericItem subclass with numeric __lt__.
    #   • Bug: pLDDT per-residue array was read from 'atom_plddts'
    #     (per-atom; wrong length). Added 'token_plddts' to the key
    #     list and a length sanity check vs. chain_ids.
    #   • Bug: contact markers were drawn as a vertical line at the
    #     centre of each off-diagonal block. Now drawn at the actual
    #     contact residue positions in both AB and BA quadrants.
    #   • Bug: changing the contact threshold did not refresh the
    #     "Best contact pair" column in the table. Fixed.
    #   • Added: TSV export of the results table.
    #   • Added: min-ipTM / max-PAEinter filter spinboxes.
    #   • Added: double-click a row to open the job folder.
    #   • Added: user-visible dialog for jobs that failed to parse.

    @staticmethod
    def _af3a_user_role():
        """Qt5/Qt6-agnostic UserRole constant."""
        try:
            return Qt.ItemDataRole.UserRole      # PyQt6
        except AttributeError:
            return Qt.UserRole                   # PyQt5

    def _create_af3_analysis_tab(self):
        """PAE heatmaps (ChimeraX style) + pLDDT + inter-chain metrics.
        Supports any number of chains. Interactive hover shows PAE value."""
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(4, 4, 4, 4)
        lay.setSpacing(4)

        if not MATPLOTLIB_AVAILABLE:
            warn = QLabel(
                "⚠  matplotlib not installed.\n"
                "   pip install matplotlib numpy")
            warn.setStyleSheet(
                "background:#FFF3CD;color:#856404;padding:8px;border-radius:4px;")
            lay.addWidget(warn)

        # ── Top toolbar ──────────────────────────────────────────
        top = QHBoxLayout()
        btn_load = QPushButton("📂 Load AF3 results folder")
        btn_load.clicked.connect(self._af3a_load_folder)
        top.addWidget(btn_load)

        btn_clear = QPushButton("🗑 Clear")
        btn_clear.clicked.connect(self._af3a_clear)
        top.addWidget(btn_clear)

        top.addWidget(QLabel("  Contact threshold:"))
        self._af3a_thresh_spin = QDoubleSpinBox()
        self._af3a_thresh_spin.setRange(1.0, 20.0)
        self._af3a_thresh_spin.setValue(5.0)
        self._af3a_thresh_spin.setSingleStep(1.0)
        self._af3a_thresh_spin.setSuffix(" Å")
        self._af3a_thresh_spin.setFixedWidth(80)
        self._af3a_thresh_spin.setToolTip("")
        self._af3a_thresh_spin.valueChanged.connect(
            lambda _: self._af3a_replot_selected())
        top.addWidget(self._af3a_thresh_spin)

        # Filters: hide rows with low ipTM / high PAE_inter
        top.addWidget(QLabel("  min ipTM:"))
        self._af3a_filter_iptm = QDoubleSpinBox()
        self._af3a_filter_iptm.setRange(0.0, 1.0)
        self._af3a_filter_iptm.setSingleStep(0.05)
        self._af3a_filter_iptm.setDecimals(2)
        self._af3a_filter_iptm.setValue(0.0)
        self._af3a_filter_iptm.setFixedWidth(60)
        self._af3a_filter_iptm.setToolTip("")
        self._af3a_filter_iptm.valueChanged.connect(
            lambda _: self._af3a_apply_filters())
        top.addWidget(self._af3a_filter_iptm)

        top.addWidget(QLabel("max PAE_inter:"))
        self._af3a_filter_paei = QDoubleSpinBox()
        self._af3a_filter_paei.setRange(0.0, 32.0)
        self._af3a_filter_paei.setSingleStep(1.0)
        self._af3a_filter_paei.setDecimals(1)
        self._af3a_filter_paei.setValue(32.0)
        self._af3a_filter_paei.setSuffix(" Å")
        self._af3a_filter_paei.setFixedWidth(80)
        self._af3a_filter_paei.setToolTip("")
        self._af3a_filter_paei.valueChanged.connect(
            lambda _: self._af3a_apply_filters())
        top.addWidget(self._af3a_filter_paei)

        btn_pdf = QPushButton("📄 Export plots PDF")
        btn_pdf.setToolTip("Export all visible PAE/pLDDT plots to a PDF file")
        btn_pdf.clicked.connect(self._af3a_export_pdf)
        top.addWidget(btn_pdf)

        btn_tsv = QPushButton("⬇ Export TSV")
        btn_tsv.setToolTip("Export the full results table (all metrics) to TSV")
        btn_tsv.clicked.connect(self._af3a_export_tsv)
        top.addWidget(btn_tsv)

        top.addStretch()
        self._af3a_status = QLabel("No results loaded")
        self._af3a_status.setStyleSheet("font-size:11px;color:#666;")
        top.addWidget(self._af3a_status)
        lay.addLayout(top)

        # ── Motif detection toolbar (v1.16) ────────────────────────
        # A second compact toolbar dedicated to the 2-D interaction
        # motif detector (see _af3a_detect_motifs). Spinbox changes
        # retrigger motif detection on the currently loaded jobs and
        # refresh the motif table + plot overlays.
        motif_bar = QHBoxLayout()
        motif_bar.setSpacing(4)

        mlabel = QLabel("🎯 Motifs")
        mlabel.setStyleSheet(
            "font-weight:bold;color:#37474F;padding:0 4px 0 2px;")
        motif_bar.addWidget(mlabel)

        motif_bar.addWidget(QLabel("core PAE<"))
        self._af3a_motif_core_spin = QDoubleSpinBox()
        self._af3a_motif_core_spin.setRange(1.0, 20.0)
        self._af3a_motif_core_spin.setValue(self._AF3_MOTIF_PAE_CORE)
        self._af3a_motif_core_spin.setSingleStep(0.5)
        self._af3a_motif_core_spin.setSuffix(" Å")
        self._af3a_motif_core_spin.setFixedWidth(70)
        self._af3a_motif_core_spin.setToolTip("")
        motif_bar.addWidget(self._af3a_motif_core_spin)

        motif_bar.addWidget(QLabel("ext PAE<"))
        self._af3a_motif_ext_spin = QDoubleSpinBox()
        self._af3a_motif_ext_spin.setRange(2.0, 30.0)
        self._af3a_motif_ext_spin.setValue(self._AF3_MOTIF_PAE_EXT)
        self._af3a_motif_ext_spin.setSingleStep(1.0)
        self._af3a_motif_ext_spin.setSuffix(" Å")
        self._af3a_motif_ext_spin.setFixedWidth(70)
        self._af3a_motif_ext_spin.setToolTip("")
        motif_bar.addWidget(self._af3a_motif_ext_spin)

        motif_bar.addWidget(QLabel("min contact≥"))
        self._af3a_motif_contact_spin = QDoubleSpinBox()
        self._af3a_motif_contact_spin.setRange(0.00, 1.00)
        self._af3a_motif_contact_spin.setValue(self._AF3_MOTIF_MIN_CONTACT)
        self._af3a_motif_contact_spin.setSingleStep(0.01)
        self._af3a_motif_contact_spin.setDecimals(2)
        self._af3a_motif_contact_spin.setFixedWidth(65)
        self._af3a_motif_contact_spin.setToolTip("")
        motif_bar.addWidget(self._af3a_motif_contact_spin)

        motif_bar.addWidget(QLabel("min size"))
        self._af3a_motif_size_spin = QSpinBox()
        self._af3a_motif_size_spin.setRange(2, 50)
        self._af3a_motif_size_spin.setValue(self._AF3_MOTIF_MIN_SIZE)
        self._af3a_motif_size_spin.setFixedWidth(50)
        self._af3a_motif_size_spin.setToolTip("")
        motif_bar.addWidget(self._af3a_motif_size_spin)

        self._af3a_motif_recip_cb = QCheckBox("reciprocal")
        self._af3a_motif_recip_cb.setChecked(True)
        self._af3a_motif_recip_cb.setToolTip("")
        motif_bar.addWidget(self._af3a_motif_recip_cb)

        btn_motifs_refresh = QPushButton("Rerun")
        btn_motifs_refresh.setToolTip(
            "Re-run motif detection on all loaded jobs with the "
            "current parameters.")
        btn_motifs_refresh.clicked.connect(self._af3a_rerun_motifs)
        motif_bar.addWidget(btn_motifs_refresh)

        btn_motifs_tsv = QPushButton("⬇ Motifs TSV")
        btn_motifs_tsv.setToolTip("Export all detected motifs to TSV.")
        btn_motifs_tsv.clicked.connect(self._af3a_export_motifs_tsv)
        motif_bar.addWidget(btn_motifs_tsv)

        motif_bar.addStretch()
        self._af3a_motif_count_lbl = QLabel("—")
        self._af3a_motif_count_lbl.setStyleSheet(
            "font-size:11px;color:#455A64;padding:0 4px;")
        motif_bar.addWidget(self._af3a_motif_count_lbl)

        lay.addLayout(motif_bar)

        # ── Explanatory info label (like AF3 mode description) ─────────
        # Shows a context-sensitive description based on which filter
        # or motif control the user last hovered over. Default text
        # summarises the confidence classification rules.
        self._af3a_info_lbl = QLabel(
            "<b>Confidence classification:</b>  "
            "<span style='color:#085041'>HIGH ★</span> = PAE_min &lt; 4 Å <i>and</i> cp_ipTM &ge; 0.50 (bold rows) &nbsp;|&nbsp; "
            "<span style='color:#633806'>MED</span> = PAE_min 4–8 Å &nbsp;|&nbsp; "
            "<span style='color:#A32D2D'>LOW</span> = PAE_min &ge; 8 Å &nbsp;|&nbsp; "
            "★ = focal metrics from <tt>chain_pair_pae_min / chain_pair_iptm</tt> "
            "(AF3 <tt>summary_confidences.json</tt>). "
            "Hover any control above for details.")
        self._af3a_info_lbl.setWordWrap(True)
        self._af3a_info_lbl.setTextFormat(
            Qt.TextFormat.RichText if QT_VERSION == 6 else Qt.RichText)
        self._af3a_info_lbl.setStyleSheet(
            "background:#fffde7;color:#555;"
            "border:1px solid #f9a825;border-radius:4px;"
            "padding:5px 8px;font-size:11px;")
        lay.addWidget(self._af3a_info_lbl)


        def _make_hover(widget, msg):
            """Install an eventFilter that shows msg in info_lbl on Enter."""
            class _HE(QObject):
                def eventFilter(self_, obj, ev):
                    try:
                        from PyQt6.QtCore import QEvent as QEv
                        enter = QEv.Type.Enter
                    except (ImportError, AttributeError):
                        from PyQt5.QtCore import QEvent as QEv
                        enter = QEv.Enter
                    try:
                        from PyQt6.QtCore import QEvent as QEv2
                        leave = QEv2.Type.Leave
                    except (ImportError, AttributeError):
                        from PyQt5.QtCore import QEvent as QEv2
                        leave = QEv2.Leave
                    if ev.type() == enter:
                        self._af3a_info_lbl.setText(msg)
                    elif ev.type() == leave:
                        self._af3a_info_lbl.setText(
                            "<b>Confidence classification:</b>  "
                            "<span style='color:#085041'>HIGH ★</span> = PAE_min &lt; 4 Å <i>and</i> cp_ipTM &ge; 0.50 &nbsp;|&nbsp; "
                            "<span style='color:#633806'>MED</span> = PAE_min 4–8 Å &nbsp;|&nbsp; "
                            "<span style='color:#A32D2D'>LOW</span> = PAE_min &ge; 8 Å &nbsp;|&nbsp; "
                            "Hover any control for details.")
                    return False
            fe = _HE(widget)
            widget.installEventFilter(fe)
            widget._hover_filter = fe   # keep reference

        _make_hover(self._af3a_thresh_spin,
            "<b>Contact threshold (Å)</b> — Residues with inter-chain PAE below this value are counted as direct contacts. "
            "Typical: <b>5 Å</b> = atomic contact; 10 Å = proximity. "
            "⚠ This does NOT change PAE_min or ipTM — those come from AF3 summary JSON directly. "
            "Affects: contact-region string, n_contacts display, and pLDDT plot markers.")

        _make_hover(self._af3a_filter_iptm,
            "<b>min ipTM filter</b> — Hide rows where ipTM &lt; this value. "
            "ipTM = interface pTM: measures how well AF3 predicts the relative positions of chains A and B (0–1). "
            "<b>&gt; 0.75</b> = high confidence &nbsp;|&nbsp; "
            "<b>0.50–0.75</b> = likely interaction (check PAE_min ★) &nbsp;|&nbsp; "
            "<b>&lt; 0.50</b> = uncertain. Set to 0.00 to show all. "
            "Reference: Evans et al. 2022 (AlphaFold-Multimer).")

        _make_hover(self._af3a_filter_paei,
            "<b>max PAE_inter filter</b> — Hide rows where global PAE_inter &gt; this value. "
            "PAE_inter = mean PAE over the ENTIRE inter-chain quadrant. "
            "⚠ <b>Warning:</b> this metric is diluted by disordered regions. "
            "A protein interacting through a single domain may show PAE_inter = 22 Å "
            "even though PAE_min = 2.1 Å (real focal contact). "
            "Prefer filtering by the <b>PAE_min ★</b> column. Set to 32 Å to show all.")

        _make_hover(self._af3a_motif_core_spin,
            "<b>core PAE (Å)</b> — Tier-1 threshold for the interaction motif detector. "
            "Pixels in the inter-chain PAE matrix below this value form the motif core (dark blue in heatmap). "
            "Typical: <b>5–8 Å</b>. Lower = stricter, fewer but higher-confidence motifs. "
            "Example: core = 8 Å detects all contacts within ~8 Å of predicted position.")

        _make_hover(self._af3a_motif_ext_spin,
            "<b>ext PAE (Å)</b> — Tier-2 (extension) threshold. "
            "After the core is found, the bounding box is expanded to include nearby pixels below this softer threshold. "
            "Must be <b>greater than core PAE</b>. Typical: <b>12–18 Å</b>. "
            "Larger ext PAE = larger motif regions reported.")

        _make_hover(self._af3a_motif_contact_spin,
            "<b>min contact ≥</b> — Minimum AF3 <tt>contact_probs</tt> score for a pixel to join the core mask. "
            "<tt>contact_probs</tt> is AF3's internal probability that two residues are in physical contact, "
            "independent of PAE. Using both together rejects 'aligned-but-far' artefacts. "
            "Set to <b>0.00</b> to disable this filter (use PAE alone).")

        _make_hover(self._af3a_motif_size_spin,
            "<b>min size</b> — Minimum bounding-box dimension (residues) for a motif to survive. "
            "Motifs smaller than this on either axis are discarded as noise. "
            "Recommended: <b>3–5</b> for domain interactions; <b>2</b> for peptide binding sites. "
            "Example: min size = 5 means the motif must span ≥ 5 residues on EACH chain.")

        _make_hover(self._af3a_motif_recip_cb,
            "<b>reciprocal</b> — Require the motif to appear in BOTH PAE quadrants: "
            "A→B (rows=A, cols=B) AND B→A (rows=B, cols=A) with ≥ 50% pixel overlap. "
            "A real interaction shows low PAE in both orientations. "
            "A ghost motif (disordered alignment artefact) appears in only one quadrant. "
            "<b>Strongly recommended: keep checked.</b> Uncheck only to investigate marginal cases.")
        # on already-loaded jobs and refresh the motif table + overlays.
        self._af3a_motif_core_spin.valueChanged.connect(
            lambda _: self._af3a_rerun_motifs())
        self._af3a_motif_ext_spin.valueChanged.connect(
            lambda _: self._af3a_rerun_motifs())
        self._af3a_motif_contact_spin.valueChanged.connect(
            lambda _: self._af3a_rerun_motifs())
        self._af3a_motif_size_spin.valueChanged.connect(
            lambda _: self._af3a_rerun_motifs())
        self._af3a_motif_recip_cb.stateChanged.connect(
            lambda _: self._af3a_rerun_motifs())

        # ── Hover value label ─────────────────────────────────────
        self._af3a_hover_lbl = QLabel("")
        self._af3a_hover_lbl.setStyleSheet(
            "font-size:10px;font-family:monospace;color:#1565C0;"
            "padding:2px 6px;background:#E3F2FD;border-radius:3px;")
        lay.addWidget(self._af3a_hover_lbl)

        # ── Splitter: job table (top) / plots scroll (bottom) ─────
        splitter = QSplitter(Vertical)

        # Job table
        tg = QGroupBox("AF3 prediction results")
        tg_l = QVBoxLayout(tg)
        tg_l.setContentsMargins(4, 4, 4, 4)
        self._af3a_table = QTableWidget()
        self._af3a_table.setColumnCount(12)
        self._af3a_table.setHorizontalHeaderLabels([
            'Job name', 'Chains (n)', 'Chain IDs',
            'ipTM', 'ptm', 'mean_pLDDT', 'ranking_score',
            'PAE_inter (Å)',
            'PAE_min ★ (Å)', 'cp_ipTM ★', 'Contact% ★',
            'Best contact pair'])

        # Standard Qt column-header tooltips — shown as native balloon when
        # hovering a header cell. Does NOT replace QHeaderView, so sorting
        # by clicking column headers is preserved.
        _HDR_TIPS = {
            0: "Job name — AF3 prediction folder identifier.",
            1: "Chains (n) — number of protein chains (usually 2).",
            2: "Chain IDs — ORF names and sizes.",
            3: ("ipTM — Interface predicted TM-score (0–1). "
                "Measures the accuracy of predicted relative chain positions. "
                "> 0.75 = high confidence  |  0.50–0.75 = likely  |  < 0.50 = uncertain. "
                "[Evans et al. 2022, AlphaFold-Multimer]"),
            4: ("ptm — Predicted TM-score for the whole complex. "
                "Reflects overall folding quality; less informative than ipTM for PPI."),
            5: ("mean_pLDDT — Mean per-residue confidence (0–100). "
                "> 70 = good structure quality  |  < 50 = likely disordered."),
            6: ("ranking_score — AF3 combined confidence metric. "
                "Used for ranking when multiple diffusion samples are available."),
            7: ("PAE_inter (Å) — GLOBAL mean PAE of the entire inter-chain quadrant. "
                "⚠ Warning: diluted by disordered regions. "
                "Use PAE_min ★ as primary metric for domain-limited interactions."),
            8: ("PAE_min ★ (Å) — Minimum PAE in the off-diagonal quadrant "
                "(chain_pair_pae_min from summary_confidences.json). "
                "< 4 Å = at least one contact with near-atomic confidence, "
                "regardless of global disorder. KEY metric for domain-limited PPI. "
                "[Evans 2022; Bryant 2022]"),
            9: ("cp_ipTM ★ — Chain-pair ipTM for the interface only "
                "(chain_pair_iptm off-diagonal). "
                "≥ 0.50 combined with PAE_min < 4 Å = HIGH confidence focal interaction."),
            10: ("Contact% ★ — Fraction of residue pairs with PAE < 5 Å. "
                 "Even 5–10% confirms a real focal interface."),
            11: ("Best contact pair — Focal residue ranges where PAE < 5 Å. "
                 "Based on per-residue minimum PAE (not mean)."),
        }
        for _col, _tip in _HDR_TIPS.items():
            _item = self._af3a_table.horizontalHeaderItem(_col)
            if _item is not None:
                _item.setToolTip(_tip)

        self._af3a_table.setSelectionBehavior(SelectRows)
        self._af3a_table.setSelectionMode(SingleSelection)
        self._af3a_table.setEditTriggers(
            QAbstractItemView.EditTrigger.NoEditTriggers
            if QT_VERSION == 6 else QAbstractItemView.NoEditTriggers)
        self._af3a_table.horizontalHeader().setStretchLastSection(True)
        self._af3a_table.setAlternatingRowColors(True)
        self._af3a_table.setSortingEnabled(True)
        for i, cw in enumerate([200, 60, 120, 55, 55, 70, 90, 85, 250]):
            self._af3a_table.setColumnWidth(i, cw)
        self._af3a_table.selectionModel().selectionChanged.connect(
            self._af3a_on_select)
        self._af3a_table.cellDoubleClicked.connect(
            self._af3a_on_double_click)
        tg_l.addWidget(self._af3a_table)
        splitter.addWidget(tg)

        # Scrollable plot area
        self._af3a_scroll = QScrollArea()
        self._af3a_scroll.setWidgetResizable(True)
        self._af3a_scroll.setFrameShape(
            QFrame.Shape.NoFrame if QT_VERSION == 6 else QFrame.NoFrame)
        self._af3a_plot_container = QWidget()
        self._af3a_plot_layout   = QVBoxLayout(self._af3a_plot_container)
        self._af3a_plot_layout.setSpacing(10)
        self._af3a_scroll.setWidget(self._af3a_plot_container)
        splitter.addWidget(self._af3a_scroll)
        splitter.setSizes([240, 600])
        lay.addWidget(splitter)

        # Keep refs to active canvases for PDF export
        self._af3a_active_canvases = []

        self._tabs.addTab(w, "📈 Interaction Results")

    # ── AF3 Analysis helpers ─────────────────────────────────────

    def _af3a_load_folder(self):
        folder = QFileDialog.getExistingDirectory(
            self, "Select AF3 results folder")
        if not folder:
            return
        self._af3_analysis_dir = folder
        self._af3a_scan_folder(folder)

    def _af3a_clear(self):
        self._af3_analysis_results = []
        self._af3a_table.setRowCount(0)
        self._af3a_clear_plots()
        self._af3a_status.setText("Cleared")

    # ── AF3 file classification (content-based) ───────────────
    #
    # v1.15: classify JSON files by their top-level keys rather than
    # their filename, because the AF3 server sometimes writes files
    # whose filenames do not clearly mark their role. Also keeps the
    # large input-data JSON (which can be 80+ MB of MSA / sequences)
    # from being loaded into memory during the scan.

    # Signatures (top-level keys) used to identify each JSON type.
    _AF3_SUMMARY_SIGNATURE_KEYS  = {'iptm', 'ptm', 'ranking_score'}
    _AF3_CONF_SIGNATURE_KEYS     = {'pae', 'atom_plddts', 'token_chain_ids'}
    _AF3_DATA_SIGNATURE_KEYS     = {'sequences', 'modelSeeds', 'dialect'}

    # Files larger than this are skipped without opening (protects us
    # from accidentally loading the input-data JSON with its full MSA).
    _AF3_MAX_JSON_SCAN_BYTES = 200 * 1024 * 1024   # 200 MB

    @classmethod
    def _af3a_classify_json(cls, path: Path) -> str:
        """Classify a JSON file by reading only its top-level keys.

        Returns one of:
          'summary'      — AF3 summary_confidences (small, has ipTM/pTM)
          'confidences'  — AF3 full confidences (PAE + pLDDT + chains)
          'data'         — AF3 input data (sequences, MSA — skip)
          'unknown'      — any other JSON

        The check is streaming-friendly: for small files (< 1 MB) we
        just json.load; for larger files we only probe with ijson if
        available, falling back to a size-based heuristic.
        """
        try:
            size = path.stat().st_size
        except OSError:
            return 'unknown'
        if size > cls._AF3_MAX_JSON_SCAN_BYTES:
            # Too large — almost certainly the input-data JSON.
            return 'data'

        # Small files: just load them
        if size < 1 * 1024 * 1024:
            try:
                with open(path, encoding='utf-8') as f:
                    obj = json.load(f)
                if not isinstance(obj, dict):
                    return 'unknown'
                keys = set(obj.keys())
                if cls._AF3_DATA_SIGNATURE_KEYS.issubset(keys):
                    return 'data'
                if cls._AF3_CONF_SIGNATURE_KEYS & keys:
                    # Summary usually small; full confidences usually
                    # huge. The confidence signature keys are also
                    # rare in summaries, so a hit here is decisive.
                    if cls._AF3_CONF_SIGNATURE_KEYS.issubset(keys):
                        return 'confidences'
                if cls._AF3_SUMMARY_SIGNATURE_KEYS.issubset(keys):
                    return 'summary'
                return 'unknown'
            except (OSError, json.JSONDecodeError):
                return 'unknown'

        # Larger files (> 1 MB): try to stream top-level keys via ijson
        try:
            import ijson
            keys_seen = set()
            with open(path, 'rb') as f:
                for prefix, event, _val in ijson.parse(f):
                    if (prefix and '.' not in prefix
                            and event in ('start_array', 'start_map',
                                          'string', 'number',
                                          'boolean', 'null')):
                        keys_seen.add(prefix)
                    # Stop once we've seen enough to classify
                    if len(keys_seen) >= 12:
                        break
            if cls._AF3_DATA_SIGNATURE_KEYS.issubset(keys_seen):
                return 'data'
            if cls._AF3_CONF_SIGNATURE_KEYS.issubset(keys_seen):
                return 'confidences'
            if cls._AF3_SUMMARY_SIGNATURE_KEYS.issubset(keys_seen):
                return 'summary'
            return 'unknown'
        except ImportError:
            # ijson not available → fall back to size heuristic: a file
            # > 1 MB in an AF3 job folder is almost always the full
            # confidences (PAE matrix is O(n²)).
            return 'confidences'
        except (OSError, ValueError):
            return 'unknown'

    def _af3a_find_job_files(self, job_dir: Path) -> dict:
        """Find AF3 output files inside a job directory.

        Returns a dict with keys:
          'summary'       : Path or None   — summary_confidences.json
          'confidences'   : Path or None   — full confidences.json
          'model'         : Path or None   — model .cif
          'ranking_csv'   : Path or None   — ranking_scores.csv
          'sample_dirs'   : list[Path]     — seed-*_sample-* subfolders

        Uses content-based classification for all .json files in the
        directory, so files whose names differ from the standard AF3
        server convention are still recognised.
        """
        result = {
            'summary': None, 'confidences': None,
            'model': None, 'ranking_csv': None,
            'input': None,         # ← NEW: ppigFinder / AF3 server input JSON
            'sample_dirs': [],
        }
        try:
            entries = list(job_dir.iterdir())
        except OSError:
            return result

        # Classify .json files by content; prefer summary over
        # confidences when multiple candidates exist (rare).
        for p in entries:
            if not p.is_file():
                continue
            name_lc = p.name.lower()
            if name_lc.endswith('.json'):
                kind = self._af3a_classify_json(p)
                if kind == 'summary' and result['summary'] is None:
                    result['summary'] = p
                elif kind == 'confidences' and result['confidences'] is None:
                    result['confidences'] = p
                elif kind == 'data' and result['input'] is None:
                    result['input'] = p   # ppigFinder-generated input JSON
                # 'unknown' → ignored
            elif name_lc.endswith('.cif'):
                # Skip known non-model CIFs (there usually aren't any)
                if result['model'] is None:
                    result['model'] = p
            elif name_lc == 'ranking_scores.csv':
                result['ranking_csv'] = p

        # Collect sample sub-dirs (detail views, not primary jobs)
        for p in entries:
            if p.is_dir() and re.match(r'seed-\d+_sample-\d+', p.name):
                result['sample_dirs'].append(p)
        result['sample_dirs'].sort()

        return result

    def _af3a_scan_folder(self, folder: str):
        """Scan a folder tree for AF3 jobs.

        v1.15 (content-based):
          • Each candidate directory is inspected for AF3 output files
            via content classification (_af3a_classify_json), not file-
            name globs. This handles AF3 server outputs whose naming
            differs from the strict _summary_confidences.json /
            _confidences.json convention.
          • Huge input-data JSONs are filtered out by size BEFORE being
            opened (they contain MSAs and would otherwise bloat memory).
          • seed-*_sample-* subfolders are treated as the per-sample
            detail of their parent job (the best model is already
            promoted to the parent), not as independent jobs.
          • Parse errors are surfaced in a single dialog at the end.
          • A progress dialog shows scan status for large folder trees.
          • The top-scoring row is auto-selected after populate.
        """
        root = Path(folder)
        try:
            top_level = sorted(root.iterdir())
        except OSError as e:
            QMessageBox.critical(
                self, "AF3 Analysis",
                f"Could not read folder:\n{folder}\n\n{e}")
            return

        # Build list of candidate job dirs:
        #   • each direct subdir whose immediate children look like an
        #     AF3 job (contain at least one summary JSON)
        #   • or the root itself if it directly contains AF3 outputs
        #   • sample subdirs (seed-*_sample-*) are explicitly skipped
        #     at the top level; the parent is scanned instead
        candidates = []
        for p in top_level:
            if not p.is_dir():
                continue
            if re.match(r'seed-\d+_sample-\d+', p.name):
                # Sample folders under the root (e.g. user pointed at a
                # single-job dir itself) → skip; the root will handle.
                continue
            candidates.append(p)
        if not candidates:
            candidates = [root]

        # Progress dialog
        try:
            if QT_VERSION == 6:
                from PyQt6.QtWidgets import QProgressDialog
            else:
                from PyQt5.QtWidgets import QProgressDialog
            progress = QProgressDialog(
                "Scanning AF3 jobs…", "Cancel",
                0, len(candidates), self)
            progress.setWindowTitle("AF3 Analysis")
            progress.setMinimumDuration(400)
            progress.setValue(0)
        except Exception:
            progress = None

        jobs = []
        errors  = []  # list of (dir_name, error_message)
        skipped = 0   # dirs that are not AF3 jobs

        for i, job_dir in enumerate(candidates):
            if progress is not None:
                progress.setValue(i)
                progress.setLabelText(
                    f"Scanning AF3 jobs… ({i}/{len(candidates)})")
                QApplication.processEvents()
                if progress.wasCanceled():
                    break

            files = self._af3a_find_job_files(job_dir)
            if files['summary'] is None:
                # Not an AF3 job folder — skip silently
                skipped += 1
                continue

            try:
                parsed = self._af3a_parse_job(
                    job_dir,
                    files['summary'],
                    files['confidences'],
                    model_cif=files['model'],
                    ranking_csv=files['ranking_csv'],
                    input_json=files['input'])
                if parsed:
                    jobs.append(parsed)
            except Exception as e:
                errors.append((job_dir.name, str(e)))
                print(f"[AF3 Analysis] {job_dir.name}: {e}")

        if progress is not None:
            progress.setValue(len(candidates))

        self._af3_analysis_results = jobs
        self._af3a_populate_table()

        status = f"{len(jobs)} job(s) loaded from {Path(folder).name}/"
        if errors:  status += f"  ⚠ {len(errors)} failed"
        if skipped: status += f"  ({skipped} non-AF3 dirs skipped)"
        self._af3a_status.setText(status)
        self._update_orfs_list()

        if self._af3a_table.rowCount() > 0:
            self._af3a_table.selectRow(0)

        if errors:
            msg = "\n".join(f"• {name}: {err}" for name, err in errors[:15])
            if len(errors) > 15:
                msg += f"\n… and {len(errors) - 15} more."
            QMessageBox.warning(
                self, "AF3 Analysis — partial load",
                f"{len(errors)} job(s) could not be parsed:\n\n{msg}")

    def _af3a_parse_job(self, job_dir: Path, sum_path: Path,
                         conf_path, model_cif: Path = None,
                         ranking_csv: Path = None,
                         input_json: Path = None) -> dict:
        """Parse one AF3 job. Handles any number of chains.

        v1.15:
          • Rich summary: uses chain_iptm, chain_pair_iptm,
            chain_pair_pae_min, chain_ptm, fraction_disordered,
            has_clash when present.
          • Per-residue pLDDT: if only 'atom_plddts' is available
            (current AF3 server output), aggregates atoms → residues
            using the atom→residue map parsed from the model .cif.
          • token_res_ids are kept so the plot can show actual residue
            numbers (1..N per chain, not a running 0..total index).
          • ranking_scores.csv is attached to the result dict when
            present, for display as a per-sample breakdown.
        """
        # ── Summary scores ─────────────────────────────────────
        with open(sum_path, encoding='utf-8') as f:
            summary = json.load(f)

        iptm          = summary.get('iptm')
        ptm           = summary.get('ptm')
        ranking_score = summary.get('ranking_score')

        # Rich per-chain + pair data (AF3 server format)
        chain_iptm         = summary.get('chain_iptm')          # list[float]
        chain_ptm          = summary.get('chain_ptm')           # list[float]
        chain_pair_iptm    = summary.get('chain_pair_iptm')     # 2D list
        chain_pair_pae_min = summary.get('chain_pair_pae_min')  # 2D list
        fraction_disordered = summary.get('fraction_disordered')
        has_clash          = summary.get('has_clash')

        # Mean pLDDT: AF3 summary does not carry it; we compute it
        # below from atom_plddts / per-residue arrays.
        mean_plddt = (summary.get('mean_plddt')
                      or summary.get('mean_pLDDT'))

        # ── PAE + pLDDT + chain info ────────────────────────────
        pae_matrix    = None
        plddt_arr     = None        # per-residue
        chain_ids     = None        # per-residue chain labels
        token_res_ids = None
        atom_plddts   = None
        contact_probs = None

        if conf_path and Path(conf_path).is_file():
            with open(conf_path, encoding='utf-8') as f:
                conf = json.load(f)
            pae_matrix    = conf.get('pae')
            contact_probs = conf.get('contact_probs')

            # Per-residue (token) chain IDs and residue IDs
            for key in ('token_chain_ids', 'chain_ids', 'asym_id'):
                v = conf.get(key)
                if v:
                    chain_ids = v
                    break
            token_res_ids = conf.get('token_res_ids')

            # Per-residue pLDDT — try direct keys first.
            for key in ('token_plddts', 'plddt',
                        'predicted_lddt', 'residue_plddt'):
                v = conf.get(key)
                if v:
                    plddt_arr = v
                    break

            # Per-atom pLDDT (always keep for fallback / CIF aggregation)
            atom_plddts    = conf.get('atom_plddts')

            # If direct per-residue pLDDT is missing or length mismatch,
            # aggregate per-atom pLDDT to per-residue using the CIF map.
            if (plddt_arr is None
                    or (chain_ids is not None
                        and len(plddt_arr) != len(chain_ids))):
                plddt_arr = self._af3a_derive_token_plddt(
                    atom_plddts, chain_ids, token_res_ids, model_cif)

        # If mean_plddt missing, compute from whatever is available.
        if mean_plddt is None:
            for arr in (plddt_arr, atom_plddts):
                if arr:
                    try:
                        mean_plddt = float(np.mean(arr))
                        break
                    except Exception:
                        pass

        # ── Chain layout ────────────────────────────────────────
        chain_order = []
        chain_lens  = {}
        if chain_ids:
            for cid in chain_ids:
                if cid not in chain_lens:
                    chain_lens[cid] = 0
                    chain_order.append(cid)
                chain_lens[cid] += 1
        elif pae_matrix:
            n = len(pae_matrix)
            chain_lens  = {'A': n}
            chain_order = ['A']

        n_chains = len(chain_order)

        # ── ORF names from directory name ────────────────────────
        orf_names = [f"ORF{m.group(1)}"
                     for m in re.finditer(r'orf(\d+)',
                                           job_dir.name, re.IGNORECASE)]

        # ── Sequence verification (v2.0) ─────────────────────────────────
        # Cross-check the sequences inside the input JSON (if present) with
        # the currently loaded ORFs in self.orfs.  Three possible outcomes:
        #
        #   seq_status = 'verified'   — folder-name ORF IDs confirmed by seq
        #   seq_status = 'corrected'  — folder name was WRONG; correct IDs
        #                               found by sequence search and applied
        #   seq_status = 'mismatch'   — sequences found in genome but indices
        #                               differ from folder name (e.g. after
        #                               re-running Pyrodigal)
        #   seq_status = 'no_input'   — input JSON absent; name-only fallback
        #   seq_status = 'no_orfs'    — no genome loaded; cannot verify
        #
        # The seq_chains list is parallel to chain_order and carries the
        # verified protein sequence for each chain (stripped of stop-codon *).
        seq_status   = 'no_input'
        seq_chains   = []        # verified sequences, one per chain
        seq_verified_names = []  # ORF names resolved by sequence (may differ)

        input_seqs = []          # sequences extracted from input JSON
        if input_json and Path(input_json).is_file():
            try:
                with open(input_json, encoding='utf-8') as _f:
                    _inp = json.load(_f)
                for _entry in _inp.get('sequences', []):
                    _pc = _entry.get('proteinChain', {})
                    _seq = _pc.get('sequence', '').strip().rstrip('*').upper()
                    if _seq:
                        input_seqs.append(_seq)
            except (OSError, json.JSONDecodeError, KeyError, TypeError):
                pass

        if input_seqs and hasattr(self, 'orfs') and self.orfs:
            seq_status = 'unmatched'
            seq_chains = input_seqs

            # Build a lookup: sequence → ORF index (exact match, case-insensitive)
            # For speed, hash by length first
            _orfs_by_len: dict = {}
            for _i, _o in enumerate(self.orfs):
                _p = _o['protein'].strip().rstrip('*').upper()
                _orfs_by_len.setdefault(len(_p), []).append((_i, _p))

            def _find_orf_by_seq(seq: str):
                """Return ORF index matching seq exactly, or -1."""
                for _idx, _p in _orfs_by_len.get(len(seq), []):
                    if _p == seq:
                        return _idx
                return -1

            verified_indices = [_find_orf_by_seq(s) for s in input_seqs]

            all_found  = all(i >= 0 for i in verified_indices)
            name_match = ([f"ORF{i+1}" for i in verified_indices]
                          == orf_names[:len(verified_indices)])

            if all_found:
                seq_verified_names = [f"ORF{i+1}" for i in verified_indices]
                if name_match:
                    seq_status = 'verified'
                else:
                    # Sequences found but folder-name ORF IDs are wrong
                    # → use verified names (handles re-numbering after re-run)
                    seq_status = 'corrected'
                    orf_names  = seq_verified_names  # ← authoritative update
            else:
                # Some sequences not found → genome may be different
                seq_status = 'mismatch'
        elif not (hasattr(self, 'orfs') and self.orfs):
            seq_status = 'no_orfs'

        # Map chain letters → ORF names (A→orf_names[0], B→orf_names[1] …)
        chain_to_orf = {}
        for ci, cid in enumerate(chain_order):
            chain_to_orf[cid] = (orf_names[ci]
                                 if ci < len(orf_names) else cid)

        # ── Inter-chain metrics ─────────────────────────────────
        thresh = self._af3a_thresh_spin.value()
        pair_metrics   = {}   # (cid_A, cid_B) → dict
        best_pae_inter = None
        best_pair      = ('?', '?')

        if pae_matrix and n_chains >= 2:
            try:
                pae_np = np.array(pae_matrix, dtype=float)
                best_pae_min = None   # for best-pair selection by PAE_min
                for i_c, ca in enumerate(chain_order):
                    for j_c, cb in enumerate(chain_order):
                        if i_c >= j_c:
                            continue
                        r0 = sum(chain_lens[chain_order[k]]
                                 for k in range(i_c))
                        r1 = r0 + chain_lens[ca]
                        c0 = sum(chain_lens[chain_order[k]]
                                 for k in range(j_c))
                        c1 = c0 + chain_lens[cb]
                        sub_AB = pae_np[r0:r1, c0:c1]
                        sub_BA = pae_np[c0:c1, r0:r1]
                        # Global mean PAE (diluted by disordered regions)
                        pi = float((sub_AB.mean() + sub_BA.mean()) / 2)
                        # Minimum PAE — best single contact point (focal metric)
                        pae_min_ab = float(sub_AB.min())
                        pae_min_ba = float(sub_BA.min())
                        pae_min_pair = min(pae_min_ab, pae_min_ba)
                        # Contact fraction: residue pairs with PAE < 5Å
                        contact_frac = float((sub_AB < 5.0).mean())
                        # Per-row MIN (not mean) to find focal residues correctly
                        minA = sub_AB.min(axis=1)   # best contact each A-res has w/ any B
                        minB = sub_BA.min(axis=1)   # best contact each B-res has w/ any A
                        mA   = sub_AB.mean(axis=1)  # kept for legacy n_contacts count
                        mB   = sub_BA.mean(axis=1)
                        contacts = int((mA < thresh).sum()
                                       + (mB < thresh).sum())
                        cr = self._af3a_contact_str(
                            minA, minB,
                            chain_to_orf.get(ca, ca),
                            chain_to_orf.get(cb, cb),
                            thresh)
                        entry = {
                            'pae_inter':      pi,
                            'pae_min':        pae_min_pair,
                            'contact_frac':   contact_frac,
                            'n_contacts':     contacts,
                            'contact_region': cr,
                        }
                        # Attach rich summary pair-data if available
                        try:
                            if chain_pair_iptm:
                                entry['pair_iptm'] = float(
                                    chain_pair_iptm[i_c][j_c])
                            if chain_pair_pae_min:
                                # Prefer summary chain_pair_pae_min (more accurate)
                                entry['pair_pae_min'] = float(
                                    chain_pair_pae_min[i_c][j_c])
                                pae_min_pair = entry['pair_pae_min']
                        except (IndexError, TypeError, ValueError):
                            pass
                        pair_metrics[(ca, cb)] = entry
                        # Select best pair by PAE_min (focal) not global mean
                        if best_pae_min is None or pae_min_pair < best_pae_min:
                            best_pae_min   = pae_min_pair
                            best_pae_inter = pi
                            best_pair      = (ca, cb)
            except Exception as e:
                print(f"[AF3 Analysis] inter-metric: {e}")

        best_cr = (pair_metrics[best_pair]['contact_region']
                   if best_pair in pair_metrics else '-')

        # ── ranking_scores.csv (per-sample) ─────────────────────
        ranking_samples = []
        if ranking_csv and Path(ranking_csv).is_file():
            try:
                with open(ranking_csv, encoding='utf-8') as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        try:
                            ranking_samples.append({
                                'seed':   int(row.get('seed', 0)),
                                'sample': int(row.get('sample', 0)),
                                'ranking_score': float(
                                    row.get('ranking_score', 'nan')),
                            })
                        except (TypeError, ValueError):
                            continue
            except OSError:
                pass

        # ── Motif detection (v1.16) ─────────────────────────────
        # Build a preliminary result dict so _af3a_detect_motifs can
        # use the standard accessors. Detection runs once at parse time
        # using the current toolbar settings; re-runs happen on demand
        # when the user changes a spinbox.
        prelim = {
            'pae_matrix':    pae_matrix,
            'contact_probs': contact_probs,
            'plddt_arr':     plddt_arr,
            'token_res_ids': token_res_ids,
            'chain_order':   chain_order,
            'chain_lens':    chain_lens,
            'chain_to_orf':  chain_to_orf,
        }
        motifs = []
        try:
            # Pull live parameters from the toolbar when already built,
            # otherwise fall back to class defaults.
            pae_core = (self._af3a_motif_core_spin.value()
                        if hasattr(self, '_af3a_motif_core_spin')
                        else self._AF3_MOTIF_PAE_CORE)
            pae_ext  = (self._af3a_motif_ext_spin.value()
                        if hasattr(self, '_af3a_motif_ext_spin')
                        else self._AF3_MOTIF_PAE_EXT)
            min_ct   = (self._af3a_motif_contact_spin.value()
                        if hasattr(self, '_af3a_motif_contact_spin')
                        else self._AF3_MOTIF_MIN_CONTACT)
            min_sz   = (self._af3a_motif_size_spin.value()
                        if hasattr(self, '_af3a_motif_size_spin')
                        else self._AF3_MOTIF_MIN_SIZE)
            require_recip = (self._af3a_motif_recip_cb.isChecked()
                             if hasattr(self, '_af3a_motif_recip_cb')
                             else True)
            motifs = self._af3a_detect_motifs(
                prelim, pae_core, pae_ext, min_ct, min_sz,
                require_reciprocal=require_recip)
        except Exception as e:
            print(f"[AF3 Analysis] motif detection failed for "
                  f"{job_dir.name}: {e}")

        return {
            'job_name':      job_dir.name,
            'job_dir':       str(job_dir),
            'orf_names':     orf_names,
            'chain_order':   chain_order,
            'chain_lens':    chain_lens,
            'chain_to_orf':  chain_to_orf,
            'n_chains':      n_chains,
            # Global metrics
            'iptm':          iptm,
            'ptm':           ptm,
            'mean_plddt':    mean_plddt,
            'ranking_score': ranking_score,
            'fraction_disordered': fraction_disordered,
            'has_clash':     has_clash,
            # Per-chain metrics (v1.15)
            'chain_iptm':    chain_iptm,
            'chain_ptm':     chain_ptm,
            # Arrays for plotting
            'pae_matrix':    pae_matrix,
            'contact_probs': contact_probs,
            'plddt_arr':     plddt_arr,
            'token_res_ids': token_res_ids,
            # Pairwise interface metrics
            'pair_metrics':  pair_metrics,
            'pae_inter':     best_pae_inter,
            'best_pair':     best_pair,
            'contact_region': best_cr,
            # Focal metrics (v2.0) — correctly identify domain-limited interactions
            'pae_min_inter':  (pair_metrics[best_pair].get('pair_pae_min') or
                               pair_metrics[best_pair].get('pae_min'))
                              if best_pair in pair_metrics else None,
            'cp_iptm_inter':  pair_metrics[best_pair].get('pair_iptm')
                              if best_pair in pair_metrics else None,
            'contact_frac':   pair_metrics[best_pair].get('contact_frac')
                              if best_pair in pair_metrics else None,
            # Motifs (v1.16)
            'motifs':        motifs,
            # Sequence verification (v2.0)
            'seq_status':    seq_status,
            'seq_chains':    seq_chains,
            # Extras
            'ranking_samples': ranking_samples,
            'partner_name':  (orf_names[1] if len(orf_names) > 1
                              else chain_order[1] if len(chain_order) > 1
                              else '-'),
        }

    def _af3a_derive_token_plddt(self, atom_plddts, token_chain_ids,
                                  token_res_ids, model_cif):
        """Aggregate per-atom pLDDT → per-residue (token) pLDDT.

        In current AF3 server outputs only 'atom_plddts' is stored. To
        obtain per-residue values we need the atom→residue mapping,
        which is available in the model .cif file (the B_iso column of
        each ATOM record equals the atom pLDDT).

        Returns a list of length len(token_chain_ids) (mean of atoms in
        that residue), or None if a reliable mapping cannot be built.
        """
        if (not atom_plddts or not token_chain_ids
                or not token_res_ids or model_cif is None):
            return None
        if not Path(model_cif).is_file():
            return None

        # Parse atom_site records to get (chain, seq_id) per atom.
        atom_key_order = []   # list of (chain, seq_id) in atom order
        try:
            in_atom_loop = False
            col_names = []
            with open(model_cif, encoding='utf-8', errors='replace') as f:
                for line in f:
                    s = line.strip()
                    if s.startswith('loop_'):
                        col_names = []
                        in_atom_loop = False
                        continue
                    if s.startswith('_atom_site.'):
                        col_names.append(s)
                        in_atom_loop = False
                        continue
                    if col_names and col_names[0] == '_atom_site.group_PDB':
                        if (line.startswith('ATOM')
                                or line.startswith('HETATM')):
                            in_atom_loop = True
                            parts = line.split()
                            # Fixed CIF layout from AF3 (verified):
                            #   parts[6] = label_asym_id (chain)
                            #   parts[8] = label_seq_id  (residue)
                            if len(parts) >= 9:
                                try:
                                    atom_key_order.append(
                                        (parts[6], int(parts[8])))
                                except ValueError:
                                    atom_key_order.append(None)
                        elif in_atom_loop and s and not s.startswith('#'):
                            # End of the atom loop
                            break
        except OSError:
            return None

        if len(atom_key_order) != len(atom_plddts):
            # Mapping couldn't be built reliably → abort
            return None

        # Group atom pLDDTs by (chain, seq_id)
        from collections import defaultdict
        per_res = defaultdict(list)
        for key, plddt in zip(atom_key_order, atom_plddts):
            if key is not None:
                per_res[key].append(plddt)

        # Build per-token array, ordered by token_chain_ids + token_res_ids
        out = []
        for c, r in zip(token_chain_ids, token_res_ids):
            vs = per_res.get((c, int(r)))
            if vs:
                out.append(sum(vs) / len(vs))
            else:
                out.append(None)
        # If too many residues lack mapping, give up (protects the plot)
        if sum(v is None for v in out) > 0.1 * len(out):
            return None
        # Replace any remaining None with 0 (rare)
        return [v if v is not None else 0.0 for v in out]

    # ────────────────────────────────────────────────────────────────
    # Motif detection (v1.16)
    # ────────────────────────────────────────────────────────────────
    #
    # Detect 2-D interaction motifs in the off-diagonal quadrants of the
    # PAE matrix using two confidence tiers + contact_probs validation:
    #
    #   core      : PAE < pae_core    (≈ 5 Å, high-confidence interface
    #               residues — the "nucleus" of the motif)
    #   extended  : PAE < pae_ext     (≈ 12 Å, exploratory peripheral
    #               residues adjacent to a core)
    #   contacts  : contact_probs >= min_contact (3-D proximity, from
    #               the *_confidences.json; filters out "aligned but far"
    #               artefacts which low-PAE alone cannot distinguish)
    #
    # The algorithm:
    #   1. Extract the A↔B off-diagonal quadrant (+ B↔A for reciprocity).
    #   2. Build core_mask = (PAE < pae_core) AND (contact_probs >= min_contact)
    #                     — contact filter is applied only if the matrix
    #                       is available, otherwise it defaults to True.
    #   3. Close small gaps with a 3×3 morphological closing.
    #   4. Label connected components. Each ≥ (min_size × min_size)
    #      component is a candidate motif.
    #   5. For each motif, extend its bounding box with neighboring
    #      pixels satisfying the extended tier.
    #   6. Reciprocity: require a B↔A component whose transposed bbox
    #      overlaps ≥ 50 % of the A↔B bbox.
    #   7. Score each motif on a 0-100 scale combining PAE, size,
    #      density, pLDDT and reciprocity.
    #
    # The result is stored on each result dict as res['motifs'], sorted
    # by score descending.

    # Default motif-detection parameters. User overrides these via
    # spinboxes in the AlphaFold Analysis toolbar.
    #
    # Defaults were calibrated against real AF3 server output — a strong
    # interface (e.g. confirmed Y2H hits) typically has core-tier PAE
    # < 5 Å over ≥ 5×5 residues; a borderline interface may only have
    # core-tier PAE < 8 Å. The contact-prob filter is intentionally soft
    # (0.05) because contact_probs are often < 0.2 even at real interfaces.
    _AF3_MOTIF_PAE_CORE    = 8.0    # Å — high-confidence tier
    _AF3_MOTIF_PAE_EXT     = 15.0   # Å — extended / peripheral tier
    _AF3_MOTIF_MIN_CONTACT = 0.05   # probability, 0 disables the filter
    _AF3_MOTIF_MIN_SIZE    = 5      # residues — min dimension of a motif

    @staticmethod
    def _af3a_label_cc(mask):
        """Label connected components in a 2-D boolean mask.

        Uses scipy.ndimage.label when available (fast C implementation),
        otherwise falls back to a pure-numpy BFS. 4-connectivity.

        Returns (labels, n) where labels is an int array of the same
        shape as mask and n is the number of components. Background
        pixels carry label 0; motifs are labelled 1..n.
        """
        if SCIPY_NDIMAGE_AVAILABLE:
            # Default structuring element is 4-connectivity
            labels, n = _scipy_ndimage.label(mask)
            return labels, n
        # Pure-numpy BFS fallback
        h, w = mask.shape
        labels = np.zeros((h, w), dtype=np.int32)
        cur = 0
        # Pre-compute neighbour offsets
        for i in range(h):
            for j in range(w):
                if not mask[i, j] or labels[i, j]:
                    continue
                cur += 1
                stack = [(i, j)]
                labels[i, j] = cur
                while stack:
                    y, x = stack.pop()
                    for dy, dx in ((-1, 0), (1, 0), (0, -1), (0, 1)):
                        ny, nx = y + dy, x + dx
                        if (0 <= ny < h and 0 <= nx < w
                                and mask[ny, nx] and not labels[ny, nx]):
                            labels[ny, nx] = cur
                            stack.append((ny, nx))
        return labels, cur

    @staticmethod
    def _af3a_binary_closing(mask, iterations=2):
        """Morphological closing (dilation → erosion) with 3×3 kernel.

        Bridges small gaps so a motif fragmented by one or two noisy
        high-PAE pixels still registers as a single component.
        Falls back to a simple numpy implementation if scipy missing.
        """
        if SCIPY_NDIMAGE_AVAILABLE:
            return _scipy_ndimage.binary_closing(mask, iterations=iterations)
        out = mask.copy()
        for _ in range(iterations):
            # dilate
            d = out.copy()
            d[1:, :]  |= out[:-1, :]
            d[:-1, :] |= out[1:, :]
            d[:, 1:]  |= out[:, :-1]
            d[:, :-1] |= out[:, 1:]
            # erode (4-connected erosion = AND of 4 shifted copies + self)
            e = d.copy()
            e[1:, :]  &= d[:-1, :]
            e[:-1, :] &= d[1:, :]
            e[:, 1:]  &= d[:, :-1]
            e[:, :-1] &= d[:, 1:]
            out = e
        return out

    @classmethod
    def _af3a_detect_motifs_one_pair(cls, pae_AB, pae_BA, contact_AB,
                                       plddt_A, plddt_B,
                                       pae_core, pae_ext, min_contact,
                                       min_size, require_reciprocal,
                                       offset_a=0, offset_b=0,
                                       token_res_A=None,
                                       token_res_B=None,
                                       res_id_to_res_num=None):
        """Detect motifs in a single A↔B quadrant pair.

        Parameters
        ----------
        pae_AB : 2-D ndarray (len_A × len_B)
            Off-diagonal PAE, rows = chain A residues, cols = chain B.
        pae_BA : 2-D ndarray (len_B × len_A)
            Transposed quadrant, for reciprocity validation.
        contact_AB : 2-D ndarray or None (len_A × len_B)
            Contact probabilities. None disables the contact filter.
        plddt_A, plddt_B : 1-D ndarray or None
            Per-residue pLDDT for chains A, B.
        pae_core, pae_ext : float
            PAE thresholds in Å for the core (strict) and extended
            (peripheral) tiers.
        min_contact : float
            Minimum contact probability for the core mask. 0 disables.
        min_size : int
            Minimum width AND height (in residues) for a surviving motif.
        require_reciprocal : bool
            Drop motifs without a matching component in the B↔A quadrant.
        offset_a, offset_b : int
            Running positional offset of this quadrant inside the full
            PAE matrix — used so plot overlays can draw rectangles in
            the correct absolute position.
        token_res_A, token_res_B : 1-D ndarray or None
            The AF3 per-residue numbers for chain A / B (from
            token_res_ids). When provided, motif positions are reported
            in the native chain-local numbering (restarting at 1 per
            chain) that the AF3 server shows.

        Returns
        -------
        list[dict]  sorted by `score` desc.
        """
        if pae_AB is None or pae_AB.size == 0:
            return []

        len_A, len_B = pae_AB.shape
        use_contact = (contact_AB is not None
                       and contact_AB.shape == pae_AB.shape
                       and min_contact > 0.0)

        # 1) Core mask — low PAE AND (optionally) high contact prob
        core_mask = (pae_AB < pae_core)
        if use_contact:
            core_mask &= (contact_AB >= min_contact)

        # 2) Close 1-2 pixel gaps
        core_mask = cls._af3a_binary_closing(core_mask, iterations=1)

        # 3) Reciprocity mask on the transposed quadrant
        if pae_BA is not None and pae_BA.shape == (len_B, len_A):
            recip_mask = (pae_BA < pae_core)
            recip_mask = cls._af3a_binary_closing(recip_mask, iterations=1)
        else:
            recip_mask = None

        # 4) Label connected components
        labels, n_cc = cls._af3a_label_cc(core_mask)
        if n_cc == 0:
            return []

        # 5) Extended mask (peripheral residues, for bbox extension)
        ext_mask = (pae_AB < pae_ext)

        motifs = []
        for k in range(1, n_cc + 1):
            comp_rows, comp_cols = np.where(labels == k)
            if len(comp_rows) == 0:
                continue
            r0, r1 = int(comp_rows.min()), int(comp_rows.max())
            c0, c1 = int(comp_cols.min()), int(comp_cols.max())
            core_h = r1 - r0 + 1
            core_w = c1 - c0 + 1
            if core_h < min_size or core_w < min_size:
                continue

            # 5a) Extend bbox with contiguous low-PAE peripheral cells
            #     (grow the bbox while any edge row/col has ≥ 1 ext pixel)
            while r0 > 0 and ext_mask[r0 - 1, c0:c1 + 1].any():
                r0 -= 1
            while r1 < len_A - 1 and ext_mask[r1 + 1, c0:c1 + 1].any():
                r1 += 1
            while c0 > 0 and ext_mask[r0:r1 + 1, c0 - 1].any():
                c0 -= 1
            while c1 < len_B - 1 and ext_mask[r0:r1 + 1, c1 + 1].any():
                c1 += 1

            block_pae = pae_AB[r0:r1 + 1, c0:c1 + 1]
            block_core_mask = (block_pae < pae_core)

            mean_pae = float(block_pae.mean())
            min_pae  = float(block_pae.min())
            density  = float(block_core_mask.mean())   # % core pixels
            h = r1 - r0 + 1
            w = c1 - c0 + 1
            area = h * w

            # Mean contact prob in the motif (only if available)
            if use_contact:
                mean_contact = float(contact_AB[r0:r1+1, c0:c1+1].mean())
                max_contact  = float(contact_AB[r0:r1+1, c0:c1+1].max())
            else:
                mean_contact = None
                max_contact  = None

            # Per-chain pLDDT for residues inside the motif
            plddt_mA = None
            plddt_mB = None
            if plddt_A is not None and len(plddt_A) >= r1 + 1:
                plddt_mA = float(np.mean(plddt_A[r0:r1 + 1]))
            if plddt_B is not None and len(plddt_B) >= c1 + 1:
                plddt_mB = float(np.mean(plddt_B[c0:c1 + 1]))
            plddt_mean = (
                None if plddt_mA is None and plddt_mB is None
                else (
                    (plddt_mA or 0.0) + (plddt_mB or 0.0)
                ) / max(1, (plddt_mA is not None) + (plddt_mB is not None))
            )

            # Reciprocity: does the transposed region in B↔A also have
            # a low-PAE block? Require ≥ 50 % overlap.
            reciprocal = False
            recip_overlap = 0.0
            if recip_mask is not None:
                recip_block = recip_mask[c0:c1 + 1, r0:r1 + 1]
                recip_overlap = float(recip_block.mean())
                reciprocal = recip_overlap >= 0.5
            if require_reciprocal and recip_mask is not None and not reciprocal:
                continue

            # Combined score (0-100). Weights chosen to match:
            #   confidence 35 %, size 20 %, density 20 %,
            #   pLDDT 15 %, reciprocity 10 %
            conf_term = max(0.0, 1.0 - mean_pae / max(pae_ext, 0.1))
            size_term = min(min(h, w) / 50.0, 1.0)   # saturates at 50
            dens_term = density
            if plddt_mean is not None:
                plddt_term = max(0.0, (plddt_mean - 50.0) / 50.0)
            else:
                plddt_term = 0.5   # neutral if unknown
            recip_term = 1.0 if reciprocal else 0.0
            score = 100.0 * (
                0.35 * conf_term
                + 0.20 * size_term
                + 0.20 * dens_term
                + 0.15 * plddt_term
                + 0.10 * recip_term)

            # Convert row/col indices to user-facing residue numbers
            def _resnum(arr, idx):
                if arr is not None and 0 <= idx < len(arr):
                    try:
                        return int(arr[idx])
                    except (TypeError, ValueError):
                        pass
                return idx + 1  # 1-based running

            a_start = _resnum(token_res_A, r0)
            a_end   = _resnum(token_res_A, r1)
            b_start = _resnum(token_res_B, c0)
            b_end   = _resnum(token_res_B, c1)

            motifs.append({
                # Bounding box in LOCAL quadrant coordinates
                'a_row0': r0, 'a_row1': r1,
                'b_col0': c0, 'b_col1': c1,
                # Residue numbering (native, per-chain)
                'a_start': a_start, 'a_end': a_end,
                'b_start': b_start, 'b_end': b_end,
                # Global (absolute) indices — for plot overlays
                'abs_row0': r0 + offset_a, 'abs_row1': r1 + offset_a,
                'abs_col0': c0 + offset_b, 'abs_col1': c1 + offset_b,
                # Metrics
                'width':  w,   'height': h,   'area': area,
                'mean_pae': mean_pae, 'min_pae': min_pae,
                'density': density,
                'mean_contact': mean_contact, 'max_contact': max_contact,
                'plddt_A': plddt_mA, 'plddt_B': plddt_mB,
                'plddt_mean': plddt_mean,
                'reciprocal': reciprocal,
                'recip_overlap': recip_overlap,
                'score': score,
            })

        motifs.sort(key=lambda m: -m['score'])
        return motifs

    @classmethod
    def _af3a_detect_motifs(cls, res: dict,
                             pae_core: float, pae_ext: float,
                             min_contact: float, min_size: int,
                             require_reciprocal: bool = True) -> list:
        """Detect motifs for every off-diagonal pair in one AF3 result.

        Returns a flat list of motifs with a 'pair' key identifying
        which chain pair the motif belongs to. Sorted by score.
        """
        pae_matrix    = res.get('pae_matrix')
        contact_probs = res.get('contact_probs')
        plddt_arr     = res.get('plddt_arr')
        token_res_ids = res.get('token_res_ids')
        chain_order   = res.get('chain_order') or []
        chain_lens    = res.get('chain_lens') or {}
        chain_to_orf  = res.get('chain_to_orf') or {}

        if (not pae_matrix or len(chain_order) < 2
                or not NUMPY_AVAILABLE):
            return []

        pae_np = np.array(pae_matrix, dtype=float)
        contact_np = (np.array(contact_probs, dtype=float)
                      if contact_probs else None)
        plddt_np = (np.array(plddt_arr, dtype=float)
                    if plddt_arr else None)
        tri_np = (np.array(token_res_ids)
                  if token_res_ids else None)

        # Running offsets for each chain (for mapping into full matrix)
        chain_offsets = {}
        acc = 0
        for cid in chain_order:
            chain_offsets[cid] = acc
            acc += chain_lens.get(cid, 0)

        all_motifs = []
        for i_c, ca in enumerate(chain_order):
            for j_c, cb in enumerate(chain_order):
                if i_c >= j_c:
                    continue  # only upper-triangle inter-chain
                r0 = chain_offsets[ca]; r1 = r0 + chain_lens[ca]
                c0 = chain_offsets[cb]; c1 = c0 + chain_lens[cb]

                pae_AB = pae_np[r0:r1, c0:c1]
                pae_BA = pae_np[c0:c1, r0:r1]
                cont_AB = (contact_np[r0:r1, c0:c1]
                           if contact_np is not None else None)
                plddt_A = (plddt_np[r0:r1]
                           if plddt_np is not None else None)
                plddt_B = (plddt_np[c0:c1]
                           if plddt_np is not None else None)
                tri_A   = (tri_np[r0:r1]
                           if tri_np is not None else None)
                tri_B   = (tri_np[c0:c1]
                           if tri_np is not None else None)

                motifs = cls._af3a_detect_motifs_one_pair(
                    pae_AB, pae_BA, cont_AB,
                    plddt_A, plddt_B,
                    pae_core=pae_core, pae_ext=pae_ext,
                    min_contact=min_contact,
                    min_size=min_size,
                    require_reciprocal=require_reciprocal,
                    offset_a=r0, offset_b=c0,
                    token_res_A=tri_A, token_res_B=tri_B)

                # Attach pair identity & friendly names
                name_A = chain_to_orf.get(ca, ca)
                name_B = chain_to_orf.get(cb, cb)
                for m in motifs:
                    m['pair']    = (ca, cb)
                    m['chain_A'] = ca
                    m['chain_B'] = cb
                    m['name_A']  = name_A
                    m['name_B']  = name_B
                    m['label']   = (f"{name_A}[{m['a_start']}-{m['a_end']}] "
                                    f"× {name_B}[{m['b_start']}-{m['b_end']}]")
                all_motifs.extend(motifs)

        all_motifs.sort(key=lambda m: -m['score'])
        # Assign stable rank indices (for plot overlay labels)
        for rank, m in enumerate(all_motifs, start=1):
            m['rank'] = rank
        return all_motifs

    @staticmethod
    def _af3a_contact_str(mean_A, mean_B, name_A, name_B, thresh):
        """Build a human-readable contact region string."""
        def _ranges(indices, prefix):
            if not len(indices):
                return ''
            ranges, s, e = [], int(indices[0]), int(indices[0])
            for idx in indices[1:]:
                idx = int(idx)
                if idx == e + 1:
                    e = idx
                else:
                    ranges.append(f"res{s+1}-{e+1}" if e > s else f"res{s+1}")
                    s = e = idx
            ranges.append(f"res{s+1}-{e+1}" if e > s else f"res{s+1}")
            return f"{prefix}: {', '.join(ranges[:3])}"

        # Use the supplied per-residue arrays directly (caller now passes min,
        # not mean, so low values correctly flag focal-contact residues)
        c5_A  = np.where(mean_A < 5.0)[0]
        c10_A = np.where(mean_A < 10.0)[0]
        c5_B  = np.where(mean_B < 5.0)[0]
        c10_B = np.where(mean_B < 10.0)[0]

        if not len(c5_A) and not len(c10_A):
            return f"no contact detected with {name_B}"

        parts = []
        if len(c5_A):
            parts.append(f"{_ranges(c5_A, name_A)} <5Å")
        elif len(c10_A):
            parts.append(f"{_ranges(c10_A, name_A)} <10Å")
        if len(c5_B):
            parts.append(f"{_ranges(c5_B, name_B)} <5Å")
        elif len(c10_B):
            parts.append(f"{_ranges(c10_B, name_B)} <10Å")
        return ' | '.join(parts) if parts else (
            f"not relevant interaction found with {name_B}")

    def _af3a_populate_table(self):
        """Populate the results table.

        BUG FIX (v1.15):
          • Disable sorting during insertion so Qt does not rearrange
            rows mid-population.
          • Store the result-list index in UserRole data on col 0 of
            each row, so `_af3a_on_select` and `_af3a_replot_selected`
            can retrieve the correct result even after the user
            re-sorts the table by any column.
          • Use a small numeric-aware QTableWidgetItem subclass for
            the float columns (ipTM, ptm, pLDDT, ranking, PAE_inter)
            so sorting is numeric, not lexicographic.
        """
        table = self._af3a_table
        user_role = self._af3a_user_role()

        # NumericItem: sorts by the stored float, displays formatted text
        class _NumericItem(QTableWidgetItem):
            def __init__(self, text: str, value):
                super().__init__(text)
                # value may be None → treat as -inf for best-first sort
                self._val = float(value) if value is not None else float('-inf')

            def __lt__(self, other):
                try:
                    return self._val < other._val
                except AttributeError:
                    return super().__lt__(other)

        # Disable sort + signals during rebuild
        was_sortable = table.isSortingEnabled()
        table.setSortingEnabled(False)
        table.blockSignals(True)
        table.setRowCount(0)

        def _c(v, thr_good, thr_ok, inv=False):
            if v is None: return None
            ok  = (v <= thr_ok)   if inv else (v >= thr_good)
            mid = (v <= thr_good) if inv else (v >= thr_ok)
            if ok:   return QColor('#C8E6C9')
            if mid:  return QColor('#FFF9C4')
            return QColor('#FFCDD2')

        for res_idx, res in enumerate(self._af3_analysis_results):
            row = table.rowCount()
            table.insertRow(row)

            chains_s = ', '.join(
                f"{res['chain_to_orf'].get(c, c)}({res['chain_lens'][c]}aa)"
                for c in res['chain_order'])
            iptm_s   = f"{res['iptm']:.3f}"       if res.get('iptm')          is not None else '-'
            ptm_s    = f"{res['ptm']:.3f}"        if res.get('ptm')           is not None else '-'
            plddt_s  = f"{res['mean_plddt']:.1f}" if res.get('mean_plddt')    is not None else '-'
            rank_s   = f"{res['ranking_score']:.4f}" if res.get('ranking_score') is not None else '-'
            pi_s     = f"{res['pae_inter']:.1f}"  if res.get('pae_inter')     is not None else '-'

            # Column 0: job name — carries the data index in UserRole
            # Add a sequence-verification badge prefix to the cell text
            _seq_st   = res.get('seq_status', 'no_input')
            _seq_icon = {
                'verified':  '✅',   # folder name confirmed by sequence
                'corrected': '🔁',   # folder name was wrong; corrected
                'mismatch':  '⚠️',   # sequences not found in current genome
                'no_input':  '📂',   # no input JSON available
                'no_orfs':   '❔',   # no genome loaded
            }.get(_seq_st, '❔')
            item0 = QTableWidgetItem(f"{_seq_icon} {res['job_name']}")
            item0.setData(user_role, res_idx)
            # QC tooltip with rich summary info
            tip_parts = []
            # ── Sequence verification status (first, most important) ──────
            _seq_msgs = {
                'verified':  '✅ Sequences verified — folder name matches genome',
                'corrected': '🔁 ORF IDs corrected by sequence match '
                             '(folder name differed from current genome numbering)',
                'mismatch':  '⚠️ SEQUENCE MISMATCH — proteins in input JSON not found '
                             'in current genome. Results may belong to a different genome.',
                'no_input':  '📂 No input JSON found — ORF IDs from folder name only '
                             '(not sequence-verified)',
                'no_orfs':   '❔ No genome loaded — cannot verify sequences',
            }
            tip_parts.append(_seq_msgs.get(_seq_st, '❔ Unknown status'))
            if res.get('has_clash') is not None:
                hc = res['has_clash']
                has = (hc if isinstance(hc, bool) else hc >= 0.5)
                tip_parts.append(
                    f"Clash: {'YES ⚠' if has else 'no'}")
            if res.get('fraction_disordered') is not None:
                tip_parts.append(
                    f"Disordered: {res['fraction_disordered']*100:.1f}%")
            # Best-pair PAE min (from chain_pair_pae_min)
            bp = res.get('best_pair')
            if bp and bp in res.get('pair_metrics', {}):
                pm = res['pair_metrics'][bp]
                if pm.get('pair_pae_min') is not None:
                    tip_parts.append(
                        f"PAE_min(best pair): {pm['pair_pae_min']:.2f} Å")
                if pm.get('pair_iptm') is not None:
                    tip_parts.append(
                        f"pair_ipTM: {pm['pair_iptm']:.2f}")
            if res.get('chain_iptm'):
                tip_parts.append(
                    "chain_iptm: ["
                    + ", ".join(f"{x:.2f}" for x in res['chain_iptm'])
                    + "]")
            if res.get('ranking_samples'):
                tip_parts.append(
                    f"{len(res['ranking_samples'])} diffusion sample(s)")
            if tip_parts:
                item0.setToolTip("\n".join(tip_parts))
            # Color the job name cell by verification status
            _seq_bg = {
                'verified':  QColor('#E8F5E9'),   # light green
                'corrected': QColor('#E3F2FD'),   # light blue
                'mismatch':  QColor('#FFEBEE'),   # light red
                'no_input':  None,                # default
                'no_orfs':   None,
            }.get(_seq_st)
            if _seq_bg:
                item0.setBackground(_seq_bg)
            table.setItem(row, 0, item0)

            # Column 1: n_chains (int; small numbers sort OK as str but still numeric)
            table.setItem(row, 1, _NumericItem(str(res['n_chains']),
                                                res['n_chains']))
            # Column 2: chain description (string)
            table.setItem(row, 2, QTableWidgetItem(chains_s))

            # Columns 3–7: numeric metrics
            it_iptm  = _NumericItem(iptm_s,   res.get('iptm'))
            it_ptm   = _NumericItem(ptm_s,    res.get('ptm'))
            it_plddt = _NumericItem(plddt_s,  res.get('mean_plddt'))
            it_rank  = _NumericItem(rank_s,   res.get('ranking_score'))
            it_paei  = _NumericItem(pi_s,     res.get('pae_inter'))

            bg_iptm = _c(res.get('iptm'),      0.75,  0.50)
            bg_paei = _c(res.get('pae_inter'), 8.0,  15.0, inv=True)
            if bg_iptm: it_iptm.setBackground(bg_iptm)
            if bg_paei: it_paei.setBackground(bg_paei)

            table.setItem(row, 3, it_iptm)
            table.setItem(row, 4, it_ptm)
            table.setItem(row, 5, it_plddt)
            table.setItem(row, 6, it_rank)
            table.setItem(row, 7, it_paei)

            # ── Columns 8–10: focal interaction metrics (v2.0) ──────────────
            pae_min_v = res.get('pae_min_inter')
            cp_iptm_v = res.get('cp_iptm_inter')
            cfrac_v   = res.get('contact_frac')

            pmin_s  = f"{pae_min_v:.2f}" if pae_min_v is not None else '-'
            cpip_s  = f"{cp_iptm_v:.2f}" if cp_iptm_v is not None else '-'
            cfrac_s = f"{cfrac_v*100:.1f}%" if cfrac_v is not None else '-'

            it_pmin = _NumericItem(pmin_s, pae_min_v)
            it_cpip = _NumericItem(cpip_s, cp_iptm_v)
            it_cfrc = _NumericItem(cfrac_s, (cfrac_v or 0) * 100)

            # PAE_min coloring: green < 4Å, yellow 4–8Å, red ≥ 8Å
            bg_pmin = _c(pae_min_v, 4.0, 8.0, inv=True)
            if bg_pmin: it_pmin.setBackground(bg_pmin)
            # cp_ipTM coloring same thresholds as global ipTM
            bg_cpip = _c(cp_iptm_v, 0.65, 0.50)
            if bg_cpip: it_cpip.setBackground(bg_cpip)

            table.setItem(row, 8,  it_pmin)
            table.setItem(row, 9,  it_cpip)
            table.setItem(row, 10, it_cfrc)

            # Column 11: focal contact region description
            table.setItem(row, 11, QTableWidgetItem(
                res.get('contact_region', '-')))

            # ── High-confidence highlight (v2.0) ───────────────────────────
            # Primary criterion: PAE_min < 4Å AND cp_ipTM ≥ 0.50 = FOCAL HIT
            # Fallback: legacy ipTM > 0.75 AND PAEinter < 8Å
            iptm_v = res.get('iptm')
            paei_v = res.get('pae_inter')
            focal_hit = (pae_min_v is not None and pae_min_v < 4.0
                         and (cp_iptm_v is None or cp_iptm_v >= 0.50))
            legacy_hit = (iptm_v is not None and paei_v is not None
                          and iptm_v > 0.75 and paei_v < 8.0)
            if focal_hit or legacy_hit:
                bold = QFont()
                bold.setBold(True)
                for c in range(table.columnCount()):
                    it = table.item(row, c)
                    if it is not None:
                        it.setFont(bold)

        table.blockSignals(False)
        if was_sortable:
            table.setSortingEnabled(True)

        # Apply any active filters
        self._af3a_apply_filters()

    def _af3a_on_select(self):
        """Plot the selected job.

        BUG FIX (v1.15): retrieve the data-list index from the
        UserRole of column 0 instead of using the view row index
        (which is scrambled whenever the user sorts the table).
        """
        rows = set(idx.row() for idx in self._af3a_table.selectedIndexes())
        if not rows:
            return
        view_row = min(rows)
        item0 = self._af3a_table.item(view_row, 0)
        if item0 is None:
            return
        res_idx = item0.data(self._af3a_user_role())
        # Defensive fallback: if UserRole is missing (e.g. legacy project
        # files loaded into this build), fall back to view_row.
        if res_idx is None:
            res_idx = view_row
        try:
            res_idx = int(res_idx)
        except (TypeError, ValueError):
            return
        if not (0 <= res_idx < len(self._af3_analysis_results)):
            return

        res = self._af3_analysis_results[res_idx]
        self._af3a_plot_job(res)

        for orf_name in res.get('orf_names', []):
            m = re.match(r'ORF(\d+)', orf_name, re.IGNORECASE)
            if m:
                idx = int(m.group(1)) - 1
                if 0 <= idx < len(self.orfs):
                    self._select_and_center_orf(idx)
                    break

    def _af3a_on_double_click(self, row: int, _col: int):
        """Open the AF3 job folder in the OS file manager."""
        item0 = self._af3a_table.item(row, 0)
        if item0 is None:
            return
        res_idx = item0.data(self._af3a_user_role())
        if res_idx is None:
            return
        try:
            res = self._af3_analysis_results[int(res_idx)]
        except (IndexError, TypeError, ValueError):
            return
        job_dir = res.get('job_dir', '')
        if not job_dir or not Path(job_dir).is_dir():
            QMessageBox.information(
                self, "Open folder",
                f"Folder not found:\n{job_dir}")
            return
        # Cross-platform folder open
        try:
            if sys.platform == 'darwin':
                subprocess.Popen(['open', job_dir])
            elif sys.platform.startswith('win'):
                os.startfile(job_dir)  # noqa: E501 - Windows-only
            else:
                subprocess.Popen(['xdg-open', job_dir])
        except Exception as e:
            QMessageBox.warning(
                self, "Open folder",
                f"Could not open folder:\n{e}")

    def _af3a_replot_selected(self):
        """Re-render plots for the currently selected row, re-using
        the latest contact threshold.

        BUG FIX (v1.15):
          • Uses UserRole-stored index (see _af3a_on_select).
          • Also refreshes the 'Best contact pair' cell in the table
            so the on-screen text matches the new threshold.
        """
        rows = set(idx.row() for idx in self._af3a_table.selectedIndexes())
        if not rows:
            return
        view_row = min(rows)
        item0 = self._af3a_table.item(view_row, 0)
        if item0 is None:
            return
        res_idx = item0.data(self._af3a_user_role())
        if res_idx is None:
            return
        try:
            res_idx = int(res_idx)
        except (TypeError, ValueError):
            return
        if not (0 <= res_idx < len(self._af3_analysis_results)):
            return

        res = self._af3_analysis_results[res_idx]
        try:
            thresh = self._af3a_thresh_spin.value()
            if res.get('pae_matrix') and res.get('n_chains', 1) >= 2:
                pae_np = np.array(res['pae_matrix'], dtype=float)
                chain_order = res['chain_order']
                chain_lens  = res['chain_lens']
                best_pair = res.get('best_pair', ('?', '?'))
                for (ca, cb) in list(res['pair_metrics'].keys()):
                    i_c = chain_order.index(ca)
                    j_c = chain_order.index(cb)
                    r0 = sum(chain_lens[chain_order[k]] for k in range(i_c))
                    r1 = r0 + chain_lens[ca]
                    c0 = sum(chain_lens[chain_order[k]] for k in range(j_c))
                    c1 = c0 + chain_lens[cb]
                    mA = pae_np[r0:r1, c0:c1].mean(axis=1)
                    mB = pae_np[c0:c1, r0:r1].mean(axis=1)
                    cr = self._af3a_contact_str(
                        mA, mB,
                        res['chain_to_orf'].get(ca, ca),
                        res['chain_to_orf'].get(cb, cb), thresh)
                    res['pair_metrics'][(ca, cb)]['contact_region'] = cr
                    res['pair_metrics'][(ca, cb)]['n_contacts'] = int(
                        (mA < thresh).sum() + (mB < thresh).sum())
                # Update 'best' contact_region for the table
                if best_pair in res['pair_metrics']:
                    res['contact_region'] = (
                        res['pair_metrics'][best_pair]['contact_region'])
                    # Refresh the table cell (col 11 = contact_region after v2.0)
                    cell = self._af3a_table.item(view_row, 11)
                    if cell is not None:
                        cell.setText(res['contact_region'])
        except Exception as e:
            print(f"[AF3 Analysis] replot threshold update: {e}")
        self._af3a_plot_job(res)

    def _af3a_clear_plots(self):
        self._af3a_active_canvases = []
        # v1.16: drop motif-highlight references since the canvas they
        # point to is about to be deleted.
        self._af3a_last_pae_ax = None
        self._af3a_last_pae_canvas = None
        self._af3a_highlight_patches = []
        while self._af3a_plot_layout.count():
            item = self._af3a_plot_layout.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()

    def _af3a_plot_job(self, res: dict):
        """Render interactive PAE heatmap(s) + pLDDT for one job."""
        if not MATPLOTLIB_AVAILABLE:
            return

        self._af3a_clear_plots()
        from matplotlib.colors import LinearSegmentedColormap

        # ── ChimeraX PAE colormap ──────────────────────────────
        pae_cmap = LinearSegmentedColormap.from_list('pae_chimerax', [
            '#1E3A8A',  # 0 Å  dark blue
            '#2563EB',  # 5 Å  blue
            '#60A5FA',  # 10 Å light blue
            '#FCD34D',  # 15 Å yellow
            '#F97316',  # 20 Å orange
            '#D1D5DB',  # 25 Å light grey
            '#F9FAFB',  # 31 Å near white
        ], N=256)

        pae   = res.get('pae_matrix')
        plddt = res.get('plddt_arr')
        chain_order = res.get('chain_order', [])
        chain_lens  = res.get('chain_lens', {})
        chain_to_orf = res.get('chain_to_orf', {})
        n_chains    = res.get('n_chains', 1)
        thresh      = self._af3a_thresh_spin.value()

        # Build ORF-labelled axis tick positions and labels
        ticks_pos  = []
        ticks_lbl  = []
        dividers   = []
        cum = 0
        for ci, cid in enumerate(chain_order):
            clen = chain_lens.get(cid, 0)
            ticks_pos.append(cum + clen / 2)
            ticks_lbl.append(chain_to_orf.get(cid, cid))
            if ci > 0:
                dividers.append(cum - 0.5)
            cum += clen

        # ── Title + metrics header ─────────────────────────────
        iptm_s  = f"{res['iptm']:.3f}"  if res.get('iptm') is not None else 'N/A'
        ptm_s   = f"{res['ptm']:.3f}"   if res.get('ptm')  is not None else 'N/A'
        pl_s    = f"{res['mean_plddt']:.1f}" if res.get('mean_plddt') is not None else 'N/A'
        pi_s    = f"{res['pae_inter']:.1f} Å" if res.get('pae_inter') is not None else 'N/A'
        nc_s    = str(res['pair_metrics'].get(
                        res.get('best_pair', ()), {}).get('n_contacts', '-'))

        hdr = QLabel(
            f"<b>{res['job_name']}</b>  —  "
            f"Chains: {n_chains}  |  "
            f"ipTM={iptm_s}  ptm={ptm_s}  "
            f"mean_pLDDT={pl_s}  "
            f"Best PAE_inter={pi_s}  contacts={nc_s} res")
        hdr.setStyleSheet(
            "font-size:11px;padding:4px 8px;background:#E8EAF6;"
            "border-radius:4px;")
        self._af3a_plot_layout.addWidget(hdr)

        # QC line: fraction_disordered + has_clash + per-chain metrics
        qc_parts = []
        fd = res.get('fraction_disordered')
        if fd is not None:
            colour = '#E65100' if fd > 0.3 else '#2E7D32'
            qc_parts.append(
                f"<span style='color:{colour}'>"
                f"disordered={fd*100:.1f}%</span>")
        hc = res.get('has_clash')
        if hc is not None:
            has = (hc if isinstance(hc, bool) else hc >= 0.5)
            if has:
                qc_parts.append(
                    "<span style='color:#C62828'>⚠ clash</span>")
            else:
                qc_parts.append(
                    "<span style='color:#2E7D32'>no clash</span>")
        # Per-chain ipTM / ptm breakdown
        c_iptm = res.get('chain_iptm') or []
        c_ptm  = res.get('chain_ptm')  or []
        for i, cid in enumerate(chain_order):
            bits = []
            if i < len(c_iptm) and c_iptm[i] is not None:
                bits.append(f"ipTM={c_iptm[i]:.2f}")
            if i < len(c_ptm) and c_ptm[i] is not None:
                bits.append(f"pTM={c_ptm[i]:.2f}")
            if bits:
                qc_parts.append(
                    f"{chain_to_orf.get(cid, cid)}({','.join(bits)})")
        if qc_parts:
            qc = QLabel("  " + "   ".join(qc_parts))
            qc.setStyleSheet(
                "font-size:10px;padding:1px 8px;color:#455A64;")
            qc.setTextFormat(
                Qt.TextFormat.RichText if QT_VERSION == 6
                else Qt.RichText)
            self._af3a_plot_layout.addWidget(qc)

        # ── Sequence verification status line ─────────────────────────────
        _seq_st   = res.get('seq_status', 'no_input')
        _seq_html = {
            'verified':  ("<span style='color:#2E7D32;font-weight:500'>"
                          "✅ Sequences verified</span>"
                          " — ORF IDs confirmed by exact sequence match "
                          "with loaded genome"),
            'corrected': ("<span style='color:#1565C0;font-weight:500'>"
                          "🔁 ORF IDs corrected</span>"
                          " — folder name differed from current genome "
                          "numbering; correct IDs found by sequence match"),
            'mismatch':  ("<span style='color:#C62828;font-weight:500'>"
                          "⚠️ SEQUENCE MISMATCH</span>"
                          " — proteins in this job were NOT found in the "
                          "current genome. Results may belong to a different "
                          "genome or a re-analysed sequence. "
                          "Interpretation may be incorrect."),
            'no_input':  ("<span style='color:#888'>"
                          "📂 Not verified</span>"
                          " — no input JSON found; ORF IDs from folder name only"),
            'no_orfs':   ("<span style='color:#888'>"
                          "❔ Cannot verify</span>"
                          " — no genome loaded"),
        }.get(_seq_st, "")
        if _seq_html:
            seq_lbl = QLabel("  " + _seq_html)
            seq_lbl.setTextFormat(
                Qt.TextFormat.RichText if QT_VERSION == 6 else Qt.RichText)
            seq_lbl.setStyleSheet(
                "font-size:10px;padding:1px 8px;"
                "background:#f8f8f8;border-left:3px solid "
                + {
                    'verified':  '#2E7D32',
                    'corrected': '#1565C0',
                    'mismatch':  '#C62828',
                }.get(_seq_st, '#ccc') + ";")
            seq_lbl.setWordWrap(True)
            self._af3a_plot_layout.addWidget(seq_lbl)

        # Contact regions + per-pair rich metrics
        if res.get('pair_metrics'):
            for (ca, cb), pm in res['pair_metrics'].items():
                extras = []
                if pm.get('pair_iptm') is not None:
                    extras.append(f"ipTM={pm['pair_iptm']:.2f}")
                if pm.get('pair_pae_min') is not None:
                    extras.append(f"PAE_min={pm['pair_pae_min']:.1f}Å")
                if pm.get('pae_inter') is not None:
                    extras.append(f"PAE_mean={pm['pae_inter']:.1f}Å")
                extras_s = (" [" + ", ".join(extras) + "]") if extras else ""
                cr_lbl = QLabel(
                    f"  {chain_to_orf.get(ca,ca)} ↔ "
                    f"{chain_to_orf.get(cb,cb)}{extras_s}: "
                    f"{pm.get('contact_region', '-')}")
                cr_lbl.setStyleSheet(
                    "font-size:10px;color:#1565C0;padding:1px 8px;")
                cr_lbl.setWordWrap(True)
                self._af3a_plot_layout.addWidget(cr_lbl)

        # ── PAE HEATMAP ────────────────────────────────────────
        if pae:
            pae_np = np.array(pae, dtype=float)
            n_total = len(pae_np)
            px_per_res = max(1.5, min(4.0, 600 / n_total))
            fig_size = max(5, n_total * px_per_res / 100)

            fig, ax = plt.subplots(
                figsize=(fig_size + 1.2, fig_size), dpi=90)
            fig.patch.set_facecolor('white')

            im = ax.imshow(pae_np, cmap=pae_cmap,
                           vmin=0, vmax=31.75,
                           aspect='equal',
                           interpolation='nearest',
                           origin='upper')

            cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
            cbar.set_label('PAE (Å)', fontsize=9)
            cbar.ax.tick_params(labelsize=8)

            # Chain dividers
            for d in dividers:
                ax.axvline(d, color='black', lw=1.2)
                ax.axhline(d, color='black', lw=1.2)

            # Axis labels: one per chain block
            ax.set_xticks(ticks_pos)
            ax.set_xticklabels(ticks_lbl, fontsize=8, rotation=30, ha='right')
            ax.set_yticks(ticks_pos)
            ax.set_yticklabels(ticks_lbl, fontsize=8)
            ax.set_xlabel('Scored residue', fontsize=9, labelpad=4)
            ax.set_ylabel('Aligned residue', fontsize=9, labelpad=4)
            ax.set_title(
                f"PAE — {res['job_name']}  ({n_chains} chain(s))",
                fontsize=9, pad=6)

            # Contact residue markers (<thresh Å in off-diag quadrants)
            #
            # v1.15 BUG FIX: the previous version drew every contact as
            # a single dot at the horizontal centre of the off-diagonal
            # block (x = cum_j + len_b/2), producing a vertical stripe
            # in the middle of each block instead of marking the actual
            # residue positions. It also only marked one quadrant per
            # pair. We now draw a small outline box at each (row_a,
            # col_b) cell whose inter-chain PAE is below the threshold,
            # for every (i_c, j_c) off-diagonal block — so both the AB
            # and BA quadrants are marked.
            marker_size = max(4, min(18, 900 / max(n_total, 1)))
            cum_i = 0
            for i_c, ca in enumerate(chain_order):
                len_a = chain_lens.get(ca, 0)
                cum_j = 0
                for j_c, cb in enumerate(chain_order):
                    len_b = chain_lens.get(cb, 0)
                    if i_c != j_c and len_a > 0 and len_b > 0:
                        sub = pae_np[cum_i:cum_i + len_a,
                                     cum_j:cum_j + len_b]
                        # Per-row mean in this quadrant → rows of chain
                        # ca that are in contact with chain cb.
                        mA = sub.mean(axis=1)
                        contact_rows = np.where(mA < thresh)[0]
                        if len(contact_rows):
                            # Draw one horizontal band per contact row,
                            # spanning the full width of the partner
                            # chain's block — this is what users expect
                            # for an "interface" annotation on a PAE map.
                            for rr in contact_rows:
                                ax.add_patch(plt.Rectangle(
                                    (cum_j - 0.5, cum_i + rr - 0.5),
                                    len_b, 1,
                                    facecolor='none',
                                    edgecolor='#00E676',
                                    linewidth=0.6,
                                    alpha=0.9,
                                    zorder=4))
                            # Emphasize the low-PAE cells themselves
                            # with small dots in each quadrant.
                            low = np.argwhere(sub < thresh)
                            if len(low):
                                # Sub-sample if too many points
                                if len(low) > 4000:
                                    step = max(1, len(low) // 4000)
                                    low = low[::step]
                                ax.scatter(
                                    cum_j + low[:, 1],
                                    cum_i + low[:, 0],
                                    s=marker_size * 0.15,
                                    color='#00E676',
                                    alpha=0.75,
                                    zorder=5, linewidths=0)
                    cum_j += len_b
                cum_i += len_a

            plt.tight_layout()

            # ── Interactive canvas with hover ───────────────────
            canvas = FigureCanvas(fig)
            canvas.setMinimumHeight(int(fig_size * 90))

            # Hover: show PAE value at cursor position
            hover_lbl = self._af3a_hover_lbl
            chain_order_ref = chain_order
            chain_lens_ref  = chain_lens
            c2o_ref         = chain_to_orf
            pae_np_ref      = pae_np
            tri_ref         = res.get('token_res_ids')

            def _on_hover(event, _pae=pae_np_ref,
                          _co=chain_order_ref, _cl=chain_lens_ref,
                          _c2o=c2o_ref, _tri=tri_ref,
                          _ax=ax, _lbl=hover_lbl):
                if event.inaxes != _ax:
                    _lbl.setText("")
                    return
                xi, yi = int(round(event.xdata)), int(round(event.ydata))
                n = len(_pae)
                if 0 <= xi < n and 0 <= yi < n:
                    val = _pae[yi, xi]
                    # Determine which chains xi and yi belong to
                    def _chain_of(pos):
                        cum = 0
                        for cid in _co:
                            cum += _cl.get(cid, 0)
                            if pos < cum:
                                return _c2o.get(cid, cid)
                        return '?'
                    cx = _chain_of(xi)
                    cy = _chain_of(yi)
                    # Use real residue numbers from token_res_ids when
                    # available (restarts at 1 per chain, matching AF3
                    # server display). Fall back to running index+1.
                    if _tri and len(_tri) == n:
                        rx = int(_tri[xi])
                        ry = int(_tri[yi])
                    else:
                        rx = xi + 1
                        ry = yi + 1
                    _lbl.setText(
                        f"  PAE  scored={cx} res{rx}  "
                        f"aligned={cy} res{ry}  →  "
                        f"{val:.2f} Å")

            # ── Motif overlays (v1.16) ────────────────────────────
            # Draw a numbered rectangle around every detected motif on
            # BOTH off-diagonal quadrants (A↔B solid, B↔A dashed mirror).
            # The number is the motif's global rank in the job. Uses a
            # warm red that is clearly visible against the blue PAE.
            motifs = res.get('motifs') or []
            for m in motifs:
                r0 = m['abs_row0']; r1 = m['abs_row1']
                c0 = m['abs_col0']; c1 = m['abs_col1']
                # Pick colour by score: green (strong) → amber → red
                if   m['score'] >= 60: colr = '#00C853'
                elif m['score'] >= 40: colr = '#FFAB00'
                else:                  colr = '#FF1744'
                # A↔B quadrant
                ax.add_patch(plt.Rectangle(
                    (c0 - 0.5, r0 - 0.5),
                    c1 - c0 + 1, r1 - r0 + 1,
                    facecolor='none', edgecolor=colr,
                    linewidth=1.6, zorder=7))
                # Rank label above upper-left corner, readable on any bg
                ax.annotate(
                    f"#{m['rank']}",
                    xy=(c0 - 0.5, r0 - 0.5),
                    xytext=(c0 + 1, r0 - 3),
                    fontsize=8, fontweight='bold',
                    color=colr,
                    bbox=dict(boxstyle='round,pad=0.15',
                              facecolor='white', edgecolor=colr,
                              linewidth=0.8, alpha=0.9),
                    zorder=8)
                # Mirrored B↔A quadrant (dashed, no label, lower alpha)
                ax.add_patch(plt.Rectangle(
                    (r0 - 0.5, c0 - 0.5),
                    r1 - r0 + 1, c1 - c0 + 1,
                    facecolor='none', edgecolor=colr,
                    linewidth=1.0, linestyle='--',
                    alpha=0.75, zorder=7))

            canvas.mpl_connect('motion_notify_event', _on_hover)

            # ── Right-click context menu on PAE canvas ────────────────────
            job_name_ref = res['job_name']
            def _pae_right_click(event, _fig=fig, _name=job_name_ref):
                if event.button != 3:   # 3 = right mouse button
                    return
                self._af3a_canvas_context_menu(_fig, _name + '_PAE')
            canvas.mpl_connect('button_press_event', _pae_right_click)

            self._af3a_active_canvases.append((fig, canvas,
                                                res['job_name'] + '_PAE'))
            # Stash refs so the motif-table click handler can draw a
            # highlight on the SAME axes without rebuilding the plot.
            self._af3a_last_pae_ax = ax
            self._af3a_last_pae_canvas = canvas
            self._af3a_highlight_patches = []
            plt.close(fig)
            self._af3a_plot_layout.addWidget(canvas)

        # ── pLDDT PLOT ─────────────────────────────────────────
        if plddt and chain_order:
            plddt_np = np.array(plddt, dtype=float)
            n_total  = len(plddt_np)
            fig_w    = max(6, n_total / 30)
            fig2, ax2 = plt.subplots(figsize=(fig_w, 2.6), dpi=90)
            fig2.patch.set_facecolor('white')

            chain_colors = ['#1565C0','#C62828','#2E7D32',
                            '#E65100','#6A1B9A','#00695C',
                            '#4A148C','#827717','#3E2723','#0D47A1']
            cum = 0
            for ci, cid in enumerate(chain_order):
                clen = chain_lens.get(cid, 0)
                if clen <= 0:
                    continue
                xs  = np.arange(cum, cum + clen)
                ys  = plddt_np[cum:cum + clen]
                col = chain_colors[ci % len(chain_colors)]
                lbl = chain_to_orf.get(cid, cid)
                ax2.plot(xs, ys, color=col, lw=0.9,
                         label=lbl, alpha=0.9)
                ax2.fill_between(xs, ys, alpha=0.12, color=col)
                cum += clen

            for d in dividers:
                ax2.axvline(d, color='black', lw=0.8)
            ax2.axhline(70, color='#2E7D32', lw=0.8, ls='--',
                        alpha=0.6, label='70 (good)')
            ax2.axhline(50, color='#E65100', lw=0.8, ls='--',
                        alpha=0.6, label='50 (low)')

            ax2.set_xticks(ticks_pos)
            ax2.set_xticklabels(ticks_lbl, fontsize=8,
                                rotation=20, ha='right')
            ax2.set_ylim(0, 100)
            ax2.set_ylabel('pLDDT', fontsize=9)
            ax2.set_title(
                f"pLDDT — {res['job_name']}", fontsize=9, pad=4)
            ax2.legend(fontsize=7, loc='lower right',
                       ncol=min(n_chains + 2, 6), framealpha=0.7)
            ax2.tick_params(labelsize=8)
            ax2.grid(axis='y', alpha=0.25)
            plt.tight_layout()

            canvas2 = FigureCanvas(fig2)
            canvas2.setMinimumHeight(220)

            # ── Right-click context menu on pLDDT canvas ─────────────────
            def _plddt_right_click(event, _fig=fig2,
                                   _name=res['job_name'] + '_pLDDT'):
                if event.button != 3:
                    return
                self._af3a_canvas_context_menu(_fig, _name)
            canvas2.mpl_connect('button_press_event', _plddt_right_click)

            self._af3a_active_canvases.append((fig2, canvas2,
                                                res['job_name'] + '_pLDDT'))
            plt.close(fig2)
            self._af3a_plot_layout.addWidget(canvas2)

        # ── Motif table (v1.16) ─────────────────────────────────────
        # Must come after both plots so it renders below them.
        try:
            self._af3a_build_motif_table(res)
        except Exception as e:
            print(f"[AF3 Analysis] motif table build failed: {e}")

        self._af3a_plot_layout.addStretch()

    def _af3a_canvas_context_menu(self, fig, default_name: str):
        """Show a right-click context menu on a matplotlib canvas.
        Allows saving/exporting the figure as PNG, PDF, SVG, or TIFF."""
        menu = QMenu(self)

        # ── PNG (high-res, lossless) ──────────────────────────────────────
        act_png = menu.addAction("🖼  Save as PNG (300 dpi)…")
        act_png.setToolTip("High-resolution raster image — best for presentations and papers.")

        act_png_screen = menu.addAction("🖼  Save as PNG (screen resolution)…")
        act_png_screen.setToolTip("Exact screen resolution — smaller file, good for quick sharing.")

        menu.addSeparator()

        # ── PDF (vector) ──────────────────────────────────────────────────
        act_pdf = menu.addAction("📄  Save as PDF (vector)…")
        act_pdf.setToolTip("Scalable vector PDF — best for publications and figure editors.")

        # ── SVG (vector) ──────────────────────────────────────────────────
        act_svg = menu.addAction("📐  Save as SVG (vector)…")
        act_svg.setToolTip("Scalable Vector Graphics — editable in Inkscape / Illustrator.")

        menu.addSeparator()

        # ── TIFF ─────────────────────────────────────────────────────────
        act_tif = menu.addAction("🔬  Save as TIFF (600 dpi, publication)…")
        act_tif.setToolTip("High-resolution TIFF at 600 dpi — standard format for journal submission.")

        menu.addSeparator()

        # ── Copy to clipboard ─────────────────────────────────────────────
        act_clip = menu.addAction("📋  Copy to clipboard (PNG)…")
        act_clip.setToolTip("Render to PNG and copy to system clipboard.")

        menu.addSeparator()

        # ── Export all canvases ───────────────────────────────────────────
        act_all_pdf = menu.addAction("📑  Export all plots to PDF…")
        act_all_pdf.setToolTip("Combine all currently visible PAE/pLDDT plots into one PDF.")

        chosen = menu.exec(
            QCursor.pos() if QT_VERSION == 6
            else QCursor.pos())
        if chosen is None:
            return

        safe = re.sub(r'[^\w\-.]', '_', default_name)

        if chosen == act_png:
            path, _ = QFileDialog.getSaveFileName(
                self, "Save as PNG", f"{safe}.png",
                "PNG Image (*.png);;All (*)")
            if path:
                try:
                    fig.savefig(path, dpi=300, bbox_inches='tight',
                                facecolor='white')
                    self._status.showMessage(f"✓ Saved: {path}")
                except Exception as e:
                    QMessageBox.critical(self, "Save PNG", str(e))

        elif chosen == act_png_screen:
            path, _ = QFileDialog.getSaveFileName(
                self, "Save as PNG (screen)", f"{safe}_screen.png",
                "PNG Image (*.png);;All (*)")
            if path:
                try:
                    fig.savefig(path, dpi=fig.dpi, bbox_inches='tight',
                                facecolor='white')
                    self._status.showMessage(f"✓ Saved: {path}")
                except Exception as e:
                    QMessageBox.critical(self, "Save PNG", str(e))

        elif chosen == act_pdf:
            path, _ = QFileDialog.getSaveFileName(
                self, "Save as PDF", f"{safe}.pdf",
                "PDF (*.pdf);;All (*)")
            if path:
                try:
                    fig.savefig(path, bbox_inches='tight',
                                facecolor='white')
                    self._status.showMessage(f"✓ Saved: {path}")
                except Exception as e:
                    QMessageBox.critical(self, "Save PDF", str(e))

        elif chosen == act_svg:
            path, _ = QFileDialog.getSaveFileName(
                self, "Save as SVG", f"{safe}.svg",
                "SVG (*.svg);;All (*)")
            if path:
                try:
                    fig.savefig(path, bbox_inches='tight',
                                facecolor='white')
                    self._status.showMessage(f"✓ Saved: {path}")
                except Exception as e:
                    QMessageBox.critical(self, "Save SVG", str(e))

        elif chosen == act_tif:
            path, _ = QFileDialog.getSaveFileName(
                self, "Save as TIFF", f"{safe}.tiff",
                "TIFF (*.tiff *.tif);;All (*)")
            if path:
                try:
                    fig.savefig(path, dpi=600, bbox_inches='tight',
                                facecolor='white')
                    self._status.showMessage(f"✓ Saved: {path}")
                except Exception as e:
                    QMessageBox.critical(self, "Save TIFF", str(e))

        elif chosen == act_clip:
            try:
                import io
                buf = io.BytesIO()
                fig.savefig(buf, format='png', dpi=150, bbox_inches='tight',
                            facecolor='white')
                buf.seek(0)
                try:
                    from PyQt6.QtGui import QImage
                except ImportError:
                    from PyQt5.QtGui import QImage
                data = buf.read()
                img  = QImage.fromData(data)
                QApplication.clipboard().setImage(img)
                self._status.showMessage("✓ PAE plot copied to clipboard")
            except Exception as e:
                QMessageBox.critical(self, "Copy to Clipboard", str(e))

        elif chosen == act_all_pdf:
            self._af3a_export_pdf()

    def _af3a_export_pdf(self):
        """Export all active PAE/pLDDT figures to a single PDF."""
        if not MATPLOTLIB_AVAILABLE:
            QMessageBox.warning(self, "Export PDF",
                "matplotlib not installed.")
            return
        canvases = getattr(self, '_af3a_active_canvases', [])
        if not canvases:
            QMessageBox.information(self, "Export PDF",
                "No plots to export.\nSelect a job first.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export plots to PDF", "",
            "PDF (*.pdf);;All (*)")
        if not path:
            return
        try:
            from matplotlib.backends.backend_pdf import PdfPages
            with PdfPages(path) as pdf:
                for fig, canvas, name in canvases:
                    pdf.savefig(fig, bbox_inches='tight',
                                facecolor='white')
            self._status.showMessage(
                f"✓ Exported {len(canvases)} plot(s) → {Path(path).name}")
        except Exception as e:
            QMessageBox.critical(self, "Export PDF",
                                 f"Error:\n{e}")

    def _af3a_apply_filters(self):
        """Hide rows that fail the min-ipTM / max-PAEinter filters.

        Rows with None values in the filtered metric are shown only when
        the filter is at its permissive default (ipTM >= 0.0 or
        PAEinter <= 32.0). This prevents surprise hiding of legacy
        jobs that lack a given metric.
        """
        if not hasattr(self, '_af3a_filter_iptm'):
            return
        min_iptm = self._af3a_filter_iptm.value()
        max_paei = self._af3a_filter_paei.value()
        user_role = self._af3a_user_role()
        n_vis = 0
        for row in range(self._af3a_table.rowCount()):
            item0 = self._af3a_table.item(row, 0)
            if item0 is None:
                continue
            idx = item0.data(user_role)
            try:
                res = self._af3_analysis_results[int(idx)]
            except (IndexError, TypeError, ValueError):
                continue
            iptm = res.get('iptm')
            paei = res.get('pae_inter')
            hide = False
            if min_iptm > 0.0:
                if iptm is None or iptm < min_iptm:
                    hide = True
            if not hide and max_paei < 32.0:
                if paei is None or paei > max_paei:
                    hide = True
            self._af3a_table.setRowHidden(row, hide)
            if not hide:
                n_vis += 1
        total = self._af3a_table.rowCount()
        if total:
            extra = (f"  ({n_vis}/{total} shown)"
                     if n_vis != total else "")
            txt = self._af3a_status.text()
            # Strip any previous "(x/y shown)" suffix
            base = re.sub(r'\s*\(\d+/\d+ shown\)\s*$', '', txt)
            self._af3a_status.setText(base + extra)

    def _af3a_export_tsv(self):
        """Export the full AF3 Analysis results table to a TSV file."""
        results = getattr(self, '_af3_analysis_results', [])
        if not results:
            QMessageBox.information(
                self, "Export TSV",
                "No results to export.\nLoad an AF3 results folder first.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export AF3 results to TSV",
            "af3_analysis.tsv",
            "TSV (*.tsv);;All (*)")
        if not path:
            return
        headers = [
            'job_name', 'job_dir', 'n_chains', 'chains',
            'iptm', 'ptm', 'mean_pLDDT', 'ranking_score',
            'pae_inter_mean', 'pae_min_inter', 'cp_iptm_inter', 'contact_frac_pct',
            'pae_min_best_pair', 'pair_iptm_best_pair',
            'best_pair', 'contact_region',
            'chain_iptm', 'chain_ptm',
            'fraction_disordered', 'has_clash',
            'n_diffusion_samples',
            'high_confidence',
        ]

        def _fmt(v, fmt='{:.4f}'):
            if v is None:
                return ''
            try:
                return fmt.format(v)
            except (TypeError, ValueError):
                return str(v)

        try:
            with open(path, 'w', encoding='utf-8') as f:
                f.write('\t'.join(headers) + '\n')
                for res in results:
                    chains = ';'.join(
                        f"{res['chain_to_orf'].get(c, c)}({res['chain_lens'][c]})"
                        for c in res.get('chain_order', []))
                    bp = res.get('best_pair', ('?', '?'))
                    bp_str = (f"{res['chain_to_orf'].get(bp[0], bp[0])}"
                              f"-{res['chain_to_orf'].get(bp[1], bp[1])}"
                              if isinstance(bp, tuple) and len(bp) == 2
                              else '')
                    iptm = res.get('iptm')
                    paei = res.get('pae_inter')
                    hc = ('yes' if (iptm is not None and paei is not None
                                    and iptm > 0.75 and paei < 8.0)
                          else 'no')
                    # Rich per-pair metrics for best pair
                    bp_pm = res.get('pair_metrics', {}).get(bp, {})
                    pae_min_s     = _fmt(bp_pm.get('pair_pae_min'), '{:.3f}')
                    cp_iptm_s     = _fmt(res.get('cp_iptm_inter'), '{:.3f}')
                    contact_frac_s= _fmt((res.get('contact_frac') or 0) * 100, '{:.1f}')
                    pair_iptm_s = _fmt(bp_pm.get('pair_iptm'), '{:.3f}')
                    # Per-chain metric stringification
                    ci_s = (';'.join(_fmt(x, '{:.3f}')
                                     for x in res['chain_iptm'])
                            if res.get('chain_iptm') else '')
                    cp_s = (';'.join(_fmt(x, '{:.3f}')
                                     for x in res['chain_ptm'])
                            if res.get('chain_ptm') else '')
                    has_clash = res.get('has_clash')
                    clash_s = ('yes' if (has_clash is True
                                         or (isinstance(has_clash,
                                                         (int, float))
                                             and has_clash >= 0.5))
                               else ('no' if has_clash is not None else ''))
                    row = [
                        res.get('job_name', ''),
                        res.get('job_dir', ''),
                        str(res.get('n_chains', '')),
                        chains,
                        _fmt(iptm),
                        _fmt(res.get('ptm')),
                        _fmt(res.get('mean_plddt'), '{:.2f}'),
                        _fmt(res.get('ranking_score')),
                        _fmt(paei, '{:.3f}'),
                        pae_min_s,
                        cp_iptm_s,
                        contact_frac_s,
                        pair_iptm_s,
                        bp_str,
                        (res.get('contact_region') or '')
                            .replace('\t', ' ').replace('\n', ' '),
                        ci_s,
                        cp_s,
                        _fmt(res.get('fraction_disordered'), '{:.4f}'),
                        clash_s,
                        str(len(res.get('ranking_samples', []))),
                        hc,
                    ]
                    f.write('\t'.join(row) + '\n')
            self._status.showMessage(
                f"✓ Exported {len(results)} row(s) → {Path(path).name}")
        except Exception as e:
            QMessageBox.critical(self, "Export TSV",
                                 f"Error:\n{e}")

    # ═══════════════════════════════════════════════════════════
    # Motif detection — UI handlers (v1.16)
    # ═══════════════════════════════════════════════════════════

    def _af3a_rerun_motifs(self):
        """Re-run motif detection on all loaded jobs and refresh UI.

        Called when any motif-control spinbox / checkbox changes. Also
        called after a job table selection to make sure the motif table
        and plot overlays are synced with the selected job.

        We do this in-place (mutating res['motifs']) so that the motif
        count label and any export reflects the current tier settings
        without requiring a folder re-scan.
        """
        results = getattr(self, '_af3_analysis_results', [])
        if not results:
            self._af3a_motif_count_lbl.setText("—")
            return
        try:
            pae_core = self._af3a_motif_core_spin.value()
            pae_ext  = self._af3a_motif_ext_spin.value()
            min_ct   = self._af3a_motif_contact_spin.value()
            min_sz   = self._af3a_motif_size_spin.value()
            recip    = self._af3a_motif_recip_cb.isChecked()
        except AttributeError:
            return

        # Enforce extended > core (otherwise extension is a no-op)
        if pae_ext <= pae_core:
            pae_ext = pae_core + 1.0

        total = 0
        with_motifs = 0
        for res in results:
            try:
                motifs = self._af3a_detect_motifs(
                    res, pae_core, pae_ext, min_ct, min_sz,
                    require_reciprocal=recip)
            except Exception as e:
                print(f"[AF3 Analysis] motif rerun failed for "
                      f"{res.get('job_name','?')}: {e}")
                motifs = []
            res['motifs'] = motifs
            total += len(motifs)
            if motifs:
                with_motifs += 1

        self._af3a_motif_count_lbl.setText(
            f"{total} motif(s) across {with_motifs}/{len(results)} job(s)")

        # Refresh the currently-selected job's plots + motif table
        self._af3a_replot_selected()

    def _af3a_export_motifs_tsv(self):
        """Export every detected motif from every loaded job to TSV."""
        results = getattr(self, '_af3_analysis_results', [])
        if not results:
            QMessageBox.information(
                self, "Export motifs TSV",
                "No results loaded.")
            return
        # Count motifs up front so we don't prompt for a path then
        # write an empty file.
        total = sum(len(r.get('motifs') or []) for r in results)
        if total == 0:
            QMessageBox.information(
                self, "Export motifs TSV",
                "No motifs detected with the current thresholds.\n"
                "Try relaxing 'core PAE' or 'min contact'.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export motifs TSV",
            "af3_motifs.tsv",
            "TSV (*.tsv);;All (*)")
        if not path:
            return

        def _fmt(v, f='{:.4f}'):
            if v is None: return ''
            try: return f.format(v)
            except (TypeError, ValueError): return str(v)

        headers = [
            'job_name', 'rank', 'score',
            'chain_A', 'orf_A', 'a_start', 'a_end', 'len_A_motif',
            'chain_B', 'orf_B', 'b_start', 'b_end', 'len_B_motif',
            'mean_PAE', 'min_PAE', 'density',
            'mean_contact', 'max_contact',
            'plddt_A', 'plddt_B', 'plddt_mean',
            'reciprocal', 'recip_overlap',
            'label',
        ]
        try:
            with open(path, 'w', encoding='utf-8') as f:
                f.write('\t'.join(headers) + '\n')
                for res in results:
                    for m in (res.get('motifs') or []):
                        row = [
                            res.get('job_name', ''),
                            str(m.get('rank', '')),
                            _fmt(m.get('score'), '{:.2f}'),
                            m.get('chain_A', ''),
                            m.get('name_A', ''),
                            str(m.get('a_start', '')),
                            str(m.get('a_end', '')),
                            str(m.get('height', '')),
                            m.get('chain_B', ''),
                            m.get('name_B', ''),
                            str(m.get('b_start', '')),
                            str(m.get('b_end', '')),
                            str(m.get('width', '')),
                            _fmt(m.get('mean_pae'), '{:.3f}'),
                            _fmt(m.get('min_pae'),  '{:.3f}'),
                            _fmt(m.get('density'),  '{:.3f}'),
                            _fmt(m.get('mean_contact')),
                            _fmt(m.get('max_contact')),
                            _fmt(m.get('plddt_A'), '{:.2f}'),
                            _fmt(m.get('plddt_B'), '{:.2f}'),
                            _fmt(m.get('plddt_mean'), '{:.2f}'),
                            'yes' if m.get('reciprocal') else 'no',
                            _fmt(m.get('recip_overlap'), '{:.2f}'),
                            m.get('label', ''),
                        ]
                        f.write('\t'.join(row) + '\n')
            self._status.showMessage(
                f"✓ Exported {total} motif(s) → {Path(path).name}")
        except OSError as e:
            QMessageBox.critical(self, "Export motifs TSV",
                                 f"Error:\n{e}")

    def _af3a_build_motif_table(self, res: dict):
        """Build a QTableWidget listing every detected motif for the
        selected job. Appended to the plot layout below the pLDDT plot.

        Clicking a row highlights that motif on the PAE heatmap.
        """
        motifs = res.get('motifs') or []
        lbl = QLabel(
            f"<b>🎯 Detected interaction motifs</b>  "
            f"— {len(motifs)} found"
            + ("" if motifs
               else "  <span style='color:#888'>(try relaxing thresholds)</span>"))
        lbl.setStyleSheet(
            "font-size:11px;padding:6px 8px;"
            "background:#ECEFF1;border-radius:4px;margin-top:8px;")
        lbl.setTextFormat(
            Qt.TextFormat.RichText if QT_VERSION == 6 else Qt.RichText)
        self._af3a_plot_layout.addWidget(lbl)

        if not motifs:
            return

        table = QTableWidget()
        # NB: no sorting — rows are already sorted by score desc and the
        # row index must match `motifs[i]` for the highlight click.
        table.setColumnCount(10)
        table.setHorizontalHeaderLabels([
            '#', 'Score', 'ORF A region', 'ORF B region',
            'Size (aa)', 'PAE mean', 'PAE min',
            'Density', 'pLDDT', 'Recip.'])
        table.setRowCount(len(motifs))
        table.setAlternatingRowColors(True)
        table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows
            if QT_VERSION == 6 else QAbstractItemView.SelectRows)
        table.setEditTriggers(
            QAbstractItemView.EditTrigger.NoEditTriggers
            if QT_VERSION == 6 else QAbstractItemView.NoEditTriggers)
        table.verticalHeader().setVisible(False)
        # Compact; it sits below the plots and is scrollable via parent
        row_h = 22
        table.setMaximumHeight(
            min(8, len(motifs)) * row_h + table.horizontalHeader().height()
            + 4)

        def _plddt_color(v):
            if v is None: return None
            if v >= 90: return QColor('#0D47A1')   # very high
            if v >= 70: return QColor('#1976D2')   # confident
            if v >= 50: return QColor('#FBC02D')   # low
            return QColor('#E64A19')               # very low

        for i, m in enumerate(motifs):
            size_s = f"{m['height']}×{m['width']}"
            pae_mean_s = f"{m['mean_pae']:.1f}"
            pae_min_s  = f"{m['min_pae']:.1f}"
            dens_s     = f"{m['density']*100:.0f}%"
            plddt_s    = (f"{m['plddt_mean']:.0f}"
                          if m.get('plddt_mean') is not None else '-')
            recip_s    = ('✓ {:.0%}'.format(m['recip_overlap'])
                          if m.get('reciprocal') else '✗')
            cells = [
                f"#{m['rank']}",
                f"{m['score']:.1f}",
                f"{m['name_A']} {m['a_start']}-{m['a_end']}",
                f"{m['name_B']} {m['b_start']}-{m['b_end']}",
                size_s,
                pae_mean_s, pae_min_s, dens_s, plddt_s, recip_s,
            ]
            for c, val in enumerate(cells):
                it = QTableWidgetItem(val)
                it.setTextAlignment(
                    Qt.AlignmentFlag.AlignCenter if QT_VERSION == 6
                    else Qt.AlignCenter)
                if c == 1:   # Score column — colour by value
                    if m['score'] >= 60:
                        it.setBackground(QColor('#C8E6C9'))
                    elif m['score'] >= 40:
                        it.setBackground(QColor('#FFF9C4'))
                    else:
                        it.setBackground(QColor('#FFCDD2'))
                if c == 8 and m.get('plddt_mean') is not None:
                    pc = _plddt_color(m['plddt_mean'])
                    if pc is not None:
                        f = it.font(); f.setBold(True); it.setFont(f)
                        it.setForeground(pc)
                table.setItem(i, c, it)
        table.resizeColumnsToContents()
        table.horizontalHeader().setStretchLastSection(True)
        self._af3a_plot_layout.addWidget(table)

        # Clicking a row highlights the motif in the PAE heatmap. We
        # do this by storing the latest PAE axes reference on self and
        # drawing a highlight rectangle that survives until another
        # motif is clicked.
        def _on_motif_select():
            rows = set(idx.row() for idx in table.selectedIndexes())
            if not rows:
                return
            i = min(rows)
            if 0 <= i < len(motifs):
                self._af3a_highlight_motif(motifs[i])
        table.selectionModel().selectionChanged.connect(
            lambda *_: _on_motif_select())

    def _af3a_highlight_motif(self, motif: dict):
        """Flash a highlight rectangle around a motif on the PAE axes.

        Uses the matplotlib axes reference stashed on self during the
        last plot, if any.
        """
        ax = getattr(self, '_af3a_last_pae_ax', None)
        canvas = getattr(self, '_af3a_last_pae_canvas', None)
        if ax is None or canvas is None:
            return
        # Remove previous highlight if any
        prev = getattr(self, '_af3a_highlight_patches', [])
        for p in prev:
            try: p.remove()
            except Exception: pass
        patches = []
        # Draw on A↔B quadrant
        r0 = motif['abs_row0']; r1 = motif['abs_row1']
        c0 = motif['abs_col0']; c1 = motif['abs_col1']
        patches.append(ax.add_patch(plt.Rectangle(
            (c0 - 0.5, r0 - 0.5), c1 - c0 + 1, r1 - r0 + 1,
            facecolor='none', edgecolor='#FF1744',
            linewidth=2.5, linestyle='-', zorder=10)))
        # And the mirrored B↔A quadrant (show the user both sides)
        patches.append(ax.add_patch(plt.Rectangle(
            (r0 - 0.5, c0 - 0.5), r1 - r0 + 1, c1 - c0 + 1,
            facecolor='none', edgecolor='#FF1744',
            linewidth=2.5, linestyle='--', alpha=0.8, zorder=10)))
        self._af3a_highlight_patches = patches
        canvas.draw_idle()

    # ═══════════════════════════════════════════════════════════
    # HPC SERVER TAB (generic)
    # ═══════════════════════════════════════════════════════════

    # ══════════════════════════════════════════════════════════════════════
    # PPI GENOMIC ARC MAP TAB  (Tab 11, v2.0)
    # Shows all AF3-predicted interactions as arcs drawn above a linear
    # genome map. Arcs are coloured by PAE_min / ipTM and height-scaled
    # by genomic distance between the two partners.
    # ══════════════════════════════════════════════════════════════════════

    def _create_ppi_arc_map_tab(self):
        w = QWidget()
        self._ppi_map_tab_widget = w
        lay = QVBoxLayout(w)
        lay.setContentsMargins(2, 2, 2, 2)
        lay.setSpacing(2)

        # ── toolbar ─────────────────────────────────────────────────────
        tb = QHBoxLayout()
        tb.setSpacing(4)

        tb.addWidget(QLabel("Show:"))
        self._ppi_arc_filter_combo = QComboBox()
        self._ppi_arc_filter_combo.addItems([
            "All results",
            "HIGH only  (PAE_min < 4 Å)",
            "HIGH + MED  (PAE_min < 8 Å)",
            "Focal hits  (PAE_min < 4 Å & cp_ipTM ≥ 0.50)",
        ])
        self._ppi_arc_filter_combo.setToolTip(
            "Filter which interaction arcs are drawn.\n"
            "HIGH = PAE_min < 4 Å (focal domain contact confirmed)\n"
            "MED  = PAE_min 4–8 Å (possible contact, check motifs)\n"
            "Focal = HIGH AND cp_ipTM ≥ 0.50 (strictest)")
        self._ppi_arc_filter_combo.currentIndexChanged.connect(
            self._ppi_arc_map_refresh)
        tb.addWidget(self._ppi_arc_filter_combo)

        tb.addWidget(QLabel("Color by:"))
        self._ppi_arc_color_combo = QComboBox()
        self._ppi_arc_color_combo.addItems(
            ["PAE_min ★", "ipTM", "cp_ipTM ★", "Contact%"])
        self._ppi_arc_color_combo.setToolTip(
            "Metric used to determine arc colour:\n"
            "  PAE_min ★ — minimum inter-chain PAE (recommended)\n"
            "  ipTM      — global interface pTM score\n"
            "  cp_ipTM ★ — chain-pair ipTM (focal, from summary JSON)\n"
            "  Contact%  — fraction of residue pairs with PAE < 5 Å")
        self._ppi_arc_color_combo.currentIndexChanged.connect(
            self._ppi_arc_map_refresh)
        tb.addWidget(self._ppi_arc_color_combo)

        tb.addWidget(QLabel("Arc height:"))
        self._ppi_arc_height_combo = QComboBox()
        self._ppi_arc_height_combo.addItems(
            ["By genomic distance", "By score (better = lower)", "Fixed"])
        self._ppi_arc_height_combo.setToolTip(
            "How arc height is determined:\n"
            "  By genomic distance — neighboring ORFs get low arcs,\n"
            "    distant ORFs get tall arcs (like published operon figures)\n"
            "  By score — high-confidence interactions drawn lower/closer\n"
            "  Fixed — all arcs same height")
        self._ppi_arc_height_combo.currentIndexChanged.connect(
            self._ppi_arc_map_refresh)
        tb.addWidget(self._ppi_arc_height_combo)

        btn_refresh = QPushButton("⟳ Refresh")
        btn_refresh.setToolTip("Rebuild the arc map from current AF3 results")
        btn_refresh.clicked.connect(self._ppi_arc_map_refresh)
        btn_refresh.setStyleSheet(
            "QPushButton{background:#2e7d32;color:white;font-weight:bold;"
            "border-radius:4px;padding:2px 8px;}"
            "QPushButton:hover{background:#388e3c;}")
        tb.addWidget(btn_refresh)

        btn_svg = QPushButton("↓ Export SVG")
        btn_svg.setToolTip("Export the arc map as a scalable vector SVG file")
        btn_svg.clicked.connect(self._ppi_arc_map_export_svg)
        tb.addWidget(btn_svg)

        btn_tsv = QPushButton("↓ Export TSV")
        btn_tsv.setToolTip("Export the displayed interactions as a TSV table")
        btn_tsv.clicked.connect(self._ppi_arc_map_export_tsv)
        tb.addWidget(btn_tsv)

        tb.addStretch()
        self._ppi_arc_count_lbl = QLabel("0 interactions")
        self._ppi_arc_count_lbl.setStyleSheet(
            "font-weight:500;color:#1D9E75;")
        tb.addWidget(self._ppi_arc_count_lbl)

        lay.addLayout(tb)

        # ── info label (shows hovered/clicked arc) ───────────────────────
        self._ppi_arc_info_lbl = QLabel(
            "Click an arc to see interaction details  |  "
            "Click an ORF to center genome map  |  "
            "Scroll = zoom  ·  Drag = pan")
        self._ppi_arc_info_lbl.setStyleSheet(
            "font-size:11px; color:#085041; background:#E1F5EE;"
            "padding:3px 7px; border-bottom:1px solid #9FE1CB;")
        self._ppi_arc_info_lbl.setWordWrap(True)
        lay.addWidget(self._ppi_arc_info_lbl)

        # ── arc map canvas ────────────────────────────────────────────────
        self._ppi_arc_widget = _PpiArcMapWidget(self)
        self._ppi_arc_widget.arc_clicked.connect(self._ppi_arc_on_click)
        self._ppi_arc_widget.orf_clicked.connect(self._ppi_arc_on_orf_click)
        self._ppi_arc_widget.arc_hovered.connect(self._ppi_arc_on_hover)
        lay.addWidget(self._ppi_arc_widget, stretch=1)

        # ── legend ────────────────────────────────────────────────────────
        leg = QHBoxLayout()
        leg.setSpacing(12)

        def _leg_line(color, dash=False):
            lbl = QLabel()
            lbl.setFixedSize(32, 12)
            lbl.setStyleSheet(
                f"border-top: {'3px solid' if not dash else '2px dashed'} {color};"
                f"margin-top:5px;")
            return lbl

        def _leg_rect(color, border):
            lbl = QLabel()
            lbl.setFixedSize(14, 10)
            lbl.setStyleSheet(
                f"background:{color};border:1px solid {border};"
                f"border-radius:2px;margin-top:1px;")
            return lbl

        for widget, text in [
            (_leg_line("#1D9E75"),          "PAE_min < 4 Å (HIGH)"),
            (_leg_line("#BA7517", dash=True),"PAE_min 4–8 Å (MED)"),
            (_leg_line("#E24B4A", dash=True),"PAE_min > 8 Å (LOW)"),
            (_leg_rect("#9FE1CB","#0F6E56"), "HMM hit ORF"),
            (_leg_rect("#FAC775","#BA7517"), "Custom / other"),
            (_leg_rect("#D3D1C7","#888780"), "No annotation"),
        ]:
            row = QHBoxLayout()
            row.setSpacing(4)
            row.addWidget(widget)
            row.addWidget(QLabel(text))
            leg.addLayout(row)
        leg.addStretch()

        leg_widget = QWidget()
        leg_widget.setLayout(leg)
        leg_widget.setStyleSheet(
            "background:var(--color-background-secondary);"
            "border-top:0.5px solid #ddd;padding:3px 0;")
        lay.addWidget(leg_widget)

        self._tabs.addTab(w, "🧬 Genomic PPI Map")

    def _ppi_arc_map_refresh(self):
        """Rebuild arc data from _af3_analysis_results and repaint."""
        results = getattr(self, '_af3_analysis_results', [])
        orfs    = getattr(self, 'orfs', [])

        filter_idx = self._ppi_arc_filter_combo.currentIndex()
        color_by   = self._ppi_arc_color_combo.currentText()
        height_by  = self._ppi_arc_height_combo.currentText()

        arcs = []
        for res in results:
            names = res.get('orf_names', [])
            if len(names) < 2:
                continue
            pae_min  = res.get('pae_min_inter')
            cp_iptm  = res.get('cp_iptm_inter')
            iptm     = res.get('iptm') or 0
            cfrac    = res.get('contact_frac') or 0

            # Apply filter
            if filter_idx == 1 and (pae_min is None or pae_min >= 4.0):
                continue
            if filter_idx == 2 and (pae_min is None or pae_min >= 8.0):
                continue
            if filter_idx == 3 and (pae_min is None or pae_min >= 4.0
                                    or (cp_iptm is not None and cp_iptm < 0.50)):
                continue

            # Resolve ORF positions — extract numeric index from name "ORF1234"
            def _find_orf(name, _orfs=orfs):
                # Primary: exact "ORF{n}" match
                if name.upper().startswith('ORF'):
                    try:
                        idx = int(name[3:]) - 1
                        if 0 <= idx < len(_orfs):
                            return idx, _orfs[idx]
                    except (ValueError, TypeError):
                        pass
                # Secondary: any digit run in the name
                try:
                    import re as _re
                    nums = _re.findall(r'\d+', name)
                    if nums:
                        idx = int(nums[-1]) - 1
                        if 0 <= idx < len(_orfs):
                            return idx, _orfs[idx]
                except (ValueError, TypeError):
                    pass
                return -1, None

            i_a, o_a = _find_orf(names[0])
            i_b, o_b = _find_orf(names[1])
            if o_a is None or o_b is None:
                continue

            # Score for colour
            if color_by.startswith("PAE_min"):
                score = pae_min
                score_inv = True   # lower = better = greener
            elif color_by.startswith("cp_ipTM"):
                score = cp_iptm
                score_inv = False
            elif color_by.startswith("Contact"):
                score = cfrac * 100
                score_inv = False
            else:
                score = iptm
                score_inv = False

            arcs.append({
                'name_a':   names[0],
                'name_b':   names[1],
                'orf_idx_a': i_a,
                'orf_idx_b': i_b,
                'start_a':  o_a['start'],
                'end_a':    o_a['end'],
                'start_b':  o_b['start'],
                'end_b':    o_b['end'],
                'pae_min':  pae_min,
                'iptm':     iptm,
                'cp_iptm':  cp_iptm,
                'contact_frac': cfrac,
                'score':    score,
                'score_inv': score_inv,
                'job_name': res.get('job_name', ''),
                'contact_region': res.get('contact_region', ''),
                'height_mode': height_by,
            })

        self._ppi_arc_widget.set_data(
            dna_length=len(getattr(self, 'dna_sequence', '') or ''),
            orfs=orfs,
            hmm_profiles=getattr(self, 'hmm_profiles', []),
            arcs=arcs,
        )
        self._ppi_arc_count_lbl.setText(
            f"{len(arcs)} interaction{'s' if len(arcs) != 1 else ''}")

    def _ppi_arc_on_click(self, arc_idx: int):
        """User clicked an arc — show details in info label."""
        arcs = self._ppi_arc_widget._arcs
        if not (0 <= arc_idx < len(arcs)):
            return
        a = arcs[arc_idx]
        pmin = f"{a['pae_min']:.2f} Å" if a['pae_min'] is not None else "—"
        cpip = f"{a['cp_iptm']:.2f}"   if a['cp_iptm'] is not None else "—"
        cfrc = f"{a['contact_frac']*100:.1f}%" if a['contact_frac'] is not None else "—"
        conf = ("HIGH ★" if a['pae_min'] is not None and a['pae_min'] < 4.0 else
                "MED"    if a['pae_min'] is not None and a['pae_min'] < 8.0 else "LOW")
        self._ppi_arc_info_lbl.setText(
            f"▶  {a['name_a']} ↔ {a['name_b']}   "
            f"ipTM={a['iptm']:.3f}  PAE_min={pmin}  cp_ipTM={cpip}  "
            f"Contact%={cfrc}  [{conf}]   "
            f"Contact zone: {a['contact_region'] or '—'}")
        # Also select in Interaction Results table
        for row in range(self._af3a_table.rowCount()):
            item = self._af3a_table.item(row, 0)
            if item and a['job_name'] in item.text():
                self._af3a_table.selectRow(row)
                break

    def _ppi_arc_on_orf_click(self, orf_idx: int):
        """User clicked an ORF node — center main genome map."""
        self._select_and_center_orf(orf_idx)

    def _ppi_arc_on_hover(self, arc_idx: int):
        """User hovered an arc — show lightweight tooltip in info label."""
        arcs = self._ppi_arc_widget._arcs
        if arc_idx < 0 or arc_idx >= len(arcs):
            self._ppi_arc_info_lbl.setText(
                "Click an arc to see interaction details  |  "
                "Click an ORF to center genome map  |  "
                "Scroll = zoom  ·  Drag = pan")
            return
        a = arcs[arc_idx]
        pmin = f"{a['pae_min']:.2f}" if a['pae_min'] is not None else "—"
        self._ppi_arc_info_lbl.setText(
            f"  {a['name_a']} ↔ {a['name_b']}  "
            f"PAE_min={pmin} Å  ipTM={a['iptm']:.2f}  "
            f"— click to select")

    def _ppi_arc_map_export_svg(self):
        """Export the arc map as SVG using QSvgGenerator."""
        try:
            from PyQt6.QtSvg import QSvgGenerator
        except ImportError:
            try:
                from PyQt5.QtSvg import QSvgGenerator
            except ImportError:
                QMessageBox.warning(self, "SVG Export",
                    "PyQt SVG module not available. Install PyQt6-Qt6 or PyQt5-sip.")
                return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export PPI Arc Map as SVG", "", "SVG (*.svg)")
        if not path:
            return
        gen = QSvgGenerator()
        gen.setFileName(path)
        gen.setSize(self._ppi_arc_widget.size())
        painter = QPainter(gen)
        self._ppi_arc_widget.render(painter)
        painter.end()
        self._status.showMessage(f"✓ PPI arc map exported: {path}")

    def _ppi_arc_map_export_tsv(self):
        """Export displayed arcs as a TSV file."""
        arcs = self._ppi_arc_widget._arcs
        if not arcs:
            QMessageBox.information(self, "Export TSV",
                "No interactions to export. Click Refresh first.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export PPI Map interactions as TSV", "", "TSV (*.tsv *.txt)")
        if not path:
            return
        headers = ["orf_a", "orf_b", "iptm", "pae_min_inter",
                   "cp_iptm_inter", "contact_frac_pct", "contact_region", "job_name"]
        rows = []
        for a in arcs:
            rows.append([
                a['name_a'], a['name_b'],
                f"{a['iptm']:.3f}" if a['iptm'] else "",
                f"{a['pae_min']:.3f}" if a['pae_min'] is not None else "",
                f"{a['cp_iptm']:.3f}" if a['cp_iptm'] is not None else "",
                f"{a['contact_frac']*100:.1f}" if a['contact_frac'] is not None else "",
                a.get('contact_region', ''),
                a.get('job_name', ''),
            ])
        try:
            with open(path, 'w', encoding='utf-8') as f:
                f.write('\t'.join(headers) + '\n')
                for r in rows:
                    f.write('\t'.join(r) + '\n')
            self._status.showMessage(f"✓ {len(rows)} interactions exported: {path}")
        except OSError as e:
            QMessageBox.warning(self, "Export TSV", str(e))

    def _create_hpc_server_tab(self):
        """Build the generic HPC Server tab with 4 sub-tabs:
        Connect | Submit | Monitor | Results."""
        w = QWidget()
        root = QVBoxLayout(w)
        root.setContentsMargins(4, 4, 4, 4)
        root.setSpacing(4)

        # Paramiko availability banner
        if not PARAMIKO_AVAILABLE:
            warn = QLabel(
                "⚠  paramiko not installed — SSH features disabled.\n"
                "   Install with:  pip install paramiko")
            warn.setStyleSheet(
                "background:#FFF3CD;color:#856404;padding:6px 10px;"
                "border-radius:4px;font-size:11px;")
            root.addWidget(warn)

        # Connection status bar (always visible)
        top_bar = QHBoxLayout()
        self._dv_status_icon = QLabel("⬤")
        self._dv_status_icon.setStyleSheet("color:#bbb;font-size:14px;")
        top_bar.addWidget(self._dv_status_icon)
        self._dv_status_lbl = QLabel("Not connected")
        self._dv_status_lbl.setStyleSheet("font-size:11px;color:#666;")
        top_bar.addWidget(self._dv_status_lbl)
        top_bar.addStretch()
        self._dv_disconnect_btn = QPushButton("Disconnect")
        self._dv_disconnect_btn.setEnabled(False)
        self._dv_disconnect_btn.clicked.connect(self._hpc_disconnect)
        self._dv_disconnect_btn.setFixedHeight(24)
        top_bar.addWidget(self._dv_disconnect_btn)
        root.addLayout(top_bar)

        # Sub-tabs (DetachableTabWidget — right-click or double-click
        # a tab to open it in its own resizable window)
        self._dv_tabs = DetachableTabWidget()
        self._dv_tabs.setTabPosition(QTabWidget.TabPosition.North
                                      if QT_VERSION == 6
                                      else QTabWidget.North)
        self._dv_tabs.setToolTip(
            "Tip: right-click any sub-tab (or double-click it) to open "
            "it in a separate window for more room.")
        root.addWidget(self._dv_tabs)

        self._dv_build_connect_tab()
        self._dv_build_submit_tab()
        self._dv_build_monitor_tab()
        self._dv_build_results_tab()

        self._tabs.addTab(w, "🖥 Submit AF3 via Server")

    # ── Sub-tab 1: Connect ──────────────────────────────────────

    def _dv_build_connect_tab(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(8, 8, 8, 8)
        lay.setSpacing(6)

        cg = QGroupBox("SSH Connection")
        cg_l = QGridLayout(cg)

        cg_l.addWidget(QLabel("User:"), 0, 0)
        self._dv_user = QLineEdit()
        self._dv_user.setPlaceholderText("e.g. my user name")
        cg_l.addWidget(self._dv_user, 0, 1)

        cg_l.addWidget(QLabel("Server:"), 1, 0)
        self._dv_host = QLineEdit()
        self._dv_host.setPlaceholderText("e.g. myserver.abc.def.gh")
        # Update groupbox title live as user types hostname
        self._dv_host.textChanged.connect(
            lambda h: cg.setTitle(f"SSH Connection — {h}" if h else "SSH Connection"))
        cg_l.addWidget(self._dv_host, 1, 1)

        cg_l.addWidget(QLabel("Port:"), 2, 0)
        self._dv_port = QSpinBox()
        self._dv_port.setRange(1, 65535)
        self._dv_port.setValue(22)
        cg_l.addWidget(self._dv_port, 2, 1)

        cg_l.addWidget(QLabel("Password:"), 3, 0)
        self._dv_pwd = QLineEdit()
        self._dv_pwd.setEchoMode(QLineEdit.EchoMode.Password
                                   if QT_VERSION == 6
                                   else QLineEdit.Password)
        self._dv_pwd.setPlaceholderText("••••••••")
        self._dv_pwd.returnPressed.connect(self._hpc_connect)
        cg_l.addWidget(self._dv_pwd, 3, 1)

        btn_row = QHBoxLayout()
        self._dv_connect_btn = QPushButton("🔗 Connect to Server")
        self._dv_connect_btn.setStyleSheet("font-weight:bold;")
        self._dv_connect_btn.clicked.connect(self._hpc_connect)
        self._dv_connect_btn.setEnabled(PARAMIKO_AVAILABLE)
        btn_row.addWidget(self._dv_connect_btn)
        btn_row.addStretch()
        cg_l.addLayout(btn_row, 4, 0, 1, 2)
        lay.addWidget(cg)

        pg = QGroupBox("Remote paths & submit command")
        pg_l = QGridLayout(pg)
        pg_l.addWidget(QLabel("Base predictions dir:"), 0, 0)
        self._dv_base_path = QLineEdit("~/af3_predictions")
        self._dv_base_path.setToolTip(
            "Remote directory where job sub-folders will be created.\n"
            "Tilde (~) is resolved automatically via SSH.")
        pg_l.addWidget(self._dv_base_path, 0, 1)
        pg_l.addWidget(QLabel("Submit command:"), 1, 0)
        self._dv_af3cmd = QLineEdit("af3_run")
        self._dv_af3cmd.setToolTip(
            "Command used to launch AF3 on the remote server.\n"
            "Examples:\n"
            "  af3_run          — wrapper that calls sbatch internally\n"
            "  sbatch run_af3.sh\n"
            "  python af3_submit.py\n"
            "  bash run.sh")
        pg_l.addWidget(self._dv_af3cmd, 1, 1)
        lay.addWidget(pg)

        # ── Scheduler selector ─────────────────────────────────
        sched_g = QGroupBox("Job scheduler")
        sched_l = QGridLayout(sched_g)
        sched_l.addWidget(QLabel("Scheduler:"), 0, 0)
        self._dv_scheduler = QComboBox()
        self._dv_scheduler.addItems([
            "SLURM  (squeue / scancel)",
            "PBS / Torque  (qstat / qdel)",
            "LSF  (bjobs / bkill)",
            "None  (direct execution)",
        ])
        self._dv_scheduler.setToolTip(
            "Select the workload manager on the remote server.\n"
            "  SLURM  — most HPC clusters (squeue, scancel, --slurm-* flags)\n"
            "  PBS    — older/government clusters (qstat, qdel)\n"
            "  LSF    — IBM clusters (bjobs, bkill)\n"
            "  None   — job runs directly, no scheduler monitoring")
        self._dv_scheduler.currentIndexChanged.connect(
            self._dv_on_scheduler_changed)
        sched_l.addWidget(self._dv_scheduler, 0, 1)
        self._dv_sched_info = QLabel(
            "Monitor: squeue -u <user>  |  Cancel: scancel <id>")
        self._dv_sched_info.setStyleSheet("font-size:10px;color:#888;")
        sched_l.addWidget(self._dv_sched_info, 1, 0, 1, 2)
        lay.addWidget(sched_g)

        # ── Environment / activation method ───────────────────
        env_g = QGroupBox("Environment activation")
        env_l = QGridLayout(env_g)
        env_l.addWidget(QLabel("Method:"), 0, 0)
        self._dv_env_method = QComboBox()
        self._dv_env_method.addItems([
            "module load  (Lmod / Env Modules)",
            "conda activate",
            "singularity exec",
            "source script",
            "None  (already on PATH)",
        ])
        self._dv_env_method.setToolTip(
            "How to activate AlphaFold3 on the remote server:\n"
            "  module load  — HPC clusters with Lmod (most common)\n"
            "  conda        — Anaconda/Miniconda environments\n"
            "  singularity  — container-based installations\n"
            "  source       — custom setup.sh / env.sh\n"
            "  None         — AF3 already on PATH, no activation needed")
        self._dv_env_method.currentIndexChanged.connect(
            self._dv_on_env_method_changed)
        env_l.addWidget(self._dv_env_method, 0, 1)
        env_l.addWidget(QLabel("Name / path:"), 1, 0)
        self._dv_module_cmd = QLineEdit("alphafold3")
        self._dv_module_cmd.setPlaceholderText(
            "module name, conda env, image path or script path")
        self._dv_module_cmd.setToolTip(
            "Argument for the activation method:\n"
            "  module load  → e.g. alphafold3  or  alphafold/3.0.0\n"
            "  conda        → e.g. af3_env\n"
            "  singularity  → e.g. /apps/af3.sif\n"
            "  source       → e.g. ~/setup_af3.sh\n"
            "  None         → leave blank")
        env_l.addWidget(self._dv_module_cmd, 1, 1)
        mod_row = QHBoxLayout()
        self._dv_mod_btn = QPushButton("⚙ Test activation")
        self._dv_mod_btn.setFixedHeight(26)
        self._dv_mod_btn.setEnabled(False)
        self._dv_mod_btn.clicked.connect(self._dv_load_module)
        self._dv_mod_btn.setToolTip(
            "Activate environment on the remote server and verify "
            "the submit command is reachable.")
        mod_row.addWidget(self._dv_mod_btn)
        self._dv_mod_status = QLabel("⬤  Environment: not checked")
        self._dv_mod_status.setStyleSheet(
            "font-size:11px; color:#888; padding-left:6px;")
        mod_row.addWidget(self._dv_mod_status)
        mod_row.addStretch()
        self._dv_mod_auto = QCheckBox("Auto-activate on connect")
        self._dv_mod_auto.setChecked(True)
        self._dv_mod_auto.setToolTip(
            "Automatically run the activation command right after SSH connection.")
        mod_row.addWidget(self._dv_mod_auto)
        env_l.addLayout(mod_row, 2, 0, 1, 2)
        lay.addWidget(env_g)
        self._dv_conn_log = QPlainTextEdit()
        self._dv_conn_log.setReadOnly(True)
        self._dv_conn_log.setFont(QFont('Courier New', 9))
        self._dv_conn_log.setMaximumHeight(160)
        self._dv_conn_log.setPlaceholderText("Connection log...")
        lay.addWidget(self._dv_conn_log)
        lay.addStretch()
        self._dv_tabs.addTab(w, "Connect")

    # ── Sub-tab 2: Submit ───────────────────────────────────────

    def _dv_build_submit_tab(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(8, 8, 8, 8)
        lay.setSpacing(6)

        # Source selector
        src_g = QGroupBox("JSON source")
        src_l = QVBoxLayout(src_g)
        self._dv_src_session = QPushButton("📋 Load jobs from current session  (one JSON per job)")
        self._dv_src_session.clicked.connect(self._dv_load_from_session)
        # ── Grouped-batch row: [Load Grouped JSON | ⚙ AF3 JSON Group Config] ──
        batch_row_w  = QWidget()
        batch_row_l  = QHBoxLayout(batch_row_w)
        batch_row_l.setContentsMargins(0, 0, 0, 0)
        batch_row_l.setSpacing(2)

        self._dv_src_batch = QPushButton("📦 Load Grouped JSON")
        self._dv_src_batch.clicked.connect(self._dv_load_from_session_batch)
        self._dv_src_batch.setToolTip(
            "Groups all AF3 jobs from the current session into partitioned\n"
            "batch JSON files (size configurable via ⚙ AF3 JSON Group Config).\n"
            "Each partition is submitted as an independent SLURM job to\n"
            "prevent GPU out-of-memory errors (anti-OOM).\n\n"
            "Current group size: configurable — click ⚙ AF3 JSON Group Config.")
        self._dv_src_batch.setFixedHeight(28)

        self._dv_batch_cfg_btn = QPushButton("⚙ AF3 JSON Group Config")
        self._dv_batch_cfg_btn.setToolTip(
            "Configure the number of AF3 jobs per partition batch.\n"
            "Default: 50 jobs / batch (anti-OOM safe).\n"
            "Increase for fast servers with ample RAM;\n"
            "decrease for large proteins or limited GPU memory.")
        self._dv_batch_cfg_btn.clicked.connect(self._dv_show_batch_group_config)
        self._dv_batch_cfg_btn.setFixedHeight(28)
        self._dv_batch_cfg_btn.setStyleSheet(
            "QPushButton { color: #0066cc; font-weight: bold; }"
            "QPushButton:hover { background: #e8f0ff; }")

        batch_row_l.addWidget(self._dv_src_batch,      stretch=1)
        batch_row_l.addWidget(self._dv_batch_cfg_btn,  stretch=1)

        self._dv_src_file = QPushButton("📂 Load JSON file(s) from disk")
        self._dv_src_file.clicked.connect(self._dv_load_from_files)
        for b in (self._dv_src_session, self._dv_src_file):
            b.setFixedHeight(28)
            src_l.addWidget(b)
        src_l.insertWidget(1, batch_row_w)   # insert after session button
        lay.addWidget(src_g)

        # Job name + remote dir
        nm_g = QGroupBox("Job naming")
        nm_l = QGridLayout(nm_g)
        nm_l.addWidget(QLabel("Job name prefix:"), 0, 0)
        self._dv_job_prefix = QLineEdit()
        self._dv_job_prefix.setPlaceholderText("e.g. Xeuvesica_contig1_XVIPCD_up_down")
        self._dv_job_prefix.textChanged.connect(self._dv_refresh_cmd_preview)
        nm_l.addWidget(self._dv_job_prefix, 0, 1)

        self._dv_ts_check = QCheckBox("Add timestamp  (Prefix_YYYYMMDD_HHMMSS)")
        self._dv_ts_check.setChecked(True)
        self._dv_ts_check.setToolTip(
            "Appends a timestamp to --job-name on each submission.\n"
            "Prevents 'Job directory already exists' errors when\n"
            "resubmitting the same experiment.")
        self._dv_ts_check.stateChanged.connect(
            lambda _: self._dv_refresh_cmd_preview())
        nm_l.addWidget(self._dv_ts_check, 1, 0, 1, 2)

        nm_l.addWidget(QLabel("Remote sub-dir:"), 2, 0)
        self._dv_remote_dir = QLineEdit()
        self._dv_remote_dir.setPlaceholderText(
            f"auto  ({datetime.now().strftime('%Y-%m-%d')})")
        nm_l.addWidget(self._dv_remote_dir, 2, 1)
        lay.addWidget(nm_g)

        # Jobs preview table
        self._dv_submit_table = QTableWidget()
        self._dv_submit_table.setColumnCount(5)
        self._dv_submit_table.setHorizontalHeaderLabels(
            ['JSON / Name', 'Chains', 'Residues', 'Size (KB)', 'Status'])
        self._dv_submit_table.setSelectionBehavior(SelectRows)
        self._dv_submit_table.setEditTriggers(
            QAbstractItemView.EditTrigger.NoEditTriggers
            if QT_VERSION == 6 else QAbstractItemView.NoEditTriggers)
        self._dv_submit_table.horizontalHeader().setStretchLastSection(True)
        self._dv_submit_table.setAlternatingRowColors(True)
        self._dv_submit_table.setMaximumHeight(160)
        for i, w_ in enumerate([200, 55, 65, 65, 70]):
            self._dv_submit_table.setColumnWidth(i, w_)
        lay.addWidget(self._dv_submit_table)

        # Summary + action buttons
        self._dv_submit_summary = QLabel("0 jobs loaded")
        self._dv_submit_summary.setStyleSheet("font-size:11px;color:#666;")
        lay.addWidget(self._dv_submit_summary)

        act_row = QHBoxLayout()
        self._dv_clear_btn = QPushButton("🗑 Clear list")
        self._dv_clear_btn.clicked.connect(self._dv_clear_submit_list)
        self._dv_upload_btn = QPushButton("⬆ Upload only (SFTP)")
        self._dv_upload_btn.clicked.connect(self._dv_upload_only)
        self._dv_run_btn = QPushButton("🚀 Upload + Submit all")
        self._dv_run_btn.setStyleSheet("font-weight:bold;")
        self._dv_run_btn.clicked.connect(self._dv_upload_and_submit)
        self._dv_dl_jsons_btn = QPushButton("📥 Download JSONs")
        self._dv_dl_jsons_btn.setToolTip(
            "Save all staged partition JSON files to a local folder.\n"
            "Useful to inspect, archive, or manually transfer the files\n"
            "before submitting to the cluster.")
        self._dv_dl_jsons_btn.clicked.connect(self._dv_download_staged_jsons)
        for b in (self._dv_clear_btn, self._dv_upload_btn,
                  self._dv_run_btn, self._dv_dl_jsons_btn):
            b.setEnabled(False)
            act_row.addWidget(b)
        lay.addLayout(act_row)

        # ── AF3 JSON modifications (modelSeeds + templates) ────────
        # Small group containing options that affect the JSON CONTENT
        # (not the command).  Kept separate so they're not duplicated
        # in every command-template profile.
        json_g = QGroupBox("⚗ AF3 JSON modifications  (modelSeeds, templates)")
        json_g.setCheckable(True)
        json_g.setChecked(False)
        json_g.setToolTip(
            "These options modify the AF3 batch JSON before upload:\n"
            "  • Model seeds  → injected into 'modelSeeds[]'\n"
            "  • Use templates → strip 'templates' from each chain\n"
            "They are independent from the command-template profile.")
        json_l = QGridLayout(json_g)
        json_l.setSpacing(5)

        json_l.addWidget(QLabel("Model seeds (n):"), 0, 0)
        self._dv_af3_seeds = QSpinBox()
        self._dv_af3_seeds.setRange(1, 10)
        self._dv_af3_seeds.setValue(1)
        self._dv_af3_seeds.setToolTip(
            "Number of random seeds injected into JSON 'modelSeeds[]'.\n"
            "More seeds = ensemble of independent models → pick best ipTM.\n"
            "1 = fastest;  5 = recommended for final confident predictions.")
        self._dv_af3_seeds.valueChanged.connect(self._dv_refresh_cmd_preview)
        json_l.addWidget(self._dv_af3_seeds, 0, 1)

        self._dv_af3_use_templates = QCheckBox("Use templates")
        self._dv_af3_use_templates.setChecked(True)
        self._dv_af3_use_templates.setToolTip(
            "If unchecked, the 'templates' array is wiped from every\n"
            "proteinChain in the batch JSON before upload — useful for\n"
            "ab-initio screening or to avoid template bias.")
        json_l.addWidget(self._dv_af3_use_templates, 0, 2, 1, 2)

        lay.addWidget(json_g)

        # ── AF3 Server Submission — Terminal & Profiles ──────────
        # Terminal-style command editor where the user types the AF3
        # command(s) for their own server.  No server-specific defaults
        # are shipped — the user saves their own commands as named
        # profiles in ~/.ppigfinder/af3_server_profiles.json
        self._dv_af3_widget = FlexibleAF3SubmitWidget(parent=w)
        self._dv_af3_widget.commandChanged.connect(self._dv_refresh_cmd_preview)
        lay.addWidget(self._dv_af3_widget)

        # Command preview
        cmd_g = QGroupBox("Command preview (editable)")
        cmd_l = QVBoxLayout(cmd_g)
        self._dv_cmd_preview = QPlainTextEdit()
        self._dv_cmd_preview.setFont(QFont('Courier New', 9))
        self._dv_cmd_preview.setStyleSheet(
            "background:#1e1e1e;color:#d4d4d4;")
        self._dv_cmd_preview.setMaximumHeight(90)
        cmd_l.addWidget(self._dv_cmd_preview)
        lay.addWidget(cmd_g)

        self._dv_submit_log = QPlainTextEdit()
        self._dv_submit_log.setReadOnly(True)
        self._dv_submit_log.setFont(QFont('Courier New', 9))
        self._dv_submit_log.setMaximumHeight(90)
        self._dv_submit_log.setPlaceholderText("Submit log...")
        lay.addWidget(self._dv_submit_log)

        self._dv_tabs.addTab(w, "Submit jobs")

    # ── Sub-tab 3: Monitor ──────────────────────────────────────

    def _dv_build_monitor_tab(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(8, 8, 8, 8)
        lay.setSpacing(6)

        # ── Row 1: auto-refresh + refresh button + timestamp ──
        ctrl = QHBoxLayout()
        ctrl.addWidget(QLabel("Auto-refresh:"))
        self._dv_poll_interval = QComboBox()
        self._dv_poll_interval.addItems(
            ["Off", "15 s", "30 s", "60 s", "2 min", "5 min"])
        self._dv_poll_interval.setCurrentText("30 s")
        self._dv_poll_interval.setToolTip(
            "How often to poll the server queue automatically.")
        self._dv_poll_interval.currentTextChanged.connect(
            self._dv_update_poll_timer)
        ctrl.addWidget(self._dv_poll_interval)

        btn_ref = QPushButton("🔄 Refresh now")
        btn_ref.setToolTip("Run squeue immediately")
        btn_ref.clicked.connect(self._hpc_poll_queue)
        ctrl.addWidget(btn_ref)

        # ── Show ALL server jobs (not only jobs from this session) ──
        self._dv_mon_show_all = QCheckBox("Show all server jobs")
        self._dv_mon_show_all.setChecked(False)
        self._dv_mon_show_all.setToolTip(
            "When checked: display every job from squeue for this user,\n"
            "not only the ones submitted in this session.\n"
            "Useful to monitor pre-existing or parallel jobs.")
        ctrl.addWidget(self._dv_mon_show_all)

        ctrl.addStretch()
        self._dv_queue_lbl = QLabel("squeue — —:—:—")
        self._dv_queue_lbl.setStyleSheet("font-size:11px;color:#666;")
        ctrl.addWidget(self._dv_queue_lbl)
        lay.addLayout(ctrl)

        # ── Row 2: filters (status + partition) + cancel button ──
        flt = QHBoxLayout()
        flt.addWidget(QLabel("Filter status:"))
        self._dv_mon_status_filter = QComboBox()
        self._dv_mon_status_filter.addItems(
            ["All", "RUNNING", "PENDING", "COMPLETED", "FAILED", "CANCELLED"])
        self._dv_mon_status_filter.setFixedWidth(110)
        self._dv_mon_status_filter.setToolTip(
            "Show only jobs in the selected SLURM state.")
        self._dv_mon_status_filter.currentTextChanged.connect(
            self._dv_refresh_monitor_table)
        flt.addWidget(self._dv_mon_status_filter)

        flt.addWidget(QLabel("Partition:"))
        self._dv_mon_partition_filter = QLineEdit()
        self._dv_mon_partition_filter.setPlaceholderText("all")
        self._dv_mon_partition_filter.setFixedWidth(80)
        self._dv_mon_partition_filter.setToolTip(
            "Filter by SLURM partition name (leave blank for all).")
        self._dv_mon_partition_filter.textChanged.connect(
            self._dv_refresh_monitor_table)
        flt.addWidget(self._dv_mon_partition_filter)

        self._dv_cancel_btn = QPushButton("⛔ Cancel selected job")
        self._dv_cancel_btn.setFixedHeight(24)
        self._dv_cancel_btn.setToolTip(
            "Run scancel on the SLURM ID of the selected row.\n"
            "Only works for RUNNING or PENDING jobs.")
        self._dv_cancel_btn.clicked.connect(self._dv_cancel_selected_job)
        flt.addWidget(self._dv_cancel_btn)

        btn_export_hist = QPushButton("💾 Export history")
        btn_export_hist.setFixedHeight(24)
        btn_export_hist.setToolTip(
            "Save the current job list (this session) to a JSON file.")
        btn_export_hist.clicked.connect(self._dv_export_job_history)
        flt.addWidget(btn_export_hist)

        flt.addStretch()
        lay.addLayout(flt)

        # ── Monitor table ──────────────────────────────────────
        self._dv_monitor_table = QTableWidget()
        self._dv_monitor_table.setColumnCount(7)
        self._dv_monitor_table.setHorizontalHeaderLabels(
            ['Job name', 'SLURM ID', 'Status', 'Partition', 'Node', 'Time', 'Action'])
        self._dv_monitor_table.setSelectionBehavior(SelectRows)
        self._dv_monitor_table.setEditTriggers(
            QAbstractItemView.EditTrigger.NoEditTriggers
            if QT_VERSION == 6 else QAbstractItemView.NoEditTriggers)
        self._dv_monitor_table.horizontalHeader().setStretchLastSection(False)
        self._dv_monitor_table.horizontalHeader().setSortIndicatorShown(True)
        self._dv_monitor_table.setSortingEnabled(True)
        self._dv_monitor_table.setAlternatingRowColors(True)
        for i, w_ in enumerate([195, 72, 90, 75, 75, 65, 100]):
            self._dv_monitor_table.setColumnWidth(i, w_)
        lay.addWidget(self._dv_monitor_table)

        # ── Raw squeue output (resizable) ─────────────────────
        squeue_g = QGroupBox("Raw squeue output")
        sq_l = QVBoxLayout(squeue_g)
        self._dv_squeue_raw = QPlainTextEdit()
        self._dv_squeue_raw.setReadOnly(True)
        self._dv_squeue_raw.setFont(QFont('Courier New', 9))
        self._dv_squeue_raw.setMinimumHeight(80)   # resizable — no max
        sq_l.addWidget(self._dv_squeue_raw)
        lay.addWidget(squeue_g)

        self._dv_tabs.addTab(w, "Monitor")

    # ── Monitor: cancel selected job ───────────────────────────
    def _dv_cancel_selected_job(self):
        """Run scancel <SLURM_ID> on the selected row."""
        if not self._ssh_client:
            QMessageBox.warning(self, "Server", "Not connected.")
            return
        rows = self._dv_monitor_table.selectedItems()
        if not rows:
            QMessageBox.information(self, "Cancel job",
                "Select a job row first.")
            return
        row = self._dv_monitor_table.currentRow()
        slurm_id = self._dv_monitor_table.item(row, 1)
        if not slurm_id or not slurm_id.text().strip().isdigit():
            QMessageBox.warning(self, "Cancel job",
                "No valid SLURM ID for the selected row.")
            return
        sid = slurm_id.text().strip()
        job_name_item = self._dv_monitor_table.item(row, 0)
        job_name = job_name_item.text() if job_name_item else sid
        ans = QMessageBox.question(self, "Cancel job",
            f"Cancel job  {job_name}  (SLURM {sid})?")
        if ans != (QMessageBox.StandardButton.Yes if QT_VERSION == 6
                   else QMessageBox.Yes):
            return

        def _do_cancel():
            out, err, rc = self._dv_ssh_exec(
                self._dv_scheduler_cmds()["cancel"].format(job_id=sid), timeout=10)
            return out, err, rc

        w = AnalysisWorker(_do_cancel)
        w.finished.connect(lambda r: (
            self._dv_log(
                f"scancel {sid} → rc={r[2]}  {(r[0]+r[1]).strip()}", 'conn'),
            self._hpc_poll_queue()
        ))
        w.error.connect(lambda e: self._dv_log(f"scancel error: {e}", 'conn'))
        self._hpc_workers.append(w)
        w.start()

    # ── Monitor: export job history ────────────────────────────
    def _dv_export_job_history(self):
        """Save the HPC job list to a JSON file."""
        path, _ = QFileDialog.getSaveFileName(
            self, "Export job history", "",
            "JSON files (*.json);;All files (*)")
        if not path:
            return
        try:
            with open(path, 'w', encoding='utf-8') as f:
                json.dump(self._hpc_jobs, f, indent=2, ensure_ascii=False)
            self._status.showMessage(
                f"✓ Job history exported — {path}")
        except Exception as ex:
            QMessageBox.critical(self, "Export error", str(ex))

    # ── Sub-tab 4: Results ──────────────────────────────────────

    def _dv_build_results_tab(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(8, 8, 8, 8)
        lay.setSpacing(6)

        # Remote path selector
        rp_g = QGroupBox("Remote output path")
        rp_l = QGridLayout(rp_g)
        rp_l.addWidget(QLabel("Remote dir:"), 0, 0)
        self._dv_remote_output = QLineEdit()
        self._dv_remote_output.setPlaceholderText(
            "~/af3_predictions/.../job_name/output/")
        self._dv_remote_output.setFont(QFont('Courier New', 9))
        rp_l.addWidget(self._dv_remote_output, 0, 1)
        btn_browse = QPushButton("Browse...")
        btn_browse.clicked.connect(self._dv_browse_remote_output)
        rp_l.addWidget(btn_browse, 0, 2)
        lay.addWidget(rp_g)

        # Remote file listing
        self._dv_remote_file_table = QTableWidget()
        self._dv_remote_file_table.setColumnCount(4)
        self._dv_remote_file_table.setHorizontalHeaderLabels(
            ['Name', 'Type', 'Size', 'Modified'])
        self._dv_remote_file_table.setSelectionBehavior(SelectRows)
        self._dv_remote_file_table.setEditTriggers(
            QAbstractItemView.EditTrigger.NoEditTriggers
            if QT_VERSION == 6 else QAbstractItemView.NoEditTriggers)
        self._dv_remote_file_table.horizontalHeader().setStretchLastSection(True)
        self._dv_remote_file_table.setAlternatingRowColors(True)
        self._dv_remote_file_table.setMaximumHeight(160)
        for i, cw in enumerate([200, 60, 80, 120]):
            self._dv_remote_file_table.setColumnWidth(i, cw)
        lay.addWidget(self._dv_remote_file_table)

        # Local destination
        lp_g = QGroupBox("Local destination")
        lp_l = QGridLayout(lp_g)
        lp_l.addWidget(QLabel("Download to:"), 0, 0)
        self._dv_local_dest = QLineEdit()
        self._dv_local_dest.setPlaceholderText("/local/path/results/")
        lp_l.addWidget(self._dv_local_dest, 0, 1)
        btn_local = QPushButton("Browse...")
        btn_local.clicked.connect(self._dv_choose_local_dest)
        lp_l.addWidget(btn_local, 0, 2)
        self._dv_auto_import = QCheckBox(
            "Auto-import ipTM/pLDDT into Ranking tab after download")
        self._dv_auto_import.setChecked(True)
        lp_l.addWidget(self._dv_auto_import, 1, 0, 1, 3)
        lay.addWidget(lp_g)

        act_row = QHBoxLayout()
        btn_list = QPushButton("📋 List remote files")
        btn_list.clicked.connect(self._dv_list_remote_output)
        self._dv_dl_btn = QPushButton("⬇ Download selected + Import")
        self._dv_dl_btn.setStyleSheet("font-weight:bold;")
        self._dv_dl_btn.clicked.connect(self._dv_download_results)
        act_row.addWidget(btn_list)
        act_row.addWidget(self._dv_dl_btn)
        act_row.addStretch()
        lay.addLayout(act_row)

        self._dv_dl_progress = QLabel("")
        self._dv_dl_progress.setStyleSheet("font-size:11px;color:#666;")
        lay.addWidget(self._dv_dl_progress)

        self._dv_tabs.addTab(w, "Results")

    # ── SSH helpers ─────────────────────────────────────────────

    def _dv_log(self, msg: str, pane: str = 'conn'):
        """Append a message to the connection or submit log pane."""
        ts = datetime.now().strftime('%H:%M:%S')
        line = f"[{ts}] {msg}"
        if pane == 'submit':
            self._dv_submit_log.appendPlainText(line)
        else:
            self._dv_conn_log.appendPlainText(line)

    def _dv_set_connected(self, connected: bool, label: str = ''):
        """Update the connection status widgets."""
        if connected:
            self._dv_status_icon.setStyleSheet("color:#2e7d32;font-size:14px;")
            self._dv_status_lbl.setText(f"Connected — {label}")
            self._dv_disconnect_btn.setEnabled(True)
            self._dv_connect_btn.setEnabled(False)
            self._dv_mod_btn.setEnabled(True)
            for b in (self._dv_upload_btn, self._dv_run_btn):
                b.setEnabled(True)
            # Start poll timer if not off
            self._dv_update_poll_timer(self._dv_poll_interval.currentText())
        else:
            self._dv_status_icon.setStyleSheet("color:#bbb;font-size:14px;")
            self._dv_status_lbl.setText("Not connected")
            self._dv_disconnect_btn.setEnabled(False)
            self._dv_connect_btn.setEnabled(PARAMIKO_AVAILABLE)
            self._dv_mod_btn.setEnabled(False)
            self._dv_mod_status.setText("⬤  Module: not checked")
            self._dv_mod_status.setStyleSheet(
                "font-size:11px; color:#888; padding-left:6px;")
            self._dv_module_loaded = False
            for b in (self._dv_upload_btn, self._dv_run_btn):
                b.setEnabled(False)
            self._hpc_poll_timer.stop()

    def _hpc_connect(self):
        """Open SSH connection to the configured HPC server."""
        if not PARAMIKO_AVAILABLE:
            QMessageBox.warning(self, "Server",
                "paramiko not installed.\n\npip install paramiko")
            return
        host = self._dv_host.text().strip()
        user = self._dv_user.text().strip()
        pwd  = self._dv_pwd.text()
        port = self._dv_port.value()
        if not host or not user:
            QMessageBox.warning(self, "Server", "Enter user and server hostname.")
            return
        self._dv_connect_btn.setEnabled(False)
        self._dv_connect_btn.setText("Connecting...")
        self._dv_log(f"Connecting to {user}@{host}:{port} ...")

        def _do_connect():
            client = paramiko.SSHClient()
            client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            client.connect(host, port=port, username=user,
                           password=pwd, timeout=15,
                           allow_agent=False, look_for_keys=False)
            sftp = client.open_sftp()
            return client, sftp

        w = AnalysisWorker(_do_connect)
        w.finished.connect(self._dv_on_connected)
        w.error.connect(self._dv_on_connect_error)
        self._hpc_workers.append(w)
        w.start()

    def _dv_on_connected(self, result):
        client, sftp = result
        self._ssh_client = client
        self._sftp_client = sftp
        # Keep TCP connection alive every 30 s to prevent idle timeout
        try:
            client.get_transport().set_keepalive(30)
        except Exception:
            pass
        label = f"{self._dv_user.text()}@{self._dv_host.text()}"
        self._dv_log(f"SSH session established — {label}")
        self._dv_set_connected(True, label)
        self._dv_connect_btn.setText("🔗 Connect to Server")
        self._status.showMessage(f"✓ Connected to server ({label})")
        # Auto-load module if checkbox is ticked
        if self._dv_mod_auto.isChecked():
            self._dv_load_module()

    def _dv_on_connect_error(self, msg):
        self._dv_log(f"ERROR: {msg}")
        self._dv_set_connected(False)
        self._dv_connect_btn.setText("🔗 Connect to Server")
        self._dv_connect_btn.setEnabled(PARAMIKO_AVAILABLE)
        QMessageBox.critical(self, "Server — connection failed", msg)

    # ── Scheduler / environment helper callbacks ───────────────

    def _dv_on_scheduler_changed(self, idx: int):
        """Update scheduler info label; disable timer when no scheduler."""
        info = {
            0: "Monitor: squeue -u <user>  |  Cancel: scancel <id>",
            1: "Monitor: qstat -u <user>   |  Cancel: qdel <id>",
            2: "Monitor: bjobs -u <user>   |  Cancel: bkill <id>",
            3: "No scheduler — direct execution (no job monitoring)",
        }
        if hasattr(self, '_dv_sched_info'):
            self._dv_sched_info.setText(info.get(idx, ''))
        if idx == 3:
            self._hpc_poll_timer.stop()
            if hasattr(self, '_dv_poll_interval'):
                self._dv_poll_interval.setEnabled(False)
        else:
            if hasattr(self, '_dv_poll_interval'):
                self._dv_poll_interval.setEnabled(True)

    def _dv_on_env_method_changed(self, idx: int):
        """Update Name/path placeholder to match the selected activation method."""
        placeholders = {
            0: "module name  e.g. alphafold3  or  alphafold/3.0.0",
            1: "conda env name  e.g. af3_env",
            2: "image path  e.g. /apps/af3.sif  or  docker://af3",
            3: "script path  e.g. ~/setup_af3.sh",
            4: "",
        }
        if hasattr(self, '_dv_module_cmd'):
            self._dv_module_cmd.setPlaceholderText(placeholders.get(idx, ''))
            self._dv_module_cmd.setEnabled(idx != 4)

    def _dv_build_activation_prefix(self) -> str:
        """Return the shell activation snippet for the selected env method.

        Used to prepend to run_cmd, e.g.:
          'module load alphafold3 2>/dev/null; '
          'conda activate af3_env && '
          ''   (None — already on PATH)
        """
        if not hasattr(self, '_dv_env_method'):
            name = self._dv_module_cmd.text().strip() if hasattr(
                self, '_dv_module_cmd') else ''
            return f"module load {name} 2>/dev/null; " if name else ""
        idx  = self._dv_env_method.currentIndex()
        name = self._dv_module_cmd.text().strip() if hasattr(
            self, '_dv_module_cmd') else ''
        if   idx == 0: return f"module load {name} 2>/dev/null; " if name else ""
        elif idx == 1: return f"conda activate {name} && "        if name else ""
        elif idx == 2: return ""   # singularity wraps the cmd itself
        elif idx == 3: return f"source {name} 2>/dev/null; "      if name else ""
        else:          return ""   # None

    def _dv_scheduler_cmds(self) -> dict:
        """Return poll/cancel commands for the selected job scheduler."""
        user = self._dv_user.text().strip() if hasattr(self, '_dv_user') else ''
        idx  = self._dv_scheduler.currentIndex() if hasattr(
            self, '_dv_scheduler') else 0
        SLURM = {
            'poll':   (f"squeue -u {user} "
                       f"--format='%.18i %.9P %.50j %.8u %.2t %.10M %.6D' "
                       f"2>/dev/null"),
            'cancel': "scancel {job_id}",
            'id_re':  r'batch job\s+(\d+)',
            'row_ok': lambda p: len(p) >= 5 and p[0].isdigit(),
        }
        PBS = {
            'poll':   f"qstat -u {user} 2>/dev/null",
            'cancel': "qdel {job_id}",
            'id_re':  r'(\d+)\.\S+',
            'row_ok': lambda p: len(p) >= 5 and p[0][0].isdigit(),
        }
        LSF = {
            'poll':   f"bjobs -u {user} 2>/dev/null",
            'cancel': "bkill {job_id}",
            'id_re':  r'Job <(\d+)>',
            'row_ok': lambda p: len(p) >= 5 and p[0].isdigit(),
        }
        NONE = {
            'poll':   'echo "No scheduler configured"',
            'cancel': 'kill {job_id}',
            'id_re':  r'(\d+)',
            'row_ok': lambda p: False,
        }
        return [SLURM, PBS, LSF, NONE][min(idx, 3)]

    def _dv_load_module(self):
        """Activate the configured environment on the remote server
        and verify the submit command is reachable.  Updates status badge."""
        if not self._ssh_client:
            return
        env_idx = self._dv_env_method.currentIndex() \
            if hasattr(self, "_dv_env_method") else 0
        if env_idx == 4:  # None — already on PATH
            self._dv_log("Environment method: None — skipping activation check.")
            return

        act_prefix = self._dv_build_activation_prefix()
        af3cmd     = self._dv_af3cmd.text().strip() or "af3_run"

        self._dv_mod_status.setText("⬤  Environment: checking…")
        self._dv_mod_status.setStyleSheet(
            "font-size:11px; color:#888; padding-left:6px;")
        self._dv_mod_btn.setEnabled(False)
        self._dv_log(f"Activating: {act_prefix.strip() or '(no prefix)'} …")

        def _do_module():
            check_cmd = (
                f'bash -lc "'
                f'{act_prefix}'
                f'echo __ACT_RC__$?; '
                f'which {af3cmd} 2>/dev/null && echo __AF3_FOUND__; '
                f'{af3cmd} --version 2>&1 | head -3'
                f'"'
            )
            out, err, rc = self._dv_ssh_exec(check_cmd, timeout=30)
            return out, err, rc

        w = AnalysisWorker(_do_module)
        w.finished.connect(self._dv_on_module_done)
        w.error.connect(self._dv_on_module_error)
        self._hpc_workers.append(w)
        w.start()

    def _dv_on_module_done(self, result):
        """Parse module load output and update status badge."""
        out, err, rc = result
        combined = out + err

        # Extract module load return code embedded in output
        m_rc = re.search(r'__ACT_RC__(\d+)', combined)
        module_rc = int(m_rc.group(1)) if m_rc else 1

        af3_found    = '__AF3_FOUND__' in combined
        # Try to extract version line (first non-empty line after module load)
        ver_lines = [ln.strip() for ln in combined.splitlines()
                     if ln.strip()
                     and not ln.startswith('__')
                     and 'module load' not in ln.lower()
                     and '__ACT_RC__' not in ln]
        version_str = ver_lines[0] if ver_lines else ''

        mod_name = self._dv_module_cmd.text().strip()
        self._dv_mod_btn.setEnabled(True)

        if module_rc == 0 and af3_found:
            self._dv_module_loaded = True
            self._dv_mod_status.setText(
                f"✅  Module loaded: {mod_name}"
                + (f"  ({version_str[:40]})" if version_str else ''))
            self._dv_mod_status.setStyleSheet(
                "font-size:11px; color:#2e7d32; "
                "font-weight:bold; padding-left:6px;")
            self._dv_log(f"✅ Module '{mod_name}' loaded — af3_run found.")
            if version_str:
                self._dv_log(f"   Version: {version_str[:80]}")
            self._status.showMessage(
                f"✓ Server: module '{mod_name}' ready")
        elif module_rc == 0 and not af3_found:
            # Module loaded but af3_run binary not on PATH
            self._dv_module_loaded = False
            self._dv_mod_status.setText(
                "⚠  Module loaded but af3_run not found in PATH")
            self._dv_mod_status.setStyleSheet(
                "font-size:11px; color:#e65100; "
                "font-weight:bold; padding-left:6px;")
            self._dv_log(
                f"⚠ Module '{mod_name}' loaded (rc=0) but "
                f"'{self._dv_af3cmd.text().strip() or 'af3_run'}' "
                f"not found in PATH.")
            self._dv_log(
                "  Check 'af3_run command' field or module name.")
        else:
            self._dv_module_loaded = False
            # Show first error line from server
            err_lines = [ln.strip() for ln in combined.splitlines()
                         if ln.strip() and '__ACT_RC__' not in ln]
            err_preview = err_lines[0][:60] if err_lines else f'rc={module_rc}'
            self._dv_mod_status.setText(
                f"❌  Module load failed: {err_preview}")
            self._dv_mod_status.setStyleSheet(
                "font-size:11px; color:#c62828; "
                "font-weight:bold; padding-left:6px;")
            self._dv_log(
                f"❌ module load {mod_name} failed (rc={module_rc}):")
            for ln in err_lines[:5]:
                self._dv_log(f"   {ln}")

    def _dv_on_module_error(self, msg):
        """SSH exec error while loading module."""
        self._dv_mod_btn.setEnabled(True)
        self._dv_module_loaded = False
        self._dv_mod_status.setText(f"❌  Module check error: {msg[:50]}")
        self._dv_mod_status.setStyleSheet(
            "font-size:11px; color:#c62828; font-weight:bold; padding-left:6px;")
        self._dv_log(f"❌ Module load error: {msg}")

    def _hpc_disconnect(self):
        """Close SSH/SFTP connections cleanly."""
        self._hpc_poll_timer.stop()
        try:
            if self._sftp_client:
                self._sftp_client.close()
            if self._ssh_client:
                self._ssh_client.close()
        except Exception:
            pass
        self._ssh_client = None
        self._sftp_client = None
        self._dv_set_connected(False)
        self._dv_log("Disconnected.")
        self._status.showMessage("Server connection closed.")

    def _dv_ssh_exec(self, cmd: str, timeout: int = 30):
        """Execute a command on the remote server and return (stdout, stderr, exit_code)."""
        if not self._ssh_client:
            raise RuntimeError("Not connected to server.")
        _, stdout, stderr = self._ssh_client.exec_command(cmd, timeout=timeout)
        exit_code = stdout.channel.recv_exit_status()
        return stdout.read().decode('utf-8', errors='replace'), \
               stderr.read().decode('utf-8', errors='replace'), exit_code

    # ── Submit helpers ──────────────────────────────────────────

    def _dv_load_from_session(self):
        """Populate the submit table from AF3 jobs in the current session."""
        if not self.af3_jobs:
            QMessageBox.information(self, "Server",
                "No AF3 jobs in session.\n"
                "Generate jobs in the AlphaFold tab first.")
            return
        self._dv_submit_table.setRowCount(0)
        for j in self.af3_jobs:
            row = self._dv_submit_table.rowCount()
            self._dv_submit_table.insertRow(row)
            n_chains = len(j.get('sequences', []))
            res = j.get('total_residues', 0)
            size_kb = res * 0.2 / 1024
            for col, val in enumerate([
                    j['name'], str(n_chains), str(res),
                    f"{size_kb:.1f}", "Pending"]):
                self._dv_submit_table.setItem(row, col, QTableWidgetItem(val))
        self._dv_submit_summary.setText(
            f"{len(self.af3_jobs)} jobs from session  (one JSON per job)")
        self._dv_pending_jobs = list(self.af3_jobs)
        for b in (self._dv_upload_btn, self._dv_run_btn,
                  self._dv_clear_btn, self._dv_dl_jsons_btn):
            b.setEnabled(True)
        self._dv_refresh_cmd_preview()

    # ── Anti-OOM partition size ────────────────────────────────────
    _AF3_PARTITION_SIZE = 50   # class-level default; overridden by instance var

    def _dv_show_batch_group_config(self):
        """Show a compact dialog to configure the AF3 batch group (partition) size."""
        current = getattr(self, '_dv_partition_size', self._AF3_PARTITION_SIZE)

        dlg = QDialog(self)
        dlg.setWindowTitle("AF3 JSON Group Config")
        dlg.setFixedSize(400, 230)
        lay = QVBoxLayout(dlg)
        lay.setSpacing(10)

        hdr = QLabel("<b>Configure AF3 batch partition size</b>")
        hdr.setAlignment(Qt.AlignmentFlag.AlignCenter
                         if QT_VERSION == 6 else Qt.AlignCenter)
        lay.addWidget(hdr)

        desc = QLabel(
            "Sets how many AF3 jobs are packed into each partition JSON.\n"
            "Each partition becomes one independent SLURM job on the server.\n"
            "Smaller values = safer for limited GPU RAM (anti-OOM).\n"
            "Larger values = fewer queue submissions, faster throughput.")
        desc.setWordWrap(True)
        desc.setStyleSheet("color:#555; font-size:11px;")
        lay.addWidget(desc)

        spin_row = QHBoxLayout()
        spin_row.addWidget(QLabel("Jobs per partition batch:"))
        spin = QSpinBox()
        spin.setRange(1, 9999)
        spin.setValue(current)
        spin.setFixedWidth(80)
        spin.setToolTip(
            "Recommended values:\n"
            "  10-20  -> very large proteins / low GPU RAM\n"
            "  50     -> default (safe for most clusters)\n"
            "  100+   -> fast clusters with ample RAM")
        spin_row.addWidget(spin)
        spin_row.addStretch()
        lay.addLayout(spin_row)

        preset_row = QHBoxLayout()
        preset_row.addWidget(QLabel("Quick presets:"))
        for label, val in [("10", 10), ("25", 25), ("50", 50),
                            ("100", 100), ("200", 200)]:
            pb = QPushButton(label)
            pb.setFixedWidth(52)
            pb.setFixedHeight(22)
            pb.clicked.connect(lambda _, v=val: spin.setValue(v))
            preset_row.addWidget(pb)
        preset_row.addStretch()
        lay.addLayout(preset_row)

        bb_row = QHBoxLayout()
        ok_btn  = QPushButton("Apply")
        ok_btn.setStyleSheet("font-weight:bold;")
        ok_btn.clicked.connect(dlg.accept)
        can_btn = QPushButton("Cancel")
        can_btn.clicked.connect(dlg.reject)
        bb_row.addStretch()
        bb_row.addWidget(ok_btn)
        bb_row.addWidget(can_btn)
        lay.addLayout(bb_row)

        accepted = ((dlg.exec() if QT_VERSION == 6 else dlg.exec_()) ==
                    (QDialog.DialogCode.Accepted if QT_VERSION == 6
                     else QDialog.Accepted))
        if accepted:
            new_size = spin.value()
            self._dv_partition_size = new_size
            self._dv_src_batch.setToolTip(
                "Groups all AF3 jobs from the current session into partitioned\n"
                "batch JSON files (size configurable via AF3 JSON Group Config).\n"
                "Each partition is submitted as an independent SLURM job to\n"
                "prevent GPU out-of-memory errors (anti-OOM).\n\n"
                f"Current group size: {new_size} jobs / partition")
            self._dv_batch_cfg_btn.setToolTip(
                f"Configure the number of AF3 jobs per partition batch.\n"
                f"Current setting: {new_size} jobs / partition\n"
                f"(Default: 50 — anti-OOM safe)")
            self._dv_log(
                f"Batch group size set to {new_size} jobs/partition.", 'submit')

    def _dv_load_from_session_batch(self):
        """Build batch JSON(s) from all session AF3 jobs and stage them for
        sequential upload/submission.

        Anti-OOM rule: if the session has more than the configured partition
        size, the list is automatically split into chunks of that size.  Each
        chunk becomes an independent JSON file submitted as a separate SLURM
        job so the server's RAM is fully released between runs.
        Configure the group size via the AF3 JSON Group Config button.
        """
        if not self.af3_jobs:
            QMessageBox.information(self, "Server",
                "No AF3 jobs in session.\n"
                "Generate jobs in the AlphaFold tab first.")
            return

        prefix     = self._dv_job_prefix.text().strip() or "af3_batch"
        n_total    = len(self.af3_jobs)
        chunk_size = getattr(self, '_dv_partition_size', self._AF3_PARTITION_SIZE)

        # ── Build AF3-format entry list ────────────────────────────
        all_entries = []
        for j in self.af3_jobs:
            all_entries.append({
                "name":       j['name'],
                "modelSeeds": [],
                "sequences":  j.get('sequences', []),
                "dialect":    "alphafoldserver",
                "version":    1,
            })

        # ── Split into partitions ──────────────────────────────────
        partitions = [all_entries[i:i + chunk_size]
                      for i in range(0, n_total, chunk_size)]
        n_parts = len(partitions)

        # ── Write each partition to a temp file ────────────────────
        self._dv_submit_table.setRowCount(0)
        self._dv_pending_jobs = []
        total_kb = 0.0

        for idx, part in enumerate(partitions):
            part_label = (f"{prefix}_part{idx+1:03d}_of_{n_parts:03d}"
                          if n_parts > 1 else f"{prefix}_all_jobs")
            fd, tmp_path = tempfile.mkstemp(
                suffix='.json',
                prefix=re.sub(r'[^\w\-]', '_', part_label) + '_')
            os.close(fd)
            with open(tmp_path, 'w', encoding='utf-8') as f:
                json.dump(part, f, indent=2, ensure_ascii=False)
            file_kb  = os.path.getsize(tmp_path) / 1024
            total_kb += file_kb
            res_part = sum(
                self.af3_jobs[i].get('total_residues', 0)
                for i in range(
                    idx * chunk_size,
                    min(idx * chunk_size + chunk_size, n_total)))

            # One table row per partition
            row = self._dv_submit_table.rowCount()
            self._dv_submit_table.insertRow(row)
            status_lbl = (f"Pending (batch {idx+1}/{n_parts})"
                          if n_parts > 1 else "Pending (batch)")
            for col, val in enumerate([
                    part_label, str(len(part)), str(res_part),
                    f"{file_kb:.0f}", status_lbl]):
                self._dv_submit_table.setItem(
                    row, col, QTableWidgetItem(val))

            self._dv_pending_jobs.append({
                'name':           part_label,
                'local_path':     tmp_path,
                'sequences':      [],       # file already written
                'total_residues': res_part,
                'status':         'pending',
                '_is_batch':      True,
                '_n_jobs':        len(part),
                '_part_index':    idx,      # 0-based
                '_n_parts':       n_parts,
            })

        # ── Summary & UI state ─────────────────────────────────────
        if n_parts > 1:
            summary = (f"{n_parts} partition files  "
                       f"({n_total} jobs total, {chunk_size}/partition, "
                       f"{total_kb:.0f} KB) — submitted sequentially")
        else:
            summary = (f"1 batch file  "
                       f"({n_total} jobs, {total_kb:.0f} KB)")

        self._dv_submit_summary.setText(summary)
        for b in (self._dv_upload_btn, self._dv_run_btn,
                  self._dv_clear_btn, self._dv_dl_jsons_btn):
            b.setEnabled(True)

        if n_parts > 1:
            self._dv_log(
                f"Auto-partitioned {n_total} jobs → {n_parts} × "
                f"{chunk_size} batch JSONs (anti-OOM).  "
                f"Will be submitted one at a time.", 'submit')
        else:
            self._dv_log(
                f"Batch JSON ready: {self._dv_pending_jobs[0]['name']}.json  "
                f"({n_total} jobs, {total_kb:.0f} KB)", 'submit')
        self._dv_refresh_cmd_preview()

    def _dv_load_from_files(self):
        """Load JSON files from disk into the submit table."""
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Select AF3 JSON files", "",
            "JSON files (*.json);;All (*)")
        if not paths:
            return
        self._dv_submit_table.setRowCount(0)
        self._dv_pending_jobs = []
        for p in paths:
            try:
                with open(p, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                # Handle both single job and batch list
                jobs = data if isinstance(data, list) else [data]
                for job_data in jobs:
                    name = job_data.get('name', Path(p).stem)
                    seqs = job_data.get('sequences', [])
                    n_chains = len(seqs)
                    res = sum(len(s.get('proteinChain', {}).get('sequence', ''))
                              * s.get('proteinChain', {}).get('count', 1)
                              for s in seqs)
                    size_kb = Path(p).stat().st_size / 1024
                    row = self._dv_submit_table.rowCount()
                    self._dv_submit_table.insertRow(row)
                    for col, val in enumerate([
                            name, str(n_chains), str(res),
                            f"{size_kb:.1f}", "Pending"]):
                        self._dv_submit_table.setItem(
                            row, col, QTableWidgetItem(val))
                    # Store minimal job info for submission
                    self._dv_pending_jobs.append({
                        'name': name, 'local_path': p,
                        'sequences': seqs, 'total_residues': res,
                        'status': 'pending'
                    })
            except Exception as e:
                self._dv_log(f"Error reading {Path(p).name}: {e}")
        n = len(self._dv_pending_jobs)
        self._dv_submit_summary.setText(f"{n} job(s) loaded from disk")
        for b in (self._dv_upload_btn, self._dv_run_btn,
                  self._dv_clear_btn, self._dv_dl_jsons_btn):
            b.setEnabled(n > 0)
        self._dv_refresh_cmd_preview()

    def _dv_clear_submit_list(self):
        self._dv_submit_table.setRowCount(0)
        self._dv_pending_jobs = []
        self._dv_submit_summary.setText("0 jobs loaded")
        for b in (self._dv_upload_btn, self._dv_run_btn,
                  self._dv_dl_jsons_btn):
            b.setEnabled(False)

    def _dv_download_staged_jsons(self):
        """Save all staged partition JSON files to a local folder.

        For single-batch sessions the one JSON is saved; for auto-partitioned
        sessions all N partition files are written to the chosen folder so the
        user can inspect, archive, or manually scp them to the cluster.
        """
        jobs = getattr(self, '_dv_pending_jobs', [])
        batch_jobs = [j for j in jobs if j.get('local_path')]
        if not batch_jobs:
            batch_jobs = jobs  # fall back to individual jobs from session

        if not batch_jobs:
            QMessageBox.information(
                self, "Download JSONs", "No staged jobs to download.")
            return

        dest_dir = QFileDialog.getExistingDirectory(
            self, "Choose folder to save JSON files", "")
        if not dest_dir:
            return

        saved  = []
        errors = []
        for job in batch_jobs:
            fname = re.sub(r'[^\w\-]', '_', job['name']) + '.json'
            dst   = os.path.join(dest_dir, fname)
            try:
                if job.get('local_path') and os.path.exists(job['local_path']):
                    shutil.copy2(job['local_path'], dst)
                else:
                    entry = {
                        "name":       job['name'],
                        "modelSeeds": [],
                        "sequences":  job.get('sequences', []),
                        "dialect":    "alphafoldserver",
                        "version":    1,
                    }
                    with open(dst, 'w', encoding='utf-8') as f:
                        json.dump(entry, f, indent=2, ensure_ascii=False)
                saved.append(fname)
            except Exception as e:
                errors.append(f"{fname}: {e}")

        msg_parts = [f"✓ {len(saved)} JSON(s) saved to:\n{dest_dir}"]
        if errors:
            msg_parts.append(f"\n⚠ {len(errors)} error(s):\n" +
                             "\n".join(errors))
        QMessageBox.information(self, "Download JSONs", "\n".join(msg_parts))
        self._dv_log(
            f"JSONs downloaded: {len(saved)} file(s) → {dest_dir}", 'submit')

    # ── AF3 preset helper ──────────────────────────────────────
    def _dv_af3_apply_preset(self, idx: int):
        """Legacy no-op stub.

        The preset selector was removed when the AF3 Advanced Options
        panel was replaced by the FlexibleAF3SubmitWidget.  This method
        is kept as a no-op for backward compatibility — older project
        files (saved before v2.0) may still try to call it.
        """
        return  # FlexibleAF3SubmitWidget handles profiles now

    def _dv_refresh_cmd_preview(self):
        """Rebuild the command preview using the FlexibleAF3SubmitWidget.

        The user-selected profile + template + parameters define the
        actual command. This function only resolves them against the
        runtime context (json_path, job_name, parent_dir, etc.) and
        formats the result for the master command-preview widget.

        For multi-partition scenarios (anti-OOM split), every partition's
        resolved command is shown so the user can see all job codes.
        """
        prefix   = self._dv_job_prefix.text().strip() or "af3_batch"
        base     = self._dv_base_path.text().strip().rstrip('/')
        rdir     = (self._dv_remote_dir.text().strip()
                    or datetime.now().strftime('%Y-%m-%d'))
        mod_name = self._dv_module_cmd.text().strip()
        jobs     = getattr(self, '_dv_pending_jobs', [])
        n_jobs   = len(jobs)
        use_ts   = (getattr(self, '_dv_ts_check', None)
                    and self._dv_ts_check.isChecked())

        ts          = datetime.now().strftime('%Y%m%d_%H%M%S') if use_ts else ''
        parent_dir  = f"{base}/{rdir}"

        # ── Build preview lines ───────────────────────────────
        lines = []
        if mod_name:
            lines.append(f"# module loaded: {mod_name}")

        # JSON-level annotations (modelSeeds + templates)
        if hasattr(self, '_dv_af3_seeds'):
            n_seeds = self._dv_af3_seeds.value()
            if n_seeds > 1:
                lines.append(
                    f"# {n_seeds} seeds will be injected into JSON modelSeeds")
        if (hasattr(self, '_dv_af3_use_templates')
                and not self._dv_af3_use_templates.isChecked()):
            lines.append("# templates will be stripped from JSON")

        # Detect multi-partition scenario
        is_partitioned = (n_jobs > 1 and jobs
                          and jobs[0].get('_is_batch')
                          and jobs[0].get('_n_parts', 1) > 1)

        # Helper: assemble runtime context for one batch
        def _ctx_for(json_fname, job_name):
            return {
                "prefix":     prefix,
                "parent_dir": parent_dir,
                "json_path":  json_fname,
                "job_name":   job_name,
                "output_dir": f"{parent_dir}/{job_name}/output",
                "date":       datetime.now().strftime('%Y-%m-%d'),
                "timestamp":  ts or datetime.now().strftime('%Y%m%d_%H%M%S'),
            }

        if is_partitioned:
            n_parts      = jobs[0]['_n_parts']
            total_afjobs = sum(j.get('_n_jobs', 0) for j in jobs)
            partition_sz = getattr(
                self, '_dv_partition_size',
                getattr(self, '_AF3_PARTITION_SIZE', 50))
            lines.append(
                f"# {total_afjobs} AF3 jobs → {n_parts} partitions "
                f"(≤{partition_sz}/batch, anti-OOM)")
            lines.append("# Submitted sequentially — one job at a time")
            lines.append("")
            for idx, j in enumerate(jobs):
                p_name      = j['name']
                json_fname  = f"{p_name}.json"
                job_name    = f"{p_name}_{ts}" if ts else p_name
                job_name    = re.sub(r'[()\[\]{}|;&!\s]+', '_',
                                      job_name).strip('_')
                lines.append(f"# ── Partition {idx+1}/{n_parts} "
                              f"({j.get('_n_jobs', '?')} jobs) ──")
                if hasattr(self, '_dv_af3_widget'):
                    lines.append(self._dv_af3_widget.build_command(
                        _ctx_for(json_fname, job_name)))
                else:
                    # Fallback: legacy bare command
                    cmd = self._dv_af3cmd.text().strip() or "af3_run"
                    lines.append(f"cd {parent_dir} && {cmd} "
                                 f"--json_path {json_fname} "
                                 f"--job-name {job_name}")
                lines.append(f"# Output → {parent_dir}/{job_name}/output/")
                lines.append("")
        else:
            # Single batch or individual jobs
            json_fname = f"{prefix}_all_jobs.json"
            job_name   = f"{prefix}_{ts}" if ts else prefix
            job_name   = re.sub(r'[()\[\]{}|;&!\s]+', '_',
                                 job_name).strip('_')
            if hasattr(self, '_dv_af3_widget'):
                lines.append(self._dv_af3_widget.build_command(
                    _ctx_for(json_fname, job_name)))
            else:
                cmd = self._dv_af3cmd.text().strip() or "af3_run"
                lines.append(f"cd {parent_dir} && {cmd} "
                             f"--json_path {json_fname} "
                             f"--job-name {job_name}")
            lines.append(f"# Output → {parent_dir}/{job_name}/output/")
            if n_jobs:
                lines.append(f"# {n_jobs} job(s) bundled in batch JSON")

        # Append profile summary
        if hasattr(self, '_dv_af3_widget'):
            summary = self._dv_af3_widget.current_profile_summary()
            if summary:
                lines.append(summary)

        self._dv_cmd_preview.setPlainText('\n'.join(lines))

    def _dv_upload_only(self):
        """Upload JSON files to the server without submitting."""
        self._dv_do_upload(submit=False)

    def _dv_upload_and_submit(self):
        """Upload JSON files and run af3_run on the server."""
        self._dv_do_upload(submit=True)

    def _dv_do_upload(self, submit: bool):
        """SFTP upload + optional af3_run — sequential multi-partition support.

        For sessions with ≤50 jobs (or loaded as a single batch), behaviour is
        unchanged: one JSON file is built, uploaded, and optionally submitted.

        For sessions that were auto-partitioned into N×50 chunks, the worker
        thread loops through every partition and submits them one at a time so
        that the server's RAM is fully released between runs (anti-OOM fix).

        Each partition is uploaded and submitted before the next one starts.
        All results are returned as a list and handled by _dv_on_upload_done.
        """
        if not self._ssh_client:
            QMessageBox.warning(self, "Server",
                "Not connected. Use the Connect tab first.")
            return
        jobs = getattr(self, '_dv_pending_jobs', [])
        if not jobs:
            QMessageBox.warning(self, "Server", "No jobs loaded.")
            return

        base   = self._dv_base_path.text().strip().rstrip('/')
        prefix = self._dv_job_prefix.text().strip() or "af3_batch"
        rdir   = (self._dv_remote_dir.text().strip()
                  or datetime.now().strftime('%Y-%m-%d'))
        cmd    = self._dv_af3cmd.text().strip() or "af3_run"

        use_ts   = getattr(self, '_dv_ts_check', None) and \
                   self._dv_ts_check.isChecked()
        ts       = datetime.now().strftime('%Y%m%d_%H%M%S') if use_ts else ''

        parent_dir = f"{base}/{rdir}"
        ssh = self._ssh_client

        # ── Collect AF3 advanced settings (read once, applied to all parts) ─
        n_seeds       = 1
        use_templates = True
        if hasattr(self, '_dv_af3_seeds'):
            n_seeds = self._dv_af3_seeds.value()
        if hasattr(self, '_dv_af3_use_templates'):
            use_templates = self._dv_af3_use_templates.isChecked()
        random.seed(42)
        model_seeds = ([random.randint(1, 2**31 - 1) for _ in range(n_seeds)]
                       if n_seeds > 1 else [])

        # Note: command-line flags / partition / extra_flags are now
        # handled by self._dv_af3_widget (FlexibleAF3SubmitWidget) — the
        # whole af3_run command is built from the user-selected profile.
        # Only the JSON-level params (model_seeds, use_templates) remain
        # as separate widgets because they modify the JSON CONTENT, not
        # the command.

        mod_prefix_str = self._dv_build_activation_prefix()

        def _build_batch_list_for_job(job):
            """Return AF3 batch list for one pending job entry."""
            batch_list = []
            if job.get('local_path'):
                try:
                    with open(job['local_path'], 'r', encoding='utf-8') as f:
                        raw = json.load(f)
                    if isinstance(raw, list):
                        for entry in raw:
                            if isinstance(entry, dict):
                                entry['modelSeeds'] = model_seeds
                                if not use_templates:
                                    for seq in entry.get('sequences', []):
                                        for chain in (seq.get(
                                                'proteinChain', []) or []):
                                            chain['templates'] = []
                        batch_list.extend(raw)
                    else:
                        raw['modelSeeds'] = model_seeds
                        batch_list.append(raw)
                except Exception:
                    batch_list.append({
                        "name":       job['name'],
                        "modelSeeds": model_seeds,
                        "sequences":  job.get('sequences', []),
                        "dialect":    "alphafoldserver",
                        "version":    1,
                    })
            else:
                batch_list.append({
                    "name":       job['name'],
                    "modelSeeds": model_seeds,
                    "sequences":  job.get('sequences', []),
                    "dialect":    "alphafoldserver",
                    "version":    1,
                })
            return batch_list

        def _upload_one_partition(job, resolved_parent):
            """Upload batch JSON for a single partition.  Returns remote path."""
            safe_name    = re.sub(r'[^\w\-]', '_', job['name'])
            batch_list   = _build_batch_list_for_job(job)
            batch_fname  = f"{job['name']}.json"

            fd, tmp_path = tempfile.mkstemp(
                suffix='.json', prefix=safe_name + '_')
            os.close(fd)
            try:
                with open(tmp_path, 'w', encoding='utf-8') as f:
                    json.dump(batch_list, f, indent=2, ensure_ascii=False)
                file_kb = os.path.getsize(tmp_path) / 1024

                sftp = ssh.open_sftp()
                try:
                    remote_json_sftp = f"{resolved_parent}/{batch_fname}"
                    sftp.put(tmp_path, remote_json_sftp)
                finally:
                    try:
                        sftp.close()
                    except Exception:
                        pass
            finally:
                if os.path.exists(tmp_path):
                    os.unlink(tmp_path)

            return batch_fname, batch_list, file_kb

        def _submit_one(batch_fname, job_name_full):
            """Run the user-defined AF3 command for a single partition JSON.

            The command is built by the FlexibleAF3SubmitWidget — i.e. the
            user-selected profile/template, with placeholders resolved
            against the runtime context (json_path, job_name, parent_dir,
            output_dir, prefix, date, timestamp).
            """
            ctx = {
                "prefix":     prefix,
                "parent_dir": parent_dir,
                "json_path":  batch_fname,
                "job_name":   job_name_full,
                "output_dir": f"{parent_dir}/{job_name_full}/output",
                "date":       datetime.now().strftime('%Y-%m-%d'),
                "timestamp":  ts or datetime.now().strftime('%Y%m%d_%H%M%S'),
            }

            if hasattr(self, '_dv_af3_widget'):
                af3_cmd_resolved = self._dv_af3_widget.build_command(ctx)
            else:
                # Defensive fallback: legacy bare af3_run command.
                af3_cmd_resolved = (
                    f"cd {parent_dir} && "
                    f"{cmd} --json_path {batch_fname} "
                    f'--job-name "{job_name_full}"'
                )

            # Strip '# comment' lines from the resolved template — they're
            # fine in shell but pollute the bash -lc string.  Also collapse
            # line continuations so multi-line templates work as one cmd.
            cmd_lines = []
            for ln in af3_cmd_resolved.splitlines():
                stripped = ln.strip()
                if not stripped or stripped.startswith('#'):
                    continue
                cmd_lines.append(ln.rstrip(' \\'))
            inner = ' '.join(cmd_lines) if cmd_lines else ''

            _inner_cmd = f"{mod_prefix_str}{inner}"
            run_cmd = f'bash -lc "{_inner_cmd}"'
            out, err, rc = self._dv_ssh_exec(run_cmd, timeout=60)

            # Match SLURM "Submitted batch job N", "batch job N",
            # PBS "<id>.<host>", LSF "Job <id> ..."
            m = re.search(
                r'(?:Submitted batch job|batch job|Job)\s+(\d+)',
                out + err, re.IGNORECASE)
            if m:
                return m.group(1), 'submitted', False
            # PBS-style "12345.cluster"
            m = re.search(r'^(\d+)\.[a-zA-Z]', out.strip())
            if m:
                return m.group(1), 'submitted', False
            combined = (out + err).lower()
            if 'already exists' in combined or 'directory already' in combined:
                return None, 'dir_exists', True
            lines_o = (out + err).strip().splitlines()
            detail  = lines_o[0] if lines_o else f'rc={rc}'
            return None, f'submit_error: {detail}', False

        def _do_upload_submit():
            # ── Resolve ~ for SFTP once ────────────────────────
            resolved_parent = parent_dir
            if parent_dir.startswith('~'):
                home_out, _, _ = self._dv_ssh_exec("echo $HOME", timeout=5)
                home = home_out.strip()
                if home:
                    resolved_parent = home + parent_dir[1:]

            # Create parent directory if needed
            sftp_chk = ssh.open_sftp()
            try:
                try:
                    sftp_chk.stat(resolved_parent)
                except FileNotFoundError:
                    self._dv_ssh_exec(f"mkdir -p {parent_dir}")
                    try:
                        sftp_chk.stat(resolved_parent)
                    except FileNotFoundError:
                        raise RuntimeError(
                            f"Could not create remote directory:\n"
                            f"{resolved_parent}")
            finally:
                try:
                    sftp_chk.close()
                except Exception:
                    pass

            all_results = []
            for job in jobs:
                # Build a per-partition job name
                j_name_clean = re.sub(r'[()\[\]{}|;&!\s]+', '_',
                                      job['name']).strip('_')
                job_name_full = (f"{j_name_clean}_{ts}"
                                 if ts else j_name_clean)

                # Upload
                batch_fname, batch_list, file_kb = _upload_one_partition(
                    job, resolved_parent)

                output_dir = f"{parent_dir}/{job_name_full}"

                status          = 'uploaded'
                slurm_id        = None
                dir_exists_error = False

                if submit:
                    slurm_id, status, dir_exists_error = _submit_one(
                        batch_fname, job_name_full)

                all_results.append({
                    'name':             job_name_full,
                    'prefix':           prefix,
                    'batch_file':       batch_fname,
                    'n_jobs':           len(batch_list),
                    'file_kb':          round(file_kb, 1),
                    'remote_json':      f"{parent_dir}/{batch_fname}",
                    'remote_dir':       output_dir,
                    'slurm_id':         slurm_id,
                    'status':           status,
                    'dir_exists_error': dir_exists_error,
                    'parent_dir':       parent_dir,
                    'rel_json':         batch_fname,
                    'mod_name':         self._dv_module_cmd.text().strip(),
                    'cmd':              cmd,
                    '_part_index':      job.get('_part_index', 0),
                    '_n_parts':         job.get('_n_parts', 1),
                })

            if not all_results:
                raise RuntimeError("No valid jobs to submit.")
            return all_results

        self._dv_run_btn.setEnabled(False)
        self._dv_upload_btn.setEnabled(False)
        n_parts = sum(1 for j in jobs if j.get('_is_batch'))
        n_label = (f"{n_parts} partition(s)" if n_parts > 1
                   else f"{len(jobs)} job(s)")
        self._dv_log(
            f"Starting {'upload+submit' if submit else 'upload only'} "
            f"({n_label})...", 'submit')

        w = AnalysisWorker(_do_upload_submit)
        w.finished.connect(self._dv_on_upload_done)
        w.error.connect(lambda e: (
            self._dv_log(f"Upload error: {e}", 'submit'),
            self._dv_run_btn.setEnabled(True),
            self._dv_upload_btn.setEnabled(True)
        ))
        self._hpc_workers.append(w)
        w.start()

    def _dv_on_upload_done(self, results):
        """Callback when upload/submit worker finishes.
        results is a list — one dict per partition (usually 1, but N for
        auto-partitioned sessions > _AF3_PARTITION_SIZE jobs).
        """
        self._dv_run_btn.setEnabled(True)
        self._dv_upload_btn.setEnabled(True)

        if not results:
            self._dv_log("Upload returned no results.", 'submit')
            return

        n_parts = len(results)
        total_jobs = sum(r.get('n_jobs', 1) for r in results)

        if n_parts > 1:
            submitted = sum(1 for r in results if r.get('slurm_id'))
            slurm_ids = [r['slurm_id'] for r in results if r.get('slurm_id')]
            self._dv_log(
                f"  {n_parts} partitions submitted — "
                f"{submitted}/{n_parts} queued", 'submit')
            for i, r in enumerate(results):
                sid = r.get('slurm_id', '')
                st  = r.get('status', '?')
                self._dv_log(
                    f"  Part {i+1:03d}: {r.get('batch_file','?')}  "
                    f"{r.get('n_jobs','?')} jobs  {r.get('file_kb',0)} KB  "
                    f"{st}" + (f"  [SLURM {sid}]" if sid else ''),
                    'submit')
                # Update the matching row in the submit table
                if i < self._dv_submit_table.rowCount():
                    self._dv_submit_table.setItem(
                        i, 4, QTableWidgetItem(
                            st + (f" [SLURM {sid}]" if sid else '')))
                # Register each in the monitor list
                self._hpc_jobs.append({
                    'name':        r.get('name', ''),
                    'slurm_id':    sid or '',
                    'remote_dir':  r.get('remote_dir', ''),
                    'remote_json': r.get('remote_json', ''),
                    'status':      st,
                    'local_output': '',
                })
            self._dv_refresh_monitor_table()
            id_summary = (", ".join(slurm_ids[:5])
                          + ("…" if len(slurm_ids) > 5 else ''))
            self._status.showMessage(
                f"✓ Server: {total_jobs} jobs in {n_parts} partitions — "
                f"SLURM [{id_summary}]" if slurm_ids else
                f"✓ Server: {total_jobs} jobs in {n_parts} partitions uploaded")
            return

        # ── Single partition / legacy path ────────────────────────
        r = results[0]
        n_jobs   = r.get('n_jobs', 1)
        file_kb  = r.get('file_kb', 0)
        status   = r.get('status', '?')
        slurm_id = r.get('slurm_id')
        batch    = r.get('batch_file', r.get('name', '?'))

        self._dv_log(
            f"  Batch file : {batch}  ({n_jobs} jobs, {file_kb} KB)",
            'submit')
        self._dv_log(
            f"  Remote dir : {r.get('remote_dir', '?')}",
            'submit')
        self._dv_log(
            f"  Status     : {status}"
            + (f"   [SLURM {slurm_id}]" if slurm_id else ''),
            'submit')

        # ── Handle "job directory already exists" ──────────────
        if r.get('dir_exists_error'):
            prefix     = r.get('prefix', '')
            remote_dir = r.get('remote_dir', '')
            parent_dir = r.get('parent_dir', '')
            rel_json   = r.get('rel_json', '')
            cmd        = r.get('cmd', 'af3_run')
            mod_name   = r.get('mod_name', '')

            dlg = QDialog(self)
            dlg.setWindowTitle("Job directory already exists")
            dlg.setFixedSize(520, 230)
            dlg_lay = QVBoxLayout(dlg)
            dlg_lay.addWidget(QLabel(
                f"<b>⚠ The job directory already exists on the server:</b><br>"
                f"<code>{remote_dir}</code><br><br>"
                f"Choose what to do:"))
            btn_row = QHBoxLayout()

            btn_delete = QPushButton("🗑 Delete old dir and resubmit")
            btn_delete.setToolTip(
                f"Runs:  rm -rf {remote_dir}  then resubmits af3_run")
            btn_rename = QPushButton("✏ Change job name prefix")
            btn_rename.setToolTip(
                "Change the Job name prefix field and resubmit manually")
            btn_cancel = QPushButton("Cancel")
            btn_delete.setFixedHeight(32)
            btn_rename.setFixedHeight(32)
            btn_cancel.setFixedHeight(32)
            btn_row.addWidget(btn_delete)
            btn_row.addWidget(btn_rename)
            btn_row.addWidget(btn_cancel)
            dlg_lay.addLayout(btn_row)

            _choice = ['']
            btn_delete.clicked.connect(
                lambda: (_choice.__setitem__(0, 'delete'), dlg.accept()))
            btn_rename.clicked.connect(
                lambda: (_choice.__setitem__(0, 'rename'), dlg.accept()))
            btn_cancel.clicked.connect(dlg.reject)

            if not (dlg.exec() if QT_VERSION == 6 else dlg.exec_()):
                return
            choice = _choice[0]

            if choice == 'rename':
                # Just focus the prefix field so user can edit it
                self._dv_job_prefix.setFocus()
                self._dv_job_prefix.selectAll()
                self._dv_log(
                    "  → Change the Job name prefix and resubmit.",
                    'submit')
                self._dv_tabs.setCurrentIndex(1)   # Submit tab
                return

            if choice == 'delete':
                self._dv_log(
                    f"  → Deleting {remote_dir} and resubmitting...",
                    'submit')
                self._dv_run_btn.setEnabled(False)

                def _do_delete_resubmit():
                    # Delete old directory
                    self._dv_ssh_exec(f"rm -rf {remote_dir}", timeout=30)
                    # Resubmit
                    mod_prefix = (f"module load {mod_name} 2>/dev/null; "
                                  if mod_name else "")
                    run_cmd = (
                        f"bash -lc '"
                        f"{mod_prefix}"
                        f"cd {parent_dir} && "
                        f"{cmd} --json_path {rel_json} "
                        f"--job-name {prefix}"
                        f"'"
                    )
                    out, err, rc = self._dv_ssh_exec(run_cmd, timeout=60)
                    m = re.search(
                        r'batch job\s+(\d+)', out + err, re.IGNORECASE)
                    if m:
                        return {'slurm_id': m.group(1), 'status': 'submitted',
                                'output': out + err}
                    else:
                        lines = (out + err).strip().splitlines()
                        detail = lines[0] if lines else f'rc={rc}'
                        return {'slurm_id': None,
                                'status': f'submit_error: {detail}',
                                'output': out + err}

                def _on_resubmit(res):
                    self._dv_run_btn.setEnabled(True)
                    sid = res.get('slurm_id')
                    st  = res.get('status', '?')
                    self._dv_log(
                        f"  Resubmit status: {st}"
                        + (f"  [SLURM {sid}]" if sid else ''),
                        'submit')
                    if sid:
                        self._hpc_jobs.append({
                            'name':        prefix,
                            'slurm_id':    sid,
                            'remote_dir':  remote_dir,
                            'remote_json': r.get('remote_json', ''),
                            'status':      'submitted',
                            'local_output': '',
                        })
                        self._dv_refresh_monitor_table()
                    self._status.showMessage(
                        f"✓ Resubmit: {st}"
                        + (f"  SLURM {sid}" if sid else ''))

                w = AnalysisWorker(_do_delete_resubmit)
                w.finished.connect(_on_resubmit)
                w.error.connect(lambda e: (
                    self._dv_log(f"  Resubmit error: {e}", 'submit'),
                    self._dv_run_btn.setEnabled(True)))
                self._hpc_workers.append(w)
                w.start()
                return

        # ── Normal success / other errors ─────────────────────
        # Update submit table status
        for row in range(self._dv_submit_table.rowCount()):
            item = self._dv_submit_table.item(row, 4)
            if item is not None:
                self._dv_submit_table.setItem(row, 4, QTableWidgetItem(status))

        # Register in monitor list
        self._hpc_jobs.append({
            'name':        r.get('name', batch),
            'slurm_id':    slurm_id or '',
            'remote_dir':  r.get('remote_dir', ''),
            'remote_json': r.get('remote_json', ''),
            'status':      status,
            'local_output': '',
        })
        self._dv_refresh_monitor_table()
        self._status.showMessage(
            f"✓ Server: {n_jobs} job(s) in batch — {status}"
            + (f"  SLURM {slurm_id}" if slurm_id else ''))

    # ── Monitor helpers ─────────────────────────────────────────

    def _dv_update_poll_timer(self, text: str):
        """Start/stop the poll timer based on combobox selection."""
        self._hpc_poll_timer.stop()
        intervals = {'15 s': 15000, '30 s': 30000,
                     '60 s': 60000, '2 min': 120000}
        ms = intervals.get(text, 0)
        if ms > 0 and self._ssh_client:
            self._hpc_poll_timer.start(ms)

    def _hpc_poll_queue(self):
        """Poll the scheduler queue on the remote server.
        Supports SLURM, PBS, LSF, or None based on the scheduler selector."""
        if not self._ssh_client:
            return
        sched = self._dv_scheduler_cmds()
        poll_cmd = sched.get("poll", "")
        if not poll_cmd or "No scheduler" in poll_cmd:
            return

        def _do_poll():
            out, err, _ = self._dv_ssh_exec(poll_cmd, timeout=15)
            return out, err

        w = AnalysisWorker(_do_poll)
        w.finished.connect(self._dv_on_poll_done)
        w.error.connect(lambda e: self._dv_squeue_raw.setPlainText(
            f"Poll error: {e}"))
        self._hpc_workers.append(w)
        w.start()

    def _dv_on_poll_done(self, result):
        """Parse squeue output, refresh monitor table, and store all-jobs list.

        Parsed fields (matches the squeue --format used in _hpc_poll_queue):
          col0=JOBID  col1=PARTITION  col2=NAME  col3=USER
          col4=ST     col5=TIME       col6=NODES
        """
        out, err = result
        self._dv_squeue_raw.setPlainText(out or err or "(no output)")
        self._dv_queue_lbl.setText(
            f"squeue — {datetime.now().strftime('%H:%M:%S')}")

        # ── Parse all running/pending jobs ────────────────────
        running = {}          # job_name → dict
        all_server_jobs = []  # for "show all" mode

        for line in out.splitlines():
            parts = line.strip().split()
            # Header or short lines
            if not parts or not parts[0].isdigit():
                continue
            if len(parts) < 5:
                continue
            slurm_id  = parts[0]
            partition = parts[1] if len(parts) > 1 else '-'
            name      = parts[2] if len(parts) > 2 else '-'
            state     = parts[4] if len(parts) > 4 else '?'
            elapsed   = parts[5] if len(parts) > 5 else '-'
            node      = parts[6] if len(parts) > 6 else '-'

            running[name] = {
                'slurm_id': slurm_id, 'status': state,
                'partition': partition, 'node': node, 'time': elapsed,
            }
            all_server_jobs.append({
                'name': name, 'slurm_id': slurm_id, 'status': state,
                'partition': partition, 'node': node, 'time': elapsed,
            })

        # Store for "show all server jobs" mode
        self._dv_squeue_parsed_all = all_server_jobs

        # ── Update internal session-job list ──────────────────
        for job in self._hpc_jobs:
            if job['name'] in running:
                info = running[job['name']]
                job['slurm_id']  = info['slurm_id']
                job['status']    = info['status']
                job['partition'] = info['partition']
                job['node']      = info['node']
                job['time']      = info['time']
            elif job.get('status') not in (
                    'COMPLETED', 'FAILED', 'CANCELLED', 'uploaded'):
                if job.get('slurm_id'):
                    job['status'] = 'COMPLETED'

        self._dv_refresh_monitor_table()

    def _dv_refresh_monitor_table(self):
        """Rebuild the monitor table applying status/partition filters.

        When "Show all server jobs" is checked the table also shows jobs
        parsed from the raw squeue output that are NOT in _hpc_jobs
        (i.e. jobs submitted outside this session).
        """
        self._dv_monitor_table.setSortingEnabled(False)  # disable during rebuild
        self._dv_monitor_table.setRowCount(0)
        status_colors = {
            'RUNNING':   '#e3f2fd', 'PENDING':   '#fff9c4',
            'COMPLETED': '#e8f5e9', 'FAILED':    '#ffebee',
            'CANCELLED': '#fce4ec', 'TIMEOUT':   '#fff3e0',
            'submitted': '#e8f5e9', 'uploaded':  '#f3e5f5',
        }

        # Active filter values
        st_filter   = getattr(self, '_dv_mon_status_filter',
                               None)
        part_filter = getattr(self, '_dv_mon_partition_filter',
                               None)
        show_all    = getattr(self, '_dv_mon_show_all', None)

        flt_st   = st_filter.currentText() if st_filter else 'All'
        flt_part = part_filter.text().strip().lower() if part_filter else ''
        do_all   = show_all.isChecked() if show_all else False

        # Build the row list ─────────────────────────────────
        # Start with session jobs
        rows_to_show = list(self._hpc_jobs)

        # If "show all" is on, merge in squeue-parsed rows
        if do_all and hasattr(self, '_dv_squeue_parsed_all'):
            session_names = {j['name'] for j in self._hpc_jobs}
            for sj in self._dv_squeue_parsed_all:
                if sj['name'] not in session_names:
                    rows_to_show.append(sj)

        for job in rows_to_show:
            st   = job.get('status', '?')
            part = job.get('partition', '-')

            # ── Apply filters ──────────────────────────────
            if flt_st != 'All' and st.upper() != flt_st.upper():
                continue
            if flt_part and flt_part not in part.lower():
                continue

            color = status_colors.get(st, '#ffffff')
            row = self._dv_monitor_table.rowCount()
            self._dv_monitor_table.insertRow(row)

            vals = [job.get('name', '-'),
                    job.get('slurm_id', '-'),
                    st,
                    part,
                    job.get('node', '-'),
                    job.get('time', '-'),
                    '']
            for col, val in enumerate(vals):
                item = QTableWidgetItem(str(val))
                item.setBackground(QColor(color))
                self._dv_monitor_table.setItem(row, col, item)

            # Action buttons in last column
            if st in ('COMPLETED', 'uploaded', 'submitted'):
                btn = QPushButton("⬇ Get results")
                btn.setFixedHeight(22)
                btn.setStyleSheet("background:#1b5e20;color:white;font-size:11px;")
                btn.clicked.connect(
                    lambda _, j=job: self._dv_autofill_results_path(j))
                self._dv_monitor_table.setCellWidget(row, 6, btn)
            elif st in ('RUNNING', 'PENDING'):
                btn = QPushButton("⛔ Cancel")
                btn.setFixedHeight(22)
                btn.setStyleSheet("background:#b71c1c;color:white;font-size:11px;")
                sid = job.get('slurm_id', '')
                btn.clicked.connect(
                    lambda _, s=sid, n=job.get('name','?'):
                        self._dv_cancel_job_by_id(s, n))
                self._dv_monitor_table.setCellWidget(row, 6, btn)

        self._dv_monitor_table.setSortingEnabled(True)

    def _dv_cancel_job_by_id(self, slurm_id: str, job_name: str):
        """Cancel a specific SLURM job by ID (used from inline table button)."""
        if not self._ssh_client or not slurm_id:
            return
        ans = QMessageBox.question(self, "Cancel job",
            f"Cancel  {job_name}  (SLURM {slurm_id})?")
        if ans != (QMessageBox.StandardButton.Yes if QT_VERSION == 6
                   else QMessageBox.Yes):
            return

        def _do():
            return self._dv_ssh_exec(self._dv_scheduler_cmds()["cancel"].format(job_id=slurm_id), timeout=10)

        worker = AnalysisWorker(_do)
        worker.finished.connect(lambda r: (
            self._dv_log(
                f"scancel {slurm_id} → rc={r[2]}  {(r[0]+r[1]).strip()}"),
            self._hpc_poll_queue()
        ))
        worker.error.connect(
            lambda e: self._dv_log(f"scancel error: {e}"))
        self._hpc_workers.append(worker)
        worker.start()

    def _dv_autofill_results_path(self, job: dict):
        """Fill the Results tab remote path from a completed job.
        Output is at: {base}/{rdir}/{prefix}/output/
        which equals: remote_dir/output/
        """
        remote_dir = job.get('remote_dir', '')
        if remote_dir:
            # af3_run creates output/ inside the job dir (= remote_dir)
            output_path = f"{remote_dir}/output"
            self._dv_remote_output.setText(output_path)
            self._dv_tabs.setCurrentIndex(3)
            self._status.showMessage(
                f"Results tab ready — {output_path}")

    # ── Results helpers ─────────────────────────────────────────

    def _dv_browse_remote_output(self):
        """Simple dialog to pick remote path by typing."""
        path, ok = (QInputDialog.getText(
            self, "Remote output path",
            "Enter full remote path to output/ directory:",
            text=self._dv_remote_output.text())
            if QT_VERSION == 6 else
            QInputDialog.getText(
            self, "Remote output path",
            "Enter full remote path to output/ directory:",
            text=self._dv_remote_output.text()))
        if ok and path:
            self._dv_remote_output.setText(path)

    def _dv_choose_local_dest(self):
        folder = QFileDialog.getExistingDirectory(
            self, "Choose local destination folder")
        if folder:
            self._dv_local_dest.setText(folder)

    def _dv_list_remote_output(self):
        """List files in the remote output directory."""
        if not self._ssh_client:
            QMessageBox.warning(self, "Server",
                "Not connected. Use the Connect tab first.")
            return
        remote_path = self._dv_remote_output.text().strip()
        if not remote_path:
            QMessageBox.warning(self, "Server",
                "Enter a remote output path.")
            return

        def _do_list():
            # Open fresh SFTP channel (avoids stale socket)
            sftp = self._ssh_client.open_sftp()
            try:
                # Expand ~ on remote
                if remote_path.startswith('~'):
                    out, _, _ = self._dv_ssh_exec("echo $HOME", timeout=5)
                    home = out.strip()
                    rpath = home + remote_path[1:]
                else:
                    rpath = remote_path
                items = []
                try:
                    for attr in sftp.listdir_attr(rpath):
                        is_dir = stat.S_ISDIR(attr.st_mode)
                        size = attr.st_size or 0
                        mtime = datetime.fromtimestamp(
                            attr.st_mtime).strftime('%Y-%m-%d %H:%M') \
                            if attr.st_mtime else '-'
                        items.append({
                            'name': attr.filename,
                            'type': 'dir' if is_dir else 'file',
                            'size': f"{size/1024:.1f} KB" if not is_dir else '-',
                            'modified': mtime,
                            'full_path': f"{rpath}/{attr.filename}",
                            'is_dir': is_dir,
                        })
                except FileNotFoundError:
                    raise RuntimeError(
                        f"Remote path not found:\n{rpath}")
                return items
            finally:
                try:
                    sftp.close()
                except Exception:
                    pass

        w = AnalysisWorker(_do_list)
        w.finished.connect(self._dv_on_list_done)
        w.error.connect(lambda e: QMessageBox.critical(
            self, "Server — list error", e))
        self._hpc_workers.append(w)
        w.start()

    def _dv_on_list_done(self, items):
        """Populate the remote file listing table."""
        self._dv_remote_file_table.setRowCount(0)
        for it in items:
            row = self._dv_remote_file_table.rowCount()
            self._dv_remote_file_table.insertRow(row)
            for col, val in enumerate([
                    it['name'], it['type'], it['size'], it['modified']]):
                self._dv_remote_file_table.setItem(
                    row, col, QTableWidgetItem(val))
        self._dv_dl_progress.setText(
            f"{len(items)} items in remote directory")

    def _dv_download_results(self):
        """SFTP download selected (or all) result files/folders."""
        if not self._ssh_client:
            QMessageBox.warning(self, "Server",
                "Not connected.")
            return
        remote_path = self._dv_remote_output.text().strip()
        local_dest  = self._dv_local_dest.text().strip()
        if not remote_path or not local_dest:
            QMessageBox.warning(self, "Server",
                "Set both remote output path and local destination.")
            return

        os.makedirs(local_dest, exist_ok=True)
        auto_import = self._dv_auto_import.isChecked()

        # Collect selected rows (or all if nothing selected)
        selected_names = set()
        for idx in self._dv_remote_file_table.selectedIndexes():
            item = self._dv_remote_file_table.item(idx.row(), 0)
            if item:
                selected_names.add(item.text())
        if not selected_names:
            # Select all
            for row in range(self._dv_remote_file_table.rowCount()):
                item = self._dv_remote_file_table.item(row, 0)
                if item:
                    selected_names.add(item.text())

        def _do_download():
            # Open fresh SFTP channel (avoids stale socket)
            sftp = self._ssh_client.open_sftp()
            try:
                # Expand ~
                if remote_path.startswith('~'):
                    out, _, _ = self._dv_ssh_exec("echo $HOME", timeout=5)
                    home = out.strip()
                    rpath = home + remote_path[1:]
                else:
                    rpath = remote_path

                downloaded = []
                for name in selected_names:
                    src = f"{rpath}/{name}"
                    dst = os.path.join(local_dest, name)
                    try:
                        attr = sftp.stat(src)
                        if stat.S_ISDIR(attr.st_mode):
                            os.makedirs(dst, exist_ok=True)
                            self._dv_sftp_get_dir(sftp, src, dst)
                        else:
                            sftp.get(src, dst)
                        downloaded.append(dst)
                    except Exception as e:
                        downloaded.append(f"ERROR:{name}:{e}")
                return downloaded
            finally:
                try:
                    sftp.close()
                except Exception:
                    pass

        self._dv_dl_btn.setEnabled(False)
        self._dv_dl_progress.setText("Downloading...")

        w = AnalysisWorker(_do_download)
        w.finished.connect(
            lambda paths: self._dv_on_download_done(paths, auto_import))
        w.error.connect(lambda e: (
            QMessageBox.critical(self, "Server — download error", e),
            self._dv_dl_btn.setEnabled(True),
            self._dv_dl_progress.setText("Download failed.")))
        self._hpc_workers.append(w)
        w.start()

    def _dv_sftp_get_dir(self, sftp, remote_dir: str, local_dir: str):
        """Recursively download a remote directory via SFTP."""
        for attr in sftp.listdir_attr(remote_dir):
            rsrc = f"{remote_dir}/{attr.filename}"
            ldst = os.path.join(local_dir, attr.filename)
            if stat.S_ISDIR(attr.st_mode):
                os.makedirs(ldst, exist_ok=True)
                self._dv_sftp_get_dir(sftp, rsrc, ldst)
            else:
                sftp.get(rsrc, ldst)

    def _dv_on_download_done(self, paths: list, auto_import: bool):
        """Callback after download: show summary, optionally import ranking."""
        self._dv_dl_btn.setEnabled(True)
        ok  = [p for p in paths if not p.startswith('ERROR:')]
        err = [p for p in paths if p.startswith('ERROR:')]
        self._dv_dl_progress.setText(
            f"✓ Downloaded {len(ok)} item(s)"
            + (f" | {len(err)} error(s)" if err else ''))
        self._status.showMessage(
            f"✓ Server: {len(ok)} result(s) downloaded to "
            f"{self._dv_local_dest.text()}")

        if auto_import and ok:
            imported = 0
            for local_path in ok:
                imported += self._dv_import_ranking_from_path(local_path)
            if imported:
                self._af3_show_ranking()
                self._status.showMessage(
                    f"✓ Imported {imported} result(s) into Ranking tab")

    def _dv_import_ranking_from_path(self, local_path: str) -> int:
        """Scan a local directory (or file) for AF3 result JSON files
        and update self.af3_jobs ipTM/pLDDT values.
        Returns the number of jobs updated."""
        count = 0
        candidates = []
        if os.path.isdir(local_path):
            for root, _, files in os.walk(local_path):
                for fn in files:
                    if fn.endswith('.json'):
                        candidates.append(os.path.join(root, fn))
        elif local_path.endswith('.json'):
            candidates.append(local_path)

        for fp in candidates:
            try:
                with open(fp, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                # AF3 result files have 'iptm' at top level or under 'summary'
                iptm = (data.get('iptm')
                        or data.get('summary', {}).get('iptm'))
                plddt = (data.get('ptm') or data.get('plddt')
                         or data.get('summary', {}).get('ptm')
                         or data.get('summary', {}).get('plddt'))
                if iptm is None:
                    continue
                name = data.get('name', Path(fp).stem)
                # Match against self.af3_jobs by name
                for job in self.af3_jobs:
                    if job['name'] == name or job['name'] in fp:
                        job['iptm']  = float(iptm)
                        job['plddt'] = float(plddt) if plddt else None
                        job['status'] = 'done'
                        count += 1
                        break
            except Exception:
                pass
        if count:
            self._af3_update_jobs_table()
        return count





# ═══════════════════════════════════════════════════════════════
# EMOJI FONT SETUP — ensures emoji rendering on all platforms
# ═══════════════════════════════════════════════════════════════

def _setup_emoji_font(app):
    """Configure font fallback chain to include emoji support."""
    try:
        if QT_VERSION == 6:
            from PyQt6.QtGui import QFontDatabase
        else:
            from PyQt5.QtGui import QFontDatabase
    except ImportError:
        _apply_text_fallback()
        return

    # Emoji fonts by platform (in order of preference)
    emoji_fonts = [
        'Segoe UI Emoji',        # Windows
        'Noto Color Emoji',      # Linux (apt install fonts-noto-color-emoji)
        'Apple Color Emoji',     # macOS
        'Noto Emoji',            # Linux fallback (B&W)
        'Symbola',               # Universal fallback
    ]

    # Detect which emoji font is available
    try:
        available = QFontDatabase.families() if QT_VERSION == 6 else QFontDatabase().families()
    except Exception:
        _apply_text_fallback()
        return

    emoji_font = None
    for ef in emoji_fonts:
        if ef in available:
            emoji_font = ef
            break

    # Configure the app font with emoji fallback
    if emoji_font:
        base_font = app.font()
        new_font = QFont(f"{base_font.family()}, {emoji_font}", base_font.pointSize())
        app.setFont(new_font)
    else:
        # No emoji font found — replace emojis with ASCII text
        _apply_text_fallback()


# Substitution map: emoji → equivalent ASCII text
_EMOJI_MAP = {
    '🧬': '[DNA]', '📂': '[Open]', '🔍': '[Find]', '🧪': '[Lab]',
    '🏷️': '[Tag]', '🖼️': '[Img]', '📎': '[Clip]', '📋': '[List]',
    '💾': '[Save]', '📊': '[Chart]', '🔬': '[Scope]', '📦': '[Box]',
    '🗑️': '[Del]', '⚙️': '[Gear]', '❓': '[?]', '❌': '[X]',
    '✅': '[OK]', '⏳': '[...]', '🔎': '[Search]', '📖': '[Book]',
    '🎓': '[Learn]', 'ℹ️': '[i]', '📥': '[In]', '📁': '[Dir]',
    '🚀': '[Run]', '🌐': '[Web]', '⚡': '[Bolt]', '🔮': '[Pred]',
    '➕': '[+]', '🗺️': '[Map]', '🎨': '[Art]', '📍': '[Pin]',
    '▶': '[>]', '⌨️': '[Cmd]', '🔄': '[Ref]', '💡': '[Tip]',
    '🇬🇧': 'EN', '🇧🇷': 'PT', '🇪🇸': 'ES', '🇫🇷': 'FR',
    '🇨🇳': 'ZH', '🇯🇵': 'JA',
}

def _apply_text_fallback():
    """Replace emojis in TRANSLATIONS with ASCII equivalents."""
    for lang_dict in TRANSLATIONS.values():
        for key, val in lang_dict.items():
            if isinstance(val, str):
                for emoji, text in _EMOJI_MAP.items():
                    val = val.replace(emoji, text)
                lang_dict[key] = val


# ═══════════════════════════════════════════════════════════════
# STARTUP DEPENDENCY CHECK
# ═══════════════════════════════════════════════════════════════
def _check_dependencies_at_startup():
    """Check for missing or outdated packages and show a console warning
    (or a Qt dialog if Qt is available) before the main window opens.
    Does NOT block startup — informational only."""
    import importlib

    CHECKS = [
        # (import_name, pip_name, min_ver_tuple, required)
        ("matplotlib", "matplotlib>=3.5",  (3, 5),  True),
        ("numpy",      "numpy>=1.21",      (1, 21), True),
        ("pyrodigal",  "pyrodigal>=2.0",   (2, 0),  False),
        ("paramiko",   "paramiko>=2.9",    (2, 9),  False),
        ("scipy",      "scipy>=1.7",       (1, 7),  False),
    ]

    missing_req  = []
    missing_opt  = []
    outdated     = []

    for imp, pip_name, min_ver, required in CHECKS:
        try:
            m = importlib.import_module(imp)
            ver_str = getattr(m, '__version__', '0')
            ver_tup = tuple(int(x) for x in ver_str.split('.')[:2]
                            if x.isdigit())
            if ver_tup and ver_tup < min_ver:
                outdated.append(
                    f"{imp} {ver_str} (need >={'.'.join(map(str,min_ver))})"
                    f" — pip install \"{pip_name}\"")
        except ImportError:
            if required:
                missing_req.append(pip_name)
            else:
                missing_opt.append(pip_name)

    # Qt is already loaded at this point — use it for the dialog
    if missing_req or outdated:
        lines = ["ppigFinder dependency warning:\n"]
        if missing_req:
            lines.append("REQUIRED (app may crash without these):")
            for p in missing_req:
                lines.append(f"  pip install {p}")
        if outdated:
            lines.append("\nOUTDATED (update recommended):")
            for p in outdated:
                lines.append(f"  {p}")
        lines.append(
            "\nRun install_ppigfinder.py to fix all issues automatically.")
        try:
            from PyQt6.QtWidgets import QMessageBox
        except ImportError:
            from PyQt5.QtWidgets import QMessageBox
        mb = QMessageBox()
        mb.setIcon(QMessageBox.Icon.Warning
                   if hasattr(QMessageBox, 'Icon')
                   else QMessageBox.Warning)
        mb.setWindowTitle("ppigFinder — Missing dependencies")
        mb.setText("\n".join(lines))
        mb.exec() if hasattr(mb, 'exec') and callable(mb.exec) else mb.exec_()

    if missing_opt:
        print(
            f"[ppigFinder] Optional packages not installed "
            f"(some features disabled): {', '.join(missing_opt)}\n"
            f"  Install with: pip install {' '.join(missing_opt)}")


if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setApplicationName('ppigFinder')
    app.setApplicationDisplayName('ppigFinder — Protein-Protein Interaction Genomic Finder')
    app.setApplicationVersion('2.0')
    app.setStyle('Fusion')
    _setup_emoji_font(app)
    _check_dependencies_at_startup()
    window = ppigFinderApp()
    window.show()
    sys.exit(app.exec() if QT_VERSION == 6 else app.exec_())
