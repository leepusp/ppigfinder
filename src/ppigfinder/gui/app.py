#!/usr/bin/env python3
"""
ppigFinder — Protein-Protein Interaction Genomic Finder
========================================================
Version  : 1.01
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
  • AF3 results analysis: PAE heatmap (ChimeraX colour scheme), pLDDT
    plots, ipTM/ptm scoring, inter-chain contact detection
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
    https://github.com/<your-org>/ppigfinder

CHANGELOG
---------
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

import sys
import os
import re
import csv
import json
import math
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
        QMenu, QToolBar, QScrollArea,
        QFrame, QSizePolicy, QAbstractItemView,
        QInputDialog,
    )
    from PyQt6.QtCore import (
        Qt, QTimer, QThread, pyqtSignal,
    )
    from PyQt6.QtGui import (
        QPainter, QPen, QBrush, QColor, QFont,
        QPolygonF,
    )
    from PyQt6.QtCore import QPointF
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
        QMenu, QToolBar, QScrollArea,
        QFrame, QSizePolicy, QAbstractItemView,
        QInputDialog,
    )
    from PyQt5.QtCore import (
        Qt, QTimer, QThread, pyqtSignal,
        QPointF,
    )
    from PyQt5.QtGui import (
        QPainter, QPen, QBrush, QColor, QFont,
        QPolygonF,
    )
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

try:
    import matplotlib
    if QT_VERSION == 6:
        matplotlib.use('Qt6Agg')
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
    else:
        matplotlib.use('Qt5Agg')
        from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except Exception:
    try:
        import matplotlib
        matplotlib.use('Agg')
        from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
        import matplotlib.pyplot as plt
        MATPLOTLIB_AVAILABLE = True
    except ImportError:
        MATPLOTLIB_AVAILABLE = False
        FigureCanvas = None

try:
    import numpy as np
    NUMPY_AVAILABLE = True
except ImportError:
    np = None  # type: ignore[assignment]
    NUMPY_AVAILABLE = False



# ═══════════════════════════════════════════════════════════════
# EXTERNAL BACKEND DETECTION
# ═══════════════════════════════════════════════════════════════

def detect_backends():
    """Detect installed external tools (including WSL)."""
    backends = {}
    for tool, cmd in [('blast+', 'blastp'), ('hmmer3', 'hmmsearch')]:
        path = shutil.which(cmd)
        if path:
            try:
                r = subprocess.run([cmd, '-version'] if tool == 'blast+' else [cmd, '-h'],
                                   capture_output=True, text=True, timeout=5)
                version = r.stdout.split('\n')[0] if r.stdout else 'detectado'
                backends[tool] = {'path': path, 'version': version, 'available': True, 'wsl': False}
            except Exception:
                backends[tool] = {'path': path, 'version': '?', 'available': True, 'wsl': False}
        else:
            # Try WSL on Windows
            wsl_found = False
            if os.name == 'nt':
                try:
                    r = subprocess.run(['wsl', 'bash', '-c', f'{cmd} -h'],
                                       capture_output=True, text=True, timeout=10)
                    if r.returncode == 0 or 'Usage' in (r.stdout + r.stderr):
                        backends[tool] = {'path': f'wsl bash -c {cmd}', 'version': 'via WSL',
                                          'available': True, 'wsl': True}
                        wsl_found = True
                except Exception:
                    pass
            if not wsl_found:
                backends[tool] = {'path': None, 'version': None, 'available': False, 'wsl': False}
    return backends

BACKENDS = detect_backends()

# Pyrodigal detection (Python module, not external tool)
BACKENDS['pyrodigal'] = {
    'path': 'pyrodigal (Python module)',
    'version': pyrodigal.__version__ if PYRODIGAL_AVAILABLE else None,
    'available': PYRODIGAL_AVAILABLE,
    'wsl': False,
}


# ═══════════════════════════════════════════════════════════════
# MODULE 0: SNAPGENE (.dna) and GENBANK (.gb/.gbk) SUPPORT
# Pure-stdlib implementation — no external dependencies
# ═══════════════════════════════════════════════════════════════

import struct as _struct
import io     as _io

# ───────────────────────────────────────────────────────────────
# 0-A  SNAPGENE READER (.dna)
# ───────────────────────────────────────────────────────────────

def _snapgene_read_packets(data: bytes):
    """Iterate TLV packets from a SnapGene .dna file."""
    buf = _io.BytesIO(data)
    cookie = buf.read(1)
    if cookie != b'\x09':
        raise ValueError(
            "File not recognised as SnapGene (.dna).\n"
            "Check whether it is corrupted or an unsupported old version."
        )
    while True:
        hdr = buf.read(5)
        if len(hdr) < 5:
            break
        pkt_type = hdr[0]
        pkt_len  = _struct.unpack('>I', hdr[1:5])[0]
        payload  = buf.read(pkt_len)
        yield pkt_type, payload


def parse_snapgene_dna(filepath: str) -> dict:
    """
    Read a SnapGene (.dna) file and return:
      sequence  (str)  : DNA sequence in uppercase
      topology  (str)  : 'circular' | 'linear'
      name      (str)  : file name without extension
      features  (list) : [{name, type, start(0-based), end, strand(+/-), color}]
      primers   (list) : [{name, start(0-based), end, strand(+/-)}]
      notes     (dict) : free metadata from SnapGene
    """
    from xml.etree import ElementTree as ET

    with open(filepath, 'rb') as fh:
        data = fh.read()

    result = {
        'sequence': '', 'topology': 'linear',
        'name': Path(filepath).stem,
        'features': [], 'primers': [], 'notes': {},
    }

    for pkt_type, payload in _snapgene_read_packets(data):

        if pkt_type == 0 and len(payload) >= 2:
            flags = payload[0]
            result['topology'] = 'circular' if (flags & 0x01) else 'linear'
            result['sequence'] = payload[1:].decode('ascii', errors='replace').upper()

        elif pkt_type == 8:
            try:
                xml = ET.fromstring(payload.decode('utf-8', errors='replace'))
                for feat in xml.findall('.//Feature'):
                    name  = feat.get('name', 'unknown')
                    ftype = feat.get('type', 'misc_feature')
                    color = feat.get('color', '#aaaaaa')
                    direc = feat.get('directionality', '1')
                    strand = '+' if direc in ('1', 'forward') else '-'
                    for seg in feat.findall('Segment'):
                        rng = seg.get('range', '')
                        if '-' in rng:
                            s, e = rng.split('-', 1)
                            try:
                                result['features'].append({
                                    'name': name, 'type': ftype, 'color': color,
                                    'start': int(s) - 1,
                                    'end':   int(e),
                                    'strand': strand,
                                })
                            except ValueError:
                                pass
            except ET.ParseError:
                pass

        elif pkt_type == 6:
            try:
                xml = ET.fromstring(payload.decode('utf-8', errors='replace'))
                for pr in xml.findall('.//Primer'):
                    pname   = pr.get('name', 'primer')
                    pstrand = '+' if pr.get('templateStrand', 'sense') == 'sense' else '-'
                    for bind in pr.findall('BindingSite'):
                        rng = bind.get('location', '')
                        if '-' in rng:
                            s, e = rng.split('-', 1)
                            try:
                                result['primers'].append({
                                    'name': pname,
                                    'start': int(s) - 1, 'end': int(e),
                                    'strand': pstrand,
                                })
                            except ValueError:
                                pass
            except ET.ParseError:
                pass

        elif pkt_type == 10:
            try:
                xml = ET.fromstring(payload.decode('utf-8', errors='replace'))
                for child in xml:
                    result['notes'][child.tag] = (child.text or '').strip()
            except ET.ParseError:
                pass

    return result


# ───────────────────────────────────────────────────────────────
# 0-B  SNAPGENE WRITER (.dna)
# ───────────────────────────────────────────────────────────────

def _sg_packet(pkt_type: int, payload: bytes) -> bytes:
    """Serialise a SnapGene TLV packet."""
    return bytes([pkt_type]) + _struct.pack('>I', len(payload)) + payload


def write_snapgene_dna(filepath: str, sequence: str, features: list,
                       primers: list = None, topology: str = 'linear',
                       name: str = '', notes: dict = None):
    """
    Write a SnapGene-compatible binary .dna file.

    Parameters
    ----------
    sequence  : str   — DNA in uppercase
    features  : list  — [{name, type, start(0-based), end, strand(+/-), color}]
    primers   : list  — [{name, start(0-based), end, strand(+/-)}]  (optional)
    topology  : 'linear' | 'circular'
    name      : str   — record name (written to Notes packet)
    notes     : dict  — additional metadata
    """
    from xml.etree.ElementTree import Element, SubElement, tostring as _tostr

    buf = _io.BytesIO()
    buf.write(b'\x09')   # magic cookie

    # ── Packet 0: sequence ──────────────────────────────────────
    flags = 0x01 if topology == 'circular' else 0x00
    seq_payload = bytes([flags]) + sequence.lower().encode('ascii')
    buf.write(_sg_packet(0, seq_payload))

    # ── Packet 8: features ──────────────────────────────────────
    if features:
        root_xml = Element('Features')
        for feat in features:
            direc = '1' if feat.get('strand', '+') == '+' else '2'
            f_el = SubElement(root_xml, 'Feature',
                              name=feat.get('name', 'feature'),
                              type=feat.get('type', 'misc_feature'),
                              directionality=direc,
                              color=feat.get('color', '#aaaaaa'),
                              swappedSegmentNumbering='0',
                              allowSegmentOverlaps='0')
            # SnapGene uses 1-based closed ranges
            s = feat['start'] + 1
            e = feat['end']
            SubElement(f_el, 'Segment', range=f"{s}-{e}",
                       name=feat.get('name', ''),
                       color=feat.get('color', '#aaaaaa'),
                       type='standard')
        xml_bytes = _tostr(root_xml, encoding='unicode').encode('utf-8')
        buf.write(_sg_packet(8, xml_bytes))

    # ── Packet 6: primers ───────────────────────────────────────
    primers = primers or []
    if primers:
        root_xml = Element('Primers')
        for pr in primers:
            tpl_strand = 'sense' if pr.get('strand', '+') == '+' else 'antisense'
            pr_el = SubElement(root_xml, 'Primer',
                               name=pr.get('name', 'primer'),
                               templateStrand=tpl_strand)
            s = pr['start'] + 1
            e = pr['end']
            SubElement(pr_el, 'BindingSite', location=f"{s}-{e}")
        xml_bytes = _tostr(root_xml, encoding='unicode').encode('utf-8')
        buf.write(_sg_packet(6, xml_bytes))

    # ── Packet 10: notes ────────────────────────────────────────
    notes = notes or {}
    if name:
        notes.setdefault('Description', name)
    notes.setdefault('CustomMapLabel', name or 'ORF Pipeline v20')
    notes.setdefault('ConfirmedExperimentally', '0')
    notes.setdefault('Created', datetime.now().strftime('%Y-%m-%d'))
    notes.setdefault('CreatedBy', 'ORF Secretion Pipeline v20')
    root_xml = Element('Notes')
    for k, v in notes.items():
        el = SubElement(root_xml, k)
        el.text = str(v)
    xml_bytes = _tostr(root_xml, encoding='unicode').encode('utf-8')
    buf.write(_sg_packet(10, xml_bytes))

    with open(filepath, 'wb') as fh:
        fh.write(buf.getvalue())


# ───────────────────────────────────────────────────────────────
# 0-C  GENBANK READER (.gb / .gbk / .genbank)
# ───────────────────────────────────────────────────────────────

def parse_genbank(filepath: str) -> dict:
    """
    Pure-stdlib parser for GenBank flat files (INSDC format).

    Returns the same dict as parse_snapgene_dna:
      sequence, topology, name, features, primers (empty), notes
    """
    result = {
        'sequence': '', 'topology': 'linear',
        'name': Path(filepath).stem,
        'features': [], 'primers': [], 'notes': {},
    }

    # Color palette by feature type (matching SnapGene defaults)
    _TYPE_COLORS = {
        'CDS':          '#99ccff',
        'gene':         '#ffcc99',
        'mRNA':         '#ff9999',
        'rRNA':         '#ccff99',
        'tRNA':         '#ffff99',
        'ncRNA':        '#ff99ff',
        'regulatory':   '#99ffff',
        'rep_origin':   '#ff9966',
        'misc_feature': '#cccccc',
        'promoter':     '#ffcc00',
        'terminator':   '#cc99ff',
        'primer_bind':  '#ff66cc',
    }

    with open(filepath, 'r', errors='replace') as fh:
        content = fh.read()

    # ── Split into records (multi-record support) ───────────────
    # Use only the first record
    record_text = content.split('//')[0]

    # ── LOCUS ───────────────────────────────────────────────────
    locus_m = re.search(r'^LOCUS\s+(\S+)\s+.*?(circular|linear)', record_text,
                        re.MULTILINE | re.IGNORECASE)
    if locus_m:
        result['name']     = locus_m.group(1)
        result['topology'] = locus_m.group(2).lower()

    # ── DEFINITION / ACCESSION for notes ───────────────────────
    for tag in ('DEFINITION', 'ACCESSION', 'VERSION', 'ORGANISM'):
        m = re.search(rf'^{tag}\s+(.+?)(?=\n[A-Z])', record_text,
                      re.MULTILINE | re.DOTALL)
        if m:
            result['notes'][tag] = ' '.join(m.group(1).split())

    # ── FEATURES ────────────────────────────────────────────────
    feat_block_m = re.search(r'^FEATURES\s+.*?\n(.*?)^(?:ORIGIN|CONTIG)',
                              record_text, re.MULTILINE | re.DOTALL)
    if feat_block_m:
        feat_text = feat_block_m.group(1)

        # Split into individual feature entries
        # Each starts at col 5 with a keyword, qualifiers at col 21
        feat_entries = re.split(r'\n(?=     \S)', feat_text)

        for entry in feat_entries:
            lines = entry.split('\n')
            if not lines:
                continue
            first = lines[0].strip()
            if not first:
                continue

            # Feature type is the first token
            parts = first.split()
            if not parts:
                continue
            ftype = parts[0]
            loc_str = parts[1] if len(parts) > 1 else ''

            # Collect remaining qualifier lines
            for ln in lines[1:]:
                s = ln.strip()
                if s and not s.startswith('/'):
                    loc_str += s   # continuation of location
                else:
                    break

            # Parse location string → start, end, strand
            strand = '-' if loc_str.startswith('complement') else '+'
            # Strip complement(...) and join(...)
            loc_clean = re.sub(r'(complement|join|order)\(', '', loc_str).rstrip(')')
            # Grab all numeric pairs
            ranges = re.findall(r'(\d+)\.\.(\d+)', loc_clean)
            if not ranges:
                # single position
                m_single = re.search(r'(\d+)', loc_clean)
                if m_single:
                    pos = int(m_single.group(1))
                    ranges = [(str(pos), str(pos))]

            if not ranges:
                continue

            # For joined features keep first + last segment
            start = int(ranges[0][0]) - 1   # 0-based
            end   = int(ranges[-1][1])

            # Collect qualifiers
            qual_text = '\n'.join(lines[1:])
            name_m = re.search(r'/(?:gene|locus_tag|product|label)="([^"]+)"', qual_text)
            feat_name = name_m.group(1) if name_m else ftype

            color = _TYPE_COLORS.get(ftype, '#cccccc')

            if ftype == 'primer_bind':
                result['primers'].append({
                    'name': feat_name,
                    'start': start, 'end': end,
                    'strand': strand,
                })
            else:
                result['features'].append({
                    'name': feat_name,
                    'type': ftype,
                    'color': color,
                    'start': start,
                    'end':   end,
                    'strand': strand,
                })

    # ── ORIGIN (sequence) ───────────────────────────────────────
    origin_m = re.search(r'^ORIGIN\s*\n(.*)', record_text,
                         re.MULTILINE | re.DOTALL)
    if origin_m:
        raw_seq = origin_m.group(1)
        result['sequence'] = re.sub(r'[^a-zA-Z]', '', raw_seq).upper()

    return result


# ───────────────────────────────────────────────────────────────
# 0-D  GENBANK WRITER (.gb)
# ───────────────────────────────────────────────────────────────

def write_genbank(filepath: str, sequence: str, orfs: list,
                  sg_features: list = None, name: str = 'sequence',
                  topology: str = 'linear'):
    """
    Write a GenBank flat-file (.gb) containing:
      - CDS features for each analysed ORF
      - misc_feature for each imported SnapGene feature (if any)
    """
    sl = len(sequence)
    now = datetime.now().strftime('%d-%b-%Y').upper()
    topo_str = 'circular' if topology == 'circular' else 'linear  '
    mol = 'DNA'

    lines = []

    # ── LOCUS ──────────────────────────────────────────────────
    lines.append(f"LOCUS       {name[:16]:<16} {sl:>9} bp    {mol}     {topo_str} {now}")
    lines.append(f"DEFINITION  {name} — exported by ORF Secretion Pipeline v20.")
    lines.append( "ACCESSION   .")
    lines.append( "VERSION     .")
    lines.append( "KEYWORDS    .")
    lines.append( "SOURCE      .")
    lines.append( "  ORGANISM  .")
    lines.append( "            .")
    lines.append( "FEATURES             Location/Qualifiers")

    def _loc(start, end, strand):
        """Format a GenBank location string (1-based, closed)."""
        s = start + 1
        e = end
        loc = f"{s}..{e}"
        return f"complement({loc})" if strand == '-' else loc

    def _wrap_qual(key, val, indent=21):
        """Wrap a qualifier value at col 80."""
        prefix = ' ' * indent
        full = f'{prefix}/{key}="{val}"'
        out = []
        while len(full) > 79:
            out.append(full[:79])
            full = prefix + full[79:]
        out.append(full)
        return '\n'.join(out)

    feat_indent = '     '

    # ── SnapGene features (preserved from original file) ────────
    for feat in (sg_features or []):
        ftype = feat.get('type', 'misc_feature')
        lines.append(f"{feat_indent}{ftype:<16}{_loc(feat['start'], feat['end'], feat.get('strand', '+'))}")
        lines.append(_wrap_qual('label', feat.get('name', ftype)))
        lines.append(_wrap_qual('note', f"Imported from SnapGene; color: {feat.get('color','#aaaaaa')}"))

    # ── ORFs as CDS ─────────────────────────────────────────────
    for i, orf in enumerate(orfs):
        loc = _loc(orf['start'], orf['end'], orf['strand'])
        lines.append(f"{feat_indent}{'CDS':<16}{loc}")
        gene_name = orf.get('gene_name') or f"orf{i+1}"
        lines.append(_wrap_qual('gene', gene_name))
        prot = orf['protein'].rstrip('*')
        lines.append(_wrap_qual('product',
            orf.get('putative_function') or orf.get('observation') or 'hypothetical protein'))
        lines.append(_wrap_qual('protein_id', f"ORF{i+1}"))
        domains = '; '.join(d['domain'] for d in orf.get('domains', []))
        if domains:
            lines.append(_wrap_qual('note', f"HMM domains: {domains}"))
        # Translate qualifier — wrap at 60 aa per line inside the qualifier
        prot_chunks = [prot[j:j+60] for j in range(0, len(prot), 60)]
        lines.append(f"{'':21}/translation=\"{prot_chunks[0]}")
        for chunk in prot_chunks[1:]:
            lines.append(f"{'':21}{chunk}")
        lines[-1] += '"'

    # ── ORIGIN ──────────────────────────────────────────────────
    lines.append("ORIGIN")
    seq_lc = sequence.lower()
    for pos in range(0, sl, 60):
        chunk = seq_lc[pos:pos+60]
        # Split into groups of 10
        groups = ' '.join(chunk[i:i+10] for i in range(0, len(chunk), 10))
        lines.append(f"{pos+1:>9} {groups}")
    lines.append("//")

    with open(filepath, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')




class AdvancedORFAnalyzer:

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
        dna_seq = dna_seq.upper()
        protein = []
        for i in range(0, len(dna_seq), 3):
            codon = dna_seq[i:i+3]
            if len(codon) == 3:
                aa = self.CODON_TABLE.get(codon, 'X')
                protein.append(aa)
                if aa == '*': break
        return ''.join(protein)

    def reverse_complement(self, seq):
        comp = {'A':'T','T':'A','G':'C','C':'G','a':'t','t':'a','g':'c','c':'g'}
        return ''.join(comp.get(b, 'N') for b in reversed(seq))

    def gc_content(self, seq):
        seq = seq.upper()
        return (seq.count('G') + seq.count('C')) / len(seq) * 100 if seq else 0


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
# MODULE I: INTERNATIONALISATION (i18n) — v1.01 — English only
# ═══════════════════════════════════════════════════════════════

TRANSLATIONS = {
    # ── English ─────────────────────────────────────────────────
    'en': {
        'app_title':        '🧬 ppigFinder v1.01 — PPI Genomic Finder',
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
        'about':            'ℹ️ About',
        # Toolbar
        'btn_open':         '📂 Load a genome file',
        'btn_translate_genome': '🧬 Translate genome',
        'btn_pyrodigal':    '🧪 Pyrodigal',
        'btn_automatic':    '⚙️ Automatic',
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
        'af3_remove':           '🗑️ Remove',
        'af3_clear_all':        '🗑️ Clear All',
        'af3_jobs_frame':       '⚡ Generate Jobs',
        'af3_neighbors':        'Neighbors:',
        'af3_mode':             'Mode:',
        'af3_generate':         '⚡ Generate',
        'af3_export_json':      '💾 Export AF3 JSON',
        'af3_export_json_single': '📄 Individual JSONs',
        'af3_export_json_batch':  '📦 Batch JSON',
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
ppigFinder v1.01  (Server Edition v2)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
User Manual

WHAT THIS APP DOES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
This application is a local desktop pipeline for exploratory analysis
of bacterial genomes. It combines ORF detection, HMM annotation,
BLAST-style homology search, genome map navigation, AlphaFold job
generation, server submission, and AlphaFold result inspection in a
single interface.

MAIN WORKFLOW
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
1. Load a genome file
   Supported inputs: FASTA, multi-FASTA, GenBank (.gb/.gbk),
   and SnapGene (.dna).

2. Translate genome
   The toolbar button "Translate genome" opens two methods:
   • Pyrodigal: model-based prokaryotic gene prediction
   • Automatic: six-frame ORF detection using start/stop codons
     and the size filters defined in Parameters

3. Inspect ORFs
   The ORF table lists genomic position, frame, strand, size, GC,
   HMM hits, AF3 partner metrics, and user notes. Clicking any ORF in
   the table selects it and centers the genome map on that ORF.

4. Annotate with HMM profiles
   Add one or multiple .hmm profiles, run the search, then use
   "Annotate HMM" to color matching ORFs directly on the genome map.

5. Run sequence similarity searches
   The BLAST tab accepts a protein query and compares it against the
   translated ORFs in the loaded genome. ORF links in the results are
   clickable and center the genome map on the selected ORF.

6. Build AlphaFold jobs
   The AlphaFold tab allows pair, trimer, all-vs-all, homodimer,
   and custom multimer job generation. Jobs can be exported as
   AlphaFold Server JSON or ColabFold FASTA.

7. Analyze AlphaFold predictions
   The AlphaFold Analysis tab loads AF3 result folders, populates a
   ranking table, and shows embedded PAE and pLDDT plots inside the
   application. Selecting a job in the table automatically centers the
   genome map on the first ORF parsed from that job name.

GENOME MAP NAVIGATION
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
• Ctrl + mouse wheel: zoom at cursor position
• Toolbar − / + buttons: zoom toward the visible center
• Shift + drag: horizontal pan
• Search box: jump directly to ORF number or protein substring
• Clicking a map arrow: select the ORF and show its details

The zoom label in the toolbar always reflects the current map zoom.

ALPHAFOLD ANALYSIS TAB
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
The AlphaFold Analysis tab is designed for in-app inspection.
Plots are embedded in the scroll area and can be exported to PDF with
"Export plots PDF". The table reports job name, chain composition,
ipTM, pTM, mean pLDDT, ranking score, best inter-chain PAE, and best
contact pair description.

HELP / DIALOGS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Manual and Tutorial windows open as modal dialogs so they remain tied
to the main application window.

KEY SHORTCUTS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Ctrl + scroll       Zoom the genome map
Shift + drag        Pan the genome map
Right-click ORF     Context menu for copy / annotation / color
Right-click AF3 job Remove selected generated jobs
""",
    },

    'tutorial': {
        'en': """\
ppigFinder v1.01  (Server Edition v2)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Step-by-Step Tutorial

STEP 1 — Load a genome file
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Use "Load a genome file" in the toolbar. FASTA, GenBank, SnapGene,
and multi-FASTA are supported. After loading, the genome information
panel is updated and the genomic map becomes available.

STEP 2 — Translate the genome
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Click "Translate genome" and choose one of the two modes:
• Pyrodigal: recommended when you want a gene-caller style prediction
• Automatic: recommended when you want to inspect all possible ORFs
  from start/stop codons under the current parameter limits

STEP 3 — Explore the ORF table
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
The ORF table is synchronized with the map. Clicking any ORF row:
• selects the ORF
• updates the DNA / Protein / Domains tabs
• centers the genomic map on that ORF

STEP 4 — Add and search HMM profiles
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Open the HMM tab, add one or multiple .hmm files, then run
"Search All ORFs". After the search, click "Annotate HMM" in the
toolbar to transfer those hits to the main ORF table and genome map.

STEP 5 — Run BLAST on genome ORFs
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Paste or load a protein sequence in the BLAST Query tab and run the
search. ORF links in the BLAST results are clickable and recenter the
map automatically.

STEP 6 — Build AlphaFold jobs
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Open the AlphaFold tab and add ORFs using the selection controls.
Generate pairwise or multimer jobs, then export them as:
• AlphaFold Server JSON
• ColabFold FASTA

STEP 7 — Load AlphaFold predictions for analysis
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Open the AlphaFold Analysis tab and click "Load AF3 results folder".
The application scans the selected folder, parses the AF3 outputs, and
creates a sortable job table.

STEP 8 — Inspect prediction quality
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Selecting a row in the AlphaFold Analysis table will:
• embed the PAE plot inside the application
• embed the pLDDT plot inside the application
• center the genome map on the first ORF found in the job

STEP 9 — Use map navigation effectively
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
• Ctrl + mouse wheel zooms at the cursor
• Toolbar − / + zoom buttons zoom around the visible center
• Shift + drag pans horizontally
• The search box jumps directly to an ORF

STEP 10 — Export figures and save the project
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Use "Export Map PDF" for the genome map and "Export plots PDF" in
AlphaFold Analysis for PAE / pLDDT figures. Save the project when you
want to preserve the current state of genome, ORFs, HMM hits, and AF3
metadata for later reopening.
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
        new_level = max(0.5, min(200.0, level))

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

class ppigFinderApp(QMainWindow):

    def __init__(self):
        super().__init__()
        self.setWindowTitle('🧬 ppigFinder v1.01 — Protein-Protein Interaction Genomic Finder')
        self.resize(1550, 980)

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

    # ═══════════════════════════════════════════════════════════
    # UI SETUP
    # ═══════════════════════════════════════════════════════════

    def _setup_ui(self):
        self._create_menus()
        self._create_toolbar()
        self._create_central()
        self._create_statusbar()

    # ─── MENUS ─────────────────────────────────────────────────

    def _create_menus(self):
        mb = self.menuBar()

        # File
        fm = mb.addMenu(t('menu_file'))
        fm.addAction(t('open_fasta'), self.load_fasta)
        fm.addAction(t('open_multifasta'), self.load_multi_fasta)
        fm.addAction(t('open_snapgene'), self.load_snapgene)
        fm.addAction(t('open_genbank'), self.load_genbank)
        fm.addAction(t('load_hmm'), self.load_hmm)
        fm.addSeparator()
        fm.addAction(t('save_project'), self.save_project)
        fm.addAction('💾 Save Project As (full copy)...', self.save_project_as)
        fm.addAction(t('open_project'), self.load_project)
        fm.addSeparator()
        fm.addAction(t('save_orfs_fasta'), self.save_fasta)
        fm.addAction(t('save_report_tsv'), self.save_report_tsv)
        fm.addSeparator()
        fm.addAction(t('export_genbank'), self.export_genbank)
        fm.addAction(t('export_snapgene'), self.export_snapgene)
        fm.addSeparator()
        fm.addAction(t('quit'), self.close)

        # Parameters
        pm = mb.addMenu(t('menu_params'))
        pm.addAction("🧬 ORF Analysis Parameters...", self._show_orf_params)
        pm.addAction(t('blast_params'), self._show_blast_params)
        pm.addAction(t('hmm_params'), self._show_hmm_params)

        # Help
        hm = mb.addMenu(t('menu_help'))
        hm.addAction(t('manual'), self._show_manual)
        hm.addAction(t('tutorial'), self._show_tutorial)
        hm.addSeparator()
        hm.addAction(t('about'), self._show_about)

    # ─── TOOLBAR ───────────────────────────────────────────────

    def _create_toolbar(self):
        tb = QToolBar("Main Toolbar")
        tb.setMovable(False)
        self.addToolBar(tb)

        self._btn_open = QPushButton(t('btn_open'))
        self._btn_open.clicked.connect(self.load_fasta)
        self._btn_open.setToolTip(t('tip_open'))
        tb.addWidget(self._btn_open)

        # Translate genome dropdown button
        self._btn_translate = QPushButton(t('btn_translate_genome'))
        translate_menu = QMenu(self._btn_translate)
        
        # Pyrodigal option with description
        pyrodigal_action = translate_menu.addAction(
            f"{t('btn_pyrodigal')} — {t('desc_pyrodigal')}")
        pyrodigal_action.triggered.connect(self.analyze_orfs_pyrodigal)
        
        # Automatic option with description  
        automatic_action = translate_menu.addAction(
            f"{t('btn_automatic')} — {t('desc_automatic')}")
        automatic_action.triggered.connect(self.analyze_orfs)
        
        self._btn_translate.setMenu(translate_menu)
        self._btn_translate.setToolTip(
            "Choose gene prediction method:\n"
            f"• {t('btn_pyrodigal')}: {t('desc_pyrodigal')}\n"
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
                'AF3','Partner','ipTM','PAE_inter','Contact_region','User_note']
        self._orf_table.setColumnCount(len(cols))
        self._orf_table.setHorizontalHeaderLabels(cols)
        self._orf_table.setSelectionBehavior(SelectRows)
        self._orf_table.setSelectionMode(SingleSelection)
        self._orf_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers
                                         if QT_VERSION == 6 else QAbstractItemView.NoEditTriggers)
        self._orf_table.horizontalHeader().setStretchLastSection(True)
        self._orf_table.setAlternatingRowColors(True)
        self._orf_table.selectionModel().selectionChanged.connect(self._on_orf_table_select)
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

        layout.addWidget(self._orf_table, stretch=1)

        return w

    # ─── RIGHT PANEL: Tabs ────────────────────────────────────

    def _create_right_panel(self):
        self._tabs = QTabWidget()
        self._tabs.setMinimumWidth(450)

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
                           (t('af3_remove'), self._af3_remove_orf),
                           (t('af3_clear_all'), self._af3_clear_all)]:
            b = QPushButton(text); b.clicked.connect(slot); sb.addWidget(b)
            # Add tooltips
            if 'add_sel' in slot.__name__:
                b.setToolTip(t('tip_af3_add_sel'))
            elif 'add_hmm' in slot.__name__:
                b.setToolTip(t('tip_af3_add_hmm'))
            elif 'remove' in slot.__name__:
                b.setToolTip(t('tip_af3_remove'))
            elif 'clear_all' in slot.__name__:
                b.setToolTip(t('tip_af3_clear_all'))
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
        self._af3_sel_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers
                                             if QT_VERSION == 6 else QAbstractItemView.NoEditTriggers)
        sel_layout.addWidget(self._af3_sel_table)
        splitter.addWidget(sel_widget)

        # ═══ SECTION 2: Job Generation + Jobs Table ═══
        job_widget = QWidget()
        jl = QVBoxLayout(job_widget)
        jl.setContentsMargins(4, 2, 4, 2)
        jl.setSpacing(2)

        jb = QHBoxLayout()
        jb.setSpacing(3)
        jb.addWidget(QLabel("Neighbors:"))
        self._af3_nb_spin = QSpinBox(); self._af3_nb_spin.setRange(1,15); self._af3_nb_spin.setValue(5)
        self._af3_nb_spin.setMaximumWidth(50)
        jb.addWidget(self._af3_nb_spin)
        jb.addWidget(QLabel("Mode:"))
        self._af3_mode_combo = QComboBox()
        self._af3_mode_combo.addItems(["Pares (hit vs vizinho)","Pares + Homodímeros",
            "Trímeros (hit + 2 vizinhos)","All vs All (neighborhood)",
            "Hits HMM entre si","Hit vs all selected","Homodímero (hit vs si mesmo)"])
        self._af3_mode_combo.setMinimumWidth(180)
        jb.addWidget(self._af3_mode_combo)
        for text, slot in [(t('af3_generate'), self._af3_generate_jobs),
                           (t('af3_export_cf'), self._af3_export_colabfold),
                           ("Ranking", self._af3_show_ranking),
                           ("Clear Jobs", self._af3_clear_jobs)]:
            b = QPushButton(text); b.clicked.connect(slot); jb.addWidget(b)
            # Add tooltips
            if 'generate' in slot.__name__:
                b.setToolTip(t('tip_af3_generate'))
            elif 'export_colabfold' in slot.__name__:
                b.setToolTip(t('tip_af3_export_cf'))
            elif 'ranking' in slot.__name__:
                b.setToolTip(t('tip_af3_ranking'))
            elif 'clear_jobs' in slot.__name__:
                b.setToolTip(t('tip_af3_clear_jobs'))
        
        # AF3 JSON export menu button (individual vs batch)
        export_btn = QPushButton(t('af3_export_json'))
        export_menu = QMenu(export_btn)
        export_menu.addAction(t('af3_export_json_single'), self._af3_export_json)
        export_menu.addAction(t('af3_export_json_batch'), self._af3_export_json_batch)
        export_btn.setMenu(export_menu)
        jb.addWidget(export_btn)
        jb.addStretch()
        jl.addLayout(jb)

        self._af3_jobs_table = QTableWidget()
        self._af3_jobs_table.setColumnCount(7)
        self._af3_jobs_table.setHorizontalHeaderLabels(['Job','Hit','Partner','Residues','Status','ipTM','pLDDT'])
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
        jl.addWidget(self._af3_jobs_table)
        splitter.addWidget(job_widget)

        # ═══ SECTION 3: Custom Jobs + Output ═══
        bottom_widget = QWidget()
        bl = QVBoxLayout(bottom_widget)
        bl.setContentsMargins(4, 2, 4, 2)
        bl.setSpacing(2)

        # ── Custom job header: subunit count control ──
        custom_header = QHBoxLayout()
        custom_header.setSpacing(4)
        custom_header.addWidget(QLabel("⚡ Custom:"))
        custom_header.addWidget(QLabel("Subunidades:"))
        self._custom_n_subunits = QSpinBox()
        self._custom_n_subunits.setRange(1, 11)   # A–K = 11 chains max
        self._custom_n_subunits.setValue(2)
        self._custom_n_subunits.setMaximumWidth(52)
        self._custom_n_subunits.setToolTip(
            "Número de subunidades (cadeias) no complexo.\n"
            "Cada subunidade recebe uma letra (A, B, C … K).\n"
            "Máximo: 11 subunidades.")
        self._custom_n_subunits.valueChanged.connect(self._af3_rebuild_custom_rows)
        custom_header.addWidget(self._custom_n_subunits)

        btn_add_custom = QPushButton("➕ Add")
        btn_add_custom.clicked.connect(self._af3_add_custom_job)
        btn_add_custom.setToolTip(t('tip_af3_add_custom'))
        custom_header.addWidget(btn_add_custom)
        custom_header.addStretch()
        bl.addLayout(custom_header)

        # ── Scroll area that holds the dynamic per-subunit rows ──
        self._custom_scroll = QScrollArea()
        self._custom_scroll.setWidgetResizable(True)
        self._custom_scroll.setFrameShape(QFrame.Shape.NoFrame
                                          if QT_VERSION == 6 else QFrame.NoFrame)
        self._custom_scroll.setMaximumHeight(110)
        self._custom_scroll.setMinimumHeight(60)

        self._custom_rows_widget = QWidget()
        self._custom_rows_layout = QVBoxLayout(self._custom_rows_widget)
        self._custom_rows_layout.setContentsMargins(0, 0, 0, 0)
        self._custom_rows_layout.setSpacing(2)
        self._custom_scroll.setWidget(self._custom_rows_widget)
        bl.addWidget(self._custom_scroll)

        # Internal list: [(QLineEdit_orf, QSpinBox_n), …]  rebuilt by _af3_rebuild_custom_rows
        self._custom_subunit_rows: list = []
        self._af3_rebuild_custom_rows(2)   # build initial 2-subunit rows

        self._af3_text = QTextEdit()
        self._af3_text.setFont(QFont('Courier New', 9))
        self._af3_text.setReadOnly(True)
        bl.addWidget(self._af3_text)
        splitter.addWidget(bottom_widget)

        # Set initial proportions: 25% selection, 50% jobs, 25% custom/output
        splitter.setSizes([150, 300, 150])
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
        if not self.dna_sequence:
            QMessageBox.warning(self, "Warning", "Load a FASTA file first!")
            return
        if not PYRODIGAL_AVAILABLE:
            QMessageBox.critical(self, "Pyrodigal", t('pyrodigal_not_avail'))
            return
        self._status.showMessage("⏳ Running Pyrodigal gene prediction...")
        min_aa = self._min_length_spin.value()

        def work():
            return self.analyzer.find_orfs_pyrodigal(self.dna_sequence, meta=True, min_aa=min_aa)

        def done(orfs):
            self.orfs = orfs
            self.filtered_orfs = orfs.copy()
            self._update_orfs_list()
            self._update_info()
            self._update_map()
            self._status.showMessage(f"✓ Pyrodigal: {len(orfs)} genes predicted")

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
            QMessageBox.warning(self, "BLAST", "Paste a protein sequence first!"); return
        qp = self._parse_fasta_query(raw)
        if not qp or len(qp) < 5:
            QMessageBox.warning(self, "BLAST", "Invalid or too short sequence!"); return
        if not self.orfs:
            QMessageBox.warning(self, "BLAST", "Run ORF analysis first!"); return

        algo = self._algo_combo.currentText()
        thresh = self._identity_spin.value()
        evalue = float(self._evalue_edit.text() or '0.05')
        self._status.showMessage(f"⏳ BLAST: {len(qp)} aa vs {len(self.orfs)} ORFs...")

        params = {'threshold': thresh, 'gap_open': -self.blast_gap_open,
                  'gap_extend': -self.blast_gap_ext, 'evalue': evalue,
                  'word_size': self.blast_word_size, 'max_targets': self.blast_max_targets,
                  'matrix': self.blast_matrix, 'low_complexity': self.blast_low_complexity}

        def work():
            hits = None; algo_used = algo
            if algo.startswith("Auto"):
                if BACKENDS.get('blast+',{}).get('available'):
                    hits = self.analyzer.run_ncbi_blast(qp, self.orfs, params)
                    algo_used = "NCBI BLAST+"
                if hits is None:
                    hits = self.analyzer.kmer_blast(qp, [o['protein'] for o in self.orfs], params)
                    algo_used = "K-mer Filter"
            elif algo.startswith("NCBI"):
                hits = self.analyzer.run_ncbi_blast(qp, self.orfs, params)
                algo_used = "NCBI BLAST+"
                if hits is None:
                    hits = self.analyzer.kmer_blast(qp, [o['protein'] for o in self.orfs], params)
                    algo_used = "K-mer (fallback)"
            elif algo.startswith("K-mer"):
                hits = self.analyzer.kmer_blast(qp, [o['protein'] for o in self.orfs], params)
                algo_used = "K-mer Filter"
            else:
                hits = self.analyzer.sw_blast(qp, [o['protein'] for o in self.orfs], params)
                algo_used = "Smith-Waterman"
            return (hits or [], algo_used, qp)

        def done(result):
            hits, algo_used, query = result
            self._show_blast_results(hits, query, algo_used)
            self._status.showMessage(f"✓ {algo_used}: {len(hits)} hits")

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
        f, _ = QFileDialog.getOpenFileName(self, "Open FASTA",
            "", "FASTA (*.fasta *.fa *.fna *.faa);;All (*)")
        if not f: return
        with open(f, 'r') as fh:
            content = fh.read()
        # Parse FASTA — get first/longest sequence
        seqs = {}; current = None
        for line in content.split('\n'):
            line = line.strip()
            if line.startswith('>'):
                current = line[1:].split()[0]
                seqs[current] = []
            elif current:
                seqs[current].append(re.sub(r'[^A-Za-z]', '', line))
        if not seqs:
            # Plain sequence
            self.dna_sequence = re.sub(r'[^A-Za-z]', '', content).upper()
            self.genome_name = Path(f).stem
        else:
            longest_name = max(seqs, key=lambda k: len(''.join(seqs[k])))
            self.dna_sequence = ''.join(seqs[longest_name]).upper()
            self.genome_name = longest_name
        self._on_sequence_loaded(f)

    def load_multi_fasta(self):
        self.load_fasta()  # Same logic — picks longest

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
        if not self.orfs: return
        f, _ = QFileDialog.getSaveFileName(self, "Save FASTA", "", "FASTA (*.fasta)")
        if not f: return
        with open(f, 'w') as fh:
            for i, orf in enumerate(self.orfs):
                ds = '|'.join(d['domain'] for d in orf.get('domains', []))
                h = f">ORF{i+1}|F{orf['frame']}{orf['strand']}|{orf['start']}-{orf['end']}"
                if ds: h += f"|{ds}"
                fh.write(h + "\n")
                p = orf['protein'].rstrip('*')
                for j in range(0, len(p), 80): fh.write(p[j:j+80] + "\n")
        self._status.showMessage(f"✓ Saved {len(self.orfs)} ORFs")

    def save_report_tsv(self):
        if not self.orfs: return
        f, _ = QFileDialog.getSaveFileName(self, "Save Report", "", "TSV (*.tsv)")
        if not f: return
        with open(f, 'w', newline='') as fh:
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
        self._status.showMessage(f"✓ Report saved: {Path(f).name}")

    # ───────────────────────────────────────────────────────────
    # PROJECT SAVE / LOAD  (directory-based, v34)
    # ───────────────────────────────────────────────────────────
    PROJECT_MANIFEST = "project.json"
    PROJECT_VERSION  = "v1.01"

    # ─────────────────────────────────────────────────────────
    # PROJECT SAVE / LOAD
    # ─────────────────────────────────────────────────────────

    def _build_manifest(self, proj_dir=None, full_af3_analysis=True):
        """Build project manifest dict.
        proj_dir: Path — if given, HMM/genome file paths are relative to it.
        full_af3_analysis: when False, AF3 analysis entries are saved in a
        lightweight form without embedded PAE/pLDDT arrays.
        """
        import base64

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
        for res in getattr(self, '_af3_analysis_results', []):
            entry = dict(res)
            if not full_af3_analysis:
                # Keep only metadata needed to repopulate the table quickly.
                # Heavy arrays are reloaded on demand by rescanning af3_predictions/.
                for k in ('pae_matrix', 'plddt_arr'):
                    entry.pop(k, None)
                entry['_lightweight'] = True
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
        }

        genome_rel = ''
        if proj_dir and self.dna_sequence:
            safe = re.sub(r'[^\w\.\-]', '_', self.genome_name or 'genome')
            genome_rel = f"genome/{safe}.fasta"

        return {
            'version':              self.PROJECT_VERSION,
            'saved_at':             datetime.now().isoformat(timespec='seconds'),
            'genome_name':          self.genome_name,
            'genome_file':          genome_rel,
            'dna_sequence':         self.dna_sequence,
            'orfs':                 self.orfs,
            'hmm_profiles':         hmm_manifest,
            'af3_jobs':             self.af3_jobs,
            'result_files':         results_manifest,
            'blast_query':          blast_query,
            'blast_results_html':   blast_html,
            'snapgene':             {'features': self.snapgene_features,
                                    'primers':  self.snapgene_primers},
            'ui_state':             ui_state,
            'hpc_server':              hpc_state,
            'hpc_jobs':         getattr(self, '_hpc_jobs', []),
            'af3_analysis_results': af3_analysis_ser,
            'af3_analysis_dir':     getattr(self, '_af3_analysis_dir', ''),
        }

    def save_project(self):
        """Save project as a single self-contained JSON file.
        Includes everything: ORFs, HMM hits, AF3 jobs, PAE matrices, pLDDT.
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
        try:
            with open(path, 'w', encoding='utf-8') as fh:
                json.dump(manifest, fh, indent=2, ensure_ascii=False)
        except OSError as e:
            QMessageBox.critical(self, "Save Project", f"Cannot write file:\n{e}")
            return

        n_af3 = len(manifest.get('af3_analysis_results', []))
        kb = Path(path).stat().st_size / 1024
        self._status.showMessage(
            f"✓ Saved: {Path(path).name}  "
            f"({len(self.orfs)} ORFs, {len(self.hmm_profiles)} HMM, "
            f"{n_af3} AF3 analysis, {kb:.0f} KB)")
        QMessageBox.information(
            self, "Project Saved",
            f"✓ Project saved!\n\n"
            f"  📄 {path}\n\n"
            f"  ORFs:         {len(self.orfs)}\n"
            f"  HMM profiles: {len(self.hmm_profiles)}\n"
            f"  AF3 jobs:     {len(self.af3_jobs)}\n"
            f"  AF3 analysis: {n_af3} job(s) with PAE/pLDDT\n"
            f"  File size:    {kb:.0f} KB\n\n"
            f"To reopen: File → Open Project → Select JSON file\n\n"
            f"⚠ Server password stored as base64 (not encrypted).")

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
        import shutil

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
                    with open(proj_dir / "genome" / genome_fname, 'w') as fh:
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
                    with open(proj_dir / "results" / jfn, 'w') as fh:
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
                with open(proj_dir / "blast" / "query.fasta", 'w') as fh:
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

        # 6. Write lightweight project.json.
        # AF3 folders were already copied above, so do not embed huge PAE/pLDDT
        # arrays again in the manifest.
        if not _step(88, "Writing project manifest..."): return
        manifest = self._build_manifest(proj_dir=proj_dir, full_af3_analysis=False)
        manifest['project_copy_mode'] = 'full_copy_light_manifest'
        # Update af3_analysis_dir and job_dir paths to relative
        manifest['af3_analysis_dir'] = 'af3_predictions'
        for entry in manifest.get('af3_analysis_results', []):
            jname = entry.get('job_name', '')
            entry['job_dir'] = f"af3_predictions/{jname}"
        try:
            with open(proj_dir / self.PROJECT_MANIFEST, 'w', encoding='utf-8') as fh:
                json.dump(manifest, fh, indent=2, ensure_ascii=False)
        except OSError as e:
            QMessageBox.critical(self, "Save Project As",
                                 f"Cannot write manifest:\n{e}")
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
        import base64

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

        # ── Core ──────────────────────────────────────────────
        self.genome_name   = data.get('genome_name', '')
        self.dna_sequence  = data.get('dna_sequence', '')
        self.orfs          = data.get('orfs', [])
        self.filtered_orfs = self.orfs.copy()
        self.af3_jobs      = data.get('af3_jobs', [])

        rel_genome = data.get('genome_file', '')
        if rel_genome:
            abs_g = proj_dir / rel_genome
            self.current_fasta_path = str(abs_g) if abs_g.is_file() else ''

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
        self.hmm_hits_all = []
        for p in self.hmm_profiles:
            for h in p.get('hits', []):
                self.hmm_hits_all.append(dict(h, profile_name=p['name']))

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
        except Exception:
            pass

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
        if needs_af3_rescan and getattr(self, '_af3_analysis_dir', ''):
            try:
                self._af3a_scan_folder(self._af3_analysis_dir)
            except Exception:
                try:
                    self._af3a_populate_table()
                except Exception:
                    pass
        else:
            # Repopulate the AF3 Analysis table
            try:
                self._af3a_populate_table()
            except Exception:
                pass

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
                except Exception: pass
            if 'af3_n_neighbors' in ui:
                try: self._af3_nb_spin.setValue(ui['af3_n_neighbors'])
                except Exception: pass
            if 'zoom_level' in ui:
                self.zoom_level = ui['zoom_level']
                self._zoom_label.setText(f"{int(self.zoom_level * 100)}%")
        except Exception:
            pass

        # ── HPC server ───────────────────────────────────────────
        dv = data.get('hpc_server', {})
        try:
            if dv.get('host'):       self._dv_host.setText(dv['host'])
            if dv.get('user'):       self._dv_user.setText(dv['user'])
            if dv.get('port'):       self._dv_port.setValue(int(dv['port']))
            if dv.get('password'):
                self._dv_pwd.setText(
                    base64.b64decode(dv['password'].encode()).decode('utf-8'))
            if dv.get('base_path'): self._dv_base_path.setText(dv['base_path'])
            if dv.get('af3cmd'):    self._dv_af3cmd.setText(dv['af3cmd'])
            if dv.get('module_cmd') is not None:
                self._dv_module_cmd.setText(dv['module_cmd'])
        except Exception:
            pass
        try:
            self._hpc_jobs = data.get('hpc_jobs', [])
            self._dv_refresh_monitor_table()
        except Exception:
            pass

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

        n_annot    = sum(1 for o in self.orfs if o.get('observation')
                         or o.get('putative_function') or o.get('gene_name')
                         or o.get('custom_color'))
        n_hmm_hits = sum(len(p.get('hits', [])) for p in self.hmm_profiles)
        n_af3a     = len(self._af3_analysis_results)

        self._status.showMessage(
            f"✓ Project loaded: {self.genome_name}  "
            f"[{len(self.orfs)} ORFs | "
            f"{len(self.hmm_profiles)} HMM ({n_hmm_hits} hits) | "
            f"{len(self.af3_jobs)} AF3 jobs | "
            f"{n_af3a} AF3 analyses | "
            f"{n_annot} annotated]  ({ver})")

        saved_at = data.get('saved_at', '')
        QMessageBox.information(
            self, "Project Loaded",
            f"✓ Project loaded!\n\n"
            f"  Genome:        {self.genome_name}\n"
            f"  ORFs:          {len(self.orfs)}\n"
            f"  HMM profiles:  {len(self.hmm_profiles)} ({n_hmm_hits} hits)\n"
            f"  AF3 jobs:      {len(self.af3_jobs)}\n"
            f"  AF3 analyses:  {n_af3a} (with PAE/pLDDT)\n"
            f"  Annotated:     {n_annot}\n"
            + (f"  Saved: {saved_at}\n" if saved_at else "")
            + f"  Format: {ver}")

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
                contact_s    = best.get('contact_region', '-')
            else:
                af3_done     = '-'
                partner      = '-'
                iptm_s       = '-'
                pae_inter_s  = '-'
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
                af3_done, partner, iptm_s, pae_inter_s, contact_s, user_note,
            ]
            for col, val in enumerate(items):
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
                # Make User_note column editable
                if col == 16:
                    item.setFlags(item.flags() |
                                  (Qt.ItemFlag.ItemIsEditable if QT_VERSION == 6
                                   else Qt.ItemIsEditable))
                self._orf_table.setItem(row, col, item)

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

    def _set_zoom(self, level):
        target = max(0.5, min(200.0, level))
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
        # Center map pan on this ORF
        if self.dna_sequence:
            orf_center_frac = (orf['start'] + orf['end']) / 2 / len(self.dna_sequence)
            w = self._genome_map.width()
            mg = 40
            bw = w - 2 * mg
            gw = int(bw * self.zoom_level)
            target_x = mg + orf_center_frac * gw
            self._genome_map.pan_offset = max(0, int(target_x - w / 2))
        self._genome_map.update()

        # Select in table
        if orf in self.filtered_orfs:
            row = self.filtered_orfs.index(orf)
            self._orf_table.selectRow(row)
            self._orf_table.scrollTo(self._orf_table.model().index(row, 0))

        self._status.showMessage(f"✓ ORF{idx+1} selected ({orf['start']:,}-{orf['end']:,})")

    # ═══════════════════════════════════════════════════════════
    # RIGHT-CLICK CONTEXT MENU (ORF Table)
    # ═══════════════════════════════════════════════════════════

    def _on_orf_right_click(self, pos):
        """Right-click on ORF → annotation/color/copy menu."""
        row = self._orf_table.rowAt(pos.y())
        if row < 0 or row >= len(self.filtered_orfs):
            return
        self._orf_table.selectRow(row)
        orf = self.filtered_orfs[row]
        orf_idx = self.orfs.index(orf) if orf in self.orfs else -1

        menu = QMenu(self)
        menu.addAction(f"📝 Annotate ORF{orf_idx+1}...",
                       lambda: self._annotate_orf(orf, orf_idx))
        menu.addAction(f"🎨 Color ORF{orf_idx+1}...",
                       lambda: self._color_orf(orf, orf_idx))
        menu.addSeparator()
        menu.addAction("📋 Copy Protein",
                       lambda: self._copy_to_clipboard(orf['protein'].rstrip('*')))
        menu.addAction("📋 Copy DNA",
                       lambda: self._copy_to_clipboard(orf['dna']))
        menu.addAction("📋 Copy FASTA",
                       lambda: self._copy_to_clipboard(
                           f">ORF{orf_idx+1}|{orf['start']}-{orf['end']}\n{orf['protein'].rstrip('*')}\n"))
        menu.addSeparator()
        menu.addAction("➕ Add to AlphaFold3",
                       lambda: self._af3_add_orf_by_index(orf_idx))
        menu.exec(self._orf_table.viewport().mapToGlobal(pos)
                   if QT_VERSION == 6 else
                   self._orf_table.viewport().mapToGlobal(pos))

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
            with open(f) as fh: self._blast_query_text.setPlainText(fh.read())

    def _copy_blast_hit(self):
        if self.selected_orf:
            self._copy_to_clipboard(self.selected_orf['protein'])

    def _copy_blast_all(self):
        self._copy_to_clipboard(self._blast_results_text.toPlainText())

    def _save_blast_results(self):
        f, _ = QFileDialog.getSaveFileName(self, "Save BLAST Results", "", "Text (*.txt)")
        if f:
            with open(f, 'w') as fh: fh.write(self._blast_results_text.toPlainText())

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

    def _af3_remove_orf(self):
        rows = sorted(set(idx.row() for idx in self._af3_sel_table.selectedIndexes()), reverse=True)
        for r in rows: self._af3_sel_table.removeRow(r)
        self._af3_sel_count.setText(f"{self._af3_sel_table.rowCount()} ORFs selected")

    def _af3_clear_all(self):
        self._af3_sel_table.setRowCount(0); self._af3_jobs_table.setRowCount(0)
        self.af3_jobs = []; self._af3_text.clear(); self._af3_sel_count.setText("0 ORFs selected")

    def _af3_generate_jobs(self):
        n_sel = self._af3_sel_table.rowCount()
        if n_sel == 0: QMessageBox.warning(self, "AF3", "Select ORFs first!"); return
        if not self.orfs: QMessageBox.warning(self, "AF3", "Run ORF analysis first!"); return
        n_nb = self._af3_nb_spin.value(); mode = self._af3_mode_combo.currentText(); self.af3_jobs = []
        sel_indices = []
        for r in range(n_sel):
            try:
                idx = int(self._af3_sel_table.item(r,0).text().replace('ORF',''))-1
                if 0 <= idx < len(self.orfs): sel_indices.append(idx)
            except: continue
        orfs_by_pos = sorted(enumerate(self.orfs), key=lambda x: x[1]['start'])
        pos_to_rank = {idx: rank for rank, (idx, _) in enumerate(orfs_by_pos)}
        for hi in sel_indices:
            ho = self.orfs[hi]; hr = pos_to_rank.get(hi,0); hp = ho['protein'].rstrip('*'); hn = f"ORF{hi+1}"
            nbs = []
            for d in range(-n_nb, n_nb+1):
                if d == 0: continue
                nr = hr + d
                if 0 <= nr < len(orfs_by_pos):
                    ni, no = orfs_by_pos[nr]; nbs.append((ni, no, d))
            if mode.startswith("Pares") and "Homodímero" not in mode:
                for ni, no, d in nbs:
                    np_s = no['protein'].rstrip('*'); tr = len(hp)+len(np_s)
                    self.af3_jobs.append({'name': f"{hn}_vs_ORF{ni+1}_{'up' if d<0 else 'down'}{abs(d)}",
                        'hit_orf_idx': hi, 'partner_orf_idx': ni, 'hit_name': hn, 'partner_name': f"ORF{ni+1}",
                        'total_residues': tr, 'status': 'pending' if tr<=self.af3_max_residues else f'>{self.af3_max_residues}!',
                        'iptm': None, 'plddt': None,
                        'sequences': [{'proteinChain':{'sequence':hp,'count':1}},{'proteinChain':{'sequence':np_s,'count':1}}]})
            elif mode.startswith("Homodímero"):
                tr = len(hp)*2
                self.af3_jobs.append({'name': f"{hn}_homodimer", 'hit_orf_idx': hi, 'partner_orf_idx': hi,
                    'hit_name': hn, 'partner_name': hn, 'total_residues': tr,
                    'status': 'pending' if tr<=self.af3_max_residues else f'>{self.af3_max_residues}!',
                    'iptm': None, 'plddt': None, 'sequences': [{'proteinChain':{'sequence':hp,'count':2}}]})
            elif mode.startswith("Pares + Homodímero"):
                for ni, no, d in nbs:
                    np_s = no['protein'].rstrip('*'); tr = len(hp)+len(np_s)
                    self.af3_jobs.append({'name': f"{hn}_vs_ORF{ni+1}_{'up' if d<0 else 'down'}{abs(d)}",
                        'hit_orf_idx': hi, 'partner_orf_idx': ni, 'hit_name': hn, 'partner_name': f"ORF{ni+1}",
                        'total_residues': tr, 'status': 'pending' if tr<=self.af3_max_residues else f'>{self.af3_max_residues}!',
                        'iptm': None, 'plddt': None,
                        'sequences': [{'proteinChain':{'sequence':hp,'count':1}},{'proteinChain':{'sequence':np_s,'count':1}}]})
                tr2 = len(hp)*2
                self.af3_jobs.append({'name': f"{hn}_homodimer", 'hit_orf_idx': hi, 'partner_orf_idx': hi,
                    'hit_name': hn, 'partner_name': hn, 'total_residues': tr2,
                    'status': 'pending' if tr2<=self.af3_max_residues else f'>{self.af3_max_residues}!',
                    'iptm': None, 'plddt': None, 'sequences': [{'proteinChain':{'sequence':hp,'count':2}}]})
        self._af3_update_jobs_table()
        self._status.showMessage(f"✓ {len(self.af3_jobs)} AF3 jobs generated")

    def _af3_update_jobs_table(self):
        self._af3_jobs_table.setRowCount(0)
        for j in self.af3_jobs:
            row = self._af3_jobs_table.rowCount(); self._af3_jobs_table.insertRow(row)
            iptm_s = f"{j['iptm']:.3f}" if j.get('iptm') is not None else '-'
            plddt_s = f"{j.get('plddt',0):.1f}" if j.get('plddt') else '-'
            for col, val in enumerate([j['name'],j['hit_name'],j['partner_name'],str(j['total_residues']),j['status'],iptm_s,plddt_s]):
                self._af3_jobs_table.setItem(row, col, QTableWidgetItem(val))

    def _af3_export_json(self):
        if not self.af3_jobs: QMessageBox.warning(self,"AF3","Generate jobs first!"); return
        folder = QFileDialog.getExistingDirectory(self, "Select output folder")
        if not folder: return
        for j in self.af3_jobs:
            data = {"name":j['name'],"modelSeeds":[],"sequences":j['sequences'],"dialect":"alphafoldserver","version":2}
            with open(os.path.join(folder, f"{j['name']}.json"), 'w') as f: json.dump(data, f, indent=2)
        self._status.showMessage(f"✓ {len(self.af3_jobs)} AF3 JSONs exported")

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
        for j in self.af3_jobs:
            chains = [s['proteinChain']['sequence'] for s in j['sequences']]
            joined = ':'.join(chains); csv_lines.append(f"{j['name']},{joined}")
            with open(os.path.join(folder, f"{j['name']}.fasta"), 'w') as f: f.write(f">{j['name']}\n{joined}\n")
        with open(os.path.join(folder, f"{self.genome_name}_batch.csv"), 'w') as f: f.write('\n'.join(csv_lines))
        self._status.showMessage(f"✓ ColabFold exported: {len(self.af3_jobs)} jobs")

    def _af3_show_ranking(self):
        done = [j for j in self.af3_jobs if j.get('iptm') is not None]
        if not done: QMessageBox.information(self,"AF3","No results imported yet."); return
        ranked = sorted(done, key=lambda x: x['iptm'], reverse=True)
        txt = f"RANKING ipTM\n{'='*60}\n>0.8=high | 0.6-0.8=likely | <0.4=unlikely\n\n"
        for i, j in enumerate(ranked):
            p = j.get('plddt',0) or 0; s = "★★★" if j['iptm']>=0.8 else "★★" if j['iptm']>=0.6 else "★" if j['iptm']>=0.4 else ""
            txt += f"{i+1:<4} {j['iptm']:.3f}  {p:>5.1f}  {j['hit_name']:<10} {j['partner_name']:<10} {j['name']} {s}\n"
        self._af3_text.setPlainText(txt)

    def _af3_clear_jobs(self):
        self._af3_jobs_table.setRowCount(0); self.af3_jobs = []; self._af3_text.clear()

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
            n_spin.setToolTip(f"Número de cópias da cadeia {letter} no complexo")
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
            QMessageBox.warning(self, "Custom Job", "Nenhuma subunidade definida.")
            return

        # ── Validate & collect all subunits ──
        subunits = []    # list of (orf_idx, n_copies, chain_letter)
        errors = []
        for i, (orf_edit, n_spin) in enumerate(self._custom_subunit_rows):
            letter = self._CHAIN_LETTERS[i] if i < len(self._CHAIN_LETTERS) else str(i + 1)
            raw = orf_edit.text().strip()
            if not raw:
                errors.append(f"Cadeia {letter}: campo ORF vazio.")
                continue
            idx = _parse_orf_idx(raw)
            if idx < 0 or idx >= len(self.orfs):
                errors.append(f"Cadeia {letter}: '{raw}' não encontrada.")
                continue
            subunits.append((idx, n_spin.value(), letter))

        if errors:
            QMessageBox.warning(self, "Custom Job", "\n".join(errors))
            return
        if not subunits:
            QMessageBox.warning(self, "Custom Job", "Nenhuma subunidade válida.")
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
            'sequences':       sequences,
        })
        self._af3_update_jobs_table()
        n_chains = len(subunits)
        self._status.showMessage(
            f"✓ Custom job: {n_chains} subunidade(s) [{chain_summary}]  "
            f"({total_residues} resíduos)")

    def _af3_jobs_right_click(self, pos):
        """Right-click on AF3 jobs table → delete selected jobs."""
        rows = sorted(set(idx.row() for idx in self._af3_jobs_table.selectedIndexes()))
        if not rows: return
        menu = QMenu(self)
        n = len(rows)
        menu.addAction(f"🗑️ Delete {n} selected job{'s' if n > 1 else ''}",
                       self._af3_delete_selected_jobs)
        menu.addSeparator()
        menu.addAction("📋 Copy job names",
                       lambda: self._copy_to_clipboard(
                           '\n'.join(self.af3_jobs[r]['name'] for r in rows if r < len(self.af3_jobs))))
        menu.exec(self._af3_jobs_table.viewport().mapToGlobal(pos)
                   if QT_VERSION == 6 else
                   self._af3_jobs_table.viewport().mapToGlobal(pos))

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
        """Dialog for ORF analysis parameters (min size, start codons)."""
        dlg = QDialog(self)
        dlg.setWindowTitle("🧬 ORF Analysis Parameters")
        dlg.setFixedSize(380, 300)
        layout = QVBoxLayout(dlg)

        # Min length
        gf = QGroupBox("ORF Detection")
        gf_l = QGridLayout(gf)
        gf_l.addWidget(QLabel("Min ORF size (aa):"), 0, 0)
        min_spin = QSpinBox(); min_spin.setRange(10, 500)
        min_spin.setValue(self._min_length_spin.value())
        gf_l.addWidget(min_spin, 0, 1)

        gf_l.addWidget(QLabel("Start codons:"), 1, 0)
        cb_atg = QCheckBox("ATG"); cb_atg.setChecked(self._cb_atg.isChecked())
        cb_gtg = QCheckBox("GTG"); cb_gtg.setChecked(self._cb_gtg.isChecked())
        cb_ttg = QCheckBox("TTG"); cb_ttg.setChecked(self._cb_ttg.isChecked())
        codon_lay = QHBoxLayout()
        codon_lay.addWidget(cb_atg); codon_lay.addWidget(cb_gtg); codon_lay.addWidget(cb_ttg)
        gf_l.addLayout(codon_lay, 1, 1)
        layout.addWidget(gf)

        # Pyrodigal options
        pf = QGroupBox("Pyrodigal")
        pf_l = QGridLayout(pf)
        pf_l.addWidget(QLabel("Mode:"), 0, 0)
        pyro_mode = QComboBox()
        pyro_mode.addItems(["Metagenomic (default)", "Single genome (train)"])
        pf_l.addWidget(pyro_mode, 0, 1)
        pyro_ok = "✅ Installed" if PYRODIGAL_AVAILABLE else "❌ Not installed"
        pf_l.addWidget(QLabel(f"Status: {pyro_ok}"), 1, 0, 1, 2)
        layout.addWidget(pf)

        btns = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
                                 if QT_VERSION == 6 else QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(dlg.accept); btns.rejected.connect(dlg.reject)
        layout.addWidget(btns)

        if dlg.exec() if QT_VERSION == 6 else dlg.exec_():
            self._min_length_spin.setValue(min_spin.value())
            self._cb_atg.setChecked(cb_atg.isChecked())
            self._cb_gtg.setChecked(cb_gtg.isChecked())
            self._cb_ttg.setChecked(cb_ttg.isChecked())

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

    def _show_help_dlg(self, title_key, content_key):
        lang = _CURRENT_LANG[0]
        content = HELP_CONTENT.get(content_key, {}).get(lang) or \
                  HELP_CONTENT.get(content_key, {}).get('en', '(no content)')
        dlg = QDialog(self)
        dlg.setWindowTitle(t(title_key))
        dlg.resize(720, 560)
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

    def _show_about(self):
        blast_ok = "✅" if BACKENDS.get('blast+',{}).get('available') else "❌"
        hmmer_ok = "✅" if BACKENDS.get('hmmer3',{}).get('available') else "❌"
        pyrod_ok = "✅" if BACKENDS.get('pyrodigal',{}).get('available') else "❌"
        mat_ok   = "✅" if MATPLOTLIB_AVAILABLE else "❌"
        QMessageBox.about(self, t('about'),
            f"🧬 ppigFinder — Protein-Protein Interaction Genomic Finder\n"
            f"Version 1.01  |  MIT License\n\n"
            f"Discovery of novel bacterial PPIs via ORF prediction,\n"
            f"HMM/BLAST annotation, genomic neighbourhood analysis\n"
            f"and AlphaFold 3 structural interaction prediction.\n\n"
            f"https://github.com/<your-org>/ppigfinder\n\n"
            f"Backends:\n"
            f"  BLAST+     {blast_ok}\n"
            f"  HMMER3     {hmmer_ok}\n"
            f"  Pyrodigal  {pyrod_ok}\n"
            f"  Matplotlib {mat_ok}\n\n"
            f"PyQt{QT_VERSION} | Python {sys.version.split()[0]}")


    # ═══════════════════════════════════════════════════════════
    # AF3 ANALYSIS TAB  (v2)
    # ═══════════════════════════════════════════════════════════

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
        self._af3a_thresh_spin.setToolTip(
            "Residues with inter-chain PAE < this value are counted as contacts.")
        self._af3a_thresh_spin.valueChanged.connect(
            lambda _: self._af3a_replot_selected())
        top.addWidget(self._af3a_thresh_spin)

        btn_pdf = QPushButton("📄 Export plots PDF")
        btn_pdf.setToolTip("Export all visible PAE/pLDDT plots to a PDF file")
        btn_pdf.clicked.connect(self._af3a_export_pdf)
        top.addWidget(btn_pdf)

        top.addStretch()
        self._af3a_status = QLabel("No results loaded")
        self._af3a_status.setStyleSheet("font-size:11px;color:#666;")
        top.addWidget(self._af3a_status)
        lay.addLayout(top)

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
        self._af3a_table.setColumnCount(9)
        self._af3a_table.setHorizontalHeaderLabels([
            'Job name', 'Chains (n)', 'Chain IDs',
            'ipTM', 'ptm', 'mean_pLDDT', 'ranking_score',
            'PAE_inter (Å)', 'Best contact pair'])
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

        self._tabs.addTab(w, "📊 AlphaFold Analysis")

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

    def _af3a_scan_folder(self, folder: str):
        """Scan folder for AF3 job sub-dirs."""
        # import re as _re — removed (unused import)
        jobs = []
        root = Path(folder)
        candidates = sorted([p for p in root.iterdir() if p.is_dir()])
        if not candidates:
            candidates = [root]

        for job_dir in candidates:
            # Root-level files: *_summary_confidences.json + *_confidences.json
            sum_files  = [f for f in job_dir.glob('*_summary_confidences.json')
                          if f.parent == job_dir]
            conf_files = [f for f in job_dir.glob('*_confidences.json')
                          if f.parent == job_dir and 'summary' not in f.name
                          and 'data' not in f.name]
            if not sum_files:
                continue
            try:
                parsed = self._af3a_parse_job(
                    job_dir, sum_files[0],
                    conf_files[0] if conf_files else None)
                if parsed:
                    jobs.append(parsed)
            except Exception as e:
                print(f"[AF3 Analysis] {job_dir.name}: {e}")

        self._af3_analysis_results = jobs
        self._af3a_populate_table()
        self._af3a_status.setText(
            f"{len(jobs)} job(s) loaded from {Path(folder).name}/")
        self._update_orfs_list()

    def _af3a_parse_job(self, job_dir: Path, sum_path: Path,
                         conf_path) -> dict:
        """Parse one AF3 job. Handles any number of chains."""
        import re as _re

        # ── Summary scores ─────────────────────────────────────
        with open(sum_path, encoding='utf-8') as f:
            summary = json.load(f)

        iptm          = summary.get('iptm')
        ptm           = summary.get('ptm')
        mean_plddt    = (summary.get('mean_plddt')
                         or summary.get('mean_pLDDT')
                         or summary.get('ptm_plddt'))  # AF3 uses lowercase
        ranking_score = summary.get('ranking_score')

        # ── PAE + pLDDT + chain info ────────────────────────────
        pae_matrix  = None
        plddt_arr   = None
        chain_ids   = None   # per-residue chain labels

        if conf_path and Path(conf_path).is_file():
            with open(conf_path, encoding='utf-8') as f:
                conf = json.load(f)
            pae_raw = conf.get('pae')
            if pae_raw:
                pae_matrix = pae_raw
            # pLDDT — try several field names used by different AF3 versions
            for key in ('plddt', 'atom_plddts', 'predicted_lddt',
                        'residue_plddt'):
                v = conf.get(key)
                if v:
                    plddt_arr = v
                    break
            # Chain IDs per residue
            for key in ('token_chain_ids', 'chain_ids', 'asym_id'):
                v = conf.get(key)
                if v:
                    chain_ids = v
                    break

        # If mean_plddt missing from summary, compute from plddt array
        if mean_plddt is None and plddt_arr:
            try:
                import numpy as np
                mean_plddt = float(np.mean(plddt_arr))
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
                     for m in _re.finditer(r'orf(\d+)',
                                           job_dir.name, _re.IGNORECASE)]

        # Map chain letters → ORF names (A→orf_names[0], B→orf_names[1] …)
        chain_to_orf = {}
        for ci, cid in enumerate(chain_order):
            chain_to_orf[cid] = (orf_names[ci]
                                 if ci < len(orf_names) else cid)

        # ── Inter-chain metrics ─────────────────────────────────
        thresh = self._af3a_thresh_spin.value()
        pair_metrics = {}   # (cid_A, cid_B) → dict
        best_pae_inter = None
        best_pair      = ('?', '?')

        if pae_matrix and n_chains >= 2:
            try:
                import numpy as np
                pae_np = np.array(pae_matrix, dtype=float)
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
                        pi = float((sub_AB.mean() + sub_BA.mean()) / 2)
                        mA = sub_AB.mean(axis=1)
                        mB = sub_BA.mean(axis=1)
                        contacts = int((mA < thresh).sum()
                                       + (mB < thresh).sum())
                        cr = self._af3a_contact_str(
                            mA, mB,
                            chain_to_orf.get(ca, ca),
                            chain_to_orf.get(cb, cb),
                            thresh)
                        pair_metrics[(ca, cb)] = {
                            'pae_inter': pi, 'n_contacts': contacts,
                            'contact_region': cr}
                        if best_pae_inter is None or pi < best_pae_inter:
                            best_pae_inter = pi
                            best_pair      = (ca, cb)
            except Exception as e:
                print(f"[AF3 Analysis] inter-metric: {e}")

        best_cr = (pair_metrics[best_pair]['contact_region']
                   if best_pair in pair_metrics else '-')

        return {
            'job_name':      job_dir.name,
            'job_dir':       str(job_dir),
            'orf_names':     orf_names,
            'chain_order':   chain_order,
            'chain_lens':    chain_lens,
            'chain_to_orf':  chain_to_orf,
            'n_chains':      n_chains,
            'iptm':          iptm,
            'ptm':           ptm,
            'mean_plddt':    mean_plddt,
            'ranking_score': ranking_score,
            'pae_matrix':    pae_matrix,
            'plddt_arr':     plddt_arr,
            'pair_metrics':  pair_metrics,
            'pae_inter':     best_pae_inter,
            'best_pair':     best_pair,
            'contact_region': best_cr,
            'partner_name':  (orf_names[1] if len(orf_names) > 1
                              else chain_order[1] if len(chain_order) > 1
                              else '-'),
        }

    @staticmethod
    def _af3a_contact_str(mean_A, mean_B, name_A, name_B, thresh):
        """Build a human-readable contact region string."""
        import numpy as np
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

        c5_A  = np.where(mean_A < 5.0)[0]
        c10_A = np.where(mean_A < 10.0)[0]
        c5_B  = np.where(mean_B < 5.0)[0]
        c10_B = np.where(mean_B < 10.0)[0]

        if not len(c5_A) and not len(c10_A):
            return f"not relevant interaction found with {name_B}"

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
        self._af3a_table.setRowCount(0)
        def _c(v, thr_good, thr_ok, inv=False):
            if v is None: return None
            ok = (v <= thr_ok) if inv else (v >= thr_good)
            mid = (v <= thr_good) if inv else (v >= thr_ok)
            if ok:   return QColor('#C8E6C9')
            if mid:  return QColor('#FFF9C4')
            return QColor('#FFCDD2')

        for res in self._af3_analysis_results:
            row = self._af3a_table.rowCount()
            self._af3a_table.insertRow(row)
            chains_s = ', '.join(
                f"{res['chain_to_orf'].get(c, c)}({res['chain_lens'][c]}aa)"
                for c in res['chain_order'])
            iptm_s   = f"{res['iptm']:.3f}"  if res.get('iptm')   is not None else '-'
            ptm_s    = f"{res['ptm']:.3f}"   if res.get('ptm')    is not None else '-'
            plddt_s  = f"{res['mean_plddt']:.1f}" if res.get('mean_plddt') is not None else '-'
            rank_s   = f"{res['ranking_score']:.4f}" if res.get('ranking_score') is not None else '-'
            pi_s     = f"{res['pae_inter']:.1f}" if res.get('pae_inter') is not None else '-'

            vals = [res['job_name'], str(res['n_chains']), chains_s,
                    iptm_s, ptm_s, plddt_s, rank_s, pi_s,
                    res.get('contact_region', '-')]
            for col, val in enumerate(vals):
                item = QTableWidgetItem(str(val))
                bg = None
                if col == 3: bg = _c(res.get('iptm'), 0.75, 0.50)
                if col == 7: bg = _c(res.get('pae_inter'), 8.0, 15.0, inv=True)
                if bg: item.setBackground(bg)
                self._af3a_table.setItem(row, col, item)

    def _af3a_on_select(self):
        rows = set(idx.row() for idx in self._af3a_table.selectedIndexes())
        if not rows:
            return
        row = min(rows)
        if row >= len(self._af3_analysis_results):
            return

        res = self._af3_analysis_results[row]
        self._af3a_plot_job(res)

        for orf_name in res.get('orf_names', []):
            m = re.match(r'ORF(\d+)', orf_name, re.IGNORECASE)
            if m:
                idx = int(m.group(1)) - 1
                if 0 <= idx < len(self.orfs):
                    self._select_and_center_orf(idx)
                    break

    def _af3a_replot_selected(self):
        rows = set(idx.row() for idx in self._af3a_table.selectedIndexes())
        if rows:
            row = min(rows)
            if row < len(self._af3_analysis_results):
                # Recompute pair metrics with new threshold
                res = self._af3_analysis_results[row]
                try:
                    import numpy as np
                    thresh = self._af3a_thresh_spin.value()
                    if res.get('pae_matrix') and res.get('n_chains', 1) >= 2:
                        pae_np = np.array(res['pae_matrix'], dtype=float)
                        for (ca, cb) in list(res['pair_metrics'].keys()):
                            chain_order = res['chain_order']
                            i_c = chain_order.index(ca)
                            j_c = chain_order.index(cb)
                            r0 = sum(res['chain_lens'][chain_order[k]] for k in range(i_c))
                            r1 = r0 + res['chain_lens'][ca]
                            c0 = sum(res['chain_lens'][chain_order[k]] for k in range(j_c))
                            c1 = c0 + res['chain_lens'][cb]
                            mA = pae_np[r0:r1, c0:c1].mean(axis=1)
                            mB = pae_np[c0:c1, r0:r1].mean(axis=1)
                            cr = self._af3a_contact_str(
                                mA, mB,
                                res['chain_to_orf'].get(ca, ca),
                                res['chain_to_orf'].get(cb, cb), thresh)
                            res['pair_metrics'][(ca, cb)]['contact_region'] = cr
                except Exception:
                    pass
                self._af3a_plot_job(res)

    def _af3a_clear_plots(self):
        self._af3a_active_canvases = []
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
        import numpy as np
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

        # Contact regions for each pair
        if res.get('pair_metrics'):
            for (ca, cb), pm in res['pair_metrics'].items():
                cr_lbl = QLabel(
                    f"  {chain_to_orf.get(ca,ca)} ↔ "
                    f"{chain_to_orf.get(cb,cb)}: "
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
            cum_i = 0
            for i_c, ca in enumerate(chain_order):
                len_a = chain_lens.get(ca, 0)
                cum_j = 0
                for j_c, cb in enumerate(chain_order):
                    len_b = chain_lens.get(cb, 0)
                    if i_c != j_c:
                        sub = pae_np[cum_i:cum_i+len_a,
                                     cum_j:cum_j+len_b]
                        mA = sub.mean(axis=1)
                        ci_contact = np.where(mA < thresh)[0]
                        if len(ci_contact):
                            ax.scatter(
                                [cum_j + len_b / 2] * len(ci_contact),
                                cum_i + ci_contact,
                                s=max(2, 6 - n_chains),
                                color='lime', alpha=0.55,
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

            def _on_hover(event, _pae=pae_np_ref,
                          _co=chain_order_ref, _cl=chain_lens_ref,
                          _c2o=c2o_ref, _ax=ax, _lbl=hover_lbl):
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
                    _lbl.setText(
                        f"  PAE  scored={cx} res{xi+1}  "
                        f"aligned={cy} res{yi+1}  →  "
                        f"{val:.2f} Å")

            canvas.mpl_connect('motion_notify_event', _on_hover)
            self._af3a_active_canvases.append((fig, canvas,
                                                res['job_name'] + '_PAE'))
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
            self._af3a_active_canvases.append((fig2, canvas2,
                                                res['job_name'] + '_pLDDT'))
            plt.close(fig2)
            self._af3a_plot_layout.addWidget(canvas2)

        self._af3a_plot_layout.addStretch()

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


    # ═══════════════════════════════════════════════════════════
    # HPC SERVER TAB (generic)
    # ═══════════════════════════════════════════════════════════

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

        # Sub-tabs
        self._dv_tabs = QTabWidget()
        self._dv_tabs.setTabPosition(QTabWidget.TabPosition.North
                                      if QT_VERSION == 6
                                      else QTabWidget.North)
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
        self._dv_src_batch = QPushButton("📦 Load session as single batch JSON  (recommended for 10+ jobs)")
        self._dv_src_batch.clicked.connect(self._dv_load_from_session_batch)
        self._dv_src_batch.setToolTip(
            "Combines all AF3 jobs from the current session into a single\n"
            "batch JSON file. Uploads one file and runs af3_run once.\n"
            "Equivalent to AlphaFold → Export AF3 JSON → Batch JSON.")
        self._dv_src_file = QPushButton("📂 Load JSON file(s) from disk")
        self._dv_src_file.clicked.connect(self._dv_load_from_files)
        for b in (self._dv_src_session, self._dv_src_batch, self._dv_src_file):
            b.setFixedHeight(28)
            src_l.addWidget(b)
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
        for b in (self._dv_clear_btn, self._dv_upload_btn, self._dv_run_btn):
            b.setEnabled(False)
            act_row.addWidget(b)
        lay.addLayout(act_row)

        # ── AF3 Advanced Options ──────────────────────────────
        af3_g = QGroupBox("⚗ AF3 Advanced Options")
        af3_g.setCheckable(True)
        af3_g.setChecked(False)   # collapsed by default
        af3_g.setToolTip(
            "Expand to control AlphaFold 3 prediction quality and speed.\n"
            "Changes here are injected into the batch JSON and/or the\n"
            "af3_run command before submission.")
        af3_l = QGridLayout(af3_g)
        af3_l.setSpacing(5)

        # Row 0: Preset selector
        af3_l.addWidget(QLabel("Preset:"), 0, 0)
        self._dv_af3_preset = QComboBox()
        self._dv_af3_preset.addItems([
            "Balanced (default)",
            "Fast  — 1 seed, no templates",
            "Accurate  — 5 seeds, templates ON",
            "Custom",
        ])
        self._dv_af3_preset.setToolTip(
            "Fast:     1 seed, templates disabled — good for screening.\n"
            "Balanced: 1 seed, templates auto — default AlphaFold3 run.\n"
            "Accurate: 5 seeds, templates ON, pick best model — slower.\n"
            "Custom:   set each parameter manually below.")
        self._dv_af3_preset.currentIndexChanged.connect(
            self._dv_af3_apply_preset)
        af3_l.addWidget(self._dv_af3_preset, 0, 1, 1, 3)

        # Row 1: Num seeds + Num models
        af3_l.addWidget(QLabel("Model seeds (n):"), 1, 0)
        self._dv_af3_seeds = QSpinBox()
        self._dv_af3_seeds.setRange(1, 10)
        self._dv_af3_seeds.setValue(1)
        self._dv_af3_seeds.setToolTip(
            "Number of random seeds for af3_run.\n"
            "More seeds = ensemble of independent models → pick best ipTM.\n"
            "1 = fastest;  5 = recommended for final confident predictions.")
        self._dv_af3_seeds.valueChanged.connect(self._dv_refresh_cmd_preview)
        af3_l.addWidget(self._dv_af3_seeds, 1, 1)

        af3_l.addWidget(QLabel("Num models:"), 1, 2)
        self._dv_af3_num_models = QSpinBox()
        self._dv_af3_num_models.setRange(1, 5)
        self._dv_af3_num_models.setValue(1)
        self._dv_af3_num_models.setToolTip(
            "Number of AF3 diffusion models to run per job\n"
            "(--num_predictions flag if supported by your af3_run).\n"
            "5 models × 1 seed = 5 structures to rank by confidence.")
        self._dv_af3_num_models.valueChanged.connect(
            self._dv_refresh_cmd_preview)
        af3_l.addWidget(self._dv_af3_num_models, 1, 3)

        # Row 2: Seeds info label (templates are not configurable via af3_run CLI)
        _note = QLabel(
            "ℹ  Seeds are injected into JSON modelSeeds[ ].  "
            "Templates are controlled by the server's AF3 installation.")
        _note.setStyleSheet("font-size:10px; color:#666;")
        _note.setWordWrap(True)
        af3_l.addWidget(_note, 2, 0, 1, 4)

        # Row 3: SLURM partition + extra flags
        af3_l.addWidget(QLabel("SLURM partition:"), 3, 0)
        self._dv_af3_partition = QLineEdit()
        self._dv_af3_partition.setPlaceholderText("e.g. gpu, batch, compute  (server default if blank)")
        self._dv_af3_partition.setToolTip(
            "SLURM partition to submit to (passed as --partition <name>).\n"
            "Leave blank to use the server default partition.")
        self._dv_af3_partition.textChanged.connect(self._dv_refresh_cmd_preview)
        af3_l.addWidget(self._dv_af3_partition, 3, 1, 1, 3)

        af3_l.addWidget(QLabel("Extra af3_run flags:"), 4, 0)
        self._dv_af3_extra_flags = QLineEdit()
        self._dv_af3_extra_flags.setPlaceholderText(
            "e.g. --slurm-time 12:00:00  --slurm-mem 64G  --slurm-gres gpu:a100:1")
        self._dv_af3_extra_flags.setToolTip(
            "Real af3_run SLURM flags (from af3_run --help):\n"
            "  --slurm-partition <p>  queue/partition name\n"
            "  --slurm-time HH:MM:SS  max wall time\n"
            "  --slurm-mem MG         RAM (e.g. 64G)\n"
            "  --slurm-gres gpu:a100:1 GPU type\n"
            "  --slurm-nodes N        number of nodes\n"
            "  --slurm-ntasks N       tasks per node\n"
            "  --workdir /path        override workdir\n"
            "  --force                overwrite existing output\n"
            "  --dry-run              preview without submitting")
        self._dv_af3_extra_flags.textChanged.connect(
            self._dv_refresh_cmd_preview)
        af3_l.addWidget(self._dv_af3_extra_flags, 4, 1, 1, 3)

        lay.addWidget(af3_g)

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
        import re as _re
        m_rc = _re.search(r'__ACT_RC__(\d+)', combined)
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
        for b in (self._dv_upload_btn, self._dv_run_btn, self._dv_clear_btn):
            b.setEnabled(True)
        self._dv_refresh_cmd_preview()

    def _dv_load_from_session_batch(self):
        """Build a single batch JSON from all session AF3 jobs and stage it
        for upload as one file — exactly like AlphaFold → Export Batch JSON
        but wired into the HPC submit workflow."""
        if not self.af3_jobs:
            QMessageBox.information(self, "Server",
                "No AF3 jobs in session.\n"
                "Generate jobs in the AlphaFold tab first.")
            return

        prefix = self._dv_job_prefix.text().strip() or "af3_batch"
        batch_name = f"{prefix}_all_jobs"

        # Build combined list in AF3 batch format
        batch_list = []
        for j in self.af3_jobs:
            batch_list.append({
                "name":       j['name'],
                "modelSeeds": [],
                "sequences":  j.get('sequences', []),
                "dialect":    "alphafoldserver",
                "version":    1,
            })

        # Write to a temp file (will be uploaded as one JSON)
        import tempfile as _tmp
        fd, tmp_path = _tmp.mkstemp(suffix='.json', prefix=batch_name + '_')
        os.close(fd)
        with open(tmp_path, 'w', encoding='utf-8') as f:
            json.dump(batch_list, f, indent=2, ensure_ascii=False)
        file_kb = os.path.getsize(tmp_path) / 1024

        # Stage as a single "virtual" job whose local_path points to the temp file
        self._dv_submit_table.setRowCount(0)
        n_jobs  = len(batch_list)
        res_tot = sum(j.get('total_residues', 0) for j in self.af3_jobs)
        row = self._dv_submit_table.rowCount()
        self._dv_submit_table.insertRow(row)
        for col, val in enumerate([
                batch_name, str(n_jobs), str(res_tot),
                f"{file_kb:.0f}", "Pending (batch)"]):
            self._dv_submit_table.setItem(row, col, QTableWidgetItem(val))

        self._dv_pending_jobs = [{
            'name':          batch_name,
            'local_path':    tmp_path,
            'sequences':     [],          # not needed — file already written
            'total_residues': res_tot,
            'status':        'pending',
            '_is_batch':     True,
            '_n_jobs':       n_jobs,
        }]

        self._dv_submit_summary.setText(
            f"1 batch file  ({n_jobs} jobs combined, {file_kb:.0f} KB)")
        for b in (self._dv_upload_btn, self._dv_run_btn, self._dv_clear_btn):
            b.setEnabled(True)
        self._dv_log(
            f"Batch JSON ready: {batch_name}.json  "
            f"({n_jobs} jobs, {file_kb:.0f} KB)", 'submit')
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
                  self._dv_clear_btn):
            b.setEnabled(n > 0)
        self._dv_refresh_cmd_preview()

    def _dv_clear_submit_list(self):
        self._dv_submit_table.setRowCount(0)
        self._dv_pending_jobs = []
        self._dv_submit_summary.setText("0 jobs loaded")
        for b in (self._dv_upload_btn, self._dv_run_btn):
            b.setEnabled(False)

    # ── AF3 preset helper ──────────────────────────────────────
    def _dv_af3_apply_preset(self, idx: int):
        """Fill AF3 parameter widgets according to the chosen preset.

        idx:  0 = Balanced (default)
              1 = Fast
              2 = Accurate
              3 = Custom  (no change — user edits manually)
        """
        if not hasattr(self, '_dv_af3_seeds'):
            return  # widgets not yet built
        # Extra flags use real af3_run SLURM flag names
        presets = {
            0: dict(seeds=1, models=1, templates=True,  tpl_date='', extra=''),
            1: dict(seeds=1, models=1, templates=False, tpl_date='', extra='--slurm-time 04:00:00'),
            2: dict(seeds=5, models=5, templates=True,  tpl_date='', extra='--slurm-time 24:00:00'),
        }
        if idx not in presets:
            return   # Custom — leave current values
        p = presets[idx]
        self._dv_af3_seeds.setValue(p['seeds'])
        self._dv_af3_num_models.setValue(p['models'])
        self._dv_af3_use_templates.setChecked(p['templates'])
        self._dv_af3_tpl_date.setText(p['tpl_date'])
        self._dv_af3_extra_flags.setText(p['extra'])
        self._dv_refresh_cmd_preview()

    def _dv_refresh_cmd_preview(self):
        """Rebuild the command preview including AF3 advanced flags.

        Correct af3_run usage:
          cd <parent_dir>
          af3_run --json_path <fname>.json --job-name <prefix>
                  [--num_seeds N] [--num_predictions N]
                  [--notemplate] [--max_template_date YYYY-MM-DD]
                  [--partition <p>] [<extra_flags>]
        """
        prefix   = self._dv_job_prefix.text().strip() or "af3_batch"
        base     = self._dv_base_path.text().strip().rstrip('/')
        rdir     = (self._dv_remote_dir.text().strip()
                    or datetime.now().strftime('%Y-%m-%d'))
        cmd      = self._dv_af3cmd.text().strip() or "af3_run"
        mod_name = self._dv_module_cmd.text().strip()
        jobs     = getattr(self, '_dv_pending_jobs', [])
        n_jobs   = len(jobs)
        use_ts   = getattr(self, '_dv_ts_check', None) and \
                   self._dv_ts_check.isChecked()

        ts          = datetime.now().strftime('%Y%m%d_%H%M%S') if use_ts else ''
        job_name    = f"{prefix}_{ts}" if ts else prefix
        # Auto-sanitize: bash special chars in job name break bash -lc '...'
        job_name    = re.sub(r'[()\[\]{}|;&!\s]+', '_', job_name).strip('_')
        parent_dir  = f"{base}/{rdir}"
        json_fname  = f"{prefix}_all_jobs.json"

        # ── Collect AF3 advanced flags ────────────────────────
        af3_flags = []
        # Seeds (n) and num_models affect the JSON (modelSeeds array),
        # not the af3_run command line — shown as comments in preview only.
        if hasattr(self, '_dv_af3_seeds'):
            n_seeds = self._dv_af3_seeds.value()
            if n_seeds > 1:
                af3_flags.append(f"# {n_seeds} seeds injected into JSON modelSeeds")

        if hasattr(self, '_dv_af3_num_models'):
            n_models = self._dv_af3_num_models.value()
            if n_models > 1:
                af3_flags.append(f"# {n_models} models via JSON numDiffusionSamples")

        # Note: --notemplate / --max_template_date / --num_seeds are NOT
        # valid af3_run flags — template/seed control is via JSON modelSeeds.
        if hasattr(self, '_dv_af3_partition'):
            part = self._dv_af3_partition.text().strip()
            if part:
                af3_flags.append(f"--slurm-partition {part}")

        if hasattr(self, '_dv_af3_extra_flags'):
            extra = self._dv_af3_extra_flags.text().strip()
            if extra:
                af3_flags.append(extra)

        # ── Build preview lines ───────────────────────────────
        lines = []
        if mod_name:
            lines.append(f"# module loaded: {mod_name}")
        lines.append(f"cd {parent_dir}")

        base_cmd = f"{cmd} --json_path {json_fname} --job-name {job_name}"
        if af3_flags:
            lines.append(base_cmd + " \\")
            for i, flag in enumerate(af3_flags):
                suffix = " \\" if i < len(af3_flags) - 1 else ""
                lines.append(f"    {flag}{suffix}")
        else:
            lines.append(base_cmd)

        lines.append(f"# Output → {parent_dir}/{job_name}/output/")
        if n_jobs:
            lines.append(f"# {n_jobs} job(s) bundled in batch JSON")
        if af3_flags:
            preset_names = ["Balanced", "Fast", "Accurate", "Custom"]
            pidx = getattr(self, '_dv_af3_preset', None)
            pname = preset_names[pidx.currentIndex()] if pidx else "Custom"
            lines.append(f"# AF3 preset: {pname}")
        self._dv_cmd_preview.setPlainText('\n'.join(lines))

    def _dv_upload_only(self):
        """Upload JSON files to the server without submitting."""
        self._dv_do_upload(submit=False)

    def _dv_upload_and_submit(self):
        """Upload JSON files and run af3_run on the server."""
        self._dv_do_upload(submit=True)

    def _dv_do_upload(self, submit: bool):
        """SFTP upload + optional af3_run.

        Always bundles all pending jobs into ONE batch JSON file and
        runs a SINGLE af3_run call, matching the server's expected usage:

            af3_run --json_path <prefix>_all_jobs.json --job-name <prefix>

        The results directory on the server will be:
            ~/af3_predictions/<rdir>/<prefix>/output/
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

        # Timestamp suffix — generated once at submit time so preview matches
        use_ts   = getattr(self, '_dv_ts_check', None) and \
                   self._dv_ts_check.isChecked()
        ts       = datetime.now().strftime('%Y%m%d_%H%M%S') if use_ts else ''
        job_name = f"{prefix}_{ts}" if ts else prefix
        # Auto-sanitize: replace shell-unsafe chars so bash -lc never breaks
        job_name = re.sub(r'[()\[\]{}|;&!\s]+', '_', job_name).strip('_')

        # Parent dir = base/rdir  (we cd here before running af3_run)
        # JSON is uploaded into parent_dir/ directly — filename only for --json_path
        parent_dir      = f"{base}/{rdir}"
        batch_json_name = f"{prefix}_all_jobs.json"   # filename only
        # remote_job_dir removed — unused variable (parent_dir used directly)
        # af3_run creates job_name/output/ inside parent_dir
        output_dir      = f"{parent_dir}/{job_name}"

        ssh = self._ssh_client

        def _do_upload_submit():
            # ── 1. Collect AF3 advanced settings ──────────────
            # Seeds: generate reproducible integer seeds based on count
            n_seeds = 1
            use_templates = True
            if hasattr(self, '_dv_af3_seeds'):
                n_seeds = self._dv_af3_seeds.value()
            if hasattr(self, '_dv_af3_use_templates'):
                use_templates = self._dv_af3_use_templates.isChecked()

            import random as _random
            _random.seed(42)
            model_seeds = [_random.randint(1, 2**31 - 1)
                           for _ in range(n_seeds)] if n_seeds > 1 else []

            # ── 2. Build batch JSON in memory ─────────────────
            batch_list = []
            for job in jobs:
                if job.get('local_path'):
                    try:
                        with open(job['local_path'], 'r', encoding='utf-8') as f:
                            raw = json.load(f)
                        if isinstance(raw, list):
                            for entry in raw:
                                if isinstance(entry, dict):
                                    entry['modelSeeds'] = model_seeds
                                    if not use_templates:
                                        # Mark sequences to skip templates
                                        for seq in entry.get('sequences', []):
                                            for chain in seq.get(
                                                    'proteinChain', []) or []:
                                                chain['templates'] = []
                            batch_list.extend(raw)
                        else:
                            raw['modelSeeds'] = model_seeds
                            batch_list.append(raw)
                    except Exception:
                        batch_list.append({
                            "name":        job['name'],
                            "modelSeeds":  model_seeds,
                            "sequences":   job.get('sequences', []),
                            "dialect":     "alphafoldserver",
                            "version":     1,
                        })
                else:
                    batch_list.append({
                        "name":       job['name'],
                        "modelSeeds": model_seeds,
                        "sequences":  job.get('sequences', []),
                        "dialect":    "alphafoldserver",
                        "version":    1,
                    })

            if not batch_list:
                raise RuntimeError("No valid jobs to submit.")

            # ── 2. Write batch JSON to temp file ───────────────
            import tempfile as _tmp
            fd, tmp_path = _tmp.mkstemp(
                suffix='.json',
                prefix=re.sub(r'[^\w\-]', '_', prefix) + '_')
            os.close(fd)
            try:
                with open(tmp_path, 'w', encoding='utf-8') as f:
                    json.dump(batch_list, f, indent=2, ensure_ascii=False)
                file_kb = os.path.getsize(tmp_path) / 1024

                # ── 3. Resolve ~ for SFTP ──────────────────────
                # SFTP doesn't expand tilde — resolve to absolute path
                resolved_parent = parent_dir
                if parent_dir.startswith('~'):
                    home_out, _, _ = self._dv_ssh_exec("echo $HOME", timeout=5)
                    home = home_out.strip()
                    if home:
                        resolved_parent = home + parent_dir[1:]

                # JSON uploaded directly into parent_dir (not a sub-dir)
                # so --json_path can be filename-only after cd
                remote_json_sftp  = f"{resolved_parent}/{batch_json_name}"
                remote_json_shell = batch_json_name   # filename only — used after cd

                # ── 4. Open fresh SFTP channel ─────────────────
                sftp = ssh.open_sftp()
                try:
                    # Create parent directory if needed
                    try:
                        sftp.stat(resolved_parent)
                    except FileNotFoundError:
                        self._dv_ssh_exec(f"mkdir -p {parent_dir}")
                        try:
                            sftp.stat(resolved_parent)
                        except FileNotFoundError:
                            raise RuntimeError(
                                f"Could not create remote directory:\n"
                                f"{resolved_parent}")

                    # ── 5. Upload batch JSON ───────────────────
                    sftp.put(tmp_path, remote_json_sftp)
                finally:
                    try:
                        sftp.close()
                    except Exception:
                        pass

            finally:
                if os.path.exists(tmp_path):
                    os.unlink(tmp_path)

            # ── 6. Run af3_run ─────────────────────────────────
            # cd to parent_dir so --json_path is filename-only and
            # --job-name creates output inside parent_dir/job_name/output/
            status   = 'uploaded'
            slurm_id = None
            dir_exists_error = False
            if submit:
                mod_prefix = self._dv_build_activation_prefix()

                # Collect AF3 advanced flags — only real af3_run flags:
                #   --slurm-partition  --slurm-nodes  --slurm-ntasks
                #   --slurm-mem  --slurm-gres  --slurm-time
                #   --workdir  --image  --model_dir  --db_dir  --force  --dry-run
                # Seeds & template control are injected into the JSON (see above).
                _af3_flags = []
                if hasattr(self, '_dv_af3_partition'):
                    _pt = self._dv_af3_partition.text().strip()
                    if _pt:
                        _af3_flags.append(f"--slurm-partition {_pt}")
                if hasattr(self, '_dv_af3_extra_flags'):
                    _ex = self._dv_af3_extra_flags.text().strip()
                    if _ex:
                        _af3_flags.append(_ex)

                _flags_str = (" " + " ".join(_af3_flags)) if _af3_flags else ""
                # Use double-quoted bash -lc to safely handle any remaining
                # special characters; job_name is already sanitized above
                _inner_cmd = (
                    f"{mod_prefix}"
                    f"cd {parent_dir} && "
                    f"{cmd} "
                    f"--json_path {remote_json_shell} "
                    f'--job-name "{job_name}"'
                    f"{_flags_str}"
                )
                run_cmd = f'bash -lc "{_inner_cmd}"'
                out, err, rc = self._dv_ssh_exec(run_cmd, timeout=60)
                import re as _re
                m = _re.search(r'batch job\s+(\d+)', out + err, _re.IGNORECASE)
                if m:
                    slurm_id = m.group(1)
                    status   = 'submitted'
                else:
                    server_msg = (out + err).strip().splitlines()
                    detail = server_msg[0] if server_msg else f'rc={rc}'
                    combined = (out + err).lower()
                    if 'already exists' in combined or 'directory already' in combined:
                        dir_exists_error = True
                        status = 'dir_exists'
                    else:
                        status = f'submit_error: {detail}'

            return [{
                'name':            job_name,
                'prefix':          prefix,
                'batch_file':      batch_json_name,
                'n_jobs':          len(batch_list),
                'file_kb':         round(file_kb, 1),
                'remote_json':     f"{parent_dir}/{batch_json_name}",
                'remote_dir':      output_dir,   # where output/ will appear
                'slurm_id':        slurm_id,
                'status':          status,
                'dir_exists_error': dir_exists_error,
                'parent_dir':      parent_dir,
                'rel_json':        remote_json_shell,  # filename only
                'mod_name':        self._dv_module_cmd.text().strip(),
                'cmd':             cmd,
            }]

        self._dv_run_btn.setEnabled(False)
        self._dv_upload_btn.setEnabled(False)
        self._dv_log(f"Starting {'upload+submit' if submit else 'upload only'} "
                     f"({len(jobs)} jobs)...", 'submit')

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
        results is always a list with ONE dict (the batch result)."""
        self._dv_run_btn.setEnabled(True)
        self._dv_upload_btn.setEnabled(True)

        if not results:
            self._dv_log("Upload returned no results.", 'submit')
            return

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
                    import re as _re
                    m = _re.search(
                        r'batch job\s+(\d+)', out + err, _re.IGNORECASE)
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
                        import stat as _stat
                        is_dir = _stat.S_ISDIR(attr.st_mode)
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
                        import stat as _stat
                        attr = sftp.stat(src)
                        if _stat.S_ISDIR(attr.st_mode):
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
        import stat as _stat
        for attr in sftp.listdir_attr(remote_dir):
            rsrc = f"{remote_dir}/{attr.filename}"
            ldst = os.path.join(local_dir, attr.filename)
            if _stat.S_ISDIR(attr.st_mode):
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
if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setApplicationName('ppigFinder')
    app.setApplicationDisplayName('ppigFinder — Protein-Protein Interaction Genomic Finder')
    app.setApplicationVersion('1.01')
    app.setStyle('Fusion')
    _setup_emoji_font(app)
    window = ppigFinderApp()
    window.show()
    sys.exit(app.exec() if QT_VERSION == 6 else app.exec_())
