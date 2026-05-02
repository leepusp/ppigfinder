#!/usr/bin/env python3
"""
Genome inspector dialog for the experimental guided UI shell.

Shows metadata, workflow options and an aligned coordinate-based sequence
preview after a genome file is loaded.
"""

from __future__ import annotations

from pathlib import Path
import re

try:
    from PyQt6.QtCore import Qt
    from PyQt6.QtGui import QFont
    from PyQt6.QtWidgets import (
        QDialog,
        QVBoxLayout,
        QHBoxLayout,
        QLabel,
        QTableWidget,
        QTableWidgetItem,
        QHeaderView,
        QSplitter,
        QWidget,
        QPushButton,
        QFrame,
        QGridLayout,
    )
    QT6 = True
except Exception:
    from PyQt5.QtCore import Qt
    from PyQt5.QtGui import QFont
    from PyQt5.QtWidgets import (
        QDialog,
        QVBoxLayout,
        QHBoxLayout,
        QLabel,
        QTableWidget,
        QTableWidgetItem,
        QHeaderView,
        QSplitter,
        QWidget,
        QPushButton,
        QFrame,
        QGridLayout,
    )
    QT6 = False

from ppigfinder.ui_shell.branding import apply_ppigfinder_branding
from ppigfinder.ui_shell.input_validation import validate_genome_input


def _window_flags():
    flags = Qt.WindowType.Window if QT6 else Qt.Window
    flags |= Qt.WindowType.WindowMinimizeButtonHint if QT6 else Qt.WindowMinimizeButtonHint
    flags |= Qt.WindowType.WindowMaximizeButtonHint if QT6 else Qt.WindowMaximizeButtonHint
    flags |= Qt.WindowType.WindowCloseButtonHint if QT6 else Qt.WindowCloseButtonHint
    return flags


def _read_text(path: str | Path) -> str:
    path = Path(path)

    for encoding in ("utf-8", "latin-1", "cp1252"):
        try:
            return path.read_text(encoding=encoding)
        except UnicodeDecodeError:
            continue

    return path.read_text(errors="ignore")


def _sequence_preview_from_fasta(text: str, max_chars: int = 12000) -> tuple[str, list[str]]:
    headers = []
    seq_chunks = []

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            headers.append(line[1:].strip())
        else:
            seq_chunks.append(re.sub(r"[^A-Za-z]", "", line))

    sequence = "".join(seq_chunks).upper()
    return sequence[:max_chars], headers[:20]


def _sequence_preview_from_genbank(text: str, max_chars: int = 12000) -> str:
    idx = text.find("ORIGIN")

    if idx == -1:
        return ""

    origin = text[idx:]
    sequence = re.sub(r"[^A-Za-z]", "", origin.replace("ORIGIN", ""))
    sequence = sequence.replace("END", "")
    return sequence.upper()[:max_chars]


def _organism_from_genbank(text: str) -> str:
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("ORGANISM"):
            return stripped.replace("ORGANISM", "", 1).strip()
        if "/organism=" in stripped:
            return stripped.split("=", 1)[1].strip().strip('"')
    return "N/A"


def _group_sequence(sequence: str, block_size: int = 10) -> str:
    return " ".join(
        sequence[i:i + block_size]
        for i in range(0, len(sequence), block_size)
    )


def inspect_genome_file(path: str) -> dict:
    summary = validate_genome_input(path)
    text = _read_text(path)
    suffix = Path(path).suffix.lower()

    organism = "N/A"
    headers = []
    preview = ""

    if suffix in {".fa", ".fasta", ".fna", ".ffn"}:
        preview, headers = _sequence_preview_from_fasta(text)
        if headers:
            organism = headers[0]
    elif suffix in {".gb", ".gbk", ".genbank"}:
        preview = _sequence_preview_from_genbank(text)
        organism = _organism_from_genbank(text)
    else:
        preview = re.sub(r"[^A-Za-z]", "", text[:12000]).upper()

    return {
        "path": str(path),
        "name": summary.name,
        "file_type": summary.file_type,
        "organism": organism,
        "sequence_count": summary.sequence_count,
        "total_length": summary.total_length,
        "longest_length": summary.longest_length,
        "gc_percent": summary.gc_percent,
        "valid": summary.valid,
        "message": summary.message,
        "headers": headers,
        "sequence_preview": preview,
    }


class SummaryCard(QFrame):
    def __init__(self, title: str, value: str, subtitle: str = "", parent=None):
        super().__init__(parent)
        self.setObjectName("InfoCard")

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 10, 12, 10)
        layout.setSpacing(5)

        title_label = QLabel(title)
        title_label.setObjectName("InfoCardTitle")

        value_label = QLabel(value)
        value_label.setObjectName("InfoCardValue")
        value_label.setWordWrap(True)

        subtitle_label = QLabel(subtitle)
        subtitle_label.setObjectName("InfoCardSubtitle")
        subtitle_label.setWordWrap(True)

        layout.addWidget(title_label)
        layout.addWidget(value_label)
        layout.addWidget(subtitle_label)


class GenomeInspectorDialog(QDialog):
    """
    Genome data inspector with metadata, coordinate-based sequence table and next workflow actions.
    """

    def __init__(self, inspection: dict, parent=None):
        super().__init__(parent)

        self.inspection = dict(inspection)
        self.workflow_parent = parent

        self.setWindowTitle("Genome data loaded")
        self.setWindowFlags(_window_flags())
        self.resize(1220, 820)
        self.setMinimumSize(980, 640)
        self.setSizeGripEnabled(True)
        apply_ppigfinder_branding(self)

        root = QVBoxLayout(self)
        root.setContentsMargins(18, 18, 18, 18)
        root.setSpacing(12)

        title = QLabel("Genome data loaded")
        title.setObjectName("SectionTitle")
        root.addWidget(title)

        subtitle = QLabel(
            "The selected genome was registered in the guided workflow. "
            "Review metadata and the coordinate-aligned sequence preview, then continue with ORF prediction or add annotation inputs."
        )
        subtitle.setWordWrap(True)
        root.addWidget(subtitle)

        cards = QGridLayout()
        cards.setHorizontalSpacing(10)
        cards.setVerticalSpacing(10)

        cards.addWidget(SummaryCard("Input type", str(self.inspection.get("file_type", "N/A")), "Detected file format"), 0, 0)
        cards.addWidget(SummaryCard("Total length", f"{int(self.inspection.get('total_length') or 0):,} bp", "Nucleotide length"), 0, 1)
        cards.addWidget(SummaryCard("GC%", str(self.inspection.get("gc_percent") if self.inspection.get("gc_percent") is not None else "N/A"), "Genome composition"), 0, 2)
        cards.addWidget(SummaryCard("Validation", "Valid" if self.inspection.get("valid") else "Problem", self.inspection.get("message", "")), 0, 3)

        root.addLayout(cards)

        splitter = QSplitter()
        splitter.setChildrenCollapsible(False)
        splitter.setOrientation(Qt.Orientation.Vertical if QT6 else Qt.Vertical)

        upper = QSplitter()
        upper.setChildrenCollapsible(False)
        upper.setOrientation(Qt.Orientation.Horizontal if QT6 else Qt.Horizontal)

        metadata_widget = QWidget()
        metadata_layout = QVBoxLayout(metadata_widget)
        metadata_layout.setContentsMargins(0, 0, 0, 0)

        metadata_title = QLabel("Metadata")
        metadata_title.setObjectName("SectionSubTitle")
        metadata_layout.addWidget(metadata_title)

        self.table = QTableWidget(0, 2)
        self.table.setHorizontalHeaderLabels(["Field", "Value"])
        header = self.table.horizontalHeader()
        try:
            header.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents if QT6 else QHeaderView.ResizeToContents)
            header.setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch if QT6 else QHeaderView.Stretch)
        except Exception:
            pass
        self.table.setAlternatingRowColors(True)
        metadata_layout.addWidget(self.table)

        actions_widget = QWidget()
        actions_layout = QVBoxLayout(actions_widget)
        actions_layout.setContentsMargins(8, 0, 0, 0)
        actions_layout.setSpacing(10)

        actions_title = QLabel("Next workflow options")
        actions_title.setObjectName("SectionSubTitle")
        actions_layout.addWidget(actions_title)

        hint = QLabel(
            "A valid genome enables ORF prediction. You may also add protein query input for BLAST "
            "or HMM profiles for domain annotation before/after ORF prediction."
        )
        hint.setWordWrap(True)
        actions_layout.addWidget(hint)

        btn_orfs = QPushButton("Go to Protein / ORFs")
        btn_orfs.clicked.connect(self._go_to_orfs)
        actions_layout.addWidget(btn_orfs)

        btn_predict = QPushButton("Predict ORFs now")
        btn_predict.clicked.connect(self._predict_orfs_now)
        actions_layout.addWidget(btn_predict)

        btn_query = QPushButton("Add protein query for BLAST")
        btn_query.clicked.connect(self._add_protein_query)
        actions_layout.addWidget(btn_query)

        btn_hmm = QPushButton("Add HMM profiles")
        btn_hmm.clicked.connect(self._add_hmm_profiles)
        actions_layout.addWidget(btn_hmm)

        btn_data = QPushButton("Return to Data / Project")
        btn_data.clicked.connect(self._go_to_data)
        actions_layout.addWidget(btn_data)

        actions_layout.addStretch(1)

        upper.addWidget(metadata_widget)
        upper.addWidget(actions_widget)
        upper.setStretchFactor(0, 3)
        upper.setStretchFactor(1, 1)

        sequence_widget = QWidget()
        sequence_layout = QVBoxLayout(sequence_widget)
        sequence_layout.setContentsMargins(0, 0, 0, 0)

        sequence_title = QLabel("Coordinate sequence preview")
        sequence_title.setObjectName("SectionSubTitle")
        sequence_layout.addWidget(sequence_title)

        self.sequence_table = QTableWidget(0, 3)
        self.sequence_table.setHorizontalHeaderLabels(["Start", "Sequence", "End"])
        self.sequence_table.setAlternatingRowColors(True)

        seq_header = self.sequence_table.horizontalHeader()
        try:
            seq_header.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents if QT6 else QHeaderView.ResizeToContents)
            seq_header.setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch if QT6 else QHeaderView.Stretch)
            seq_header.setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents if QT6 else QHeaderView.ResizeToContents)
        except Exception:
            pass

        font = QFont("Courier New")
        font.setPointSize(10)
        self.sequence_table.setFont(font)

        sequence_layout.addWidget(self.sequence_table)

        splitter.addWidget(upper)
        splitter.addWidget(sequence_widget)
        splitter.setStretchFactor(0, 2)
        splitter.setStretchFactor(1, 3)

        root.addWidget(splitter, 1)

        self._populate_metadata()
        self._populate_sequence_table()

    def _populate_metadata(self) -> None:
        rows = [
            ("Path", self.inspection.get("path", "")),
            ("Name / locus", self.inspection.get("name", "N/A")),
            ("File type", self.inspection.get("file_type", "N/A")),
            ("Organism / first header", self.inspection.get("organism", "N/A")),
            ("Valid", "Yes" if self.inspection.get("valid") else "No"),
            ("Sequence count", str(self.inspection.get("sequence_count") or "N/A")),
            ("Total length", f"{int(self.inspection.get('total_length') or 0):,} bp"),
            ("Longest sequence", f"{int(self.inspection.get('longest_length') or 0):,} bp"),
            ("GC%", str(self.inspection.get("gc_percent") if self.inspection.get("gc_percent") is not None else "N/A")),
            ("Message", self.inspection.get("message", "")),
        ]

        headers = self.inspection.get("headers") or []
        if headers:
            rows.append(("FASTA headers", "\n".join(headers[:10])))

        self.table.setRowCount(len(rows))

        for row, (field, value) in enumerate(rows):
            self.table.setItem(row, 0, QTableWidgetItem(field))
            self.table.setItem(row, 1, QTableWidgetItem(str(value)))

        self.table.resizeRowsToContents()

    def _populate_sequence_table(self) -> None:
        sequence = re.sub(r"[^A-Za-z]", "", self.inspection.get("sequence_preview", "")).upper()
        width = 60
        rows = []

        for i in range(0, len(sequence), width):
            chunk = sequence[i:i + width]
            rows.append((i + 1, _group_sequence(chunk), i + len(chunk)))

        self.sequence_table.setRowCount(len(rows))

        for row, (start, grouped, end) in enumerate(rows):
            self.sequence_table.setItem(row, 0, QTableWidgetItem(str(start)))
            self.sequence_table.setItem(row, 1, QTableWidgetItem(grouped))
            self.sequence_table.setItem(row, 2, QTableWidgetItem(str(end)))

        self.sequence_table.resizeRowsToContents()

    def _go_to_orfs(self) -> None:
        if hasattr(self.workflow_parent, "show_route"):
            self.workflow_parent.show_route("orfs")

    def _go_to_data(self) -> None:
        if hasattr(self.workflow_parent, "show_route"):
            self.workflow_parent.show_route("data")

    def _predict_orfs_now(self) -> None:
        if hasattr(self.workflow_parent, "_predict_guided_orfs"):
            self.workflow_parent._predict_guided_orfs()

    def _add_protein_query(self) -> None:
        if hasattr(self.workflow_parent, "_select_file"):
            self.workflow_parent._select_file(
                "protein_query_file",
                "Open protein query FASTA",
                "Protein FASTA (*.faa *.fa *.fasta *.txt);;All files (*)",
            )

    def _add_hmm_profiles(self) -> None:
        if hasattr(self.workflow_parent, "_select_file"):
            self.workflow_parent._select_file(
                "hmm_profile_file",
                "Open HMM profile file",
                "HMM profiles (*.hmm);;All files (*)",
            )


def show_genome_inspector(path: str, parent=None) -> GenomeInspectorDialog:
    inspection = inspect_genome_file(path)
    dialog = GenomeInspectorDialog(inspection, parent=parent)
    dialog.show()
    return dialog
