#!/usr/bin/env python3
"""
Genome inspector dialog for the experimental guided UI shell.

This window opens automatically after selecting a genome file. It allows the
user to keep seeing what was loaded while the workflow advances to ORF options.
"""

from __future__ import annotations

from pathlib import Path
import re

try:
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import (
        QDialog,
        QVBoxLayout,
        QLabel,
        QTableWidget,
        QTableWidgetItem,
        QHeaderView,
        QTextEdit,
        QSplitter,
        QWidget,
    )
    QT6 = True
except Exception:
    from PyQt5.QtCore import Qt
    from PyQt5.QtWidgets import (
        QDialog,
        QVBoxLayout,
        QLabel,
        QTableWidget,
        QTableWidgetItem,
        QHeaderView,
        QTextEdit,
        QSplitter,
        QWidget,
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


def _sequence_preview_from_fasta(text: str, max_chars: int = 6000) -> tuple[str, list[str]]:
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


def _sequence_preview_from_genbank(text: str, max_chars: int = 6000) -> str:
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


def _format_sequence(sequence: str, width: int = 80) -> str:
    if not sequence:
        return "No sequence preview available."

    lines = []
    for i in range(0, len(sequence), width):
        lines.append(sequence[i:i + width])
    return "\n".join(lines)


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
        preview = text[:6000]

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


class GenomeInspectorDialog(QDialog):
    """
    Genome data inspector.
    """

    def __init__(self, inspection: dict, parent=None):
        super().__init__(parent)

        self.inspection = dict(inspection)

        self.setWindowTitle("Genome data loaded")
        self.setWindowFlags(_window_flags())
        self.resize(1100, 760)
        self.setMinimumSize(860, 580)
        self.setSizeGripEnabled(True)
        apply_ppigfinder_branding(self)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(18, 18, 18, 18)
        layout.setSpacing(12)

        title = QLabel("Genome data loaded")
        title.setObjectName("SectionTitle")
        layout.addWidget(title)

        subtitle = QLabel(
            "The selected genome was registered in the guided workflow. "
            "The workspace will move to ORF Discovery while this window keeps the loaded data visible."
        )
        subtitle.setWordWrap(True)
        layout.addWidget(subtitle)

        splitter = QSplitter()
        splitter.setOrientation(Qt.Orientation.Vertical if QT6 else Qt.Vertical)

        meta_widget = QWidget()
        meta_layout = QVBoxLayout(meta_widget)

        self.table = QTableWidget(0, 2)
        self.table.setHorizontalHeaderLabels(["Field", "Value"])
        header = self.table.horizontalHeader()
        try:
            mode = QHeaderView.ResizeMode.Stretch if QT6 else QHeaderView.Stretch
            header.setSectionResizeMode(mode)
        except Exception:
            pass
        self.table.setAlternatingRowColors(True)
        meta_layout.addWidget(self.table)

        self.preview = QTextEdit()
        self.preview.setReadOnly(True)
        self.preview.setPlainText(_format_sequence(self.inspection.get("sequence_preview", "")))

        splitter.addWidget(meta_widget)
        splitter.addWidget(self.preview)
        splitter.setSizes([300, 460])

        layout.addWidget(splitter, 1)

        self._populate_metadata()

    def _populate_metadata(self) -> None:
        rows = [
            ("Path", self.inspection.get("path", "")),
            ("Name / locus", self.inspection.get("name", "N/A")),
            ("File type", self.inspection.get("file_type", "N/A")),
            ("Organism / first header", self.inspection.get("organism", "N/A")),
            ("Valid", "Yes" if self.inspection.get("valid") else "No"),
            ("Sequence count", str(self.inspection.get("sequence_count") or "N/A")),
            ("Total length", str(self.inspection.get("total_length") or "N/A")),
            ("Longest sequence", str(self.inspection.get("longest_length") or "N/A")),
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


def show_genome_inspector(path: str, parent=None) -> GenomeInspectorDialog:
    inspection = inspect_genome_file(path)
    dialog = GenomeInspectorDialog(inspection, parent=parent)
    dialog.show()
    return dialog
