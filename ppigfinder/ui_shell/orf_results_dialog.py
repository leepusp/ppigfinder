#!/usr/bin/env python3
"""
ORF Discovery results window for the experimental guided UI shell.

This dialog turns ORF prediction into a visual workflow result:
- summary cards
- ORF table
- strand/frame statistics
- compact genome/ORF map preview
- next workflow actions
"""

from __future__ import annotations

from collections import Counter
from statistics import mean

try:
    from PyQt6.QtCore import Qt, QRectF
    from PyQt6.QtGui import QPainter, QColor, QPen, QBrush, QPolygonF, QFont
    from PyQt6.QtCore import QPointF
    from PyQt6.QtWidgets import (
        QDialog,
        QVBoxLayout,
        QHBoxLayout,
        QGridLayout,
        QLabel,
        QPushButton,
        QFrame,
        QSplitter,
        QWidget,
        QTableWidget,
        QTableWidgetItem,
        QHeaderView,
    )
    QT6 = True
except Exception:
    from PyQt5.QtCore import Qt, QRectF, QPointF
    from PyQt5.QtGui import QPainter, QColor, QPen, QBrush, QPolygonF, QFont
    from PyQt5.QtWidgets import (
        QDialog,
        QVBoxLayout,
        QHBoxLayout,
        QGridLayout,
        QLabel,
        QPushButton,
        QFrame,
        QSplitter,
        QWidget,
        QTableWidget,
        QTableWidgetItem,
        QHeaderView,
    )
    QT6 = False

from ppigfinder.ui_shell.branding import apply_ppigfinder_branding


def _window_flags():
    flags = Qt.WindowType.Window if QT6 else Qt.Window
    flags |= Qt.WindowType.WindowMinimizeButtonHint if QT6 else Qt.WindowMinimizeButtonHint
    flags |= Qt.WindowType.WindowMaximizeButtonHint if QT6 else Qt.WindowMaximizeButtonHint
    flags |= Qt.WindowType.WindowCloseButtonHint if QT6 else Qt.WindowCloseButtonHint
    return flags


def _alignment_center():
    return Qt.AlignmentFlag.AlignCenter if QT6 else Qt.AlignCenter


def _safe_int(value, default=0):
    try:
        return int(value)
    except Exception:
        return default


def _orf_length(orf) -> int:
    if hasattr(orf, "aa_length"):
        return _safe_int(orf.aa_length)
    return len(getattr(orf, "protein_sequence", "") or "")


def _protein_preview(orf, length: int = 48) -> str:
    seq = getattr(orf, "protein_sequence", "") or ""
    if len(seq) <= length:
        return seq
    return seq[:length] + "..."


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


class ORFMapWidget(QWidget):
    """
    Compact ORF map preview.

    Plus strand is drawn above the axis, minus strand below the axis.
    """

    def __init__(self, orfs, genome_length: int = 0, parent=None):
        super().__init__(parent)
        self.orfs = list(orfs or [])
        self.genome_length = int(genome_length or self._infer_length())
        self.setMinimumHeight(210)

    def _infer_length(self) -> int:
        if not self.orfs:
            return 0
        return max(_safe_int(getattr(orf, "end", 0)) for orf in self.orfs)

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing if QT6 else QPainter.Antialiasing)

        width = self.width()
        height = self.height()
        margin_x = 36
        axis_y = height // 2
        usable = max(1, width - 2 * margin_x)

        painter.fillRect(self.rect(), QColor("#f8fbfc"))

        painter.setPen(QPen(QColor("#607d8b"), 2))
        painter.drawLine(margin_x, axis_y, width - margin_x, axis_y)

        if not self.orfs or self.genome_length <= 0:
            painter.setPen(QColor("#607d8b"))
            painter.drawText(self.rect(), _alignment_center(), "No ORF map available.")
            painter.end()
            return

        visible = self.orfs[:160]

        plus_color = QColor("#43a047")
        minus_color = QColor("#1e88e5")
        border = QColor("#263238")

        for orf in visible:
            start = _safe_int(getattr(orf, "start", 0))
            end = _safe_int(getattr(orf, "end", 0))
            strand = getattr(orf, "strand", "+")

            if end < start:
                start, end = end, start

            x1 = margin_x + (start / self.genome_length) * usable
            x2 = margin_x + (end / self.genome_length) * usable

            if x2 - x1 < 5:
                x2 = x1 + 5

            y = axis_y - 42 if strand == "+" else axis_y + 24
            h = 18
            arrow = min(13, max(5, (x2 - x1) * 0.35))

            if strand == "+":
                points = [
                    QPointF(x1, y),
                    QPointF(x2 - arrow, y),
                    QPointF(x2, y + h / 2),
                    QPointF(x2 - arrow, y + h),
                    QPointF(x1, y + h),
                ]
                color = plus_color
            else:
                points = [
                    QPointF(x2, y),
                    QPointF(x1 + arrow, y),
                    QPointF(x1, y + h / 2),
                    QPointF(x1 + arrow, y + h),
                    QPointF(x2, y + h),
                ]
                color = minus_color

            painter.setPen(QPen(border, 0.5))
            painter.setBrush(QBrush(color))
            painter.drawPolygon(QPolygonF(points))

        painter.setPen(QColor("#455a64"))
        painter.drawText(36, 24, f"Genome length: {self.genome_length:,} nt")
        painter.drawText(36, height - 16, f"Showing {len(visible)} of {len(self.orfs)} ORFs")

        painter.setPen(QColor("#43a047"))
        painter.drawText(width - 190, 24, "+ strand")
        painter.setPen(QColor("#1e88e5"))
        painter.drawText(width - 105, 24, "- strand")

        painter.end()


class ORFStatsWidget(QWidget):
    """
    Simple text-based visual summary for strand/frame counts.
    """

    def __init__(self, orfs, parent=None):
        super().__init__(parent)
        self.orfs = list(orfs or [])

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)

        title = QLabel("Prediction summary")
        title.setObjectName("SectionSubTitle")
        layout.addWidget(title)

        self.table = QTableWidget(0, 2)
        self.table.setHorizontalHeaderLabels(["Metric", "Value"])
        self.table.verticalHeader().setVisible(False)
        self.table.setAlternatingRowColors(True)

        header = self.table.horizontalHeader()
        try:
            header.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents if QT6 else QHeaderView.ResizeToContents)
            header.setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch if QT6 else QHeaderView.Stretch)
        except Exception:
            pass

        layout.addWidget(self.table)
        self._populate()

    def _populate(self):
        strands = Counter(getattr(orf, "strand", "?") for orf in self.orfs)
        frames = Counter(str(getattr(orf, "frame", "?")) for orf in self.orfs)
        lengths = [_orf_length(orf) for orf in self.orfs]

        rows = [
            ("ORFs predicted", str(len(self.orfs))),
            ("Plus strand", str(strands.get("+", 0))),
            ("Minus strand", str(strands.get("-", 0))),
            ("Frames", ", ".join(f"{k}: {v}" for k, v in sorted(frames.items())) or "N/A"),
            ("Mean length", f"{mean(lengths):.1f} aa" if lengths else "N/A"),
            ("Longest", f"{max(lengths)} aa" if lengths else "N/A"),
            ("Shortest", f"{min(lengths)} aa" if lengths else "N/A"),
        ]

        self.table.setRowCount(len(rows))

        for row, (key, value) in enumerate(rows):
            self.table.setItem(row, 0, QTableWidgetItem(key))
            self.table.setItem(row, 1, QTableWidgetItem(value))

        self.table.resizeRowsToContents()


class GuidedORFResultsDialog(QDialog):
    """
    Rich ORF Discovery results window.
    """

    def __init__(self, orfs, parent=None):
        super().__init__(parent)

        self.orfs = list(orfs or [])
        self.workflow_parent = parent

        lengths = [_orf_length(orf) for orf in self.orfs]
        plus = sum(1 for orf in self.orfs if getattr(orf, "strand", "") == "+")
        minus = sum(1 for orf in self.orfs if getattr(orf, "strand", "") == "-")
        genome_length = self._infer_genome_length()

        self.setWindowTitle("ORF Discovery Results")
        self.setWindowFlags(_window_flags())
        self.resize(1280, 820)
        self.setMinimumSize(980, 640)
        self.setSizeGripEnabled(True)
        apply_ppigfinder_branding(self)

        root = QVBoxLayout(self)
        root.setContentsMargins(18, 18, 18, 18)
        root.setSpacing(12)

        title = QLabel("ORF Discovery Results")
        title.setObjectName("SectionTitle")
        root.addWidget(title)

        subtitle = QLabel(
            "Predicted ORFs were generated from the loaded genome and are now available "
            "for BLAST, HMM/domain annotation, neighbourhood analysis and AlphaFold/PPI candidate selection."
        )
        subtitle.setWordWrap(True)
        root.addWidget(subtitle)

        cards = QGridLayout()
        cards.setHorizontalSpacing(10)
        cards.setVerticalSpacing(10)

        cards.addWidget(SummaryCard("ORFs predicted", f"{len(self.orfs):,}", "Protein-coding candidates"), 0, 0)
        cards.addWidget(SummaryCard("Mean length", f"{mean(lengths):.1f} aa" if lengths else "N/A", "Average protein length"), 0, 1)
        cards.addWidget(SummaryCard("Longest ORF", f"{max(lengths)} aa" if lengths else "N/A", "Maximum protein length"), 0, 2)
        cards.addWidget(SummaryCard("Strands", f"+ {plus:,} / - {minus:,}", "Predicted orientation"), 0, 3)

        root.addLayout(cards)

        splitter = QSplitter()
        splitter.setChildrenCollapsible(False)
        splitter.setOrientation(Qt.Orientation.Horizontal if QT6 else Qt.Horizontal)

        # Left: table
        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 8, 0)
        left_layout.setSpacing(8)

        table_title = QLabel("ORF table")
        table_title.setObjectName("SectionSubTitle")
        left_layout.addWidget(table_title)

        self.table = QTableWidget(0, 8, self)
        self.table.setHorizontalHeaderLabels(
            [
                "ORF ID",
                "Start",
                "End",
                "Strand",
                "Frame",
                "AA length",
                "NT length",
                "Protein preview",
            ]
        )
        self.table.verticalHeader().setVisible(False)
        self.table.setAlternatingRowColors(True)

        header = self.table.horizontalHeader()
        try:
            header.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents if QT6 else QHeaderView.ResizeToContents)
            header.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents if QT6 else QHeaderView.ResizeToContents)
            header.setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents if QT6 else QHeaderView.ResizeToContents)
            header.setSectionResizeMode(3, QHeaderView.ResizeMode.ResizeToContents if QT6 else QHeaderView.ResizeToContents)
            header.setSectionResizeMode(4, QHeaderView.ResizeMode.ResizeToContents if QT6 else QHeaderView.ResizeToContents)
            header.setSectionResizeMode(5, QHeaderView.ResizeMode.ResizeToContents if QT6 else QHeaderView.ResizeToContents)
            header.setSectionResizeMode(6, QHeaderView.ResizeMode.ResizeToContents if QT6 else QHeaderView.ResizeToContents)
            header.setSectionResizeMode(7, QHeaderView.ResizeMode.Stretch if QT6 else QHeaderView.Stretch)
        except Exception:
            pass

        font = QFont("Courier New")
        font.setPointSize(9)
        self.table.setFont(font)

        left_layout.addWidget(self.table, 1)

        # Right: map + stats + next actions
        right = QWidget()
        right_layout = QVBoxLayout(right)
        right_layout.setContentsMargins(8, 0, 0, 0)
        right_layout.setSpacing(10)

        map_title = QLabel("ORF map preview")
        map_title.setObjectName("SectionSubTitle")
        right_layout.addWidget(map_title)

        self.map_widget = ORFMapWidget(self.orfs, genome_length=genome_length)
        right_layout.addWidget(self.map_widget, 2)

        self.stats_widget = ORFStatsWidget(self.orfs)
        right_layout.addWidget(self.stats_widget, 2)

        next_title = QLabel("Next workflow options")
        next_title.setObjectName("SectionSubTitle")
        right_layout.addWidget(next_title)

        btn_annotation = QPushButton("Go to Annotation")
        btn_annotation.clicked.connect(self._go_to_annotation)
        right_layout.addWidget(btn_annotation)

        btn_candidates = QPushButton("Review candidate ORFs")
        btn_candidates.clicked.connect(self._review_candidates)
        right_layout.addWidget(btn_candidates)

        btn_export = QPushButton("Export ORF FASTA")
        btn_export.clicked.connect(self._export_fasta)
        right_layout.addWidget(btn_export)

        btn_af3 = QPushButton("Continue to AlphaFold / PPI")
        btn_af3.clicked.connect(self._go_to_alphafold)
        right_layout.addWidget(btn_af3)

        splitter.addWidget(left)
        splitter.addWidget(right)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 2)

        root.addWidget(splitter, 1)

        self._populate()

    def _infer_genome_length(self) -> int:
        if not self.orfs:
            return 0
        return max(_safe_int(getattr(orf, "end", 0)) for orf in self.orfs)

    def _populate(self) -> None:
        self.table.setRowCount(len(self.orfs))

        for row, orf in enumerate(self.orfs):
            values = [
                getattr(orf, "id", f"orf_{row + 1}"),
                str(getattr(orf, "start", "")),
                str(getattr(orf, "end", "")),
                str(getattr(orf, "strand", "")),
                str(getattr(orf, "frame", "")),
                str(_orf_length(orf)),
                str(len(getattr(orf, "nt_sequence", "") or "")),
                _protein_preview(orf),
            ]

            for col, value in enumerate(values):
                item = QTableWidgetItem(value)
                if col in {1, 2, 4, 5, 6}:
                    item.setTextAlignment(_alignment_center())
                self.table.setItem(row, col, item)

        self.table.resizeRowsToContents()

    def _go_to_annotation(self):
        if hasattr(self.workflow_parent, "show_route"):
            self.workflow_parent.show_route("annotation")

    def _go_to_alphafold(self):
        if hasattr(self.workflow_parent, "show_route"):
            self.workflow_parent.show_route("alphafold")

    def _review_candidates(self):
        if hasattr(self.workflow_parent, "_show_annotation_candidates"):
            self.workflow_parent._show_annotation_candidates()

    def _export_fasta(self):
        if hasattr(self.workflow_parent, "_export_guided_orfs_fasta"):
            self.workflow_parent._export_guided_orfs_fasta()


def show_guided_orf_results(orfs, parent=None) -> None:
    dialog = GuidedORFResultsDialog(orfs, parent=parent)
    dialog.showMaximized()
    dialog.exec() if hasattr(dialog, "exec") else dialog.exec_()
