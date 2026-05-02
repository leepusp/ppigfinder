#!/usr/bin/env python3
"""
Annotation candidate table for the experimental guided UI shell.

This dialog uses guided ORF predictions as candidate entries for annotation
and downstream AlphaFold/PPI prioritization.
"""

from __future__ import annotations

try:
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import (
        QDialog,
        QVBoxLayout,
        QLabel,
        QTableWidget,
        QTableWidgetItem,
        QHeaderView,
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
    )
    QT6 = False

from ppigfinder.ui_shell.branding import apply_ppigfinder_branding


def _window_flags():
    flags = Qt.WindowType.Window if QT6 else Qt.Window
    flags |= Qt.WindowType.WindowMinimizeButtonHint if QT6 else Qt.WindowMinimizeButtonHint
    flags |= Qt.WindowType.WindowMaximizeButtonHint if QT6 else Qt.WindowMaximizeButtonHint
    flags |= Qt.WindowType.WindowCloseButtonHint if QT6 else Qt.WindowCloseButtonHint
    return flags


class AnnotationCandidatesDialog(QDialog):
    """
    Candidate ORF table for annotation/PPI prioritization.
    """

    def __init__(self, orfs, parent=None):
        super().__init__(parent)

        self.orfs = list(orfs or [])

        self.setWindowTitle("Guided Annotation Candidates")
        self.setWindowFlags(_window_flags())
        self.resize(1150, 720)
        self.setMinimumSize(860, 560)
        self.setSizeGripEnabled(True)
        apply_ppigfinder_branding(self)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(18, 18, 18, 18)
        layout.setSpacing(12)

        title = QLabel(f"Candidate ORFs for annotation: {len(self.orfs)}")
        title.setObjectName("SectionTitle")
        layout.addWidget(title)

        subtitle = QLabel(
            "This table currently lists guided ORF predictions. BLAST, HMM/domain "
            "and neighbourhood evidence will be connected progressively to refine "
            "candidate prioritization."
        )
        subtitle.setWordWrap(True)
        layout.addWidget(subtitle)

        self.table = QTableWidget(0, 9, self)
        self.table.setHorizontalHeaderLabels(
            [
                "Candidate",
                "ORF ID",
                "Start",
                "End",
                "Strand",
                "Frame",
                "AA length",
                "Annotation status",
                "Suggested next step",
            ]
        )

        header = self.table.horizontalHeader()
        try:
            mode = QHeaderView.ResizeMode.Stretch if QT6 else QHeaderView.Stretch
            header.setSectionResizeMode(mode)
        except Exception:
            pass

        self.table.setAlternatingRowColors(True)
        self._populate()
        layout.addWidget(self.table)

    def _populate(self) -> None:
        self.table.setRowCount(len(self.orfs))

        for row, orf in enumerate(self.orfs):
            values = [
                str(row + 1),
                orf.id,
                str(orf.start),
                str(orf.end),
                orf.strand,
                str(orf.frame),
                str(orf.aa_length),
                "Pending BLAST/HMM",
                "Review neighbourhood / select for AF3",
            ]

            for col, value in enumerate(values):
                item = QTableWidgetItem(value)
                self.table.setItem(row, col, item)


def show_annotation_candidates(orfs, parent=None) -> None:
    dialog = AnnotationCandidatesDialog(orfs, parent=parent)
    dialog.showMaximized()
    dialog.exec() if hasattr(dialog, "exec") else dialog.exec_()
