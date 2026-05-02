#!/usr/bin/env python3
from __future__ import annotations

try:
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import (
        QDialog,
        QVBoxLayout,
        QHBoxLayout,
        QGridLayout,
        QLabel,
        QPushButton,
        QFrame,
    )
    QT6 = True
except Exception:
    from PyQt5.QtCore import Qt
    from PyQt5.QtWidgets import (
        QDialog,
        QVBoxLayout,
        QHBoxLayout,
        QGridLayout,
        QLabel,
        QPushButton,
        QFrame,
    )
    QT6 = False

from ppigfinder.ui_shell.branding import apply_ppigfinder_branding


def _window_flags():
    flags = Qt.WindowType.Window if QT6 else Qt.Window
    flags |= Qt.WindowType.WindowCloseButtonHint if QT6 else Qt.WindowCloseButtonHint
    flags |= Qt.WindowType.WindowMinimizeButtonHint if QT6 else Qt.WindowMinimizeButtonHint
    return flags


class DataStartDialog(QDialog):
    """
    First workflow decision: which data type should enter the analysis?
    """

    def __init__(self, parent=None):
        super().__init__(parent)

        self.choice = ""

        self.setWindowTitle("Start ppigFinder workflow")
        self.setWindowFlags(_window_flags())
        self.resize(980, 560)
        self.setMinimumSize(820, 460)
        apply_ppigfinder_branding(self)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(22, 22, 22, 22)
        layout.setSpacing(16)

        title = QLabel("Start with input data")
        title.setObjectName("SectionTitle")
        layout.addWidget(title)

        subtitle = QLabel(
            "Choose the data type. ppigFinder will enable the next workflow options according "
            "to the input: genome loading, ORF prediction, BLAST query, HMM annotation, AF3 import or project restoration."
        )
        subtitle.setWordWrap(True)
        layout.addWidget(subtitle)

        grid = QGridLayout()
        grid.setHorizontalSpacing(14)
        grid.setVerticalSpacing(14)

        grid.addWidget(
            self._card(
                "Genome file",
                "FASTA, multi-FASTA, GenBank or SnapGene DNA file. Recommended for a new full analysis.",
                "Open genome",
                "genome",
            ),
            0,
            0,
        )

        grid.addWidget(
            self._card(
                "Protein query",
                "Protein FASTA or amino-acid sequence file used later as a BLAST query against predicted ORFs.",
                "Open protein",
                "protein",
            ),
            0,
            1,
        )

        grid.addWidget(
            self._card(
                "HMM profiles",
                "HMM profile database or custom profiles for domain annotation of predicted ORFs.",
                "Open HMM",
                "hmm",
            ),
            0,
            2,
        )

        grid.addWidget(
            self._card(
                "Project",
                "Resume a saved ppigFinder project with previous data and analysis state.",
                "Open project",
                "project",
            ),
            1,
            0,
        )

        grid.addWidget(
            self._card(
                "Snapshot",
                "Import a portable Project Snapshot JSON for reproducibility or continuation.",
                "Import snapshot",
                "snapshot",
            ),
            1,
            1,
        )

        grid.addWidget(
            self._card(
                "AF3 results",
                "Import AlphaFold 3 output folders for metrics and interaction review.",
                "Import results",
                "af3_results",
            ),
            1,
            2,
        )

        layout.addLayout(grid, 1)

        hint = QLabel(
            "A genome starts the complete workflow. Protein and HMM inputs are stored as annotation inputs "
            "and become active after ORF prediction."
        )
        hint.setWordWrap(True)
        hint.setObjectName("InfoFooter")
        layout.addWidget(hint)

    def _card(self, title: str, description: str, button_text: str, choice: str) -> QFrame:
        frame = QFrame()
        frame.setObjectName("InfoCard")

        layout = QVBoxLayout(frame)
        layout.setContentsMargins(14, 12, 14, 12)
        layout.setSpacing(8)

        title_label = QLabel(title)
        title_label.setObjectName("InfoCardTitle")

        desc_label = QLabel(description)
        desc_label.setObjectName("InfoCardSubtitle")
        desc_label.setWordWrap(True)

        button = QPushButton(button_text)
        button.clicked.connect(lambda: self._accept_choice(choice))

        layout.addWidget(title_label)
        layout.addWidget(desc_label)
        layout.addStretch(1)
        layout.addWidget(button)

        return frame

    def _accept_choice(self, choice: str) -> None:
        self.choice = choice
        self.accept()


def choose_initial_data(parent=None) -> str:
    dialog = DataStartDialog(parent=parent)
    dialog.exec() if hasattr(dialog, "exec") else dialog.exec_()
    return dialog.choice
