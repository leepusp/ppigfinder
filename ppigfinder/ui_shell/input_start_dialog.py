#!/usr/bin/env python3
from __future__ import annotations

try:
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import (
        QDialog,
        QVBoxLayout,
        QHBoxLayout,
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
    Unified first input dialog.

    The user selects one input file and ppigFinder detects whether it is a
    genome, protein query, HMM profile, project or snapshot.
    """

    def __init__(self, parent=None):
        super().__init__(parent)

        self.choice = ""

        self.setWindowTitle("Start ppigFinder workflow")
        self.setWindowFlags(_window_flags())
        self.resize(760, 420)
        self.setMinimumSize(660, 360)
        apply_ppigfinder_branding(self)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(22, 22, 22, 22)
        layout.setSpacing(16)

        title = QLabel("Add input data")
        title.setObjectName("SectionTitle")
        layout.addWidget(title)

        subtitle = QLabel(
            "Select an input file and ppigFinder will detect how it fits in the workflow. "
            "Genome files start the full analysis; protein queries and HMM profiles are stored "
            "for annotation after ORF prediction; project/snapshot files restore previous work."
        )
        subtitle.setWordWrap(True)
        layout.addWidget(subtitle)

        row = QHBoxLayout()
        row.setSpacing(14)

        row.addWidget(
            self._card(
                "Select input file",
                "Accepted examples: FASTA, GenBank, SnapGene, protein FASTA, HMM profiles, project JSON or snapshot JSON.",
                "Choose file",
                "file",
            ),
            2,
        )

        row.addWidget(
            self._card(
                "Import AF3 results",
                "Use this when you already have an AlphaFold 3 output folder and want to review metrics later in the workflow.",
                "Choose folder",
                "af3_results",
            ),
            1,
        )

        layout.addLayout(row, 1)

        hint = QLabel(
            "The guided workflow will update available steps after the input is recognized."
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
