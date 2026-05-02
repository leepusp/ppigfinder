#!/usr/bin/env python3
"""
AlphaFold 3 candidate pair builder dialog for the guided UI shell.
"""

from __future__ import annotations

try:
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import (
        QDialog,
        QVBoxLayout,
        QHBoxLayout,
        QLabel,
        QPushButton,
        QTableWidget,
        QTableWidgetItem,
        QHeaderView,
        QFileDialog,
        QMessageBox,
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
        QTableWidget,
        QTableWidgetItem,
        QHeaderView,
        QFileDialog,
        QMessageBox,
    )
    QT6 = False

from ppigfinder.ui_shell.branding import apply_ppigfinder_branding
from ppigfinder.alphafold.guided_pair_builder import (
    make_adjacent_orf_pairs,
    write_af3_json,
)


def _window_flags():
    flags = Qt.WindowType.Window if QT6 else Qt.Window
    flags |= Qt.WindowType.WindowMinimizeButtonHint if QT6 else Qt.WindowMinimizeButtonHint
    flags |= Qt.WindowType.WindowMaximizeButtonHint if QT6 else Qt.WindowMaximizeButtonHint
    flags |= Qt.WindowType.WindowCloseButtonHint if QT6 else Qt.WindowCloseButtonHint
    return flags


class AF3PairBuilderDialog(QDialog):
    """
    Build and export candidate ORF pairs for AlphaFold 3.
    """

    def __init__(self, orfs, parent=None):
        super().__init__(parent)

        self.orfs = list(orfs or [])
        self.pairs = make_adjacent_orf_pairs(self.orfs, max_pairs=50)
        self.exported_json_path = ""

        self.setWindowTitle("AlphaFold / PPI Candidate Pair Builder")
        self.setWindowFlags(_window_flags())
        self.resize(1180, 760)
        self.setMinimumSize(900, 580)
        self.setSizeGripEnabled(True)
        apply_ppigfinder_branding(self)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(18, 18, 18, 18)
        layout.setSpacing(12)

        title = QLabel(f"Candidate AlphaFold 3 pairs: {len(self.pairs)}")
        title.setObjectName("SectionTitle")
        layout.addWidget(title)

        description = QLabel(
            "This first guided strategy builds adjacent ORF pairs from the predicted ORF list. "
            "Later strategies can include query-vs-neighbours, HMM-positive ORFs, selected ORFs "
            "all-vs-all, homodimers and custom stoichiometry."
        )
        description.setWordWrap(True)
        layout.addWidget(description)

        self.table = QTableWidget(0, 7, self)
        self.table.setHorizontalHeaderLabels(
            [
                "Job name",
                "ORF A",
                "ORF B",
                "A length",
                "B length",
                "Distance nt",
                "Strategy",
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
        layout.addWidget(self.table, 1)

        buttons = QHBoxLayout()

        export_button = QPushButton("Export AF3 Server JSON")
        export_button.clicked.connect(self._export_json)
        buttons.addWidget(export_button)

        buttons.addStretch(1)
        layout.addLayout(buttons)

    def _populate(self) -> None:
        self.table.setRowCount(len(self.pairs))

        for row, pair in enumerate(self.pairs):
            values = [
                pair.name,
                pair.orf_a_id,
                pair.orf_b_id,
                str(len(pair.sequence_a)),
                str(len(pair.sequence_b)),
                str(pair.distance_nt),
                "Adjacent ORFs",
            ]

            for col, value in enumerate(values):
                self.table.setItem(row, col, QTableWidgetItem(value))

    def _export_json(self) -> None:
        if not self.pairs:
            QMessageBox.warning(
                self,
                "AF3 JSON export",
                "No candidate pairs are available. Predict ORFs first.",
            )
            return

        path, _ = QFileDialog.getSaveFileName(
            self,
            "Export AF3 Server JSON",
            "ppigfinder_guided_af3_pairs.json",
            "JSON (*.json);;All files (*)",
        )

        if not path:
            return

        write_af3_json(path, self.pairs, seeds=[1])
        self.exported_json_path = path

        QMessageBox.information(
            self,
            "AF3 JSON export",
            "AlphaFold Server JSON exported:\n\n" + path,
        )


def open_af3_pair_builder_dialog(orfs, parent=None) -> dict:
    dialog = AF3PairBuilderDialog(orfs, parent=parent)
    dialog.showMaximized()
    dialog.exec() if hasattr(dialog, "exec") else dialog.exec_()

    return {
        "pair_count": len(dialog.pairs),
        "json_path": dialog.exported_json_path,
    }
