#!/usr/bin/env python3
"""
GUI helpers for AlphaFold Server JSON export.
"""

from __future__ import annotations

from itertools import combinations
from pathlib import Path

try:
    from PyQt6.QtWidgets import (
        QFileDialog,
        QMessageBox,
        QDialog,
        QVBoxLayout,
        QHBoxLayout,
        QLabel,
        QComboBox,
        QSpinBox,
        QCheckBox,
        QLineEdit,
        QDialogButtonBox,
    )
except Exception:
    from PyQt5.QtWidgets import (
        QFileDialog,
        QMessageBox,
        QDialog,
        QVBoxLayout,
        QHBoxLayout,
        QLabel,
        QComboBox,
        QSpinBox,
        QCheckBox,
        QLineEdit,
        QDialogButtonBox,
    )

from ppigfinder.services.alphafold_service import AlphaFoldService


class AF3ServerJsonExportDialog(QDialog):
    """
    Small configuration dialog for AlphaFold Server JSON export.
    """

    MODE_SELECTED_ALL_VS_ALL = "Selected ORFs all-vs-all"
    MODE_ALL_ORFS_NEIGHBORS = "All ORFs neighbor pairs"
    MODE_SELECTED_VS_NEIGHBORS = "Selected ORFs vs neighbors"

    def __init__(self, parent=None):
        super().__init__(parent)

        self.setWindowTitle("Export AlphaFold Server JSON")
        self.setMinimumWidth(460)

        layout = QVBoxLayout(self)

        layout.addWidget(QLabel("Export mode:"))
        self.mode_combo = QComboBox()
        self.mode_combo.addItems([
            self.MODE_SELECTED_ALL_VS_ALL,
            self.MODE_ALL_ORFS_NEIGHBORS,
            self.MODE_SELECTED_VS_NEIGHBORS,
        ])
        layout.addWidget(self.mode_combo)

        row = QHBoxLayout()
        row.addWidget(QLabel("Neighbor window:"))
        self.window_spin = QSpinBox()
        self.window_spin.setRange(1, 20)
        self.window_spin.setValue(2)
        row.addWidget(self.window_spin)
        layout.addLayout(row)

        self.templates_check = QCheckBox("Use structure templates when supported")
        self.templates_check.setChecked(False)
        layout.addWidget(self.templates_check)

        layout.addWidget(QLabel("Model seeds, comma-separated. Empty = automatic:"))
        self.seeds_edit = QLineEdit()
        self.seeds_edit.setPlaceholderText("e.g. 1,2,3")
        layout.addWidget(self.seeds_edit)

        self.summary_label = QLabel(
            "The exported JSON follows the AlphaFold Server/web format."
        )
        self.summary_label.setWordWrap(True)
        layout.addWidget(self.summary_label)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok
            | QDialogButtonBox.StandardButton.Cancel
            if hasattr(QDialogButtonBox, "StandardButton")
            else QDialogButtonBox.Ok | QDialogButtonBox.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def model_seeds(self) -> list[str]:
        raw = self.seeds_edit.text().strip()

        if not raw:
            return []

        return [item.strip() for item in raw.split(",") if item.strip()]

    def use_structure_template(self) -> bool | None:
        return True if self.templates_check.isChecked() else None

    def mode(self) -> str:
        return self.mode_combo.currentText()

    def neighbor_window(self) -> int:
        return self.window_spin.value()


def _selected_orf_indexes_from_table(table) -> list[int]:
    """
    Try to recover selected ORF indexes from a QTableWidget-like object.
    """
    indexes: list[int] = []

    try:
        rows = sorted({item.row() for item in table.selectedItems()})
    except Exception:
        return indexes

    for row in rows:
        indexes.append(row)

    return indexes


def _find_orf_table(window):
    for attr in [
        "_orf_table",
        "orf_table",
        "_orfs_table",
        "orfs_table",
        "table_orfs",
        "orfTable",
    ]:
        candidate = getattr(window, attr, None)
        if candidate is not None:
            return candidate

    return None


def _valid_orf_indexes(indexes: list[int], n_orfs: int) -> list[int]:
    return sorted({idx for idx in indexes if 0 <= idx < n_orfs})


def _all_vs_all_pairs(indexes: list[int]) -> list[tuple[int, int]]:
    return list(combinations(indexes, 2))


def _neighbor_pairs(indexes: list[int], window: int) -> list[tuple[int, int]]:
    wanted = set(indexes)
    pairs = []

    for i in indexes:
        for j in range(i + 1, i + window + 1):
            if j in wanted:
                pairs.append((i, j))

    return pairs


def _selected_vs_neighbors(selected: list[int], n_orfs: int, window: int) -> list[tuple[int, int]]:
    pairs = set()

    for i in selected:
        start = max(0, i - window)
        end = min(n_orfs, i + window + 1)

        for j in range(start, end):
            if j == i:
                continue

            pairs.add(tuple(sorted((i, j))))

    return sorted(pairs)


def _build_pairs(window, dialog: AF3ServerJsonExportDialog) -> list[tuple[int, int]]:
    orfs = getattr(window, "orfs", []) or []
    n_orfs = len(orfs)

    table = _find_orf_table(window)
    selected = []

    if table is not None:
        selected = _valid_orf_indexes(_selected_orf_indexes_from_table(table), n_orfs)

    mode = dialog.mode()
    neighbor_window = dialog.neighbor_window()

    if mode == dialog.MODE_SELECTED_ALL_VS_ALL:
        if len(selected) < 2:
            raise ValueError("Select at least two ORFs in the ORF table.")
        return _all_vs_all_pairs(selected)

    if mode == dialog.MODE_ALL_ORFS_NEIGHBORS:
        indexes = list(range(n_orfs))
        return _neighbor_pairs(indexes, neighbor_window)

    if mode == dialog.MODE_SELECTED_VS_NEIGHBORS:
        if not selected:
            raise ValueError("Select at least one ORF in the ORF table.")
        return _selected_vs_neighbors(selected, n_orfs, neighbor_window)

    raise ValueError(f"Unsupported export mode: {mode}")


def export_selected_orfs_as_server_json(window) -> bool:
    """
    Export ORF pairs as AlphaFold Server JSON.
    """
    orfs = getattr(window, "orfs", []) or []

    if not orfs:
        QMessageBox.warning(
            window,
            "Export AF3 Server JSON",
            "No ORFs available. Run ORF analysis first.",
        )
        return False

    dialog = AF3ServerJsonExportDialog(window)

    if hasattr(dialog, "exec"):
        accepted = dialog.exec() == (
            QDialog.DialogCode.Accepted
            if hasattr(QDialog, "DialogCode")
            else QDialog.Accepted
        )
    else:
        accepted = dialog.exec_() == QDialog.Accepted

    if not accepted:
        return False

    try:
        pairs = _build_pairs(window, dialog)
    except Exception as exc:
        QMessageBox.warning(
            window,
            "Export AF3 Server JSON",
            str(exc),
        )
        return False

    if not pairs:
        QMessageBox.information(
            window,
            "Export AF3 Server JSON",
            "No ORF pairs were generated for this mode.",
        )
        return False

    f, _ = QFileDialog.getSaveFileName(
        window,
        "Export AlphaFold Server JSON",
        "af3_server_jobs.json",
        "JSON (*.json)",
    )

    if not f:
        return False

    try:
        AlphaFoldService().export_server_json_from_legacy_orf_pairs(
            f,
            orfs=orfs,
            pairs=pairs,
            model_seeds=dialog.model_seeds(),
            use_structure_template=dialog.use_structure_template(),
        )
    except Exception as exc:
        QMessageBox.critical(
            window,
            "Export AF3 Server JSON",
            "Could not export AlphaFold Server JSON.\n\n"
            "A common cause is an ORF protein sequence containing unsupported "
            "residues such as X, B, Z, U or internal stop codons.\n\n"
            f"Details:\n{exc}",
        )
        return False

    try:
        window._status.showMessage(
            f"Exported {len(pairs)} AF3 Server jobs to {Path(f).name}"
        )
    except Exception:
        pass

    QMessageBox.information(
        window,
        "Export AF3 Server JSON",
        f"Exported {len(pairs)} pairwise jobs:\n{f}",
    )

    return True
