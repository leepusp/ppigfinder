#!/usr/bin/env python3
"""
GUI helpers for AlphaFold Server JSON export.
"""

from __future__ import annotations

from itertools import combinations
from pathlib import Path

try:
    from PyQt6.QtWidgets import QFileDialog, QMessageBox
except Exception:
    from PyQt5.QtWidgets import QFileDialog, QMessageBox

from ppigfinder.services.alphafold_service import AlphaFoldService


def _selected_orf_indexes_from_table(table) -> list[int]:
    """
    Try to recover selected ORF indexes from a QTableWidget-like object.

    This is intentionally permissive because the legacy table may store
    indexes in different columns or only use row numbers.
    """
    indexes: list[int] = []

    try:
        rows = sorted({item.row() for item in table.selectedItems()})
    except Exception:
        return indexes

    for row in rows:
        indexes.append(row)

    return indexes


def export_selected_orfs_as_server_json(window) -> bool:
    """
    Export selected ORFs as AlphaFold Server JSON all-vs-all pairs.
    """
    orfs = getattr(window, "orfs", []) or []

    if not orfs:
        QMessageBox.warning(
            window,
            "Export AF3 Server JSON",
            "No ORFs available. Run ORF analysis first.",
        )
        return False

    table = None

    for attr in [
        "_orf_table",
        "orf_table",
        "_orfs_table",
        "orfs_table",
        "table_orfs",
    ]:
        candidate = getattr(window, attr, None)
        if candidate is not None:
            table = candidate
            break

    if table is None:
        QMessageBox.warning(
            window,
            "Export AF3 Server JSON",
            "Could not find the ORF table in the current interface.",
        )
        return False

    selected = _selected_orf_indexes_from_table(table)

    if len(selected) < 2:
        QMessageBox.information(
            window,
            "Export AF3 Server JSON",
            "Select at least two ORFs in the ORF table to generate pairwise AF3 jobs.",
        )
        return False

    pairs = list(combinations(selected, 2))

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
            model_seeds=[],
            use_structure_template=None,
        )
    except Exception as exc:
        QMessageBox.critical(
            window,
            "Export AF3 Server JSON",
            f"Could not export AlphaFold Server JSON:\n{exc}",
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
