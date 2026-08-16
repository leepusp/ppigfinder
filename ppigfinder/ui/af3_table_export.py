#!/usr/bin/env python3
"""
GUI bridge for exporting AF3 results tables.
"""

from __future__ import annotations

try:
    from PyQt6.QtGui import QAction
    from PyQt6.QtWidgets import QFileDialog, QMessageBox
except Exception:
    from PyQt5.QtWidgets import QAction, QFileDialog, QMessageBox

from ppigfinder.io.af3_table_export import write_af3_results_table


def export_af3_results_table(window) -> bool:
    """
    Export current AF3 results to TSV.
    """
    results = getattr(window, "af3_results", []) or []

    if not results:
        QMessageBox.information(
            window,
            "Export AF3 results",
            "No AF3 results available. Import or analyze AF3 results first.",
        )
        return False

    f, _ = QFileDialog.getSaveFileName(
        window,
        "Export AF3 Results Table",
        "af3_results.tsv",
        "TSV (*.tsv);;CSV (*.csv)",
    )

    if not f:
        return False

    delimiter = "," if f.lower().endswith(".csv") else "\t"

    try:
        write_af3_results_table(f, results, delimiter=delimiter)
    except Exception as exc:
        QMessageBox.critical(
            window,
            "Export AF3 results",
            f"Could not export AF3 results table:\n{exc}",
        )
        return False

    try:
        window._status.showMessage(f"Exported {len(results)} AF3 result(s): {f}")
    except Exception:
        pass

    return True


def install_af3_results_export_action(window) -> None:
    """
    Add AF3 results table export to the File menu without modifying legacy_v20.py.
    """
    if getattr(window, "_ppig_af3_results_export_action_installed", False):
        return

    try:
        menu_bar = window.menuBar()
    except Exception:
        return

    file_menu = None

    for action in menu_bar.actions():
        text = action.text().replace("&", "").strip()
        if "File" in text:
            file_menu = action.menu()
            break

    if file_menu is None:
        file_menu = menu_bar.addMenu("File")

    action = QAction("Export AF3 Results Table", window)
    action.triggered.connect(lambda checked=False: export_af3_results_table(window))
    file_menu.addAction(action)

    window._ppig_af3_results_export_action_installed = True
