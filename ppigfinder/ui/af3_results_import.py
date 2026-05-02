#!/usr/bin/env python3
"""
GUI bridge for importing AlphaFold/AF3 result folders.
"""

from __future__ import annotations

from pathlib import Path

try:
    from PyQt6.QtGui import QAction
    from PyQt6.QtWidgets import QFileDialog, QMessageBox
except Exception:
    from PyQt5.QtWidgets import QAction, QFileDialog, QMessageBox

from ppigfinder.services.alphafold_results_service import AlphaFoldResultsService
from ppigfinder.ui.workers import run_in_thread


def import_af3_results_folder(window) -> bool:
    """
    Import AF3 result folders into the current GUI state using the modular parser.
    """
    folder = QFileDialog.getExistingDirectory(
        window,
        "Select AlphaFold/AF3 results folder",
        "",
    )

    if not folder:
        return False

    try:
        window._status.showMessage(f"Parsing AF3 results: {folder}")
    except Exception:
        pass

    service = AlphaFoldResultsService()

    def work():
        cache_dir = Path(folder) / ".ppigfinder_cache"
        return service.parse_root_parallel(
            folder,
            workers=None,
            cache_dir=cache_dir,
        )

    def done(results):
        existing = getattr(window, "af3_results", []) or []
        window.af3_results = existing + list(results)

        try:
            window._status.showMessage(
                f"Imported {len(results)} AF3 result(s) from {folder}"
            )
        except Exception:
            pass

        for method_name in [
            "_populate_af3_results_table",
            "populate_af3_results_table",
            "_update_af3_results_table",
            "update_af3_results_table",
            "_refresh_af3_results",
            "refresh_af3_results",
        ]:
            method = getattr(window, method_name, None)
            if callable(method):
                try:
                    method()
                    break
                except Exception:
                    pass

        QMessageBox.information(
            window,
            "Import AF3 results",
            f"Imported {len(results)} AF3 result(s).",
        )

    def failed(traceback_text):
        try:
            window._status.showMessage("AF3 result import failed")
        except Exception:
            pass

        QMessageBox.critical(
            window,
            "Import AF3 results",
            f"Could not import AF3 results:\n\n{traceback_text}",
        )

    run_in_thread(
        window,
        work,
        on_finished=done,
        on_failed=failed,
    )

    return True


def install_af3_results_import_action(window) -> None:
    """
    Add 'Import AF3 Results Folder' to the File menu without modifying legacy_v20.py.
    """
    if getattr(window, "_ppig_af3_results_import_action_installed", False):
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

    file_menu.addSeparator()

    action = QAction("Import AF3 Results Folder", window)
    action.triggered.connect(lambda checked=False: import_af3_results_folder(window))
    file_menu.addAction(action)

    window._ppig_af3_results_import_action_installed = True
