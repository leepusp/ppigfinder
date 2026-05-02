#!/usr/bin/env python3
"""
GUI bridge for ppigFinder project snapshots.
"""

from __future__ import annotations

try:
    from PyQt6.QtWidgets import QFileDialog, QMessageBox
except Exception:
    from PyQt5.QtWidgets import QFileDialog, QMessageBox

from ppigfinder.services.project_service import ProjectService


def export_project_snapshot_from_window(window) -> bool:
    """
    Export current GUI state to the new versioned project snapshot format.
    """
    f, _ = QFileDialog.getSaveFileName(
        window,
        "Export ppigFinder Project Snapshot",
        "",
        "ppigFinder snapshot (*.ppigfinder.json);;JSON (*.json)",
    )

    if not f:
        return False

    try:
        project = ProjectService().build_snapshot_from_legacy_window(window)
        ProjectService().save(f, project)
    except Exception as exc:
        QMessageBox.critical(
            window,
            "Export project snapshot",
            f"Could not export project snapshot:\n{exc}",
        )
        return False

    try:
        window._status.showMessage(f"Project snapshot exported: {f}")
    except Exception:
        pass

    return True


def import_project_snapshot_into_window(window) -> bool:
    """
    Import a new versioned project snapshot into the current GUI.
    """
    f, _ = QFileDialog.getOpenFileName(
        window,
        "Import ppigFinder Project Snapshot",
        "",
        "ppigFinder snapshot (*.ppigfinder.json *.json);;All (*)",
    )

    if not f:
        return False

    try:
        project = ProjectService().load(f)
        ProjectService().apply_snapshot_to_legacy_window(window, project)
    except Exception as exc:
        QMessageBox.critical(
            window,
            "Import project snapshot",
            f"Could not import project snapshot:\n{exc}",
        )
        return False

    try:
        window._status.showMessage(f"Project snapshot imported: {f}")
    except Exception:
        pass

    return True
