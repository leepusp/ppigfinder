#!/usr/bin/env python3
"""
Backend status dialog and GUI installer.
"""

from __future__ import annotations

try:
    from PyQt6.QtGui import QAction
    from PyQt6.QtWidgets import QDialog, QVBoxLayout, QPushButton, QHBoxLayout
except Exception:
    from PyQt5.QtWidgets import QAction, QDialog, QVBoxLayout, QPushButton, QHBoxLayout

from ppigfinder.infrastructure.backends import refresh_backends, BACKENDS
from ppigfinder.ui.widgets.backend_status_widget import BackendStatusWidget


def _backend_rows(detailed: bool = False) -> list[dict]:
    """
    Build rows from backend detection.
    """
    backends = refresh_backends(detailed=True) if detailed else BACKENDS

    rows = []

    for key, label in [
        ("blast+", "BLAST+"),
        ("hmmer3", "HMMER3"),
        ("pyrodigal", "Pyrodigal"),
        ("paramiko", "Paramiko/SSH"),
    ]:
        item = backends.get(key, {}) or {}

        rows.append(
            {
                "backend": label,
                "available": item.get("available", False),
                "version": item.get("version", ""),
                "details": item.get("path", ""),
            }
        )

    return rows


class BackendStatusDialog(QDialog):
    """
    Dialog showing backend availability.
    """

    def __init__(self, parent=None):
        super().__init__(parent)

        self.setWindowTitle("Backend Status")
        self.resize(720, 320)

        layout = QVBoxLayout(self)

        self.widget = BackendStatusWidget(self)
        self.widget.set_backends(_backend_rows(detailed=False))
        layout.addWidget(self.widget)

        button_row = QHBoxLayout()

        refresh_button = QPushButton("Refresh detailed status")
        refresh_button.clicked.connect(self.refresh_detailed)
        button_row.addWidget(refresh_button)

        close_button = QPushButton("Close")
        close_button.clicked.connect(self.accept)
        button_row.addWidget(close_button)

        layout.addLayout(button_row)

    def refresh_detailed(self) -> None:
        self.widget.set_backends(_backend_rows(detailed=True))


def show_backend_status_dialog(window) -> None:
    dialog = BackendStatusDialog(window)
    dialog.exec() if hasattr(dialog, "exec") else dialog.exec_()


def install_backend_status_action(window) -> None:
    """
    Add Backend Status action without modifying legacy_v20.py.
    """
    if getattr(window, "_ppig_backend_status_action_installed", False):
        return

    menu_bar = window.menuBar()

    system_menu = None
    for action in menu_bar.actions():
        text = action.text().replace("&", "").strip().lower()
        if text in {"system", "tools", "window"}:
            system_menu = action.menu()
            break

    if system_menu is None:
        system_menu = menu_bar.addMenu("System")

    action = QAction("Backend Status", window)
    action.triggered.connect(lambda checked=False: show_backend_status_dialog(window))
    system_menu.addAction(action)

    window._ppig_backend_status_action_installed = True
