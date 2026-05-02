#!/usr/bin/env python3
"""
Recent file support for ppigFinder.

This module adds a portable Recent Files menu using QSettings.
"""

from __future__ import annotations

from pathlib import Path


try:
    from PyQt6.QtGui import QAction
    from PyQt6.QtCore import QSettings
    from PyQt6.QtWidgets import QMessageBox
except Exception:
    from PyQt5.QtWidgets import QAction, QMessageBox
    from PyQt5.QtCore import QSettings


_SETTINGS_ORG = "LEEPBioinfo"
_SETTINGS_APP = "ppigFinder"
_SETTINGS_KEY = "recent_files"
_MAX_RECENT = 10


def _settings() -> QSettings:
    return QSettings(_SETTINGS_ORG, _SETTINGS_APP)


def get_recent_files() -> list[str]:
    """
    Return existing recent files.
    """
    settings = _settings()
    value = settings.value(_SETTINGS_KEY, [])

    if isinstance(value, str):
        value = [value]

    files = []
    for item in value or []:
        path = str(item)
        if path and Path(path).exists() and path not in files:
            files.append(path)

    return files[:_MAX_RECENT]


def add_recent_file(path: str | Path) -> None:
    """
    Add a file to the recent files list.
    """
    path = str(Path(path).expanduser().resolve())
    files = get_recent_files()

    if path in files:
        files.remove(path)

    files.insert(0, path)
    files = files[:_MAX_RECENT]

    _settings().setValue(_SETTINGS_KEY, files)


def clear_recent_files() -> None:
    """
    Clear recent file list.
    """
    _settings().setValue(_SETTINGS_KEY, [])


def install_recent_files_menu(window, open_callback) -> None:
    """
    Install a Recent Files menu into the main window.

    open_callback must accept one path argument.
    """
    try:
        menu_bar = window.menuBar()
    except Exception:
        return

    recent_menu = menu_bar.addMenu("Recent Files")

    def rebuild_menu():
        recent_menu.clear()
        files = get_recent_files()

        if not files:
            empty = QAction("No recent files", window)
            empty.setEnabled(False)
            recent_menu.addAction(empty)
            return

        for path in files:
            label = Path(path).name
            action = QAction(label, window)
            action.setToolTip(path)
            action.triggered.connect(lambda checked=False, p=path: open_callback(p))
            recent_menu.addAction(action)

        recent_menu.addSeparator()
        clear_action = QAction("Clear recent files", window)

        def _clear():
            clear_recent_files()
            rebuild_menu()

        clear_action.triggered.connect(_clear)
        recent_menu.addAction(clear_action)

    recent_menu.aboutToShow.connect(rebuild_menu)
    rebuild_menu()
