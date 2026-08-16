#!/usr/bin/env python3
"""
Toolbar polish helpers for ppigFinder.
"""

from __future__ import annotations

try:
    from PyQt6.QtCore import QSize, Qt
    from PyQt6.QtWidgets import QToolBar
    QT6 = True
except Exception:
    from PyQt5.QtCore import QSize, Qt
    from PyQt5.QtWidgets import QToolBar
    QT6 = False


def _text_beside_icon_style():
    return Qt.ToolButtonStyle.ToolButtonTextBesideIcon if QT6 else Qt.ToolButtonTextBesideIcon


def polish_toolbars(window, icon_size: int = 20) -> None:
    """
    Make toolbars more readable and stable across environments.
    """
    try:
        toolbars = window.findChildren(QToolBar)
    except Exception:
        return

    for index, toolbar in enumerate(toolbars, start=1):
        try:
            if not toolbar.objectName():
                toolbar.setObjectName(f"toolbar_{index}")

            toolbar.setMovable(False)
            toolbar.setIconSize(QSize(icon_size, icon_size))
            toolbar.setToolButtonStyle(_text_beside_icon_style())
        except Exception:
            pass
