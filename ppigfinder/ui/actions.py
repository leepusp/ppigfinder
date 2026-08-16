#!/usr/bin/env python3
"""
Standard QAction factory for ppigFinder.

Centralizing QAction creation makes icons, shortcuts and text fallback
consistent across menus and toolbars.
"""

from __future__ import annotations

try:
    from PyQt6.QtGui import QAction
except Exception:
    from PyQt5.QtWidgets import QAction

from ppigfinder.ui.icon_provider import make_icon
from ppigfinder.ui.text_fallback import clean_ui_text


def create_action(
    parent,
    text: str,
    callback=None,
    icon_name: str | None = None,
    shortcut: str | None = None,
    tooltip: str | None = None,
    checkable: bool = False,
) -> QAction:
    """
    Create a portable QAction with generated icon and cleaned text.
    """
    action = QAction(clean_ui_text(text), parent)

    if icon_name:
        action.setIcon(make_icon(icon_name))

    if shortcut:
        action.setShortcut(shortcut)

    if tooltip:
        action.setToolTip(clean_ui_text(tooltip))
        action.setStatusTip(clean_ui_text(tooltip))

    if checkable:
        action.setCheckable(True)

    if callback:
        action.triggered.connect(callback)

    return action
