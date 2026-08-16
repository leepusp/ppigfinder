#!/usr/bin/env python3
"""
Text fallback utilities for portable UI rendering.

Some Linux/X11/VNC environments do not render emoji reliably. This module
normalizes labels without relying on emoji fonts and removes textual icon
fallback prefixes such as [Open], [Dir], [DNA], [Gear], etc.
"""

from __future__ import annotations

import re

try:
    from PyQt6.QtGui import QAction
    from PyQt6.QtWidgets import (
        QPushButton,
        QToolButton,
        QCheckBox,
        QGroupBox,
        QLabel,
        QTabWidget,
        QComboBox,
    )
except Exception:
    from PyQt5.QtWidgets import (
        QAction,
        QPushButton,
        QToolButton,
        QCheckBox,
        QGroupBox,
        QLabel,
        QTabWidget,
        QComboBox,
    )


_EMOJI_RE = re.compile(
    "["
    "\U0001F300-\U0001FAFF"
    "\U00002700-\U000027BF"
    "\U00002600-\U000026FF"
    "]+",
    flags=re.UNICODE,
)

_PREFIX_RE = re.compile(r"^\s*\[[^\]]{1,16}\]\s*")


_REPLACEMENTS = {
    "📂": "",
    "💾": "",
    "🧬": "",
    "🔬": "",
    "🧪": "",
    "🖥": "",
    "⚙": "",
    "❓": "",
    "▶": "",
    "■": "",
    "✓": "",
    "✗": "",
    "★": "*",
    "→": "->",
    "↔": "<->",
}


def strip_icon_prefix(text: str) -> str:
    """
    Remove textual icon prefixes like [Open], [Dir], [DNA], [Chart].
    """
    return _PREFIX_RE.sub("", text or "")


def clean_ui_text(text: str) -> str:
    """
    Remove or normalize symbols that commonly break in remote Qt sessions.
    """
    if not text:
        return text

    text = strip_icon_prefix(text)

    for old, new in _REPLACEMENTS.items():
        text = text.replace(old, new)

    text = _EMOJI_RE.sub("", text)
    text = re.sub(r"\s+", " ", text).strip()

    return text


def _clean_text_object(obj) -> None:
    if not hasattr(obj, "text") or not hasattr(obj, "setText"):
        return

    try:
        original = obj.text()
        cleaned = clean_ui_text(original)
        if cleaned and cleaned != original:
            obj.setText(cleaned)
    except Exception:
        pass


def apply_text_fallback_to_window(window) -> None:
    """
    Apply text cleanup to menus, actions, buttons, labels, groups, tabs and
    combo boxes already created in a window.
    """
    for cls in (QAction, QPushButton, QToolButton, QCheckBox, QGroupBox, QLabel):
        try:
            for obj in window.findChildren(cls):
                _clean_text_object(obj)
        except Exception:
            pass

    try:
        for tab_widget in window.findChildren(QTabWidget):
            for index in range(tab_widget.count()):
                original = tab_widget.tabText(index)
                cleaned = clean_ui_text(original)
                if cleaned and cleaned != original:
                    tab_widget.setTabText(index, cleaned)
        # Some main windows may themselves be a QTabWidget subclass.
        if isinstance(window, QTabWidget):
            for index in range(window.count()):
                original = window.tabText(index)
                cleaned = clean_ui_text(original)
                if cleaned and cleaned != original:
                    window.setTabText(index, cleaned)
    except Exception:
        pass

    try:
        for combo in window.findChildren(QComboBox):
            for index in range(combo.count()):
                original = combo.itemText(index)
                cleaned = clean_ui_text(original)
                if cleaned and cleaned != original:
                    combo.setItemText(index, cleaned)
    except Exception:
        pass
