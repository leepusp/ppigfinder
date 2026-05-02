#!/usr/bin/env python3
"""
Text fallback utilities for portable UI rendering.

Some Linux/X11/VNC environments do not render emoji reliably. This module
normalizes labels without relying on emoji fonts.
"""

from __future__ import annotations

import re


_EMOJI_RE = re.compile(
    "["
    "\U0001F300-\U0001FAFF"
    "\U00002700-\U000027BF"
    "\U00002600-\U000026FF"
    "]+",
    flags=re.UNICODE,
)


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


def clean_ui_text(text: str) -> str:
    """
    Remove or normalize symbols that commonly break in remote Qt sessions.
    """
    if not text:
        return text

    for old, new in _REPLACEMENTS.items():
        text = text.replace(old, new)

    text = _EMOJI_RE.sub("", text)
    text = re.sub(r"\s+", " ", text).strip()

    return text


def apply_text_fallback_to_window(window) -> None:
    """
    Apply text cleanup to common UI elements already created in a window.
    """
    try:
        objects = window.findChildren(object)
    except Exception:
        return

    allowed_class_names = {
        "QAction",
        "QPushButton",
        "QToolButton",
        "QCheckBox",
        "QGroupBox",
    }

    for obj in objects:
        cls_name = obj.__class__.__name__

        if cls_name not in allowed_class_names:
            continue

        if not hasattr(obj, "text") or not hasattr(obj, "setText"):
            continue

        try:
            original = obj.text()
            cleaned = clean_ui_text(original)
            if cleaned and cleaned != original:
                obj.setText(cleaned)
        except Exception:
            pass
