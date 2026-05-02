#!/usr/bin/env python3
"""
Portable icon provider for ppigFinder.

The icons are generated with QPainter, avoiding dependency on emoji fonts,
external icon themes, or system-specific resources.
"""

from __future__ import annotations


try:
    from PyQt6.QtGui import QIcon, QPixmap, QPainter, QColor, QFont
    from PyQt6.QtCore import Qt
    QT6 = True
except Exception:
    from PyQt5.QtGui import QIcon, QPixmap, QPainter, QColor, QFont
    from PyQt5.QtCore import Qt
    QT6 = False


def _align_center():
    return Qt.AlignmentFlag.AlignCenter if QT6 else Qt.AlignCenter


def _transparent():
    return Qt.GlobalColor.transparent if QT6 else Qt.transparent


def _no_pen():
    return Qt.PenStyle.NoPen if QT6 else Qt.NoPen


_ICON_DEFINITIONS = {
    "app": ("P", "#263238", "#ffffff"),
    "open": ("OP", "#1565c0", "#ffffff"),
    "save": ("SV", "#2e7d32", "#ffffff"),
    "export": ("EX", "#00695c", "#ffffff"),
    "orf": ("ORF", "#6a1b9a", "#ffffff"),
    "blast": ("BL", "#ef6c00", "#ffffff"),
    "hmm": ("HM", "#ad1457", "#ffffff"),
    "af3": ("AF3", "#283593", "#ffffff"),
    "hpc": ("HPC", "#37474f", "#ffffff"),
    "run": ("▶", "#2e7d32", "#ffffff"),
    "stop": ("■", "#b71c1c", "#ffffff"),
    "settings": ("⚙", "#455a64", "#ffffff"),
    "help": ("?", "#0277bd", "#ffffff"),
}


def make_icon(name: str, size: int = 32) -> QIcon:
    """
    Create a portable QIcon by drawing a simple labeled tile.
    """
    label, bg, fg = _ICON_DEFINITIONS.get(
        name,
        (name[:3].upper() if name else "?", "#455a64", "#ffffff"),
    )

    pixmap = QPixmap(size, size)
    pixmap.fill(_transparent())

    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing if QT6 else QPainter.Antialiasing)
    painter.setPen(_no_pen())
    painter.setBrush(QColor(bg))
    painter.drawRoundedRect(1, 1, size - 2, size - 2, 7, 7)

    font = QFont()
    font.setBold(True)

    if len(label) <= 1:
        font.setPointSize(max(10, int(size * 0.55)))
    elif len(label) == 2:
        font.setPointSize(max(8, int(size * 0.34)))
    else:
        font.setPointSize(max(7, int(size * 0.25)))

    painter.setFont(font)
    painter.setPen(QColor(fg))
    painter.drawText(pixmap.rect(), _align_center(), label)
    painter.end()

    return QIcon(pixmap)


def set_window_icon(window) -> None:
    """
    Apply the ppigFinder application icon to a window.
    """
    try:
        window.setWindowIcon(make_icon("app", 48))
    except Exception:
        pass
