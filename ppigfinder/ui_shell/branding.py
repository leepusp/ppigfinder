#!/usr/bin/env python3
"""
Branding helpers for the future ppigFinder UI shell.

The current icon is generated programmatically so the interface does not
depend on external image assets during development.
"""

from __future__ import annotations

try:
    from PyQt6.QtCore import Qt, QRectF
    from PyQt6.QtGui import QIcon, QPixmap, QPainter, QColor, QFont, QBrush
    QT6 = True
except Exception:
    from PyQt5.QtCore import Qt, QRectF
    from PyQt5.QtGui import QIcon, QPixmap, QPainter, QColor, QFont, QBrush
    QT6 = False


def _no_pen():
    return Qt.PenStyle.NoPen if QT6 else Qt.NoPen


def _align_center():
    return Qt.AlignmentFlag.AlignCenter if QT6 else Qt.AlignCenter


def _antialiasing():
    return QPainter.RenderHint.Antialiasing if QT6 else QPainter.Antialiasing


def create_ppigfinder_pixmap(size: int = 96) -> QPixmap:
    """
    Create a simple ppigFinder icon inspired by the current P logo.
    """
    pixmap = QPixmap(size, size)
    pixmap.fill(QColor(0, 0, 0, 0))

    painter = QPainter(pixmap)
    painter.setRenderHint(_antialiasing())

    # Background
    painter.setPen(_no_pen())
    painter.setBrush(QBrush(QColor("#1b3a4b")))
    radius = max(8, int(size * 0.16))
    painter.drawRoundedRect(QRectF(0, 0, size, size), radius, radius)

    # Subtle left accent
    painter.setBrush(QBrush(QColor("#66bb6a")))
    painter.drawRoundedRect(
        QRectF(size * 0.08, size * 0.08, size * 0.12, size * 0.84),
        radius * 0.55,
        radius * 0.55,
    )

    # Letter
    font = QFont("Arial")
    font.setBold(True)
    font.setPointSize(max(22, int(size * 0.52)))
    painter.setFont(font)
    painter.setPen(QColor("#ffffff"))
    painter.drawText(QRectF(0, 0, size, size), _align_center(), "P")

    # Small interaction dot
    painter.setPen(_no_pen())
    painter.setBrush(QBrush(QColor("#90caf9")))
    dot = size * 0.12
    painter.drawEllipse(QRectF(size * 0.70, size * 0.20, dot, dot))

    painter.end()
    return pixmap


def create_ppigfinder_icon(size: int = 96) -> QIcon:
    return QIcon(create_ppigfinder_pixmap(size))


def apply_ppigfinder_branding(window) -> None:
    """
    Apply generated ppigFinder icon to any Qt window.
    """
    try:
        window.setWindowIcon(create_ppigfinder_icon())
    except Exception:
        pass
