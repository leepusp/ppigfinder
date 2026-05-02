#!/usr/bin/env python3
"""
Simple ORF map preview for the experimental guided UI shell.

This is intentionally Qt-only for now. Later it can be replaced or extended
with HTML/D3-based visualization.
"""

from __future__ import annotations

try:
    from PyQt6.QtCore import Qt, QRectF
    from PyQt6.QtGui import QPainter, QColor, QPen, QBrush, QPolygonF
    from PyQt6.QtCore import QPointF
    from PyQt6.QtWidgets import QWidget
    QT6 = True
except Exception:
    from PyQt5.QtCore import Qt, QRectF, QPointF
    from PyQt5.QtGui import QPainter, QColor, QPen, QBrush, QPolygonF
    from PyQt5.QtWidgets import QWidget
    QT6 = False


def _antialiasing():
    return QPainter.RenderHint.Antialiasing if QT6 else QPainter.Antialiasing


class ORFMapPreviewWidget(QWidget):
    """
    Minimal genome/ORF map preview.

    Input format:
        [
            {"start": 1, "end": 900, "strand": "+", "id": "..."},
            ...
        ]
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.orfs: list[dict] = []
        self.genome_length: int = 0
        self.setMinimumHeight(180)

    def set_data(self, orfs: list[dict], genome_length: int) -> None:
        self.orfs = list(orfs or [])
        self.genome_length = int(genome_length or 0)
        self.update()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(_antialiasing())

        width = self.width()
        height = self.height()

        bg = QColor("#f8fbfc")
        painter.fillRect(self.rect(), bg)

        margin_x = 28
        axis_y = height // 2
        usable_w = max(1, width - 2 * margin_x)

        # Axis
        painter.setPen(QPen(QColor("#607d8b"), 2))
        painter.drawLine(margin_x, axis_y, width - margin_x, axis_y)

        if not self.orfs or self.genome_length <= 0:
            painter.setPen(QColor("#607d8b"))
            painter.drawText(
                self.rect(),
                Qt.AlignmentFlag.AlignCenter if QT6 else Qt.AlignCenter,
                "No ORF map available yet. Run Predict ORFs.",
            )
            painter.end()
            return

        # Limit display to avoid visual clutter.
        visible_orfs = self.orfs[:120]

        plus_color = QColor("#43a047")
        minus_color = QColor("#1e88e5")

        for orf in visible_orfs:
            start = int(orf.get("start", 0) or 0)
            end = int(orf.get("end", 0) or 0)
            strand = str(orf.get("strand", "+"))

            if end < start:
                start, end = end, start

            x1 = margin_x + (start / self.genome_length) * usable_w
            x2 = margin_x + (end / self.genome_length) * usable_w

            if x2 - x1 < 4:
                x2 = x1 + 4

            y = axis_y - 34 if strand == "+" else axis_y + 18
            h = 16
            arrow = min(12, max(5, (x2 - x1) * 0.35))

            if strand == "+":
                points = [
                    QPointF(x1, y),
                    QPointF(x2 - arrow, y),
                    QPointF(x2, y + h / 2),
                    QPointF(x2 - arrow, y + h),
                    QPointF(x1, y + h),
                ]
                color = plus_color
            else:
                points = [
                    QPointF(x2, y),
                    QPointF(x1 + arrow, y),
                    QPointF(x1, y + h / 2),
                    QPointF(x1 + arrow, y + h),
                    QPointF(x2, y + h),
                ]
                color = minus_color

            painter.setPen(QPen(QColor("#263238"), 0.5))
            painter.setBrush(QBrush(color))
            painter.drawPolygon(QPolygonF(points))

        painter.setPen(QColor("#455a64"))
        painter.drawText(28, 22, f"Genome length: {self.genome_length:,} nt")
        painter.drawText(28, height - 16, f"Showing {len(visible_orfs)} of {len(self.orfs)} ORFs")

        painter.setPen(QColor("#43a047"))
        painter.drawText(width - 170, 22, "+ strand")
        painter.setPen(QColor("#1e88e5"))
        painter.drawText(width - 90, 22, "- strand")

        painter.end()
