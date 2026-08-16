#!/usr/bin/env python3
"""
Backend status table widget.

Shows availability of external and optional ppigFinder backends without
depending on emoji or system icon fonts.
"""

from __future__ import annotations

try:
    from PyQt6.QtWidgets import (
        QWidget,
        QVBoxLayout,
        QTableWidget,
        QTableWidgetItem,
        QHeaderView,
    )
    from PyQt6.QtGui import QColor
    QT6 = True
except Exception:
    from PyQt5.QtWidgets import (
        QWidget,
        QVBoxLayout,
        QTableWidget,
        QTableWidgetItem,
        QHeaderView,
    )
    from PyQt5.QtGui import QColor
    QT6 = False


class BackendStatusWidget(QWidget):
    """
    Small table showing backend availability.
    """

    def __init__(self, parent=None):
        super().__init__(parent)

        layout = QVBoxLayout(self)

        self.table = QTableWidget(0, 4, self)
        self.table.setHorizontalHeaderLabels(
            ["Backend", "Status", "Version", "Details"]
        )

        header = self.table.horizontalHeader()
        try:
            header.setSectionResizeMode(QHeaderView.ResizeMode.Stretch if QT6 else QHeaderView.Stretch)
        except Exception:
            pass

        self.table.setAlternatingRowColors(True)

        try:
            self.table.setEditTriggers(
                QTableWidget.EditTrigger.NoEditTriggers
                if QT6
                else QTableWidget.NoEditTriggers
            )
        except Exception:
            pass

        layout.addWidget(self.table)

    def set_backends(self, rows: list[dict]) -> None:
        """
        Populate table from backend rows.
        """
        self.table.setRowCount(len(rows))

        for row_index, row in enumerate(rows):
            backend = row.get("backend", "")
            available = bool(row.get("available", False))
            version = row.get("version", "")
            details = row.get("details", "")

            status_text = "Available" if available else "Missing"

            items = [
                QTableWidgetItem(str(backend)),
                QTableWidgetItem(status_text),
                QTableWidgetItem(str(version or "")),
                QTableWidgetItem(str(details or "")),
            ]

            for col_index, item in enumerate(items):
                self.table.setItem(row_index, col_index, item)

            color = QColor("#2e7d32") if available else QColor("#b71c1c")
            items[1].setForeground(color)
