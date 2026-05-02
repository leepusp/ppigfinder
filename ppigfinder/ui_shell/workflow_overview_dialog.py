#!/usr/bin/env python3
"""
Workflow overview dialog for the experimental guided UI shell.
"""

from __future__ import annotations

try:
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import (
        QDialog,
        QVBoxLayout,
        QLabel,
        QTableWidget,
        QTableWidgetItem,
        QHeaderView,
    )
    QT6 = True
except Exception:
    from PyQt5.QtCore import Qt
    from PyQt5.QtWidgets import (
        QDialog,
        QVBoxLayout,
        QLabel,
        QTableWidget,
        QTableWidgetItem,
        QHeaderView,
    )
    QT6 = False

from ppigfinder.ui_shell.branding import apply_ppigfinder_branding
from ppigfinder.ui_shell.workflow_model import ordered_steps


def _window_flags():
    flags = Qt.WindowType.Window if QT6 else Qt.Window
    flags |= Qt.WindowType.WindowMinimizeButtonHint if QT6 else Qt.WindowMinimizeButtonHint
    flags |= Qt.WindowType.WindowMaximizeButtonHint if QT6 else Qt.WindowMaximizeButtonHint
    flags |= Qt.WindowType.WindowCloseButtonHint if QT6 else Qt.WindowCloseButtonHint
    return flags


class WorkflowOverviewDialog(QDialog):
    """
    Dialog showing the complete guided workflow with dependencies and outputs.
    """

    def __init__(self, parent=None):
        super().__init__(parent)

        self.setWindowTitle("ppigFinder Guided Workflow Map")
        self.setWindowFlags(_window_flags())
        self.resize(1280, 760)
        self.setMinimumSize(960, 600)
        self.setSizeGripEnabled(True)
        apply_ppigfinder_branding(self)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(18, 18, 18, 18)
        layout.setSpacing(12)

        title = QLabel("ppigFinder guided workflow map")
        title.setObjectName("SectionTitle")
        layout.addWidget(title)

        subtitle = QLabel(
            "The guided interface is organized around data dependencies, operations, "
            "generated outputs and visualizations. Each step unlocks the next analysis decision."
        )
        subtitle.setWordWrap(True)
        layout.addWidget(subtitle)

        self.table = QTableWidget(0, 7, self)
        self.table.setHorizontalHeaderLabels(
            [
                "Step",
                "Objective",
                "Required input",
                "Operations",
                "Outputs",
                "Visualizations",
                "Next",
            ]
        )

        header = self.table.horizontalHeader()
        try:
            mode = QHeaderView.ResizeMode.Stretch if QT6 else QHeaderView.Stretch
            header.setSectionResizeMode(mode)
        except Exception:
            pass

        self.table.setAlternatingRowColors(True)
        self._populate()
        layout.addWidget(self.table, 1)

    def _populate(self) -> None:
        steps = list(ordered_steps())
        self.table.setRowCount(len(steps))

        for row, step in enumerate(steps):
            values = [
                step.title + (" (optional)" if step.optional else ""),
                step.objective,
                "\n".join(step.required_inputs),
                "\n".join(step.operations),
                "\n".join(step.outputs),
                "\n".join(step.visualizations),
                ", ".join(step.next_steps) if step.next_steps else "End",
            ]

            for col, value in enumerate(values):
                self.table.setItem(row, col, QTableWidgetItem(value))


def show_workflow_overview(parent=None) -> None:
    dialog = WorkflowOverviewDialog(parent=parent)
    dialog.showMaximized()
    dialog.exec() if hasattr(dialog, "exec") else dialog.exec_()
