#!/usr/bin/env python3
"""
Reusable UI shell components for ppigFinder.
"""

from __future__ import annotations

from ppigfinder.ui_shell.qt import (
    QFrame,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QHBoxLayout,
    QWidget,
)


class InfoCard(QFrame):
    """
    Simple information card used in the guided workspace.
    """

    def __init__(self, title: str, value: str, description: str = "", parent=None):
        super().__init__(parent)
        self.setObjectName("Card")

        layout = QVBoxLayout(self)
        layout.setContentsMargins(14, 12, 14, 12)
        layout.setSpacing(6)

        title_label = QLabel(title)
        title_label.setObjectName("CardTitle")
        layout.addWidget(title_label)

        value_label = QLabel(value)
        value_label.setObjectName("MetricValue")
        value_label.setWordWrap(True)
        layout.addWidget(value_label)

        if description:
            desc_label = QLabel(description)
            desc_label.setObjectName("CardDescription")
            desc_label.setWordWrap(True)
            layout.addWidget(desc_label)


class ActionCard(QFrame):
    """
    Action card with description and a button.
    """

    def __init__(self, title: str, description: str, button_text: str, callback, parent=None):
        super().__init__(parent)
        self.setObjectName("Card")

        layout = QVBoxLayout(self)
        layout.setContentsMargins(14, 12, 14, 12)
        layout.setSpacing(8)

        title_label = QLabel(title)
        title_label.setObjectName("CardTitle")
        layout.addWidget(title_label)

        desc_label = QLabel(description)
        desc_label.setObjectName("CardDescription")
        desc_label.setWordWrap(True)
        layout.addWidget(desc_label)

        layout.addStretch(1)

        button = QPushButton(button_text)
        button.clicked.connect(callback)
        layout.addWidget(button)


class FlowStrip(QWidget):
    """
    Horizontal workflow strip.
    """

    def __init__(self, steps: list[str], active_index: int = 0, parent=None):
        super().__init__(parent)

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)

        for index, step in enumerate(steps):
            label = QLabel(step)
            label.setObjectName("FlowStepActive" if index == active_index else "FlowStep")
            label.setWordWrap(True)
            layout.addWidget(label)

            if index < len(steps) - 1:
                arrow = QLabel("→")
                arrow.setObjectName("FlowArrow")
                layout.addWidget(arrow)
