#!/usr/bin/env python3
"""
Reusable module pages for the workspace shell.
"""

from __future__ import annotations

from ppigfinder.ui_shell.components import InfoCard, ActionCard, FlowStrip
from ppigfinder.ui_shell.docs_content import docs_for
from ppigfinder.ui_shell.docs_dialog import show_module_documentation
from ppigfinder.ui_shell.module_visuals import ModuleVisualizationPanel
from ppigfinder.ui_shell.qt import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QGridLayout,
    QLabel,
    QFrame,
    QPushButton,
)


WORKFLOW_STEPS = [
    "Data",
    "Genome",
    "ORFs",
    "Annotation",
    "AlphaFold",
    "HPC",
    "Reports",
]


MODULE_STEP_INDEX = {
    "overview": 0,
    "data": 0,
    "genome": 1,
    "orfs": 2,
    "annotation": 3,
    "alphafold": 4,
    "hpc": 5,
    "reports": 6,
}


class ModulePage(QWidget):
    """
    Generic page for one module of the guided workspace.
    """

    def __init__(self, route, actions=None, parent=None):
        super().__init__(parent)
        self.route = route
        self.actions = actions or []

        docs = docs_for(route.id)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(24, 22, 24, 22)
        layout.setSpacing(16)

        title = QLabel(route.title)
        title.setObjectName("HeroTitle")
        layout.addWidget(title)

        subtitle = QLabel(route.description)
        subtitle.setObjectName("HeroSubtitle")
        subtitle.setWordWrap(True)
        layout.addWidget(subtitle)

        layout.addWidget(
            FlowStrip(
                WORKFLOW_STEPS,
                active_index=MODULE_STEP_INDEX.get(route.id, 0),
            )
        )

        cards = QGridLayout()
        cards.setSpacing(12)
        layout.addLayout(cards)

        cards.addWidget(
            InfoCard(
                "Input data",
                docs["input"],
                "Main data required by this module.",
            ),
            0,
            0,
        )
        cards.addWidget(
            InfoCard(
                "Expected output",
                docs["output"],
                "Result produced after completing this step.",
            ),
            0,
            1,
        )
        cards.addWidget(
            InfoCard(
                "Status",
                route.status,
                "Current readiness of this workflow module.",
            ),
            0,
            2,
        )

        purpose = QLabel(
            f"<b>Purpose</b><br>{docs['purpose']}<br><br>"
            f"<b>Next</b><br>{docs['next']}"
        )
        purpose.setWordWrap(True)
        layout.addWidget(purpose)

        docs_button = QPushButton("Read module guide")
        docs_button.clicked.connect(
            lambda checked=False: show_module_documentation(
                self.route.id,
                self.route.title,
                self,
            )
        )
        layout.addWidget(docs_button)

        body = QHBoxLayout()
        body.setSpacing(14)
        layout.addLayout(body, 1)

        actions_box = QFrame()
        actions_box.setObjectName("Card")
        actions_layout = QVBoxLayout(actions_box)
        actions_layout.setContentsMargins(16, 16, 16, 16)
        actions_layout.setSpacing(10)

        actions_title = QLabel("Actions")
        actions_title.setObjectName("SectionTitle")
        actions_layout.addWidget(actions_title)

        if self.actions:
            for action in self.actions:
                label = action.get("label", "Action")
                description = action.get("description", "")
                callback = action.get("callback")
                button_text = action.get("button_text", "Run")
                card = ActionCard(
                    label,
                    description,
                    button_text,
                    callback,
                )
                actions_layout.addWidget(card)
        else:
            actions_layout.addWidget(QLabel("No actions connected yet for this module."))

        actions_layout.addStretch(1)
        body.addWidget(actions_box, 1)

        self.visualization_panel = ModuleVisualizationPanel(route.id)
        body.addWidget(self.visualization_panel, 1)

    def _visualization_hint(self, module_id: str) -> str:
        hints = {
            "overview": (
                "Future dashboard with project status, workflow progress, recently opened projects, "
                "loaded genome summary and next-step suggestions."
            ),
            "data": (
                "Future data-entry panel showing supported input formats, current loaded file, "
                "project restoration status and validation messages before moving to ORF prediction."
            ),
            "genome": (
                "Future interactive genome map inspired by tools such as lovis4u, "
                "showing ORFs, coordinates, GC context and genomic neighborhoods."
            ),
            "orfs": (
                "Future ORF/protein summary with length distribution, strand/frame distribution "
                "and sequence export status."
            ),
            "annotation": (
                "Future annotation dashboard with BLAST hit summaries, HMM/domain composition "
                "and neighborhood context."
            ),
            "alphafold": (
                "Future AlphaFold/PPI dashboard with ipTM, cp_ipTM, PAE_min, PAE_inter, "
                "contact percentage and interaction classification cards."
            ),
            "reports": (
                "Future report/export dashboard with HTML, JSON snapshot and TSV/CSV export status."
            ),
        }

        return hints.get(module_id, "Future dynamic visualization panel.")
