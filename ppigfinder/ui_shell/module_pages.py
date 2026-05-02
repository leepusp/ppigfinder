#!/usr/bin/env python3
from __future__ import annotations

try:
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import (
        QWidget,
        QVBoxLayout,
        QHBoxLayout,
        QLabel,
        QPushButton,
        QFrame,
        QSplitter,
        QScrollArea,
        QMessageBox,
    )
    QT6 = True
except Exception:
    from PyQt5.QtCore import Qt
    from PyQt5.QtWidgets import (
        QWidget,
        QVBoxLayout,
        QHBoxLayout,
        QLabel,
        QPushButton,
        QFrame,
        QSplitter,
        QScrollArea,
        QMessageBox,
    )
    QT6 = False

from ppigfinder.ui_shell.visual_panels import FlowRibbonWidget, InfoCard, ModuleVisualizationPanel


class ActionItem(QFrame):
    def __init__(self, label: str, description: str, button_text: str, callback, parent=None):
        super().__init__(parent)
        self.setObjectName("ActionItem")

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 10, 12, 10)
        layout.setSpacing(6)

        title = QLabel(label)
        title.setObjectName("ActionTitle")

        desc = QLabel(description)
        desc.setWordWrap(True)
        desc.setObjectName("ActionDescription")

        btn_row = QHBoxLayout()
        btn_row.addStretch(1)

        button = QPushButton(button_text)
        button.clicked.connect(callback)

        btn_row.addWidget(button)

        layout.addWidget(title)
        layout.addWidget(desc)
        layout.addLayout(btn_row)


class ModulePage(QWidget):
    def __init__(self, route, actions: list[dict] | None = None, parent=None, on_route_selected=None):
        super().__init__(parent)

        self.route = route
        self.actions = actions or []

        root = QVBoxLayout(self)
        root.setContentsMargins(6, 6, 6, 6)
        root.setSpacing(12)

        self.page_title = QLabel(route.title)
        self.page_title.setObjectName("SectionTitle")

        self.page_subtitle = QLabel(route.description)
        self.page_subtitle.setWordWrap(True)
        self.page_subtitle.setObjectName("SectionSubtitle")

        self.ribbon = FlowRibbonWidget(on_route_selected=on_route_selected)

        root.addWidget(self.page_title)
        root.addWidget(self.page_subtitle)
        root.addWidget(self.ribbon)

        summary_row = QHBoxLayout()
        summary_row.setSpacing(10)

        self.input_card = InfoCard("Input data", "", "Main data required by this module.")
        self.output_card = InfoCard("Expected output", "", "Result produced after completing this step.")
        self.status_card = InfoCard("Status", route.status, "Current readiness of this workflow module.")

        summary_row.addWidget(self.input_card, 3)
        summary_row.addWidget(self.output_card, 2)
        summary_row.addWidget(self.status_card, 2)

        root.addLayout(summary_row)

        self.purpose_label = QLabel("")
        self.purpose_label.setWordWrap(True)
        root.addWidget(self.purpose_label)

        self.next_label = QLabel("")
        self.next_label.setWordWrap(True)
        root.addWidget(self.next_label)

        self.guide_button = QPushButton("Read module guide")
        self.guide_button.clicked.connect(self._show_guide)
        root.addWidget(self.guide_button)

        splitter = QSplitter()
        splitter.setChildrenCollapsible(False)
        try:
            splitter.setOrientation(Qt.Orientation.Horizontal if QT6 else Qt.Horizontal)
        except Exception:
            pass

        # Actions side
        action_container = QWidget()
        action_layout = QVBoxLayout(action_container)
        action_layout.setContentsMargins(0, 0, 0, 0)
        action_layout.setSpacing(10)

        action_title = QLabel("Actions")
        action_title.setObjectName("SectionSubTitle")
        action_layout.addWidget(action_title)

        self.action_scroll = QScrollArea()
        self.action_scroll.setWidgetResizable(True)

        self.action_body = QWidget()
        self.action_body_layout = QVBoxLayout(self.action_body)
        self.action_body_layout.setContentsMargins(0, 0, 0, 0)
        self.action_body_layout.setSpacing(10)

        for action in self.actions:
            item = ActionItem(
                label=action.get("label", "Action"),
                description=action.get("description", ""),
                button_text=action.get("button_text", "Run"),
                callback=action.get("callback", lambda: None),
            )
            self.action_body_layout.addWidget(item)

        self.action_body_layout.addStretch(1)
        self.action_scroll.setWidget(self.action_body)
        action_layout.addWidget(self.action_scroll, 1)

        splitter.addWidget(action_container)

        # Visualization side
        self.visualization_panel = ModuleVisualizationPanel(route.id)
        splitter.addWidget(self.visualization_panel)

        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 2)

        root.addWidget(splitter, 1)

        self._set_static_texts()

    def _show_guide(self):
        try:
            from ppigfinder.ui_shell.docs_dialog import open_module_guide_dialog
            open_module_guide_dialog(self.route.id, parent=self)
            return
        except Exception:
            pass

        QMessageBox.information(
            self,
            f"{self.route.title} guide",
            f"Module: {self.route.title}\n\n"
            f"Data type: {self.route.data_type}\n\n"
            f"Description:\n{self.route.description}\n\n"
            f"Documentation dialog is not fully connected yet.",
        )

    def _set_static_texts(self):
        text_map = {
            "overview": (
                "Purpose: Visual overview of the guided workflow and its dependencies.",
                "Next: Start with Data / Project.",
                "Workflow status overview",
                "Progress and next recommended action",
            ),
            "data": (
                "Purpose: Start from input data. Load a genome, project or snapshot and ppigFinder will enable the next analysis steps according to the data type.",
                "Next: A valid genome automatically enables ORF prediction. A project or snapshot restores workflow state for review/export.",
                "Genome, protein query, HMM profile, project/snapshot JSON, or AF3 result folder",
                "Recognized input routed to the appropriate workflow step",
            ),
            "genome": (
                "Purpose: Inspect the loaded genome and validate sequence-scale information before prediction.",
                "Next: Continue to Protein / ORFs.",
                "Loaded genome sequence and metadata",
                "Genome inspection state and coordinate-ready genome context",
            ),
            "orfs": (
                "Purpose: Identify ORFs and generate protein sequences for downstream annotation.",
                "Next: Continue to Annotation.",
                "Loaded genome sequence",
                "Predicted ORFs, protein sequences and ORF tables",
            ),
            "annotation": (
                "Purpose: Assign similarity hits, domains and genomic context.",
                "Next: Select candidates for AlphaFold / PPI analysis.",
                "Predicted protein sequences and ORF coordinates",
                "BLAST/HMM/neighbourhood annotations",
            ),
            "alphafold": (
                "Purpose: Build AF3 candidate jobs and organize structural interaction analysis.",
                "Next: Optionally continue to HPC or directly to Reports.",
                "Annotated candidate proteins / ORFs",
                "AF3 JSON, candidate pairs and imported AF3 results",
            ),
            "hpc": (
                "Purpose: Configure optional server execution and DaVinci/HPC connectivity.",
                "Next: Continue to Reports.",
                "AF3 jobs and execution configuration",
                "HPC configuration, connection state and submission preparation",
            ),
            "reports": (
                "Purpose: Export summaries, reports and reproducible workflow outputs.",
                "Next: Finalize and share results.",
                "Current project/workflow state",
                "Markdown/HTML/JSON/TSV outputs",
            ),
        }

        purpose, next_step, input_text, output_text = text_map.get(
            self.route.id,
            ("Purpose: N/A", "Next: N/A", "N/A", "N/A"),
        )

        self.purpose_label.setText(f"<b>{purpose.split(':',1)[0]}:</b> {purpose.split(':',1)[1].strip()}")
        self.next_label.setText(f"<b>{next_step.split(':',1)[0]}:</b> {next_step.split(':',1)[1].strip()}")

        self.input_card.set_content("Input data", input_text, "Main data required by this module.")
        self.output_card.set_content("Expected output", output_text, "Result produced after completing this step.")

    def update_state(self, workflow_state) -> None:
        self.ribbon.update_state(workflow_state)
        self.visualization_panel.update_state(workflow_state)
