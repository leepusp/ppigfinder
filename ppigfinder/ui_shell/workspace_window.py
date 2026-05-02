#!/usr/bin/env python3
"""
Workspace shell for the future ppigFinder guided interface.
"""

from __future__ import annotations

from ppigfinder.ui_shell.module_pages import ModulePage
from ppigfinder.ui_shell.navigation import ModuleRoute
from ppigfinder.ui_shell.qt import (
    QMainWindow,
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QStackedWidget,
)
from ppigfinder.ui_shell.theme import APP_TITLE, shell_stylesheet
from ppigfinder.ui_shell.branding import apply_ppigfinder_branding


DEFAULT_ROUTES = [
    ModuleRoute(
        id="overview",
        title="Overview",
        description="General overview of the ppigFinder analysis workflow.",
        data_type="Workflow",
        status="Ready",
    ),
    ModuleRoute(
        id="genome",
        title="DNA / Genome",
        description="Genome loading, translation and sequence inspection.",
        data_type="DNA / Genome",
        status="Ready",
    ),
    ModuleRoute(
        id="orfs",
        title="Protein / ORFs",
        description="ORF prediction, ORF review and protein export.",
        data_type="Protein / ORFs",
        status="Ready",
    ),
    ModuleRoute(
        id="annotation",
        title="Annotation",
        description="BLAST, HMM/domain and neighborhood analyses.",
        data_type="Functional annotation",
        status="Ready",
    ),
    ModuleRoute(
        id="alphafold",
        title="AlphaFold / PPI",
        description="AF3 export, result import and interaction interpretation.",
        data_type="Protein interaction",
        status="Ready",
    ),
    ModuleRoute(
        id="reports",
        title="Reports",
        description="Generate HTML, JSON and tabular outputs.",
        data_type="Reporting",
        status="Ready",
    ),
]


class WorkspaceWindow(QMainWindow):
    """
    Future stepwise workspace window.
    """

    def __init__(self, bridge=None, routes=None):
        super().__init__()

        self.bridge = bridge
        self.routes = routes or DEFAULT_ROUTES

        self.setWindowTitle(f"{APP_TITLE} — Workspace")
        self.resize(1320, 820)
        self.setStyleSheet(shell_stylesheet())
        apply_ppigfinder_branding(self)

        central = QWidget()
        self.setCentralWidget(central)

        root = QHBoxLayout(central)
        root.setContentsMargins(18, 18, 18, 18)
        root.setSpacing(18)

        self.navigation = QListWidget()
        self.navigation.setMaximumWidth(280)
        root.addWidget(self.navigation)

        main_area = QVBoxLayout()
        root.addLayout(main_area, 1)

        self.breadcrumb = QLabel("Workspace > Overview")
        self.breadcrumb.setObjectName("SectionTitle")
        main_area.addWidget(self.breadcrumb)

        self.pages = QStackedWidget()
        main_area.addWidget(self.pages, 1)

        self._build_pages()
        self.navigation.currentRowChanged.connect(self._on_route_changed)

        if self.navigation.count():
            self.navigation.setCurrentRow(0)

    def _action(self, label: str, description: str, action_name: str) -> dict:
        def run():
            if self.bridge:
                self.bridge.call(action_name)

        return {
            "label": label,
            "description": description,
            "callback": run,
        }

    def _build_pages(self):
        for route in self.routes:
            item = QListWidgetItem(f"{route.title}\n{route.status}")
            item.setToolTip(
                f"{route.title}\n\n"
                f"Data type: {route.data_type}\n"
                f"Description: {route.description}\n"
                f"Status: {route.status}"
            )
            self.navigation.addItem(item)

            actions = []

            if route.id == "overview":
                actions = [
                    self._action(
                        "Open current interface",
                        "Open the current full ppigFinder interface while this modular workspace evolves.",
                        "open_legacy_interface",
                    ),
                ]
            elif route.id == "genome":
                actions = [
                    self._action(
                        "Open genome",
                        "Load a FASTA, GenBank or SnapGene nucleotide file.",
                        "load_fasta",
                    ),
                    self._action(
                        "Translate genome",
                        "Translate genome sequence or selected regions into protein sequence.",
                        "translate_genome",
                    ),
                    self._action(
                        "Open current genome tools",
                        "Open the current interface with all genome tools available.",
                        "open_legacy_interface",
                    ),
                ]
            elif route.id == "orfs":
                actions = [
                    self._action(
                        "Predict ORFs",
                        "Detect protein-coding open reading frames from the loaded genome.",
                        "analyze_orfs",
                    ),
                    self._action(
                        "Open current ORF tools",
                        "Open the current interface for ORF table inspection and filters.",
                        "open_legacy_interface",
                    ),
                ]
            elif route.id == "annotation":
                actions = [
                    self._action(
                        "Run BLAST",
                        "Search a protein query against predicted ORFs.",
                        "run_blast",
                    ),
                    self._action(
                        "Annotate HMM",
                        "Annotate conserved domains with HMMER or internal fallback scanner.",
                        "annotate_hmm",
                    ),
                    self._action(
                        "Open current annotation tools",
                        "Open the current annotation tabs and tools.",
                        "open_legacy_interface",
                    ),
                ]
            elif route.id == "alphafold":
                actions = [
                    self._action(
                        "Export AF3 Server JSON",
                        "Generate AlphaFold Server JSON for selected ORF protein pairs.",
                        "export_af3_server_json",
                    ),
                    self._action(
                        "Import AF3 results",
                        "Import AlphaFold 3 output folders and extract interaction metrics.",
                        "import_af3_results_folder",
                    ),
                    self._action(
                        "Open current AlphaFold tools",
                        "Open the current interface with AlphaFold-related tabs.",
                        "open_legacy_interface",
                    ),
                ]
            elif route.id == "reports":
                actions = [
                    self._action(
                        "Export HTML report",
                        "Generate a standalone HTML report from the current project snapshot.",
                        "export_html_report",
                    ),
                    self._action(
                        "Open current export tools",
                        "Open the current interface for additional export options.",
                        "open_legacy_interface",
                    ),
                ]

            page = ModulePage(route, actions=actions)
            self.pages.addWidget(page)

    def _on_route_changed(self, index):
        if index < 0 or index >= len(self.routes):
            return

        route = self.routes[index]
        self.breadcrumb.setText(f"Workspace > {route.title}")
        self.pages.setCurrentIndex(index)
