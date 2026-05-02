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
    QFileDialog,
    QMessageBox,
)
from ppigfinder.ui_shell.theme import APP_TITLE, shell_stylesheet
from ppigfinder.ui_shell.branding import apply_ppigfinder_branding
from ppigfinder.ui_shell.input_validation import validate_genome_input, summary_to_state


DEFAULT_ROUTES = [
    ModuleRoute(
        id="overview",
        title="Overview",
        description="General overview of the ppigFinder analysis workflow.",
        data_type="Workflow",
        status="Ready",
    ),
    ModuleRoute(
        id="data",
        title="Data / Project",
        description="Start a project, open genome data or restore previous analyses.",
        data_type="Project / Input data",
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
        description="BLAST, HMM/domain and neighbourhood analyses.",
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

    At this stage, the shell keeps actions inside the guided workflow.
    Full backend execution will be connected progressively.
    """

    def __init__(self, bridge=None, routes=None):
        super().__init__()

        self.bridge = bridge
        self.routes = routes or DEFAULT_ROUTES
        self.route_index_by_id = {}
        self.selected_inputs = {}
        self.module_pages_by_id = {}

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

    def _button_text_for_action(self, label: str) -> str:
        label_lower = label.lower()

        if label_lower.startswith("open"):
            return "Open"
        if label_lower.startswith("import"):
            return "Import"
        if label_lower.startswith("export"):
            return "Export"
        if label_lower.startswith("predict"):
            return "Predict"
        if label_lower.startswith("annotate"):
            return "Annotate"
        if label_lower.startswith("run"):
            return "Run"
        if label_lower.startswith("continue"):
            return "Continue"

        return "Continue"

    def _show_message(self, title: str, message: str) -> None:
        QMessageBox.information(self, title, message)

    def _select_file(self, key: str, title: str, file_filter: str) -> bool:
        path, _ = QFileDialog.getOpenFileName(self, title, "", file_filter)

        if not path:
            return False

        self.selected_inputs[key] = path

        if key == "genome_file":
            summary = validate_genome_input(path)
            self.selected_inputs.update(summary_to_state(summary))

            message = (
                "Selected genome file:\n\n"
                + path
                + "\n\nValidation: "
                + ("OK" if summary.valid else "Problem detected")
                + "\nType: "
                + summary.file_type
                + "\nLength: "
                + str(summary.total_length or "N/A")
                + "\nGC%: "
                + str(summary.gc_percent if summary.gc_percent is not None else "N/A")
                + "\n\n"
                + summary.message
            )
        else:
            message = (
                "Selected file:\n\n"
                + path
                + "\n\nThis file is now registered in the guided shell. "
                + "Full backend parsing will be connected progressively."
            )

        self._show_message("Input selected", message)
        self._refresh_visualization_panels()
        return True

    def _select_folder(self, key: str, title: str) -> bool:
        path = QFileDialog.getExistingDirectory(self, title, "")

        if not path:
            return False

        self.selected_inputs[key] = path

        message = (
            "Selected folder:\n\n"
            + path
            + "\n\nThis folder is now registered in the guided shell. "
            + "Full backend parsing will be connected progressively."
        )

        self._show_message("Folder selected", message)
        self._refresh_visualization_panels()
        return True

    def _handle_action(self, action_name: str) -> None:
        if action_name.startswith("route:"):
            route_id = action_name.split(":", 1)[1]
            self.show_route(route_id)
            return

        if action_name == "select_file:genome":
            self._select_file(
                "genome_file",
                "Open genome file",
                "Genome files (*.fasta *.fa *.fna *.gb *.gbk *.dna);;All files (*)",
            )
            return

        if action_name == "select_file:project":
            self._select_file(
                "project_file",
                "Open ppigFinder project",
                "ppigFinder project (*.json *.ppigfinder.json);;All files (*)",
            )
            return

        if action_name == "select_file:snapshot":
            self._select_file(
                "snapshot_file",
                "Import Project Snapshot v3",
                "Project Snapshot (*.json *.ppigfinder.json);;All files (*)",
            )
            return

        if action_name == "select_folder:af3_results":
            self._select_folder(
                "af3_results_folder",
                "Import AlphaFold/AF3 results folder",
            )
            return

        if action_name.startswith("info:"):
            message = action_name.split(":", 1)[1]
            self._show_message("Guided workflow", message)
            return

        if action_name == "open_legacy_interface":
            if self.bridge:
                self.bridge.call(action_name)
            return

        message = (
            "This action is not connected in the guided shell yet:\n\n"
            + action_name
        )
        self._show_message("Guided workflow", message)

    def _action(self, label: str, description: str, action_name: str) -> dict:
        def run():
            self._handle_action(action_name)

        return {
            "label": label,
            "description": description,
            "button_text": self._button_text_for_action(label),
            "callback": run,
        }

    def _build_pages(self):
        for index, route in enumerate(self.routes):
            self.route_index_by_id[route.id] = index

            item = QListWidgetItem(f"{route.title}\n{route.status}")
            item.setToolTip(
                f"{route.title}\n\n"
                f"Data type: {route.data_type}\n"
                f"Description: {route.description}\n"
                f"Status: {route.status}"
            )
            self.navigation.addItem(item)

            actions = self._actions_for_route(route.id)
            page = ModulePage(route, actions=actions)
            self.module_pages_by_id[route.id] = page
            self.pages.addWidget(page)

    def _actions_for_route(self, route_id: str) -> list[dict]:
        if route_id == "overview":
            return [
                self._action(
                    "Start with Data / Project",
                    "Begin by adding genome data, opening a project or importing a project snapshot.",
                    "route:data",
                ),
            ]

        if route_id == "data":
            return [
                self._action(
                    "Open genome file",
                    "Insert the main nucleotide dataset for a new bacterial genome analysis.",
                    "select_file:genome",
                ),
                self._action(
                    "Open project",
                    "Resume a previous ppigFinder session with genome data, ORFs, annotations and analysis state.",
                    "select_file:project",
                ),
                self._action(
                    "Import Project Snapshot v3",
                    "Load a portable JSON snapshot for reproducible analysis, reporting or continuation.",
                    "select_file:snapshot",
                ),
                self._action(
                    "Continue to Protein / ORFs",
                    "Move to ORF prediction after adding or restoring input data.",
                    "route:orfs",
                ),
            ]

        if route_id == "genome":
            return [
                self._action(
                    "Review genome data",
                    "Inspect genome information and prepare for downstream ORF prediction.",
                    "info:Genome inspection will be connected to guided state after input parsing is wired.",
                ),
                self._action(
                    "Continue to Protein / ORFs",
                    "Proceed to ORF prediction and protein sequence generation.",
                    "route:orfs",
                ),
            ]

        if route_id == "orfs":
            return [
                self._action(
                    "Predict ORFs",
                    "Detect protein-coding open reading frames from the loaded genome.",
                    "info:ORF prediction will be connected directly to this guided module after backend state binding is completed.",
                ),
                self._action(
                    "Review ORF table",
                    "Inspect coordinates, strand, frame, size, GC content and protein sequence output.",
                    "info:ORF table review will be embedded in this module as the guided interface evolves.",
                ),
                self._action(
                    "Continue to Annotation",
                    "Move to BLAST, HMM/domain and neighbourhood analysis.",
                    "route:annotation",
                ),
            ]

        if route_id == "annotation":
            return [
                self._action(
                    "Run BLAST",
                    "Search a protein query against predicted ORFs.",
                    "info:BLAST search will be connected inside the Annotation module after guided query input is added.",
                ),
                self._action(
                    "Annotate HMM",
                    "Annotate conserved domains with HMMER or internal fallback scanner.",
                    "info:HMM/domain annotation will be connected after guided profile selection is added.",
                ),
                self._action(
                    "Continue to AlphaFold / PPI",
                    "Move from annotation evidence to structural interaction analysis.",
                    "route:alphafold",
                ),
            ]

        if route_id == "alphafold":
            return [
                self._action(
                    "Export AF3 Server JSON",
                    "Generate AlphaFold Server JSON for selected ORF protein pairs.",
                    "info:AF3 Server JSON export will be connected after guided ORF/protein pair selection is implemented.",
                ),
                self._action(
                    "Import AF3 results",
                    "Import AlphaFold 3 output folders and extract interaction metrics.",
                    "select_folder:af3_results",
                ),
                self._action(
                    "Continue to Reports",
                    "Move to result export and reporting.",
                    "route:reports",
                ),
            ]

        if route_id == "reports":
            return [
                self._action(
                    "Export HTML report",
                    "Generate a standalone HTML report from the current project snapshot.",
                    "info:HTML report export will be connected after guided project state is synchronized.",
                ),
                self._action(
                    "Export Project Snapshot v3",
                    "Export a versioned JSON snapshot for reproducibility.",
                    "info:Project Snapshot export will be connected after guided state synchronization.",
                ),
                self._action(
                    "Open full current interface",
                    "Advanced option: open the complete current interface after reviewing the guided workflow.",
                    "open_legacy_interface",
                ),
            ]

        return []


    def _refresh_visualization_panels(self) -> None:
        """
        Refresh visualization/status panels using the current guided state.
        """
        for page in self.module_pages_by_id.values():
            panel = getattr(page, "visualization_panel", None)
            if panel is not None and hasattr(panel, "update_state"):
                panel.update_state(self.selected_inputs)

    def show_route(self, route_id: str) -> bool:
        """
        Switch workspace to a route by ID.
        """
        index = self.route_index_by_id.get(route_id)

        if index is None:
            return False

        self.navigation.setCurrentRow(index)
        return True

    def _on_route_changed(self, index):
        if index < 0 or index >= len(self.routes):
            return

        route = self.routes[index]
        self.breadcrumb.setText(f"Workspace > {route.title}")
        self.pages.setCurrentIndex(index)
