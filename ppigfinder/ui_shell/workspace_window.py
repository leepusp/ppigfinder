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
from ppigfinder.ui_shell.genome_inspector_dialog import show_genome_inspector
from ppigfinder.ui_shell.orf_results_dialog import GuidedORFResultsDialog
from ppigfinder.ui_shell.annotation_candidates_dialog import AnnotationCandidatesDialog
from ppigfinder.ui_shell.hpc_dialog import open_hpc_connection_dialog
from ppigfinder.ui_shell.workflow_overview_dialog import show_workflow_overview
from ppigfinder.ui_shell.af3_pair_builder_dialog import open_af3_pair_builder_dialog
from ppigfinder.ui_shell.guided_backend import (
    predict_orfs_from_file,
    write_orfs_fasta,
    write_guided_summary,
)


DEFAULT_ROUTES = [
    ModuleRoute("overview", "Overview", "General overview of the ppigFinder analysis workflow.", "Workflow", "Ready"),
    ModuleRoute("data", "Data / Project", "Start a project, open genome data or restore previous analyses.", "Project / Input data", "Ready"),
    ModuleRoute("genome", "DNA / Genome", "Genome loading, translation and sequence inspection.", "DNA / Genome", "Ready"),
    ModuleRoute("orfs", "Protein / ORFs", "ORF prediction, ORF review and protein export.", "Protein / ORFs", "Ready"),
    ModuleRoute("annotation", "Annotation", "BLAST, HMM/domain and neighbourhood analyses.", "Functional annotation", "Ready"),
    ModuleRoute("alphafold", "AlphaFold / PPI", "AF3 export, result import and interaction interpretation.", "Protein interaction", "Ready"),
    ModuleRoute("hpc", "DaVinci / HPC", "Optional server/HPC execution and workflow preparation.", "HPC / Remote execution", "Optional"),
    ModuleRoute("reports", "Reports", "Generate HTML, JSON and tabular outputs.", "Reporting", "Ready"),
]


class WorkspaceWindow(QMainWindow):
    """
    Data-driven guided workflow window.
    """

    def __init__(self, bridge=None, routes=None):
        super().__init__()

        self.bridge = bridge
        self.routes = routes or DEFAULT_ROUTES
        self.route_index_by_id = {}
        self.selected_inputs = {}
        self.module_pages_by_id = {}
        self.guided_orfs = []
        self._floating_windows = []

        self.setWindowTitle(f"{APP_TITLE} — Guided Workflow")
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

        self.statusBar().showMessage("Ready. Start by loading genome data or opening a project.")

    def _button_text_for_action(self, label: str) -> str:
        lower = label.lower()
        if lower.startswith("open"):
            return "Open"
        if lower.startswith("import"):
            return "Import"
        if lower.startswith("export"):
            return "Export"
        if lower.startswith("predict"):
            return "Predict"
        if lower.startswith("annotate"):
            return "Annotate"
        if lower.startswith("run"):
            return "Run"
        if lower.startswith("continue"):
            return "Continue"
        if lower.startswith("show"):
            return "Show"
        if lower.startswith("build"):
            return "Build"
        return "Continue"

    def _show_message(self, title: str, message: str) -> None:
        QMessageBox.information(self, title, message)

    def _remember_window(self, dialog) -> None:
        self._floating_windows.append(dialog)
        try:
            dialog.destroyed.connect(lambda *_: self._cleanup_windows())
        except Exception:
            pass

    def _cleanup_windows(self) -> None:
        self._floating_windows = [w for w in self._floating_windows if w is not None]

    def _select_file(self, key: str, title: str, file_filter: str) -> bool:
        path, _ = QFileDialog.getOpenFileName(self, title, "", file_filter)

        if not path:
            return False

        self.selected_inputs[key] = path

        if key == "genome_file":
            summary = validate_genome_input(path)
            self.selected_inputs.update(summary_to_state(summary))

            dialog = show_genome_inspector(path, parent=self)
            self._remember_window(dialog)

            self._refresh_visualization_panels()

            if summary.valid:
                self.statusBar().showMessage(
                    "Genome loaded and validated. Moving to Protein / ORFs.",
                    10000,
                )
                self.show_route("orfs")
            else:
                self.statusBar().showMessage(
                    "Genome selected, but validation found a problem.",
                    10000,
                )
                self.show_route("data")

            return True

        self.statusBar().showMessage(f"Selected file: {path}", 10000)
        self._refresh_visualization_panels()
        return True

    def _select_folder(self, key: str, title: str) -> bool:
        path = QFileDialog.getExistingDirectory(self, title, "")

        if not path:
            return False

        self.selected_inputs[key] = path
        self.statusBar().showMessage(f"Selected folder: {path}", 10000)
        self._refresh_visualization_panels()
        return True

    def _predict_guided_orfs(self) -> None:
        genome_file = self.selected_inputs.get("genome_file", "")

        if not genome_file:
            self._show_message(
                "Predict ORFs",
                "Select a genome file in Data / Project before predicting ORFs.",
            )
            self.show_route("data")
            return

        summary = predict_orfs_from_file(genome_file, min_aa=30)
        self.guided_orfs = summary.orfs

        self.selected_inputs.update(
            {
                "guided_orf_count": summary.orf_count,
                "guided_longest_orf_aa": summary.longest_orf_aa,
                "guided_shortest_orf_aa": summary.shortest_orf_aa,
                "guided_orf_min_aa": summary.min_aa,
                "guided_orf_source": summary.source_file,
                "guided_orf_map": [
                    {
                        "id": orf.id,
                        "start": orf.start,
                        "end": orf.end,
                        "strand": orf.strand,
                        "frame": orf.frame,
                        "aa_length": orf.aa_length,
                    }
                    for orf in summary.orfs
                ],
            }
        )

        self._refresh_visualization_panels()

        dialog = GuidedORFResultsDialog(self.guided_orfs, parent=self)
        dialog.show()
        self._remember_window(dialog)

        self.statusBar().showMessage(
            f"ORF prediction completed: {summary.orf_count} ORFs. Moving to Annotation.",
            12000,
        )
        self.show_route("annotation")

    def _show_guided_orfs(self) -> None:
        if not self.guided_orfs:
            self._show_message(
                "Review ORF table",
                "No guided ORFs are available yet. Run Predict ORFs first.",
            )
            return

        dialog = GuidedORFResultsDialog(self.guided_orfs, parent=self)
        dialog.show()
        self._remember_window(dialog)

    def _export_guided_orfs_fasta(self) -> None:
        if not self.guided_orfs:
            self._show_message(
                "Export ORF FASTA",
                "No guided ORFs are available yet. Run Predict ORFs first.",
            )
            return

        path, _ = QFileDialog.getSaveFileName(
            self,
            "Export guided ORF proteins",
            "guided_orfs.faa",
            "Protein FASTA (*.faa *.fasta *.fa);;All files (*)",
        )

        if not path:
            return

        write_orfs_fasta(path, self.guided_orfs)
        self.statusBar().showMessage(f"ORF protein FASTA exported: {path}", 12000)

    def _show_annotation_candidates(self) -> None:
        if not self.guided_orfs:
            self._show_message(
                "Annotation candidates",
                "No ORFs are available yet. Run Predict ORFs first.",
            )
            self.show_route("orfs")
            return

        self.selected_inputs["guided_annotation_candidates_count"] = len(self.guided_orfs)
        self._refresh_visualization_panels()

        dialog = AnnotationCandidatesDialog(self.guided_orfs, parent=self)
        dialog.show()
        self._remember_window(dialog)

    def _mark_annotation_step(self, key: str, label: str) -> None:
        self.selected_inputs[key] = True
        self._refresh_visualization_panels()
        self.statusBar().showMessage(f"{label} selected in guided workflow.", 10000)

    def _open_af3_pair_builder(self) -> None:
        if not self.guided_orfs:
            self._show_message(
                "AlphaFold / PPI",
                "No ORFs are available yet. Run Predict ORFs first.",
            )
            self.show_route("orfs")
            return

        result = open_af3_pair_builder_dialog(self.guided_orfs, parent=self)

        self.selected_inputs["af3_pair_count"] = result.get("pair_count", 0)

        if result.get("json_path"):
            self.selected_inputs["af3_json_path"] = result["json_path"]
            self.selected_inputs["af3_json_exported"] = True

        self._refresh_visualization_panels()

    def _open_hpc_dialog(self) -> None:
        status, config = open_hpc_connection_dialog(parent=self)

        self.selected_inputs["hpc_profile"] = config.profile
        self.selected_inputs["hpc_host"] = config.host
        self.selected_inputs["hpc_user"] = config.user
        self.selected_inputs["hpc_port"] = config.port

        if status is not None:
            self.selected_inputs["hpc_status"] = "OK" if status.connection_ok else "Problem"
            self.selected_inputs["hpc_mode"] = "Local cluster" if status.running_on_cluster else "SSH"
            self.selected_inputs["hpc_message"] = status.message
        else:
            self.selected_inputs["hpc_status"] = "Configured"
            self.selected_inputs["hpc_mode"] = "Not tested"

        self._refresh_visualization_panels()

    def _export_guided_summary(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Export guided workflow summary",
            "ppigfinder_guided_summary.md",
            "Markdown (*.md);;Text files (*.txt);;All files (*)",
        )

        if not path:
            return

        write_guided_summary(path, self.selected_inputs)
        self.statusBar().showMessage(f"Guided workflow summary exported: {path}", 12000)

    def _handle_action(self, action_name: str) -> None:
        if action_name.startswith("route:"):
            self.show_route(action_name.split(":", 1)[1])
            return

        if action_name == "guided:workflow_map":
            show_workflow_overview(parent=self)
            return

        if action_name == "select_file:genome":
            self._select_file("genome_file", "Open genome file", "Genome files (*.fasta *.fa *.fna *.gb *.gbk *.dna);;All files (*)")
            return

        if action_name == "select_file:project":
            self._select_file("project_file", "Open ppigFinder project", "ppigFinder project (*.json *.ppigfinder.json);;All files (*)")
            return

        if action_name == "select_file:snapshot":
            self._select_file("snapshot_file", "Import Project Snapshot v3", "Project Snapshot (*.json *.ppigfinder.json);;All files (*)")
            return

        if action_name == "select_folder:af3_results":
            self._select_folder("af3_results_folder", "Import AlphaFold/AF3 results folder")
            return

        if action_name == "guided:predict_orfs":
            self._predict_guided_orfs()
            return

        if action_name == "guided:review_orfs":
            self._show_guided_orfs()
            return

        if action_name == "guided:export_orfs_fasta":
            self._export_guided_orfs_fasta()
            return

        if action_name == "guided:annotation_candidates":
            self._show_annotation_candidates()
            return

        if action_name == "guided:blast":
            self._mark_annotation_step("guided_blast_planned", "BLAST")
            return

        if action_name == "guided:hmm":
            self._mark_annotation_step("guided_hmm_planned", "HMM/domain annotation")
            return

        if action_name == "guided:neighborhood":
            self._mark_annotation_step("guided_neighborhood_planned", "Neighbourhood analysis")
            return

        if action_name == "guided:af3_pair_builder":
            self._open_af3_pair_builder()
            return

        if action_name == "guided:hpc_connection":
            self._open_hpc_dialog()
            return

        if action_name == "guided:summary":
            self._export_guided_summary()
            return

        if action_name.startswith("info:"):
            self._show_message("Guided workflow", action_name.split(":", 1)[1])
            return

        if action_name == "open_legacy_interface":
            if self.bridge:
                self.bridge.call(action_name)
            return

        self._show_message("Guided workflow", "Action not connected yet:\n\n" + action_name)

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

            page = ModulePage(route, actions=self._actions_for_route(route.id))
            self.module_pages_by_id[route.id] = page
            self.pages.addWidget(page)

    def _actions_for_route(self, route_id: str) -> list[dict]:
        if route_id == "overview":
            return [
                self._action("Show workflow map", "Review the complete data-driven ppigFinder workflow.", "guided:workflow_map"),
                self._action("Start with Data / Project", "Begin by adding genome data, opening a project or importing a project snapshot.", "route:data"),
            ]

        if route_id == "data":
            return [
                self._action("Open genome file", "Load and validate the main nucleotide dataset. The workflow will automatically move to ORF Discovery if valid.", "select_file:genome"),
                self._action("Open project", "Resume a previous ppigFinder session.", "select_file:project"),
                self._action("Import Project Snapshot v3", "Load a portable JSON snapshot for reproducible analysis.", "select_file:snapshot"),
            ]

        if route_id == "genome":
            return [
                self._action("Open genome file", "Replace or load a genome file and inspect its metadata.", "select_file:genome"),
                self._action("Continue to Protein / ORFs", "Proceed to ORF prediction and protein sequence generation.", "route:orfs"),
            ]

        if route_id == "orfs":
            return [
                self._action("Predict ORFs", "Run a lightweight six-frame ORF scan inside the guided shell.", "guided:predict_orfs"),
                self._action("Review guided ORF table", "Inspect guided ORF coordinates, strand, frame and protein previews.", "guided:review_orfs"),
                self._action("Export ORF FASTA", "Export guided ORF protein sequences.", "guided:export_orfs_fasta"),
                self._action("Continue to Annotation", "Move to BLAST, HMM/domain and neighbourhood analysis.", "route:annotation"),
            ]

        if route_id == "annotation":
            return [
                self._action("Review candidate ORFs", "Inspect guided ORFs as candidates for annotation and PPI analysis.", "guided:annotation_candidates"),
                self._action("Run BLAST", "Mark BLAST as selected in the guided flow.", "guided:blast"),
                self._action("Annotate HMM", "Mark HMM/domain annotation as selected in the guided flow.", "guided:hmm"),
                self._action("Neighbourhood analysis", "Mark neighbourhood analysis as selected in the guided flow.", "guided:neighborhood"),
                self._action("Continue to AlphaFold / PPI", "Move from annotation evidence to structural interaction analysis.", "route:alphafold"),
            ]

        if route_id == "alphafold":
            return [
                self._action("Build AF3 candidate pairs", "Generate adjacent ORF candidate pairs and optionally export AlphaFold Server JSON.", "guided:af3_pair_builder"),
                self._action("Import AF3 results", "Import AlphaFold 3 output folders.", "select_folder:af3_results"),
                self._action("Continue to DaVinci / HPC", "Optionally prepare server execution before reporting.", "route:hpc"),
            ]

        if route_id == "hpc":
            return [
                self._action("Configure DaVinci / HPC", "Open graphical configuration and connection testing for DaVinci/HPC workflows.", "guided:hpc_connection"),
                self._action("Prepare AF3 Slurm template", "Use the HPC dialog to inspect a basic AF3 Slurm template.", "guided:hpc_connection"),
                self._action("Continue to Reports", "Proceed to result export and reporting.", "route:reports"),
            ]

        if route_id == "reports":
            return [
                self._action("Export guided summary", "Export the current guided shell state as a Markdown summary.", "guided:summary"),
                self._action("Export Project Snapshot v3", "Project Snapshot export will be connected after full guided state synchronization.", "info:Project Snapshot export will be connected after guided state synchronization."),
                self._action("Open full current interface", "Advanced option: open the complete current interface.", "open_legacy_interface"),
            ]

        return []

    def _refresh_visualization_panels(self) -> None:
        for page in self.module_pages_by_id.values():
            panel = getattr(page, "visualization_panel", None)
            if panel is not None and hasattr(panel, "update_state"):
                panel.update_state(self.selected_inputs)

    def show_route(self, route_id: str) -> bool:
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
        self._refresh_visualization_panels()
