#!/usr/bin/env python3
from __future__ import annotations

try:
    from PyQt6.QtCore import QTimer
    from PyQt6.QtWidgets import (
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
except Exception:
    from PyQt5.QtCore import QTimer
    from PyQt5.QtWidgets import (
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

from ppigfinder.ui_shell.module_pages import ModulePage
from ppigfinder.ui_shell.navigation import ModuleRoute
from ppigfinder.ui_shell.theme import APP_TITLE, shell_stylesheet
from ppigfinder.ui_shell.branding import apply_ppigfinder_branding
from ppigfinder.ui_shell.workflow_state import WorkflowState
from ppigfinder.ui_shell.input_start_dialog import choose_initial_data


DEFAULT_ROUTES = [
    ModuleRoute("overview", "Overview", "General overview of the ppigFinder analysis workflow.", "Workflow", "Ready"),
    ModuleRoute("data", "Data / Project", "Start a project, open genome data or restore previous analyses.", "Project / Input data", "Start here"),
    ModuleRoute("genome", "DNA / Genome", "Genome loading, translation and sequence inspection.", "DNA / Genome", "Waiting input"),
    ModuleRoute("orfs", "Protein / ORFs", "ORF prediction, ORF review and protein export.", "Protein / ORFs", "Waiting genome"),
    ModuleRoute("annotation", "Annotation", "BLAST, HMM/domain and neighbourhood analyses.", "Functional annotation", "Waiting ORFs"),
    ModuleRoute("alphafold", "AlphaFold / PPI", "AF3 export, result import and interaction interpretation.", "Protein interaction", "Waiting candidates"),
    ModuleRoute("hpc", "DaVinci / HPC", "Optional server/HPC execution and workflow preparation.", "HPC / Remote execution", "Optional"),
    ModuleRoute("reports", "Reports", "Generate HTML, JSON and tabular outputs.", "Reporting", "Waiting outputs"),
]


class WorkspaceWindow(QMainWindow):
    """
    Data-first guided workflow.

    The workspace opens on Data / Project. If no input is loaded, it prompts
    the user to choose a genome/project/snapshot and then enables downstream
    steps according to the generated state.
    """

    def __init__(self, bridge=None, routes=None):
        super().__init__()

        self.bridge = bridge
        self.routes = routes or DEFAULT_ROUTES
        self.route_index_by_id = {}
        self.navigation_items_by_id = {}
        self.module_pages_by_id = {}
        self.workflow_state = WorkflowState()
        self.guided_orfs = []
        self._floating_windows = []
        self._auto_data_prompt_done = False

        self.setWindowTitle(f"{APP_TITLE} — Guided Workflow")
        self.resize(1440, 900)
        self.setMinimumSize(1200, 760)
        self.setStyleSheet(shell_stylesheet())
        apply_ppigfinder_branding(self)

        central = QWidget()
        self.setCentralWidget(central)

        root = QHBoxLayout(central)
        root.setContentsMargins(18, 18, 18, 18)
        root.setSpacing(18)

        self.navigation = QListWidget()
        self.navigation.setMaximumWidth(260)
        root.addWidget(self.navigation)

        main_area = QVBoxLayout()
        root.addLayout(main_area, 1)

        self.breadcrumb = QLabel("Workspace > Data / Project")
        self.breadcrumb.setObjectName("SectionTitle")
        main_area.addWidget(self.breadcrumb)

        self.pages = QStackedWidget()
        main_area.addWidget(self.pages, 1)

        self._build_pages()
        self.navigation.currentRowChanged.connect(self._on_route_changed)

        # Data-first behavior.
        self.show_route("data")
        self.statusBar().showMessage("Start by loading genome data, opening a project, or importing a snapshot.")

        QTimer.singleShot(450, self._maybe_prompt_for_initial_data)

    # --------------------------------------------------------
    # State helpers
    # --------------------------------------------------------

    def _has_any_input(self) -> bool:
        return bool(
            self.workflow_state.get("genome_file")
            or self.workflow_state.get("project_file")
            or self.workflow_state.get("snapshot_file")
        )

    def _show_message(self, title: str, message: str) -> None:
        QMessageBox.information(self, title, message)

    def _remember_window(self, dialog) -> None:
        if dialog is not None:
            self._floating_windows.append(dialog)

    def _safe_import(self, module_name: str, object_name: str):
        module = __import__(module_name, fromlist=[object_name])
        return getattr(module, object_name)

    def _set_current_route(self, route_id: str) -> None:
        self.workflow_state.set_current_route(route_id)
        self._refresh_pages()

    def _refresh_pages(self) -> None:
        self._update_navigation_status()

        for page in self.module_pages_by_id.values():
            page.update_state(self.workflow_state)

    def _update_navigation_status(self) -> None:
        completed = self.workflow_state.completed_steps()
        genome_loaded = bool(self.workflow_state.get("genome_file"))
        orfs_ready = bool(self.workflow_state.get("guided_orf_count"))
        annotation_selected = bool(
            self.workflow_state.get("guided_blast_planned")
            or self.workflow_state.get("guided_hmm_planned")
            or self.workflow_state.get("guided_neighborhood_planned")
        )
        af3_ready = bool(self.workflow_state.get("af3_pair_count") or self.workflow_state.get("af3_json_path"))
        hpc_status = self.workflow_state.get("hpc_status")

        status_by_route = {
            "overview": "Ready",
            "data": "Loaded" if self._has_any_input() else "Start here",
            "genome": "Ready" if genome_loaded else "Waiting input",
            "orfs": "Completed" if orfs_ready else ("Ready" if genome_loaded else "Waiting genome"),
            "annotation": "Selected" if annotation_selected else ("Ready" if orfs_ready else "Waiting ORFs"),
            "alphafold": "Ready" if orfs_ready else "Waiting candidates",
            "hpc": str(hpc_status) if hpc_status else "Optional",
            "reports": "Ready" if self._has_any_input() else "Waiting outputs",
        }

        for route in self.routes:
            item = self.navigation_items_by_id.get(route.id)
            if item is None:
                continue

            status = status_by_route.get(route.id, route.status)
            item.setText(f"{route.title}\n{status}")

    # --------------------------------------------------------
    # Initial data behavior
    # --------------------------------------------------------

    def _maybe_prompt_for_initial_data(self) -> None:
        if self._auto_data_prompt_done:
            return

        if self.workflow_state.current_route != "data":
            return

        if self._has_any_input():
            return

        self._auto_data_prompt_done = True
        self._open_initial_data_dialog()

    def _detect_input_key_from_path(self, path: str) -> str:
        """
        Detect input role from file extension/name.
        """
        lower = path.lower()

        if lower.endswith((".gb", ".gbk", ".genbank", ".dna", ".fna", ".ffn")):
            return "genome_file"

        if lower.endswith((".hmm",)):
            return "hmm_profile_file"

        if lower.endswith((".ppigfinder.json",)):
            return "project_file"

        if lower.endswith((".json",)):
            # JSON can be project/snapshot/AF3, but at this stage a file JSON is
            # treated as a project/snapshot input. AF3 results are normally folders.
            return "snapshot_file"

        if lower.endswith((".faa", ".pep", ".protein", ".prot")):
            return "protein_query_file"

        if lower.endswith((".fa", ".fasta")):
            # Ambiguous. Prefer genome for complete workflow; protein FASTA can
            # also be loaded through the same selector if named .faa/.pep.
            return "genome_file"

        return "generic_input_file"

    def _select_any_input_file(self) -> bool:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Add input data",
            "",
            (
                "Supported inputs (*.fa *.fasta *.fna *.ffn *.gb *.gbk *.genbank *.dna "
                "*.faa *.pep *.protein *.prot *.hmm *.json *.ppigfinder.json);;"
                "Genome files (*.fa *.fasta *.fna *.ffn *.gb *.gbk *.genbank *.dna);;"
                "Protein FASTA (*.faa *.pep *.protein *.prot *.fa *.fasta);;"
                "HMM profiles (*.hmm);;"
                "Project/Snapshot JSON (*.json *.ppigfinder.json);;"
                "All files (*)"
            ),
        )

        if not path:
            return False

        key = self._detect_input_key_from_path(path)
        return self._register_selected_file(key, path)

    def _register_selected_file(self, key: str, path: str) -> bool:
        self.workflow_state.set_input(key, path)
        self.workflow_state.add_event(self.workflow_state.current_route, "select_file", path)

        if key == "genome_file":
            self._load_genome_file(path)
            return True

        if key == "project_file":
            self.workflow_state.set_flag("project_loaded", True)
            self.workflow_state.add_event("data", "open_project", path)
            self.statusBar().showMessage("Project file selected. Project restore service will be connected progressively.", 10000)
            self._refresh_pages()
            self.show_route("reports")
            return True

        if key == "snapshot_file":
            self.workflow_state.set_flag("snapshot_loaded", True)
            self.workflow_state.add_event("data", "import_snapshot", path)
            self.statusBar().showMessage("Snapshot/JSON input selected. Snapshot import service will be connected progressively.", 10000)
            self._refresh_pages()
            self.show_route("reports")
            return True

        if key == "protein_query_file":
            self.workflow_state.set_flag("protein_query_loaded", True)
            self.workflow_state.add_event("data", "load_protein_query", path)
            self.statusBar().showMessage("Protein query loaded. It will be available for BLAST after ORF prediction.", 10000)
            self._refresh_pages()
            self.show_route("annotation")
            return True

        if key == "hmm_profile_file":
            self.workflow_state.set_flag("hmm_profile_loaded", True)
            self.workflow_state.add_event("data", "load_hmm_profile", path)
            self.statusBar().showMessage("HMM profile loaded. It will be available for domain annotation after ORF prediction.", 10000)
            self._refresh_pages()
            self.show_route("annotation")
            return True

        self.workflow_state.set_flag("generic_input_loaded", True)
        self.statusBar().showMessage("Input file selected, but no specialized handler was detected yet.", 10000)
        self._refresh_pages()
        return True

    def _open_initial_data_dialog(self) -> None:
        choice = choose_initial_data(parent=self)

        if choice == "file":
            self._select_any_input_file()
        elif choice == "af3_results":
            self._select_folder(
                "af3_results_folder",
                "Import AF3 results folder",
            )

    # --------------------------------------------------------
    # File/folder loading
    # --------------------------------------------------------

    def _select_file(self, key: str, title: str, file_filter: str) -> bool:
        path, _ = QFileDialog.getOpenFileName(self, title, "", file_filter)

        if not path:
            self._refresh_pages()
            return False

        # If caller passes a generic key, auto-detect. Otherwise preserve the
        # explicit key.
        if key in {"generic_input_file", "auto"}:
            key = self._detect_input_key_from_path(path)

        return self._register_selected_file(key, path)

    def _load_genome_file(self, path: str) -> None:
        try:
            validate_genome_input = self._safe_import(
                "ppigfinder.ui_shell.input_validation",
                "validate_genome_input",
            )
            summary = validate_genome_input(path)

            self.workflow_state.set_metric("sequence_count", getattr(summary, "sequence_count", None))
            self.workflow_state.set_metric("genome_sequence_count", getattr(summary, "sequence_count", None))
            self.workflow_state.set_metric("total_length", getattr(summary, "total_length", None))
            self.workflow_state.set_metric("genome_total_length", getattr(summary, "total_length", None))
            self.workflow_state.set_metric("longest_length", getattr(summary, "longest_length", None))
            self.workflow_state.set_metric("gc_percent", getattr(summary, "gc_percent", None))
            self.workflow_state.set_metric("genome_gc_percent", getattr(summary, "gc_percent", None))
            self.workflow_state.set_input("genome_file_type", getattr(summary, "file_type", "Unknown"))
            self.workflow_state.set_input("genome_name", getattr(summary, "name", "Genome"))
            self.workflow_state.set_flag("genome_valid", getattr(summary, "valid", False))
            self.workflow_state.add_event("data", "load_genome", path)

            try:
                show_genome_inspector = self._safe_import(
                    "ppigfinder.ui_shell.genome_inspector_dialog",
                    "show_genome_inspector",
                )
                dialog = show_genome_inspector(path, parent=self)
                self._remember_window(dialog)
            except Exception as exc:
                self.statusBar().showMessage(f"Genome inspector unavailable: {exc}", 10000)

            self._refresh_pages()

            if getattr(summary, "valid", False):
                self.statusBar().showMessage(
                    "Genome loaded and validated. ORF prediction is now available.",
                    12000,
                )
                self.show_route("orfs")
            else:
                self.statusBar().showMessage(
                    "Genome selected, but validation detected a problem.",
                    12000,
                )
                self.show_route("data")

        except Exception as exc:
            self._show_message("Genome loading", f"Could not validate genome file:\n\n{exc}")
            self._refresh_pages()

    def _select_folder(self, key: str, title: str) -> bool:
        path = QFileDialog.getExistingDirectory(self, title, "")

        if not path:
            return False

        self.workflow_state.set_input(key, path)
        self.workflow_state.add_event(self.workflow_state.current_route, "select_folder", path)
        self._refresh_pages()
        return True

    # --------------------------------------------------------
    # Workflow operations
    # --------------------------------------------------------

    def _predict_guided_orfs(self) -> None:
        genome_file = self.workflow_state.get("genome_file")

        if not genome_file:
            self._show_message("Predict ORFs", "Load a genome file first in Data / Project.")
            self.show_route("data")
            QTimer.singleShot(300, self._open_initial_data_dialog)
            return

        try:
            predict_orfs_from_file = self._safe_import(
                "ppigfinder.ui_shell.guided_backend",
                "predict_orfs_from_file",
            )
            summary = predict_orfs_from_file(genome_file, min_aa=30)
            self.guided_orfs = list(getattr(summary, "orfs", []))

            self.workflow_state.set_metric("guided_orf_count", getattr(summary, "orf_count", len(self.guided_orfs)))
            self.workflow_state.set_metric("guided_longest_orf_aa", getattr(summary, "longest_orf_aa", None))
            self.workflow_state.set_metric("guided_shortest_orf_aa", getattr(summary, "shortest_orf_aa", None))
            self.workflow_state.set_metric("guided_orf_min_aa", getattr(summary, "min_aa", 30))

            lengths = [len(getattr(orf, "protein_sequence", "") or "") for orf in self.guided_orfs]
            plus_count = sum(1 for orf in self.guided_orfs if getattr(orf, "strand", "") == "+")
            minus_count = sum(1 for orf in self.guided_orfs if getattr(orf, "strand", "") == "-")

            self.workflow_state.set_metric("guided_orf_mean_aa", round(sum(lengths) / len(lengths), 2) if lengths else 0)
            self.workflow_state.set_metric("guided_orf_plus_count", plus_count)
            self.workflow_state.set_metric("guided_orf_minus_count", minus_count)
            self.workflow_state.add_event("orfs", "predict_orfs", f"{len(self.guided_orfs)} ORFs")

            try:
                GuidedORFResultsDialog = self._safe_import(
                    "ppigfinder.ui_shell.orf_results_dialog",
                    "GuidedORFResultsDialog",
                )
                dialog = GuidedORFResultsDialog(self.guided_orfs, parent=self)
                dialog.show()
                self._remember_window(dialog)
            except Exception as exc:
                self.statusBar().showMessage(f"ORF table unavailable: {exc}", 10000)

            self._refresh_pages()
            self.statusBar().showMessage(
                f"ORF prediction completed: {len(self.guided_orfs)} ORFs. Moving to Annotation.",
                12000,
            )
            self.show_route("annotation")

        except Exception as exc:
            self._show_message("Predict ORFs", f"Could not run guided ORF prediction:\n\n{exc}")

    def _show_guided_orfs(self) -> None:
        if not self.guided_orfs:
            self._show_message("ORFs", "No guided ORFs available yet.")
            return

        try:
            GuidedORFResultsDialog = self._safe_import(
                "ppigfinder.ui_shell.orf_results_dialog",
                "GuidedORFResultsDialog",
            )
            dialog = GuidedORFResultsDialog(self.guided_orfs, parent=self)
            dialog.show()
            self._remember_window(dialog)
        except Exception as exc:
            self._show_message("ORF table", f"Could not open ORF viewer:\n\n{exc}")

    def _export_guided_orfs_fasta(self) -> None:
        if not self.guided_orfs:
            self._show_message("Export ORF FASTA", "No guided ORFs available yet.")
            return

        path, _ = QFileDialog.getSaveFileName(
            self,
            "Export guided ORFs",
            "guided_orfs.faa",
            "Protein FASTA (*.faa *.fasta *.fa);;All files (*)",
        )

        if not path:
            return

        try:
            write_orfs_fasta = self._safe_import(
                "ppigfinder.ui_shell.guided_backend",
                "write_orfs_fasta",
            )
            write_orfs_fasta(path, self.guided_orfs)
            self.workflow_state.add_event("orfs", "export_orfs_fasta", path)
            self._refresh_pages()
            self.statusBar().showMessage(f"ORF FASTA exported: {path}", 10000)
        except Exception as exc:
            self._show_message("Export ORF FASTA", f"Could not export FASTA:\n\n{exc}")

    def _show_annotation_candidates(self) -> None:
        if not self.guided_orfs:
            self._show_message("Annotation", "No ORFs available yet. Run Predict ORFs first.")
            self.show_route("orfs")
            return

        try:
            AnnotationCandidatesDialog = self._safe_import(
                "ppigfinder.ui_shell.annotation_candidates_dialog",
                "AnnotationCandidatesDialog",
            )
            dialog = AnnotationCandidatesDialog(self.guided_orfs, parent=self)
            dialog.show()
            self._remember_window(dialog)
        except Exception as exc:
            self.statusBar().showMessage(f"Candidate table unavailable: {exc}", 10000)

        self.workflow_state.set_metric("guided_annotation_candidates_count", len(self.guided_orfs))
        self.workflow_state.add_event("annotation", "review_candidates", str(len(self.guided_orfs)))
        self._refresh_pages()

    def _mark_annotation_step(self, key: str, label: str) -> None:
        self.workflow_state.set_flag(key, True)
        self.workflow_state.add_event("annotation", key, label)
        self._refresh_pages()
        self.statusBar().showMessage(f"{label} selected in guided workflow.", 8000)

    def _open_af3_pair_builder(self) -> None:
        if not self.guided_orfs:
            self._show_message("AlphaFold / PPI", "No ORFs available yet. Run Predict ORFs first.")
            self.show_route("orfs")
            return

        try:
            open_af3_pair_builder_dialog = self._safe_import(
                "ppigfinder.ui_shell.af3_pair_builder_dialog",
                "open_af3_pair_builder_dialog",
            )
            result = open_af3_pair_builder_dialog(self.guided_orfs, parent=self)

            if isinstance(result, dict):
                self.workflow_state.set_metric("af3_pair_count", result.get("pair_count", 0))
                if result.get("json_path"):
                    self.workflow_state.set_input("af3_json_path", result["json_path"])

            self.workflow_state.add_event("alphafold", "build_af3_pairs", "")
            self._refresh_pages()

        except Exception as exc:
            self._show_message("AlphaFold / PPI", f"Could not open AF3 builder:\n\n{exc}")

    def _open_hpc_dialog(self) -> None:
        try:
            open_hpc_connection_dialog = self._safe_import(
                "ppigfinder.ui_shell.hpc_dialog",
                "open_hpc_connection_dialog",
            )
            status, config = open_hpc_connection_dialog(parent=self)

            if config is not None:
                self.workflow_state.set_input("hpc_profile", getattr(config, "profile", ""))
                self.workflow_state.set_input("hpc_host", getattr(config, "host", ""))
                self.workflow_state.set_input("hpc_user", getattr(config, "user", ""))
                self.workflow_state.set_input("hpc_port", getattr(config, "port", ""))

            if status is not None:
                self.workflow_state.set_input("hpc_status", "OK" if getattr(status, "connection_ok", False) else "Problem")
                self.workflow_state.set_input("hpc_mode", "Local cluster" if getattr(status, "running_on_cluster", False) else "SSH")
                self.workflow_state.set_flag("hpc_connected", getattr(status, "connection_ok", False))

            self.workflow_state.add_event("hpc", "configure_hpc", "")
            self._refresh_pages()

        except Exception as exc:
            self._show_message("DaVinci / HPC", f"Could not open HPC dialog:\n\n{exc}")

    def _export_guided_summary(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Export guided summary",
            "ppigfinder_guided_summary.md",
            "Markdown (*.md);;Text files (*.txt);;All files (*)",
        )

        if not path:
            return

        try:
            write_guided_summary = self._safe_import(
                "ppigfinder.ui_shell.guided_backend",
                "write_guided_summary",
            )
            payload = {}
            payload.update(self.workflow_state.loaded_inputs)
            payload.update(self.workflow_state.metrics)
            payload.update(self.workflow_state.flags)

            write_guided_summary(path, payload)
            self.workflow_state.set_flag("guided_summary_exported", True)
            self.workflow_state.add_event("reports", "export_guided_summary", path)
            self._refresh_pages()
            self.statusBar().showMessage(f"Guided summary exported: {path}", 10000)

        except Exception as exc:
            self._show_message("Reports", f"Could not export guided summary:\n\n{exc}")

    def _open_workflow_map(self) -> None:
        try:
            show_workflow_overview = self._safe_import(
                "ppigfinder.ui_shell.workflow_overview_dialog",
                "show_workflow_overview",
            )
            show_workflow_overview(parent=self)
        except Exception as exc:
            self._show_message("Workflow map", f"Could not open workflow map:\n\n{exc}")

    def _open_legacy_interface(self) -> None:
        if self.bridge:
            self.bridge.call("open_legacy_interface")
            return

        self._show_message("Legacy interface", "Bridge not available yet.")

    # --------------------------------------------------------
    # Actions
    # --------------------------------------------------------

    def _action(self, label: str, description: str, callback):
        lower = label.lower()
        button_text = "Run"

        if lower.startswith("add"):
            button_text = "Add"
        elif lower.startswith("open"):
            button_text = "Open"
        elif lower.startswith("import"):
            button_text = "Import"
        elif lower.startswith("export"):
            button_text = "Export"
        elif lower.startswith("predict"):
            button_text = "Predict"
        elif lower.startswith("annotate"):
            button_text = "Annotate"
        elif lower.startswith("continue"):
            button_text = "Continue"
        elif lower.startswith("show"):
            button_text = "Show"
        elif lower.startswith("configure"):
            button_text = "Configure"
        elif lower.startswith("build"):
            button_text = "Build"

        return {
            "label": label,
            "description": description,
            "button_text": button_text,
            "callback": callback,
        }

    def _actions_for_route(self, route_id: str) -> list[dict]:
        if route_id == "overview":
            return [
                self._action("Add input data", "Choose an input file or AF3 results folder and let ppigFinder route it into the workflow.", self._open_initial_data_dialog),
                self._action("Show workflow map", "Review the complete guided workflow and step dependencies.", self._open_workflow_map),
                self._action("Start with Data / Project", "Go to the input/status page.", lambda: self.show_route("data")),
            ]

        if route_id == "data":
            return [
                self._action("Add input data", "Choose an input file or AF3 results folder. ppigFinder will detect the data type and update the workflow.", self._open_initial_data_dialog),
            ]

        if route_id == "genome":
            return [
                self._action("Add input data", "Load or replace the current genome/project input.", self._open_initial_data_dialog),
                self._action("Continue to Protein / ORFs", "Go to ORF prediction.", lambda: self.show_route("orfs")),
            ]

        if route_id == "orfs":
            return [
                self._action("Predict ORFs", "Run guided ORF prediction and automatically advance to Annotation.", self._predict_guided_orfs),
                self._action("Review guided ORF table", "Inspect ORFs already generated in the guided workflow.", self._show_guided_orfs),
                self._action("Export ORF FASTA", "Export protein sequences for predicted ORFs.", self._export_guided_orfs_fasta),
                self._action("Continue to Annotation", "Go to annotation stage.", lambda: self.show_route("annotation")),
            ]

        if route_id == "annotation":
            return [
                self._action("Review candidate ORFs", "Inspect guided ORFs as candidates for annotation and PPI analysis.", self._show_annotation_candidates),
                self._action("Run BLAST", "Mark BLAST as part of the current annotation plan.", lambda: self._mark_annotation_step("guided_blast_planned", "BLAST")),
                self._action("Annotate HMM", "Mark HMM/domain annotation as part of the current plan.", lambda: self._mark_annotation_step("guided_hmm_planned", "HMM/domain annotation")),
                self._action("Neighbourhood analysis", "Mark genomic neighbourhood analysis as part of the current plan.", lambda: self._mark_annotation_step("guided_neighborhood_planned", "Neighbourhood analysis")),
                self._action("Continue to AlphaFold / PPI", "Proceed to structural candidate construction.", lambda: self.show_route("alphafold")),
            ]

        if route_id == "alphafold":
            return [
                self._action("Build AF3 candidate pairs", "Build candidate pairs and optionally export AlphaFold Server JSON.", self._open_af3_pair_builder),
                self._action("Import AF3 results", "Select a results folder from AlphaFold / AF3 runs.", lambda: self._select_folder("af3_results_folder", "Import AF3 results folder")),
                self._action("Continue to DaVinci / HPC", "Proceed to optional execution on DaVinci/HPC.", lambda: self.show_route("hpc")),
            ]

        if route_id == "hpc":
            return [
                self._action("Configure DaVinci / HPC", "Open the DaVinci/HPC connection dialog.", self._open_hpc_dialog),
                self._action("Continue to Reports", "Proceed to summary and export steps.", lambda: self.show_route("reports")),
            ]

        if route_id == "reports":
            return [
                self._action("Export guided summary", "Export the current guided workflow state to Markdown.", self._export_guided_summary),
                self._action("Open full current interface", "Advanced option: open the complete current ppigFinder interface.", self._open_legacy_interface),
            ]

        return []

    def _build_pages(self):
        for index, route in enumerate(self.routes):
            self.route_index_by_id[route.id] = index

            item = QListWidgetItem(f"{route.title}\n{route.status}")
            self.navigation.addItem(item)
            self.navigation_items_by_id[route.id] = item

            page = ModulePage(route, actions=self._actions_for_route(route.id))
            self.module_pages_by_id[route.id] = page
            self.pages.addWidget(page)

    def show_route(self, route_id: str) -> bool:
        index = self.route_index_by_id.get(route_id)

        if index is None:
            return False

        self.navigation.setCurrentRow(index)
        return True

    def _on_route_changed(self, index: int):
        if index < 0 or index >= len(self.routes):
            return

        route = self.routes[index]
        self.breadcrumb.setText(f"Workspace > {route.title}")
        self.pages.setCurrentIndex(index)
        self._set_current_route(route.id)

        if route.id == "data":
            QTimer.singleShot(350, self._maybe_prompt_for_initial_data)
