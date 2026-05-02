#!/usr/bin/env python3
from __future__ import annotations

try:
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import (
        QWidget,
        QVBoxLayout,
        QHBoxLayout,
        QLabel,
        QFrame,
        QGridLayout,
        QScrollArea,
        QPushButton,
        QTableWidget,
        QTableWidgetItem,
        QHeaderView,
        QLineEdit,
        QComboBox,
    )
    QT6 = True
except Exception:
    from PyQt5.QtCore import Qt
    from PyQt5.QtWidgets import (
        QWidget,
        QVBoxLayout,
        QHBoxLayout,
        QLabel,
        QFrame,
        QGridLayout,
        QScrollArea,
        QPushButton,
        QTableWidget,
        QTableWidgetItem,
        QHeaderView,
        QLineEdit,
        QComboBox,
    )
    QT6 = False

from ppigfinder.ui_shell.workflow_state import WORKFLOW_ORDER, WorkflowState


def _set_alignment(widget, alignment):
    try:
        widget.setAlignment(alignment)
    except Exception:
        pass


class InfoCard(QFrame):
    def __init__(self, title: str = "", value: str = "", subtitle: str = "", parent=None):
        super().__init__(parent)
        self.setObjectName("InfoCard")
        self.setFrameShape(QFrame.Shape.StyledPanel if QT6 else QFrame.StyledPanel)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(14, 12, 14, 12)
        layout.setSpacing(6)

        self.title_label = QLabel(title)
        self.title_label.setObjectName("InfoCardTitle")

        self.value_label = QLabel(value)
        self.value_label.setObjectName("InfoCardValue")
        self.value_label.setWordWrap(True)

        self.subtitle_label = QLabel(subtitle)
        self.subtitle_label.setObjectName("InfoCardSubtitle")
        self.subtitle_label.setWordWrap(True)

        layout.addWidget(self.title_label)
        layout.addWidget(self.value_label)
        layout.addWidget(self.subtitle_label)

    def set_content(self, title: str, value: str, subtitle: str = "") -> None:
        self.title_label.setText(title)
        self.value_label.setText(value)
        self.subtitle_label.setText(subtitle)


class FlowRibbonWidget(QWidget):
    """
    Clickable workflow ribbon.

    Each chip behaves like a navigation control. This allows the user to click
    Data, Genome, ORFs, Annotation, AlphaFold, HPC or Reports and jump directly
    to the corresponding page.
    """

    STEP_LABELS = {
        "data": "Data",
        "genome": "Genome",
        "orfs": "ORFs",
        "annotation": "Annotation",
        "alphafold": "AlphaFold",
        "hpc": "HPC",
        "reports": "Reports",
    }

    def __init__(self, parent=None, on_route_selected=None):
        super().__init__(parent)
        self._buttons = {}
        self._on_route_selected = on_route_selected

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)

        for step in WORKFLOW_ORDER[1:]:
            button = QPushButton(self.STEP_LABELS.get(step, step.capitalize()))
            button.setFlat(True)
            button.setCursor(Qt.CursorShape.PointingHandCursor if QT6 else Qt.PointingHandCursor)
            button.clicked.connect(lambda _checked=False, route_id=step: self._route_clicked(route_id))
            button.setMinimumHeight(34)

            layout.addWidget(button)
            self._buttons[step] = button

            if step != WORKFLOW_ORDER[-1]:
                arrow = QLabel("→")
                arrow.setObjectName("FlowArrow")
                _set_alignment(arrow, Qt.AlignmentFlag.AlignCenter if QT6 else Qt.AlignCenter)
                layout.addWidget(arrow)

        layout.addStretch(1)

    def _route_clicked(self, route_id: str) -> None:
        if callable(self._on_route_selected):
            self._on_route_selected(route_id)

    def set_route_callback(self, callback) -> None:
        self._on_route_selected = callback

    def _build_annotation_candidate_panel(self) -> None:
        """
        Embedded candidate table for the Annotation page.

        This avoids opening a separate dialog for the normal workflow. A separate
        full-screen table can still be used later for deep inspection.
        """
        title = QLabel("Candidate ORFs")
        title.setObjectName("SectionSubTitle")
        self.layout_main.addWidget(title)

        controls = QHBoxLayout()

        self.annotation_search = QLineEdit()
        self.annotation_search.setPlaceholderText("Filter by ORF ID")
        self.annotation_search.textChanged.connect(self._refresh_annotation_candidate_table)
        controls.addWidget(self.annotation_search, 2)

        self.annotation_strand = QComboBox()
        self.annotation_strand.addItems(["All strands", "+", "-"])
        self.annotation_strand.currentIndexChanged.connect(self._refresh_annotation_candidate_table)
        controls.addWidget(self.annotation_strand, 1)

        self.layout_main.addLayout(controls)

        self.annotation_candidate_table = QTableWidget(0, 8)
        self.annotation_candidate_table.setHorizontalHeaderLabels(
            [
                "Candidate",
                "ORF ID",
                "Start",
                "End",
                "Strand",
                "Frame",
                "AA length",
                "Status",
            ]
        )
        self.annotation_candidate_table.setAlternatingRowColors(True)

        try:
            self.annotation_candidate_table.verticalHeader().setVisible(False)
            header = self.annotation_candidate_table.horizontalHeader()
            for col in range(7):
                header.setSectionResizeMode(
                    col,
                    QHeaderView.ResizeMode.ResizeToContents if QT6 else QHeaderView.ResizeToContents,
                )
            header.setSectionResizeMode(
                7,
                QHeaderView.ResizeMode.Stretch if QT6 else QHeaderView.Stretch,
            )
        except Exception:
            pass

        self.layout_main.addWidget(self.annotation_candidate_table, 2)
        self._refresh_annotation_candidate_table()

    def set_annotation_candidates(self, orfs) -> None:
        self.annotation_candidates = list(orfs or [])
        self._refresh_annotation_candidate_table()

    def _refresh_annotation_candidate_table(self) -> None:
        if self.annotation_candidate_table is None:
            return

        query = ""
        strand_filter = "All strands"

        if self.annotation_search is not None:
            query = self.annotation_search.text().strip().upper()

        if self.annotation_strand is not None:
            strand_filter = self.annotation_strand.currentText()

        filtered = []

        for orf in self.annotation_candidates:
            oid = _candidate_orf_id(orf).upper()
            strand = _candidate_strand(orf)

            if query and query not in oid:
                continue

            if strand_filter in {"+", "-"} and strand != strand_filter:
                continue

            filtered.append(orf)

        max_rows = 1200
        visible = filtered[:max_rows]

        self.annotation_candidate_table.setRowCount(len(visible))

        for row, orf in enumerate(visible, start=0):
            values = [
                str(row + 1),
                _candidate_orf_id(orf, f"orf_{row + 1}"),
                str(_candidate_start(orf)),
                str(_candidate_end(orf)),
                _candidate_strand(orf),
                _candidate_frame(orf),
                str(_candidate_aa_length(orf)),
                "Pending BLAST/HMM/neighbourhood",
            ]

            for col, value in enumerate(values):
                self.annotation_candidate_table.setItem(row, col, QTableWidgetItem(value))

        self.annotation_candidate_table.resizeRowsToContents()

        if self.footer is not None and self.route_id == "annotation":
            if not self.annotation_candidates:
                self.footer.setText(
                    "No candidate ORFs embedded yet. Run Predict ORFs, then use Review candidate ORFs."
                )
            else:
                self.footer.setText(
                    f"Embedded candidate table: showing {len(visible):,} of {len(filtered):,} filtered "
                    f"ORFs from {len(self.annotation_candidates):,} total candidates. "
                    "Use the filters above to narrow the table."
                )

    def update_state(self, state: WorkflowState) -> None:
        completed = state.completed_steps()
        current = state.current_route

        for step, button in self._buttons.items():
            if step == current:
                button.setStyleSheet(
                    "QPushButton {background:#17384d;color:white;border-radius:10px;padding:6px 12px;font-weight:700;}"
                    "QPushButton:hover {background:#1f4a63;}"
                )
            elif step in completed:
                button.setStyleSheet(
                    "QPushButton {background:#cfe3ef;color:#17384d;border-radius:10px;padding:6px 12px;font-weight:600;}"
                    "QPushButton:hover {background:#b8d8e8;}"
                )
            else:
                button.setStyleSheet(
                    "QPushButton {background:#edf3f7;color:#60717f;border-radius:10px;padding:6px 12px;}"
                    "QPushButton:hover {background:#ddeaf0;color:#17384d;}"
                )




def _candidate_orf_id(orf, fallback=""):
    return str(getattr(orf, "id", fallback) or fallback)


def _candidate_start(orf):
    try:
        return int(getattr(orf, "start", 0))
    except Exception:
        return 0


def _candidate_end(orf):
    try:
        return int(getattr(orf, "end", 0))
    except Exception:
        return 0


def _candidate_strand(orf):
    return str(getattr(orf, "strand", "?") or "?")


def _candidate_frame(orf):
    return str(getattr(orf, "frame", "?") or "?")


def _candidate_aa_length(orf):
    if hasattr(orf, "aa_length"):
        try:
            return int(getattr(orf, "aa_length"))
        except Exception:
            pass
    return len(getattr(orf, "protein_sequence", "") or "")

class ModuleVisualizationPanel(QScrollArea):
    def __init__(self, route_id: str, parent=None):
        super().__init__(parent)
        self.route_id = route_id
        self.setWidgetResizable(True)

        self.body = QWidget()
        self.setWidget(self.body)

        self.layout_main = QVBoxLayout(self.body)
        self.layout_main.setContentsMargins(10, 10, 10, 10)
        self.layout_main.setSpacing(10)

        self.title = QLabel("Module status")
        self.title.setObjectName("SectionSubTitle")
        self.layout_main.addWidget(self.title)

        self.grid = QGridLayout()
        self.grid.setHorizontalSpacing(10)
        self.grid.setVerticalSpacing(10)
        self.layout_main.addLayout(self.grid)

        self.card1 = InfoCard()
        self.card2 = InfoCard()
        self.card3 = InfoCard()
        self.card4 = InfoCard()

        self.grid.addWidget(self.card1, 0, 0)
        self.grid.addWidget(self.card2, 0, 1)
        self.grid.addWidget(self.card3, 1, 0)
        self.grid.addWidget(self.card4, 1, 1)

        self.footer = QLabel("")
        self.footer.setWordWrap(True)
        self.footer.setObjectName("InfoFooter")
        self.layout_main.addWidget(self.footer)

        self.annotation_candidates = []
        self.annotation_candidate_table = None
        self.annotation_search = None
        self.annotation_strand = None

        if self.route_id == "annotation":
            self._build_annotation_candidate_panel()

        self.layout_main.addStretch(1)

    def update_state(self, state: WorkflowState) -> None:
        route = self.route_id

        if route == "data":
            genome_file = state.get("genome_file", "No genome file loaded")
            protein_file = state.get("protein_query_file", "No protein query loaded")
            hmm_file = state.get("hmm_profile_file", "No HMM profile loaded")

            self.card1.set_content("Genome file", str(genome_file), "Main nucleotide input")
            self.card2.set_content("Protein query", str(protein_file), "BLAST query input")
            self.card3.set_content("HMM profiles", str(hmm_file), "Domain annotation input")
            self.card4.set_content(
                "Next step",
                state.next_recommended_step().capitalize(),
                "Recommended route after data loading",
            )
            self.footer.setText(
                "This panel should evolve into a richer data-entry preview with file validation, "
                "organism metadata and input summaries."
            )
            return

        if route == "genome":
            self.card1.set_content(
                "Genome length",
                str(state.get("total_length", state.get("genome_total_length", "N/A"))),
                "Loaded nucleotide length",
            )
            self.card2.set_content(
                "GC%",
                str(state.get("gc_percent", "N/A")),
                "Genome composition",
            )
            self.card3.set_content(
                "Sequence count",
                str(state.get("sequence_count", "N/A")),
                "Number of sequences/contigs",
            )
            self.card4.set_content(
                "Next step",
                "Protein / ORFs",
                "Proceed to ORF prediction and protein generation",
            )
            self.footer.setText(
                "Future visual direction: genome strip preview, contig overview and coordinate-aware "
                "mini-map inspired by bacterial genome browsers."
            )
            return

        if route == "orfs":
            self.card1.set_content(
                "Predicted ORFs",
                str(state.get("guided_orf_count", 0)),
                "ORFs currently generated in guided flow",
            )
            self.card2.set_content(
                "Mean / longest",
                f"{state.get('guided_orf_mean_aa', 'N/A')} / {state.get('guided_longest_orf_aa', 'N/A')} aa",
                "Average and maximum protein length",
            )
            self.card3.set_content(
                "Strands",
                f"+ {state.get('guided_orf_plus_count', 0)} / - {state.get('guided_orf_minus_count', 0)}",
                "Predicted ORF orientation",
            )
            self.card4.set_content(
                "Next step",
                "Annotation",
                "Run BLAST/HMM and contextual filtering",
            )
            self.footer.setText(
                "ORF Discovery generates the protein set used by BLAST, HMM/domain annotation, "
                "neighbourhood analysis and AlphaFold/PPI candidate construction."
            )
            return

        if route == "annotation":
            candidate_count = state.get("guided_orf_count", 0)
            blast = "Selected" if state.get("guided_blast_planned") else ("Query loaded" if state.get("protein_query_file") else "Pending")
            hmm = "Selected" if state.get("guided_hmm_planned") else ("Profiles loaded" if state.get("hmm_profile_file") else "Pending")
            neigh = "Selected" if state.get("guided_neighborhood_planned") else "Pending"

            self.card1.set_content("BLAST", blast, "Protein similarity search")
            self.card2.set_content("HMM domains", hmm, "Conserved domain profile scanning")
            self.card3.set_content("Neighbourhood", neigh, "Candidate selection by genomic context")
            self.card4.set_content(
                "Next step",
                f"{candidate_count} candidates for AF3/PPI",
                "Prepare candidate pairs",
            )
            self.footer.setText(
                "Functional annotation reduces the search space and selects biologically plausible "
                "candidates for structural interaction analysis."
            )
            return

        if route == "alphafold":
            self.card1.set_content(
                "AF3 pairs",
                str(state.get("af3_pair_count", 0)),
                "Candidate pairs created from ORFs",
            )
            self.card2.set_content(
                "AF3 JSON",
                str(state.get("af3_json_path", "Not exported")),
                "Server-compatible AlphaFold JSON",
            )
            self.card3.set_content(
                "Results folder",
                str(state.get("af3_results_folder", "Not imported")),
                "Imported AF3 result directory",
            )
            self.card4.set_content(
                "Next step",
                "DaVinci / HPC or Reports",
                "Depending on execution mode",
            )
            self.footer.setText(
                "Future visual direction: pair network preview, ipTM/PAE overview and miniature "
                "interaction confidence dashboard."
            )
            return

        if route == "hpc":
            self.card1.set_content(
                "Profile",
                str(state.get("hpc_profile", "Not configured")),
                "Selected HPC/server profile",
            )
            self.card2.set_content(
                "Host",
                str(state.get("hpc_host", "N/A")),
                "Target server/cluster",
            )
            self.card3.set_content(
                "Status",
                str(state.get("hpc_status", "Not tested")),
                "Connection or execution state",
            )
            self.card4.set_content(
                "Execution mode",
                str(state.get("hpc_mode", "Not defined")),
                "Local, SSH or cluster mode",
            )
            self.footer.setText(
                "Future visual direction: queue cards, job submission progress and a more graphical "
                "DaVinci cluster connector."
            )
            return

        if route == "reports":
            self.card1.set_content(
                "Guided summary",
                "Exported" if state.get("guided_summary_exported") else "Pending",
                "Markdown summary of the guided workflow",
            )
            self.card2.set_content(
                "Project snapshot",
                "Available later",
                "Guided shell will be synchronized with project snapshots",
            )
            self.card3.set_content(
                "HTML report",
                "Pending",
                "HTML reporting integration",
            )
            self.card4.set_content(
                "Workflow next step",
                "Review / export",
                "Finalize and share results",
            )
            self.footer.setText(
                "Future visual direction: report templates, figure previews, export checklist and "
                "provenance tracking."
            )
            return

        # overview
        completed = state.completed_steps()

        self.card1.set_content("Current step", state.current_route.capitalize(), "Current guided route")
        self.card2.set_content("Completed", ", ".join(sorted(completed)) if completed else "None", "Workflow steps with data")
        self.card3.set_content("Next", state.next_recommended_step().capitalize(), "Recommended next step")
        self.card4.set_content("Events", str(len(state.events)), "Actions recorded in this session")
        self.footer.setText(
            "This overview should evolve into a fully visual workflow dashboard with process tracking, "
            "data dependencies and direct access to generated figures and tables."
        )
