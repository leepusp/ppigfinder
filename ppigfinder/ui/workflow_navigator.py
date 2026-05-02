#!/usr/bin/env python3
"""
Workflow navigator for ppigFinder.

This module adds a data-type-oriented navigation panel to the legacy GUI
without modifying legacy_v20.py.
"""

from __future__ import annotations

from dataclasses import dataclass

try:
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import (
        QDockWidget,
        QWidget,
        QVBoxLayout,
        QLabel,
        QTreeWidget,
        QTreeWidgetItem,
    )
    QT6 = True
except Exception:
    from PyQt5.QtCore import Qt
    from PyQt5.QtWidgets import (
        QDockWidget,
        QWidget,
        QVBoxLayout,
        QLabel,
        QTreeWidget,
        QTreeWidgetItem,
    )
    QT6 = False


@dataclass(slots=True)
class WorkflowItem:
    group: str
    label: str
    data_type: str
    description: str
    output: str = ""
    method_name: str | None = None
    tab_hint: str | None = None


WORKFLOW_ITEMS = [
    WorkflowItem(
        group="Project / Data",
        label="Open genome",
        data_type="FASTA, GenBank or SnapGene DNA",
        description="Load a genome or nucleotide sequence into ppigFinder.",
        output="Genome sequence and metadata",
        method_name="load_fasta",
        tab_hint="Genome",
    ),
    WorkflowItem(
        group="Project / Data",
        label="Save project",
        data_type="Current ppigFinder session",
        description="Save the current project using the legacy project format.",
        output="Project file",
        method_name="save_project",
    ),
    WorkflowItem(
        group="Project / Data",
        label="Export Project Snapshot v3",
        data_type="Genome, ORFs, BLAST, HMM, AF3 and UI state",
        description="Export a versioned JSON snapshot for reproducible reporting and future workflows.",
        output="Versioned .ppigfinder.json snapshot",
        method_name="export_project_snapshot",
    ),

    WorkflowItem(
        group="DNA / Genome",
        label="Genome overview",
        data_type="DNA/genome sequence",
        description="Show genome length, GC content, ORF count and backend availability.",
        output="Genome summary",
        tab_hint="Genome",
    ),
    WorkflowItem(
        group="DNA / Genome",
        label="Translate genome",
        data_type="DNA sequence",
        description="Translate genomic sequence or selected regions into protein sequence.",
        output="Protein translation",
        method_name="translate_genome",
        tab_hint="DNA",
    ),
    WorkflowItem(
        group="DNA / Genome",
        label="Export genome map",
        data_type="Genome coordinates and ORF annotations",
        description="Export the current genome/ORF map as a PDF figure.",
        output="PDF genome map",
        method_name="export_map_pdf",
    ),

    WorkflowItem(
        group="Protein / ORFs",
        label="Predict ORFs",
        data_type="DNA/genome sequence",
        description="Predict protein-coding open reading frames using Pyrodigal, six-frame scanning or hybrid mode.",
        output="ORF table and protein sequences",
        method_name="analyze_orfs",
    ),
    WorkflowItem(
        group="Protein / ORFs",
        label="ORF table",
        data_type="Predicted ORFs",
        description="Inspect predicted ORFs, genomic coordinates, strand, frame, size, GC and annotation columns.",
        output="ORF table",
        tab_hint="Genome",
    ),
    WorkflowItem(
        group="Protein / ORFs",
        label="Export ORF FASTA",
        data_type="Predicted protein sequences",
        description="Export predicted ORF protein sequences in FASTA format for downstream analyses.",
        output="Protein FASTA",
        method_name="save_fasta",
    ),

    WorkflowItem(
        group="Annotation",
        label="BLAST query",
        data_type="Protein sequence vs predicted ORFs",
        description="Search a protein query against predicted ORFs using BLAST+, k-mer fallback or Smith-Waterman fallback.",
        output="Similarity hits mapped to ORFs",
        method_name="run_blast",
        tab_hint="BLAST Query",
    ),
    WorkflowItem(
        group="Annotation",
        label="BLAST results",
        data_type="BLAST/K-mer/Smith-Waterman hits",
        description="Review protein similarity hits, identity, score, e-value and candidate ORFs.",
        output="BLAST result table",
        tab_hint="BLAST Results",
    ),
    WorkflowItem(
        group="Annotation",
        label="HMM domains",
        data_type="Predicted protein sequences",
        description="Annotate conserved domains using HMMER3 or built-in domain scanner fallback.",
        output="Domain annotations per ORF",
        method_name="annotate_hmm",
        tab_hint="HMM",
    ),
    WorkflowItem(
        group="Annotation",
        label="Neighborhood",
        data_type="Genomic context around ORFs",
        description="Inspect ORF genomic neighborhoods to support operon/context interpretation.",
        output="Neighborhood table/map",
        tab_hint="Neighborhood",
    ),

    WorkflowItem(
        group="AlphaFold / PPI",
        label="Export AF3 Server JSON",
        data_type="Selected ORF protein pairs",
        description="Generate AlphaFold Server/web JSON for pairwise protein-protein interaction predictions.",
        output="AlphaFold Server JSON",
        method_name="export_af3_server_json",
        tab_hint="Submit AF3 via Server",
    ),
    WorkflowItem(
        group="AlphaFold / PPI",
        label="Import AF3 results",
        data_type="AlphaFold 3 result folders",
        description="Parse AF3 output folders, extract ipTM, pTM, PAE, cp_ipTM and interaction classification.",
        output="AF3 result table",
        method_name="import_af3_results_folder",
        tab_hint="AlphaFold",
    ),
    WorkflowItem(
        group="AlphaFold / PPI",
        label="Export AF3 results table",
        data_type="Parsed AF3 results",
        description="Export parsed AF3 interaction metrics as TSV/CSV for Excel, R, Python or reporting.",
        output="TSV/CSV table",
        method_name="export_af3_results_table",
        tab_hint="AlphaFold",
    ),

    WorkflowItem(
        group="Reports / System",
        label="HTML report",
        data_type="Current project snapshot",
        description="Generate a standalone HTML report summarizing genome, ORFs, BLAST, HMM and AF3 results.",
        output="HTML report",
        method_name="export_html_report",
    ),
    WorkflowItem(
        group="Reports / System",
        label="Backend status",
        data_type="External tools and Python packages",
        description="Check availability of BLAST+, HMMER3, Pyrodigal and Paramiko/SSH.",
        output="Backend availability table",
        method_name="show_backend_status",
    ),
]


def _left_dock_area():
    return Qt.DockWidgetArea.LeftDockWidgetArea if QT6 else Qt.LeftDockWidgetArea


def _switch_to_tab(window, hint: str | None) -> bool:
    if not hint:
        return False

    try:
        from PyQt6.QtWidgets import QTabWidget
    except Exception:
        from PyQt5.QtWidgets import QTabWidget

    hint_lower = hint.lower()

    for tabs in window.findChildren(QTabWidget):
        for index in range(tabs.count()):
            text = tabs.tabText(index).lower()
            tooltip = tabs.tabToolTip(index).lower()
            if hint_lower in text or hint_lower in tooltip:
                tabs.setCurrentIndex(index)
                return True

    return False


def _call_method(window, method_name: str | None) -> bool:
    if not method_name:
        return False

    method = getattr(window, method_name, None)

    if callable(method):
        method()
        return True

    return False


class WorkflowNavigator(QWidget):
    """
    Data-type-oriented analysis navigator.
    """

    def __init__(self, window):
        super().__init__(window)
        self.window = window

        layout = QVBoxLayout(self)

        title = QLabel("Analysis Navigator")
        title.setStyleSheet("font-weight: bold;")
        layout.addWidget(title)

        subtitle = QLabel("Choose by data type or analysis step.")
        subtitle.setWordWrap(True)
        layout.addWidget(subtitle)

        self.tree = QTreeWidget()
        self.tree.setColumnCount(3)
        self.tree.setHeaderLabels(["Step", "Data", "Output"])
        self.tree.itemDoubleClicked.connect(self._activate_item)
        self.tree.itemClicked.connect(self._preview_item)

        layout.addWidget(self.tree)

        self.description = QLabel("Select an item to see what it does.")
        self.description.setWordWrap(True)
        self.description.setMinimumHeight(80)
        layout.addWidget(self.description)

        self._populate()

    def _populate(self) -> None:
        groups: dict[str, QTreeWidgetItem] = {}

        for item in WORKFLOW_ITEMS:
            if item.group not in groups:
                group_item = QTreeWidgetItem([item.group, "", ""])
                group_item.setExpanded(True)
                self.tree.addTopLevelItem(group_item)
                groups[item.group] = group_item

            tooltip = (
                f"{item.label}\n\n"
                f"Input data: {item.data_type}\n"
                f"Purpose: {item.description}\n"
                f"Output: {item.output or 'N/A'}"
            )

            child = QTreeWidgetItem([
                item.label,
                item.data_type,
                item.output,
            ])
            child.setToolTip(0, tooltip)
            child.setToolTip(1, tooltip)
            child.setToolTip(2, tooltip)
            child.setData(0, 32, item)
            groups[item.group].addChild(child)

        for col in range(3):
            self.tree.resizeColumnToContents(col)

    def _preview_item(self, tree_item, column) -> None:
        item = tree_item.data(0, 32)
        if not isinstance(item, WorkflowItem):
            self.description.setText("Choose an analysis step.")
            return

        self.description.setText(
            f"<b>{item.label}</b><br>"
            f"<b>Input:</b> {item.data_type}<br>"
            f"<b>Purpose:</b> {item.description}<br>"
            f"<b>Output:</b> {item.output or 'N/A'}"
        )

    def _activate_item(self, tree_item, column) -> None:
        item = tree_item.data(0, 32)
        if not isinstance(item, WorkflowItem):
            return

        switched = _switch_to_tab(self.window, item.tab_hint)
        called = _call_method(self.window, item.method_name)

        if not switched and not called:
            try:
                self.window._status.showMessage(f"No direct action available for: {item.label}")
            except Exception:
                pass


def install_workflow_navigator(window) -> None:
    """
    Install left-side workflow navigator.
    """
    if getattr(window, "_ppig_workflow_navigator_installed", False):
        return

    dock = QDockWidget("Analysis Navigator", window)
    dock.setObjectName("analysis_navigator_dock")
    dock.setWidget(WorkflowNavigator(window))
    dock.setMinimumWidth(360)

    window.addDockWidget(_left_dock_area(), dock)
    window._ppig_workflow_navigator_installed = True
