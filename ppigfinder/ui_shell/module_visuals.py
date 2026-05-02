#!/usr/bin/env python3
"""
Module-specific visualization/status panels for the guided UI shell.
"""

from __future__ import annotations

from pathlib import Path

from ppigfinder.ui_shell.qt import (
    QFrame,
    QLabel,
    QVBoxLayout,
    QGridLayout,
)


def _short_path(path: str, max_len: int = 62) -> str:
    if not path:
        return "Not selected"

    p = str(path)
    if len(p) <= max_len:
        return p

    return "..." + p[-max_len:]


class SmallInfoBox(QFrame):
    """
    Compact information box used inside module visualization panels.
    """

    def __init__(self, title: str, value: str, description: str = "", parent=None):
        super().__init__(parent)
        self.setObjectName("Card")

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 10, 12, 10)
        layout.setSpacing(5)

        title_label = QLabel(title)
        title_label.setObjectName("CardTitle")
        layout.addWidget(title_label)

        self.value_label = QLabel(value)
        self.value_label.setWordWrap(True)
        layout.addWidget(self.value_label)

        if description:
            desc = QLabel(description)
            desc.setObjectName("CardDescription")
            desc.setWordWrap(True)
            layout.addWidget(desc)

    def set_value(self, value: str) -> None:
        self.value_label.setText(value)


class ModuleVisualizationPanel(QFrame):
    """
    Context-aware visualization/status panel for each guided module.

    These panels are designed as stepping stones toward future dynamic
    visualizations based on genome maps, charts and AF3 dashboards.
    """

    def __init__(self, module_id: str, parent=None):
        super().__init__(parent)

        self.module_id = module_id
        self.setObjectName("VisualizationPlaceholder")

        self.layout = QVBoxLayout(self)
        self.layout.setContentsMargins(18, 18, 18, 18)
        self.layout.setSpacing(12)

        title = QLabel("Module status")
        title.setObjectName("SectionTitle")
        self.layout.addWidget(title)

        self.summary = QLabel("")
        self.summary.setWordWrap(True)
        self.layout.addWidget(self.summary)

        self.grid = QGridLayout()
        self.grid.setSpacing(10)
        self.layout.addLayout(self.grid)

        self.boxes: dict[str, SmallInfoBox] = {}

        self._build_static_panel()
        self.layout.addStretch(1)

    def _add_box(self, key: str, title: str, value: str, description: str, row: int, col: int) -> None:
        box = SmallInfoBox(title, value, description)
        self.boxes[key] = box
        self.grid.addWidget(box, row, col)

    def _build_static_panel(self) -> None:
        module_id = self.module_id

        if module_id == "data":
            self.summary.setText(
                "Entrada e validação inicial do projeto. Selecione dados antes de avançar para ORFs, anotação e AF3."
            )
            self._add_box("genome_file", "Genome file", "Not selected", "FASTA, multi-FASTA, GenBank or SnapGene.", 0, 0)
            self._add_box("project_file", "Project file", "Not selected", "Saved ppigFinder project.", 0, 1)
            self._add_box("snapshot_file", "Snapshot", "Not selected", "Project Snapshot v3 JSON.", 1, 0)
            self._add_box("next", "Next step", "Protein / ORFs", "Proceed after loading or restoring data.", 1, 1)

        elif module_id == "genome":
            self.summary.setText(
                "Resumo futuro do genoma: tamanho, GC, contigs, coordenadas, tradução e preparo para mapa genômico."
            )
            self._add_box("input", "Input state", "Waiting for genome", "Depends on Data / Project.", 0, 0)
            self._add_box("summary", "Genome summary", "Pending", "Length, GC content and metadata.", 0, 1)
            self._add_box("map", "Genome map", "Pending", "Future interactive genome visualization.", 1, 0)
            self._add_box("next", "Next step", "Protein / ORFs", "Predict protein-coding regions.", 1, 1)

        elif module_id == "orfs":
            self.summary.setText(
                "Predição e revisão de ORFs que alimentarão BLAST, HMM, neighbourhood e AF3."
            )
            self._add_box("mode", "Prediction mode", "Pyrodigal / six-frame / hybrid", "Configurable ORF prediction logic.", 0, 0)
            self._add_box("table", "ORF table", "Pending", "Coordinates, strand, frame, size and sequence.", 0, 1)
            self._add_box("export", "Export", "Protein FASTA", "Downstream sequence export.", 1, 0)
            self._add_box("next", "Next step", "Annotation", "BLAST, HMM/domain and neighbourhood.", 1, 1)

        elif module_id == "annotation":
            self.summary.setText(
                "Anotação funcional para reduzir o espaço de busca e selecionar candidatos biologicamente plausíveis."
            )
            self._add_box("blast", "BLAST", "Pending", "Protein similarity search.", 0, 0)
            self._add_box("hmm", "HMM domains", "Pending", "Conserved domain profile scanning.", 0, 1)
            self._add_box("neighborhood", "Neighbourhood", "Pending", "Candidate selection by genomic context.", 1, 0)
            self._add_box("next", "Next step", "AlphaFold / PPI", "Prepare candidate pairs.", 1, 1)

        elif module_id == "alphafold":
            self.summary.setText(
                "Preparação e análise de predições estruturais AF3 para pares ou complexos proteicos candidatos."
            )
            self._add_box("json", "AF3 JSON", "Pending", "AlphaFold Server-compatible job files.", 0, 0)
            self._add_box("af3_results_folder", "AF3 results", "Not imported", "Select AF3 output folders.", 0, 1)
            self._add_box("metrics", "Metrics", "ipTM / PAE / contacts", "Interaction confidence.", 1, 0)
            self._add_box("next", "Next step", "Reports", "Export tables and HTML summaries.", 1, 1)

        elif module_id == "reports":
            self.summary.setText(
                "Exportação e reprodutibilidade: relatórios HTML, snapshots versionados e tabelas."
            )
            self._add_box("html", "HTML report", "Pending", "Standalone summary report.", 0, 0)
            self._add_box("snapshot", "Project Snapshot", "Pending", "Versioned JSON state.", 0, 1)
            self._add_box("tables", "Tables", "Pending", "TSV/CSV exports.", 1, 0)
            self._add_box("advanced", "Advanced", "Full interface available", "Access legacy tools if needed.", 1, 1)

        else:
            self.summary.setText(
                "Visão geral do fluxo guiado do ppigFinder."
            )
            self._add_box("flow", "Workflow", "Data → Genome → ORFs → Annotation → AlphaFold → Reports", "Recommended path.", 0, 0)
            self._add_box("status", "Project status", "Ready", "Start by adding data.", 0, 1)

    def update_state(self, state: dict) -> None:
        """
        Update panel from the guided workspace state.
        """
        if self.module_id == "data":
            for key in ("genome_file", "project_file", "snapshot_file"):
                if key in self.boxes:
                    self.boxes[key].set_value(_short_path(state.get(key, "")))

        elif self.module_id == "genome":
            if "input" in self.boxes:
                if state.get("genome_file"):
                    self.boxes["input"].set_value("Genome selected")
                elif state.get("project_file") or state.get("snapshot_file"):
                    self.boxes["input"].set_value("Project restored")
                else:
                    self.boxes["input"].set_value("Waiting for genome")

        elif self.module_id == "orfs":
            if "table" in self.boxes:
                self.boxes["table"].set_value(
                    "Ready to run" if state.get("genome_file") or state.get("project_file") or state.get("snapshot_file") else "Waiting for genome"
                )

        elif self.module_id == "annotation":
            if "blast" in self.boxes:
                self.boxes["blast"].set_value("Ready after ORFs")
            if "hmm" in self.boxes:
                self.boxes["hmm"].set_value("Ready after ORFs")

        elif self.module_id == "alphafold":
            if "af3_results_folder" in self.boxes:
                self.boxes["af3_results_folder"].set_value(_short_path(state.get("af3_results_folder", "")))

        elif self.module_id == "reports":
            if "snapshot" in self.boxes:
                if state.get("project_file") or state.get("snapshot_file") or state.get("genome_file"):
                    self.boxes["snapshot"].set_value("Project state available")
                else:
                    self.boxes["snapshot"].set_value("Pending")
