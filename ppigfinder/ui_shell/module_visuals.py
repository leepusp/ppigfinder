#!/usr/bin/env python3
"""
Module-specific visualization/status panels for the guided UI shell.
"""

from __future__ import annotations

from ppigfinder.ui_shell.orf_map_preview import ORFMapPreviewWidget
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


def _yes_no(value: object) -> str:
    return "Yes" if bool(value) else "No"


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
        self.map_preview = None

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
            self._add_box("genome_valid", "Validation", "Waiting", "Lightweight validation in guided shell.", 0, 1)
            self._add_box("genome_length", "Length", "N/A", "Total nucleotide length detected.", 1, 0)
            self._add_box("genome_gc", "GC%", "N/A", "GC content from parsed sequence.", 1, 1)
            self._add_box("project_file", "Project file", "Not selected", "Saved ppigFinder project.", 2, 0)
            self._add_box("snapshot_file", "Snapshot", "Not selected", "Project Snapshot v3 JSON.", 2, 1)

        elif module_id == "genome":
            self.summary.setText(
                "Resumo do genoma selecionado para inspeção, tradução e visualização."
            )
            self._add_box("genome_name", "Genome name", "Waiting for genome", "File name or GenBank locus.", 0, 0)
            self._add_box("genome_type", "Input type", "N/A", "FASTA, GenBank or SnapGene.", 0, 1)
            self._add_box("genome_length", "Length", "N/A", "Total nucleotide length.", 1, 0)
            self._add_box("map", "Genome map", "Pending", "Future interactive genome visualization.", 1, 1)

        elif module_id == "orfs":
            self.summary.setText(
                "Predição e revisão de ORFs que alimentarão BLAST, HMM, neighbourhood e AF3."
            )
            self._add_box("input", "Input state", "Waiting for genome", "Requires loaded genome/project.", 0, 0)
            self._add_box("mode", "Prediction mode", "Pyrodigal / six-frame / hybrid", "Configurable ORF prediction logic.", 0, 1)
            self._add_box("table", "ORF table", "Pending", "Coordinates, strand, frame, size and sequence.", 1, 0)
            self._add_box("next", "Next step", "Annotation", "BLAST, HMM/domain and neighbourhood.", 1, 1)

            map_title = QLabel("ORF map preview")
            map_title.setObjectName("SectionTitle")
            self.layout.addWidget(map_title)

            self.map_preview = ORFMapPreviewWidget()
            self.layout.addWidget(self.map_preview)

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
            self.summary.setText("Visão geral do fluxo guiado do ppigFinder.")
            self._add_box("flow", "Workflow", "Data → Genome → ORFs → Annotation → AlphaFold → Reports", "Recommended path.", 0, 0)
            self._add_box("status", "Project status", "Ready", "Start by adding data.", 0, 1)

    def update_state(self, state: dict) -> None:
        """
        Update panel from the guided workspace state.
        """
        if self.module_id == "data":
            if "genome_file" in self.boxes:
                self.boxes["genome_file"].set_value(_short_path(state.get("genome_file", "")))
            if "project_file" in self.boxes:
                self.boxes["project_file"].set_value(_short_path(state.get("project_file", "")))
            if "snapshot_file" in self.boxes:
                self.boxes["snapshot_file"].set_value(_short_path(state.get("snapshot_file", "")))
            if "genome_valid" in self.boxes:
                if state.get("genome_file"):
                    self.boxes["genome_valid"].set_value("OK" if state.get("genome_valid") else "Problem")
                else:
                    self.boxes["genome_valid"].set_value("Waiting")
            if "genome_length" in self.boxes:
                self.boxes["genome_length"].set_value(str(state.get("genome_total_length") or "N/A"))
            if "genome_gc" in self.boxes:
                gc = state.get("genome_gc_percent")
                self.boxes["genome_gc"].set_value("N/A" if gc is None else f"{gc}%")

        elif self.module_id == "genome":
            if "genome_name" in self.boxes:
                self.boxes["genome_name"].set_value(state.get("genome_name") or "Waiting for genome")
            if "genome_type" in self.boxes:
                self.boxes["genome_type"].set_value(state.get("genome_file_type") or "N/A")
            if "genome_length" in self.boxes:
                self.boxes["genome_length"].set_value(str(state.get("genome_total_length") or "N/A"))

        elif self.module_id == "orfs":
            if "input" in self.boxes:
                if state.get("genome_file") or state.get("project_file") or state.get("snapshot_file"):
                    self.boxes["input"].set_value("Input available")
                else:
                    self.boxes["input"].set_value("Waiting for genome")

            if "table" in self.boxes:
                count = state.get("guided_orf_count") or 0
                self.boxes["table"].set_value(f"{count} ORFs" if count else "Pending")

            if "mode" in self.boxes:
                if state.get("guided_orf_count"):
                    self.boxes["mode"].set_value("Guided six-frame scan")
                else:
                    self.boxes["mode"].set_value("Pyrodigal / six-frame / hybrid")

            if "next" in self.boxes:
                if state.get("guided_orf_count"):
                    longest = state.get("guided_longest_orf_aa") or 0
                    shortest = state.get("guided_shortest_orf_aa") or 0
                    self.boxes["next"].set_value(f"Longest {longest} aa / shortest {shortest} aa")
                else:
                    self.boxes["next"].set_value("Annotation")

            if self.map_preview is not None:
                self.map_preview.set_data(
                    state.get("guided_orf_map", []),
                    int(state.get("genome_total_length") or 0),
                )

        elif self.module_id == "annotation":
            if "blast" in self.boxes:
                if state.get("guided_blast_planned"):
                    self.boxes["blast"].set_value("Selected")
                elif state.get("guided_orf_count"):
                    self.boxes["blast"].set_value("Ready")
                else:
                    self.boxes["blast"].set_value("Waiting for ORFs")

            if "hmm" in self.boxes:
                if state.get("guided_hmm_planned"):
                    self.boxes["hmm"].set_value("Selected")
                elif state.get("guided_orf_count"):
                    self.boxes["hmm"].set_value("Ready")
                else:
                    self.boxes["hmm"].set_value("Waiting for ORFs")

            if "neighborhood" in self.boxes:
                if state.get("guided_neighborhood_planned"):
                    self.boxes["neighborhood"].set_value("Selected")
                elif state.get("guided_orf_count"):
                    self.boxes["neighborhood"].set_value("Ready")
                else:
                    self.boxes["neighborhood"].set_value("Waiting for ORFs")

        elif self.module_id == "alphafold":
            if "af3_results_folder" in self.boxes:
                self.boxes["af3_results_folder"].set_value(_short_path(state.get("af3_results_folder", "")))

        elif self.module_id == "reports":
            if "snapshot" in self.boxes:
                if state.get("project_file") or state.get("snapshot_file") or state.get("genome_file"):
                    self.boxes["snapshot"].set_value("Project state available")
                else:
                    self.boxes["snapshot"].set_value("Pending")
