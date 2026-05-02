#!/usr/bin/env python3
"""
State container for the experimental guided UI shell.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class GuidedShellState:
    """
    Lightweight state used by the guided shell.

    This is intentionally independent from the legacy window. Backend parsing
    and full project synchronization will be connected progressively.
    """

    genome_file: str = ""
    project_file: str = ""
    snapshot_file: str = ""
    af3_results_folder: str = ""

    genome_loaded: bool = False
    project_loaded: bool = False
    snapshot_loaded: bool = False
    orfs_predicted: bool = False
    annotation_ready: bool = False
    af3_results_imported: bool = False
    report_ready: bool = False

    notes: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {
            "genome_file": self.genome_file,
            "project_file": self.project_file,
            "snapshot_file": self.snapshot_file,
            "af3_results_folder": self.af3_results_folder,
            "genome_loaded": self.genome_loaded,
            "project_loaded": self.project_loaded,
            "snapshot_loaded": self.snapshot_loaded,
            "orfs_predicted": self.orfs_predicted,
            "annotation_ready": self.annotation_ready,
            "af3_results_imported": self.af3_results_imported,
            "report_ready": self.report_ready,
            "notes": list(self.notes),
        }

    def set_input(self, key: str, value: str) -> None:
        setattr(self, key, value)

        if key == "genome_file":
            self.genome_loaded = bool(value)
            self.notes.append(f"Genome file selected: {value}")
        elif key == "project_file":
            self.project_loaded = bool(value)
            self.notes.append(f"Project file selected: {value}")
        elif key == "snapshot_file":
            self.snapshot_loaded = bool(value)
            self.notes.append(f"Project snapshot selected: {value}")
        elif key == "af3_results_folder":
            self.af3_results_imported = bool(value)
            self.notes.append(f"AF3 results folder selected: {value}")

    def mark_orf_prediction_planned(self) -> None:
        self.notes.append("ORF prediction step selected in guided shell.")

    def mark_annotation_planned(self) -> None:
        self.notes.append("Annotation step selected in guided shell.")

    def mark_report_planned(self) -> None:
        self.notes.append("Report/export step selected in guided shell.")
