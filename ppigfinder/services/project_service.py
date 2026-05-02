#!/usr/bin/env python3
"""
Project snapshot service for ppigFinder.

This service does not replace the legacy save/load yet. It provides a new
versioned project snapshot format that can be used by the GUI, CLI, batch
mode and future HTML reporting.
"""

from __future__ import annotations

from pathlib import Path

from ppigfinder.domain.project import ProjectState
from ppigfinder.io.project_json import read_project_json, write_project_json, make_json_safe


class ProjectService:
    """
    High-level project persistence service.
    """

    def save(self, path: str | Path, project: ProjectState) -> None:
        write_project_json(path, project)

    def load(self, path: str | Path) -> ProjectState:
        return read_project_json(path)

    def build_snapshot_from_legacy_window(self, window) -> ProjectState:
        """
        Build ProjectState from the current legacy GUI state.
        """
        genome = {
            "name": getattr(window, "genome_name", ""),
            "sequence": getattr(window, "dna_sequence", "") or "",
            "topology": getattr(window, "topology", "linear"),
            "features": getattr(window, "sg_features", []) or [],
            "primers": getattr(window, "sg_primers", []) or [],
            "notes": getattr(window, "notes", {}) or {},
        }

        blast = {
            "results": getattr(window, "blast_results", []) or [],
            "query": getattr(window, "blast_query", ""),
        }

        hmm = {
            "hits": getattr(window, "hmm_hits_all", []) or [],
        }

        alphafold = {
            "jobs": getattr(window, "af3_jobs", []) or [],
            "results": getattr(window, "af3_results", []) or [],
            "analysis_dir": getattr(window, "af3_analysis_dir", ""),
        }

        hpc = {
            "jobs": getattr(window, "hpc_jobs", []) or [],
            "server": getattr(window, "hpc_server", ""),
        }

        ui = {
            "current_file": getattr(window, "current_file", ""),
        }

        return ProjectState(
            genome=make_json_safe(genome),
            orfs=make_json_safe(getattr(window, "orfs", []) or []),
            blast=make_json_safe(blast),
            hmm=make_json_safe(hmm),
            alphafold=make_json_safe(alphafold),
            hpc=make_json_safe(hpc),
            ui=make_json_safe(ui),
        )

    def apply_snapshot_to_legacy_window(self, window, project: ProjectState) -> None:
        """
        Apply ProjectState to the current legacy GUI state.

        This is conservative: it restores core state and then tries known UI
        refresh methods if they exist.
        """
        genome = project.genome or {}

        window.genome_name = genome.get("name", "")
        window.dna_sequence = genome.get("sequence", "") or ""
        window.topology = genome.get("topology", "linear")
        window.sg_features = genome.get("features", []) or []
        window.sg_primers = genome.get("primers", []) or []
        window.notes = genome.get("notes", {}) or {}

        window.orfs = project.orfs or []

        window.blast_results = (project.blast or {}).get("results", []) or []
        window.hmm_hits_all = (project.hmm or {}).get("hits", []) or []
        window.af3_jobs = (project.alphafold or {}).get("jobs", []) or []
        window.af3_results = (project.alphafold or {}).get("results", []) or []

        for method_name in [
            "_on_sequence_loaded",
            "_update_orf_table",
            "update_orf_table",
            "_populate_orf_table",
            "populate_orf_table",
            "_refresh_after_project_load",
        ]:
            method = getattr(window, method_name, None)
            if not callable(method):
                continue

            try:
                if method_name == "_on_sequence_loaded":
                    method("")
                else:
                    method()
            except TypeError:
                try:
                    method(None)
                except Exception:
                    pass
            except Exception:
                pass

        try:
            window._status.showMessage("Project snapshot loaded")
        except Exception:
            pass
