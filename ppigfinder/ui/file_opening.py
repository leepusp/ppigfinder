#!/usr/bin/env python3
"""
GUI-level file opening helpers.

This module connects backend genome loading services to the current GUI state.
"""

from __future__ import annotations

from pathlib import Path

try:
    from PyQt6.QtWidgets import QMessageBox
except Exception:
    from PyQt5.QtWidgets import QMessageBox

from ppigfinder.services.genome_service import GenomeService
from ppigfinder.ui.recent_files import add_recent_file


def open_genome_file_into_window(window, path: str | Path) -> bool:
    """
    Load a genome/sequence file and apply it to the current ppigFinder window.

    Returns True on success.
    """
    path = Path(path)

    try:
        genome = GenomeService().load_by_extension(path)
    except Exception as exc:
        QMessageBox.critical(
            window,
            "Open file",
            f"Could not open file:\n{path}\n\n{exc}",
        )
        return False

    if not genome.sequence:
        QMessageBox.warning(
            window,
            "Open file",
            f"No sequence was found in:\n{path}",
        )
        return False

    window.dna_sequence = genome.sequence.upper()
    window.genome_name = genome.name or path.stem

    # Preserve extra metadata when the legacy GUI supports it.
    if hasattr(window, "topology"):
        window.topology = genome.topology

    if hasattr(window, "sg_features"):
        window.sg_features = genome.features

    if hasattr(window, "sg_primers"):
        window.sg_primers = genome.primers

    if hasattr(window, "notes"):
        window.notes = genome.notes

    if hasattr(window, "_on_sequence_loaded"):
        window._on_sequence_loaded(str(path))

    add_recent_file(path)
    return True
