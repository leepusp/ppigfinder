#!/usr/bin/env python3
"""
Documentation loader for the experimental UI shell.
"""

from __future__ import annotations

from pathlib import Path

from ppigfinder.ui_shell.docs_content import docs_for


DOC_MAP = {
    "data": "data_project.md",
    "genome": "genome.md",
    "orfs": "orf_prediction.md",
    "annotation": "annotation.md",
    "alphafold": "alphafold.md",
    "hpc": "davinci_hpc.md",
    "reports": "reports.md",
}


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def user_manual_dir() -> Path:
    return project_root() / "docs" / "user_manual"


def load_module_markdown(module_id: str) -> str:
    """
    Load module documentation from docs/user_manual.

    Falls back to embedded short documentation if the file is missing.
    """
    filename = DOC_MAP.get(module_id)

    if filename:
        path = user_manual_dir() / filename
        if path.exists():
            return path.read_text(encoding="utf-8")

    docs = docs_for(module_id)
    return (
        f"# {module_id.title()}\n\n"
        f"## Purpose\n\n{docs['purpose']}\n\n"
        f"## Input\n\n{docs['input']}\n\n"
        f"## Output\n\n{docs['output']}\n\n"
        f"## Next\n\n{docs['next']}\n"
    )
