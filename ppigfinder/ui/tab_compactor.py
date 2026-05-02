#!/usr/bin/env python3
"""
Compact long tab labels while preserving full labels as tooltips.
"""

from __future__ import annotations

try:
    from PyQt6.QtWidgets import QTabWidget
except Exception:
    from PyQt5.QtWidgets import QTabWidget


TAB_LABELS = {
    "BLAST Query": "BLAST-Q",
    "BLAST Results": "BLAST-R",
    "Neighborhood": "Neighbor",
    "AlphaFold": "AF",
    "Submit AF3 via Server": "AF3 Server",
    "AlphaFold Results": "AF Results",
    "Protein": "Prot",
    "Domains": "Dom",
}


def compact_tab_labels(window) -> None:
    """
    Compact tab labels for every QTabWidget in the window.
    """
    try:
        tab_widgets = window.findChildren(QTabWidget)
    except Exception:
        return

    for tabs in tab_widgets:
        for index in range(tabs.count()):
            original = tabs.tabText(index).strip()
            compact = TAB_LABELS.get(original)

            if not compact:
                continue

            tabs.setTabToolTip(index, original)
            tabs.setTabText(index, compact)
