#!/usr/bin/env python3
"""
AlphaFold interaction classification helpers.
"""

from __future__ import annotations

from ppigfinder.alphafold.metrics import classify_interaction, InteractionMetrics


def classify_metrics(metrics: InteractionMetrics) -> str:
    """
    Classify an InteractionMetrics object.
    """
    return classify_interaction(
        pae_min=metrics.pae_min,
        cp_iptm=metrics.cp_iptm,
        pae_inter=metrics.pae_inter,
    )


def classification_color(label: str) -> str:
    """
    Return a simple semantic color for report/GUI usage.
    """
    label = (label or "").upper()

    if label == "HIGH":
        return "#2e7d32"

    if label == "MED":
        return "#ef6c00"

    if label == "LOW":
        return "#b71c1c"

    return "#607d8b"
