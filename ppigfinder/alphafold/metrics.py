#!/usr/bin/env python3
"""
AlphaFold/AF3 interaction metrics.

This module contains backend-only logic for scoring predicted protein-protein
interactions from PAE/contact probability data.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isnan
from typing import Any


@dataclass
class InteractionMetrics:
    """
    Summary metrics for one predicted interaction.
    """

    iptm: float | None = None
    ptm: float | None = None
    ranking_score: float | None = None

    pae_inter: float | None = None
    pae_min: float | None = None
    cp_iptm: float | None = None
    contact_percent: float | None = None

    chain_a: str = "A"
    chain_b: str = "B"
    classification: str = "UNKNOWN"


def _safe_float(value) -> float | None:
    try:
        if value is None:
            return None
        result = float(value)
        if isnan(result):
            return None
        return result
    except Exception:
        return None


# Public alias used by newer parser code.
safe_float = _safe_float


def mean(values: list[float]) -> float | None:
    """
    Return mean of numeric values or None.
    """
    values = [float(v) for v in values if v is not None]

    if not values:
        return None

    return sum(values) / len(values)


def min_value(values: list[float]) -> float | None:
    """
    Return minimum numeric value or None.
    """
    values = [float(v) for v in values if v is not None]

    if not values:
        return None

    return min(values)


def flatten_matrix(matrix: Any) -> list[float]:
    """
    Flatten a list-like matrix into numeric values.
    """
    if matrix is None:
        return []

    values: list[float] = []

    try:
        for row in matrix:
            if isinstance(row, (list, tuple)):
                for item in row:
                    value = _safe_float(item)
                    if value is not None:
                        values.append(value)
            else:
                value = _safe_float(row)
                if value is not None:
                    values.append(value)
    except Exception:
        return []

    return values


def off_diagonal_blocks(
    matrix: list[list[float]],
    chain_a_length: int,
    chain_b_length: int,
) -> tuple[list[list[float]], list[list[float]]]:
    """
    Extract AB and BA off-diagonal blocks from a full PAE/contact matrix.

    Matrix layout expected:
        [ A x A | A x B ]
        [ B x A | B x B ]
    """
    if not matrix:
        return [], []

    a = int(chain_a_length)
    b = int(chain_b_length)

    ab = []
    ba = []

    for i in range(0, a):
        row = matrix[i]
        ab.append(row[a:a + b])

    for i in range(a, a + b):
        row = matrix[i]
        ba.append(row[0:a])

    return ab, ba


def pae_inter_from_blocks(ab_block, ba_block=None) -> float | None:
    """
    Mean inter-chain PAE from AB and optionally BA blocks.
    """
    values = flatten_matrix(ab_block)

    if ba_block is not None:
        values.extend(flatten_matrix(ba_block))

    return mean(values)


def pae_min_from_blocks(ab_block, ba_block=None) -> float | None:
    """
    Minimum inter-chain PAE from AB and optionally BA blocks.
    """
    values = flatten_matrix(ab_block)

    if ba_block is not None:
        values.extend(flatten_matrix(ba_block))

    return min_value(values)


def contact_percent_from_pae(
    ab_block,
    ba_block=None,
    threshold: float = 5.0,
) -> float | None:
    """
    Percentage of inter-chain PAE cells below threshold.
    """
    values = flatten_matrix(ab_block)

    if ba_block is not None:
        values.extend(flatten_matrix(ba_block))

    if not values:
        return None

    n_contact = sum(1 for value in values if value < threshold)
    return 100.0 * n_contact / len(values)


def classify_interaction(
    pae_min: float | None,
    cp_iptm: float | None = None,
    pae_inter: float | None = None,
) -> str:
    """
    Classify interaction confidence.

    Conservative focal rule:
    - HIGH: PAE_min < 4 and cp_ipTM >= 0.50
    - MED:  PAE_min < 8
    - LOW:  otherwise
    """
    if pae_min is None:
        return "UNKNOWN"

    if cp_iptm is not None and pae_min < 4.0 and cp_iptm >= 0.50:
        return "HIGH"

    if pae_min < 8.0:
        return "MED"

    if pae_inter is not None and pae_inter < 8.0:
        return "MED"

    return "LOW"


def extract_pair_metric(matrix, i: int = 0, j: int = 1) -> float | None:
    """
    Extract off-diagonal value from a chain-pair matrix such as
    chain_pair_iptm or chain_pair_pae_min.
    """
    try:
        return _safe_float(matrix[i][j])
    except Exception:
        return None


def metrics_from_af3_summary(summary: dict) -> InteractionMetrics:
    """
    Build InteractionMetrics from an AF3 summary-like dictionary.

    This accepts flexible key variants used by local AF3, AF3 Server and
    older ppigFinder parsers.
    """
    iptm = _safe_float(summary.get("iptm", summary.get("ipTM")))
    ptm = _safe_float(summary.get("ptm", summary.get("pTM")))
    ranking_score = _safe_float(summary.get("ranking_score"))

    cp_iptm = safe_float(
        summary.get("cp_iptm", summary.get("chain_pair_iptm_value"))
    )
    pae_min = safe_float(
        summary.get("pae_min", summary.get("PAE_min", summary.get("chain_pair_pae_min_value")))
    )

    if "chain_pair_iptm" in summary:
        cp_iptm = extract_pair_metric(summary.get("chain_pair_iptm"), 0, 1)

    if "chain_pair_pae_min" in summary:
        pae_min = extract_pair_metric(summary.get("chain_pair_pae_min"), 0, 1)

    pae_inter = _safe_float(summary.get("pae_inter", summary.get("PAE_inter")))

    classification = classify_interaction(
        pae_min=pae_min,
        cp_iptm=cp_iptm,
        pae_inter=pae_inter,
    )

    return InteractionMetrics(
        iptm=iptm,
        ptm=ptm,
        ranking_score=ranking_score,
        pae_inter=pae_inter,
        pae_min=pae_min,
        cp_iptm=cp_iptm,
        classification=classification,
    )
