#!/usr/bin/env python3
"""
AlphaFold/AF3 result parser.

This parser is intentionally conservative and stdlib-only at import time.
It reads common AF3 output files and returns serializable dictionaries that
can be stored in ProjectState and used by HTML reports.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any
import csv
import json

from ppigfinder.alphafold.metrics import (
    InteractionMetrics,
    contact_percent_from_pae,
    metrics_from_af3_summary,
    off_diagonal_blocks,
    pae_inter_from_blocks,
    pae_min_from_blocks,
)


@dataclass(slots=True)
class AF3ParsedResult:
    """
    Serializable AF3 result summary.
    """

    name: str
    result_dir: str

    summary_file: str | None = None
    confidence_file: str | None = None
    ranking_file: str | None = None
    model_file: str | None = None

    iptm: float | None = None
    ptm: float | None = None
    ranking_score: float | None = None

    pae_inter: float | None = None
    pae_min: float | None = None
    cp_iptm: float | None = None
    contact_percent: float | None = None

    classification: str = "UNKNOWN"
    notes: list[str] | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _read_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _read_first_csv_row(path: Path) -> dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            return dict(row)
    return {}


def _find_first(root: Path, patterns: list[str]) -> Path | None:
    for pattern in patterns:
        matches = sorted(root.rglob(pattern))
        if matches:
            return matches[0]
    return None


def _find_all(root: Path, patterns: list[str]) -> list[Path]:
    result: list[Path] = []

    for pattern in patterns:
        result.extend(root.rglob(pattern))

    return sorted(set(result))


def _merge_dicts(*items: dict | None) -> dict:
    merged = {}

    for item in items:
        if isinstance(item, dict):
            merged.update(item)

    return merged


def _extract_chain_lengths(data: dict) -> tuple[int | None, int | None]:
    """
    Try to recover two-chain lengths from confidence JSON variants.
    """
    for key in ["chain_lengths", "token_chain_lengths"]:
        value = data.get(key)
        if isinstance(value, list) and len(value) >= 2:
            try:
                return int(value[0]), int(value[1])
            except Exception:
                pass

    chains = data.get("chains")
    if isinstance(chains, list) and len(chains) >= 2:
        lengths = []
        for chain in chains[:2]:
            if isinstance(chain, dict):
                for key in ["length", "num_residues", "n_res"]:
                    if key in chain:
                        lengths.append(chain[key])
                        break
        if len(lengths) >= 2:
            try:
                return int(lengths[0]), int(lengths[1])
            except Exception:
                pass

    return None, None


def _extract_matrix(data: dict, keys: list[str]):
    for key in keys:
        value = data.get(key)
        if value is not None:
            return value
    return None


def _metrics_from_confidence_json(confidence: dict, current: InteractionMetrics) -> InteractionMetrics:
    """
    Improve metrics using confidence JSON matrices when available.

    This function must preserve metrics already obtained from summary.json,
    especially chain-pair metrics such as cp_ipTM and PAE_min.
    """
    chain_a_len, chain_b_len = _extract_chain_lengths(confidence)

    pae_matrix = _extract_matrix(
        confidence,
        ["pae", "predicted_aligned_error", "pae_matrix"],
    )

    pae_inter = current.pae_inter
    pae_min = current.pae_min
    contact_percent = current.contact_percent

    if pae_matrix is not None and chain_a_len and chain_b_len:
        ab, ba = off_diagonal_blocks(pae_matrix, chain_a_len, chain_b_len)

        computed_pae_inter = pae_inter_from_blocks(ab, ba)
        computed_pae_min = pae_min_from_blocks(ab, ba)
        computed_contact_percent = contact_percent_from_pae(ab, ba, threshold=5.0)

        if pae_inter is None:
            pae_inter = computed_pae_inter

        # Prefer the focal minimum from the PAE matrix when available.
        if computed_pae_min is not None:
            pae_min = computed_pae_min

        if contact_percent is None:
            contact_percent = computed_contact_percent

    classification = "UNKNOWN"

    try:
        from ppigfinder.alphafold.metrics import classify_interaction
        classification = classify_interaction(
            pae_min=pae_min,
            cp_iptm=current.cp_iptm,
            pae_inter=pae_inter,
        )
    except Exception:
        classification = current.classification

    return InteractionMetrics(
        iptm=current.iptm,
        ptm=current.ptm,
        ranking_score=current.ranking_score,
        pae_inter=pae_inter,
        pae_min=pae_min,
        cp_iptm=current.cp_iptm,
        contact_percent=contact_percent,
        chain_a=current.chain_a,
        chain_b=current.chain_b,
        classification=classification,
    )


def parse_af3_result_directory(path: str | Path) -> dict:
    """
    Parse one AF3 result directory into a serializable dictionary.
    """
    root = Path(path)
    notes: list[str] = []

    if not root.exists():
        raise FileNotFoundError(f"AF3 result directory not found: {root}")

    summary_file = _find_first(
        root,
        [
            "*summary*.json",
            "*ranking*.json",
            "summary.json",
        ],
    )

    confidence_file = _find_first(
        root,
        [
            "*confidences*.json",
            "*confidence*.json",
            "confidences.json",
        ],
    )

    ranking_file = _find_first(
        root,
        [
            "*ranking_scores*.csv",
            "*ranking*.csv",
        ],
    )

    model_file = _find_first(
        root,
        [
            "*.cif",
            "*.pdb",
        ],
    )

    summary_data = {}
    confidence_data = {}
    ranking_data = {}

    if summary_file:
        try:
            data = _read_json(summary_file)
            if isinstance(data, dict):
                summary_data = data
        except Exception as exc:
            notes.append(f"Could not read summary JSON: {exc}")

    if confidence_file:
        try:
            data = _read_json(confidence_file)
            if isinstance(data, dict):
                confidence_data = data
        except Exception as exc:
            notes.append(f"Could not read confidence JSON: {exc}")

    if ranking_file:
        try:
            ranking_data = _read_first_csv_row(ranking_file)
        except Exception as exc:
            notes.append(f"Could not read ranking CSV: {exc}")

    merged = _merge_dicts(summary_data, confidence_data, ranking_data)
    metrics = metrics_from_af3_summary(merged)

    if confidence_data:
        metrics = _metrics_from_confidence_json(confidence_data, metrics)

    result = AF3ParsedResult(
        name=root.name,
        result_dir=str(root),
        summary_file=str(summary_file) if summary_file else None,
        confidence_file=str(confidence_file) if confidence_file else None,
        ranking_file=str(ranking_file) if ranking_file else None,
        model_file=str(model_file) if model_file else None,
        iptm=metrics.iptm,
        ptm=metrics.ptm,
        ranking_score=metrics.ranking_score,
        pae_inter=metrics.pae_inter,
        pae_min=metrics.pae_min,
        cp_iptm=metrics.cp_iptm,
        contact_percent=metrics.contact_percent,
        classification=metrics.classification,
        notes=notes,
    )

    return result.to_dict()


def parse_af3_results_root(root: str | Path) -> list[dict]:
    """
    Parse all candidate result directories below a root directory.
    """
    root = Path(root)

    if not root.exists():
        raise FileNotFoundError(f"AF3 results root not found: {root}")

    candidates = []

    for directory in sorted([root] + [p for p in root.rglob("*") if p.is_dir()]):
        has_result_file = any(
            directory.glob(pattern)
            for pattern in [
                "*summary*.json",
                "*confidences*.json",
                "*confidence*.json",
                "*ranking_scores*.csv",
                "*.cif",
                "*.pdb",
            ]
        )
        if has_result_file:
            candidates.append(directory)

    seen = set()
    unique_candidates = []

    for item in candidates:
        resolved = str(item.resolve())
        if resolved not in seen:
            seen.add(resolved)
            unique_candidates.append(item)

    return [parse_af3_result_directory(path) for path in unique_candidates]
