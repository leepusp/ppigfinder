#!/usr/bin/env python3
"""
Guided AlphaFold pair builder.

This module is UI-independent and converts guided ORF objects into
AlphaFold Server-compatible candidate pair jobs.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import re


@dataclass
class GuidedPair:
    name: str
    orf_a_id: str
    orf_b_id: str
    sequence_a: str
    sequence_b: str
    start_a: int = 0
    start_b: int = 0
    distance_nt: int = 0


def _clean_job_name(name: str) -> str:
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name.strip())
    name = re.sub(r"_+", "_", name).strip("_")
    return name[:80] or "ppigfinder_pair"


def _getattr_safe(obj, name: str, default=None):
    if isinstance(obj, dict):
        return obj.get(name, default)
    return getattr(obj, name, default)


def make_adjacent_orf_pairs(orfs, max_pairs: int = 50) -> list[GuidedPair]:
    """
    Build adjacent ORF pairs from guided ORF predictions.

    This is the first guided strategy. Later strategies can add:
    - query versus neighbours
    - HMM-positive ORFs
    - selected ORFs all-vs-all
    - custom stoichiometry
    """
    sorted_orfs = sorted(
        list(orfs or []),
        key=lambda item: int(_getattr_safe(item, "start", 0) or 0),
    )

    pairs: list[GuidedPair] = []

    for index in range(max(0, len(sorted_orfs) - 1)):
        if len(pairs) >= max_pairs:
            break

        a = sorted_orfs[index]
        b = sorted_orfs[index + 1]

        a_id = str(_getattr_safe(a, "id", f"orf_{index + 1}"))
        b_id = str(_getattr_safe(b, "id", f"orf_{index + 2}"))

        seq_a = str(_getattr_safe(a, "protein_sequence", "") or "")
        seq_b = str(_getattr_safe(b, "protein_sequence", "") or "")

        if not seq_a or not seq_b:
            continue

        start_a = int(_getattr_safe(a, "start", 0) or 0)
        end_a = int(_getattr_safe(a, "end", start_a) or start_a)
        start_b = int(_getattr_safe(b, "start", 0) or 0)

        distance = max(0, start_b - end_a)

        name = _clean_job_name(f"{a_id}_vs_{b_id}")

        pairs.append(
            GuidedPair(
                name=name,
                orf_a_id=a_id,
                orf_b_id=b_id,
                sequence_a=seq_a,
                sequence_b=seq_b,
                start_a=start_a,
                start_b=start_b,
                distance_nt=distance,
            )
        )

    return pairs


def pair_to_af3_job(pair: GuidedPair, seeds: list[int] | None = None) -> dict:
    """
    Convert one guided pair to AlphaFold Server-style JSON.
    """
    seeds = seeds or [1]

    return {
        "name": pair.name,
        "dialect": "alphafoldserver",
        "version": 1,
        "modelSeeds": seeds,
        "sequences": [
            {
                "proteinChain": {
                    "sequence": pair.sequence_a,
                    "count": 1,
                }
            },
            {
                "proteinChain": {
                    "sequence": pair.sequence_b,
                    "count": 1,
                }
            },
        ],
        "metadata": {
            "source": "ppigFinder guided shell",
            "orf_a": pair.orf_a_id,
            "orf_b": pair.orf_b_id,
            "distance_nt": pair.distance_nt,
        },
    }


def pairs_to_af3_payload(pairs: list[GuidedPair], seeds: list[int] | None = None) -> list[dict]:
    return [pair_to_af3_job(pair, seeds=seeds) for pair in pairs]


def write_af3_json(path: str, pairs: list[GuidedPair], seeds: list[int] | None = None) -> None:
    payload = pairs_to_af3_payload(pairs, seeds=seeds)

    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")
