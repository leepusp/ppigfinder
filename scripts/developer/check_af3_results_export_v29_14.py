#!/usr/bin/env python3
"""Self-check for ppigFinder v29.14 AF3 TSV export."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ppigfinder.alphafold.results_export import (
    AF3_RESULTS_TSV_HEADERS,
    write_af3_results_tsv_file,
)


def main() -> int:
    result = {
        "job_name": "ORF0001_vs_ORF0002",
        "job_dir": "/tmp/example",
        "n_chains": 2,
        "chain_order": ["A", "B"],
        "chain_to_orf": {"A": "ORF0001", "B": "ORF0002"},
        "chain_lens": {"A": 100, "B": 120},
        "iptm": 0.81,
        "ptm": 0.61,
        "mean_plddt": 87.65,
        "ranking_score": 0.45,
        "pae_inter": 4.32,
        "contact_frac": 0.256,
        "best_pair": ("A", "B"),
        "pair_metrics": {("A", "B"): {"pair_pae_min": 1.23456}},
        "hotspot": {
            "hotspot_score": 0.901,
            "hotspot_mean_pae": 2.97,
            "hotspot_size_A": 5,
            "hotspot_size_B": 6,
        },
        "contact_region": "A:10-20\tB:30-40",
        "chain_iptm": [0.8, 0.7],
        "chain_ptm": [0.6, 0.5],
        "ranking_samples": [1, 2, 3],
    }

    out = Path(tempfile.mkdtemp()) / "af3.tsv"
    write_af3_results_tsv_file(
        out,
        [result],
        anchor_sequence_provider=lambda _: ("AAA\tBBB\nCCC", "DDD\nEEE"),
    )

    lines = out.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    row = lines[1].split("\t")
    idx = {name: i for i, name in enumerate(header)}

    assert header == AF3_RESULTS_TSV_HEADERS
    assert "cp_ipTM" not in header
    assert len(row) == len(header)
    assert row[idx["pae_min_inter"]] == "1.235"
    assert row[idx["pae_min_best_pair"]] == "1.235"
    assert row[idx["hotspot_score"]] == "0.901"
    assert row[idx["hotspot_mean_pae"]] == "2.97"
    assert row[idx["hotspot_size"]] == "5x6"
    assert row[idx["anchor_seq_A"]] == "AAA BBB CCC"
    assert row[idx["anchor_seq_B"]] == "DDD EEE"

    print("OK: AF3 results TSV export v29.14 self-check passed.")
    print("columns:", len(header))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
