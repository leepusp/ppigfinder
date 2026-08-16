#!/usr/bin/env python3
"""Self-check for v29.14 AF3 parsed-result builder."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ppigfinder.alphafold.parsed_result import build_af3_job_result


def main() -> int:
    best_pair = ("A", "B")
    ctx = {
        "job_dir": Path("/tmp/ORF0001_vs_ORF0002"),
        "sum_path": Path("/tmp/summary.json"),
        "conf_path": Path("/tmp/confidences.json"),
        "model_cif": Path("/tmp/model.cif"),
        "ranking_csv": None,
        "input_json": None,
        "lazy": False,
        "orf_names": ["ORF0001", "ORF0002"],
        "chain_order": ["A", "B"],
        "chain_lens": {"A": 100, "B": 120},
        "chain_to_orf": {"A": "ORF0001", "B": "ORF0002"},
        "n_chains": 2,
        "iptm": 0.81,
        "ptm": 0.62,
        "mean_plddt": 87.5,
        "ranking_score": 0.45,
        "fraction_disordered": 0.01,
        "has_clash": False,
        "chain_iptm": [0.8, 0.7],
        "chain_ptm": [0.6, 0.5],
        "pae_matrix": None,
        "contact_probs": None,
        "plddt_arr": None,
        "token_res_ids": None,
        "pair_metrics": {
            best_pair: {
                "pair_pae_min": 1.234,
                "pair_iptm": 0.81,
                "contact_frac": 0.25,
                "hotspot": {"hotspot_score": 0.9},
            }
        },
        "best_pae_inter": 4.2,
        "best_pair": best_pair,
        "best_cr": "A:10-20 x B:30-40",
        "motifs": [],
        "seq_status": "ok",
        "seq_status_legacy": "ok",
        "seq_chains": [],
        "chain_layout_source": "summary_chain_iptm",
        "seq_fingerprint": "abc",
        "seq_seed_fingerprint": "abc:1",
        "model_seeds": [1],
        "completeness": "complete",
        "truncation_info": [],
        "ambiguity_info": [],
        "validation_warnings": [],
        "ranking_samples": [],
    }

    res = build_af3_job_result(ctx)

    assert res["job_name"] == "ORF0001_vs_ORF0002"
    assert res["heavy_loaded"] is True
    assert res["pae_min_inter"] == 1.234
    assert res["cp_iptm_inter"] == 0.81
    assert res["contact_frac"] == 0.25
    assert res["hotspot"]["hotspot_score"] == 0.9
    assert res["partner_name"] == "ORF0002"
    assert len(res) == 50, len(res)

    print("OK: AF3 parsed-result builder v29.14 self-check passed.")
    print("keys:", len(res))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
