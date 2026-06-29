#!/usr/bin/env python3
from __future__ import annotations

import tempfile
from pathlib import Path

import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from ppigfinder.services.structure_prediction_service import (
    StructuralPredictionRequest,
    make_sequence_target,
    make_sequence_targets_from_records,
    prepare_structural_prediction_batch,
)


def main() -> int:
    root = Path(tempfile.mkdtemp(prefix="ppigfinder_structure_service_check_"))

    records = [
        {"target_id": "ORF00001", "sequence": "M" * 90, "chain_id": "A"},
        {"target_id": "ORF2601", "sequence": "M" * 120, "chain_id": "B"},
        {"target_id": "ORF00002", "sequence": "M" * 70, "chain_id": "C"},
    ]

    targets = make_sequence_targets_from_records(records)

    if len(targets) != 3:
        raise SystemExit(f"Expected 3 targets, observed {len(targets)}.")

    extra = make_sequence_target("ORF_EXTRA", "M" * 10, chain_id="D")
    if extra.target_id != "ORF_EXTRA" or extra.chain_id != "D":
        raise SystemExit("make_sequence_target returned unexpected values.")

    request = StructuralPredictionRequest(
        base_job_id="ORF00001_ORF2601_vs_ORF00002",
        root_dir=root,
        targets=targets,
        backend_ids=("af3", "boltz2", "foldcp"),
        run_id="dry_run",
        cluster="davinci",
    )

    prepared = prepare_structural_prediction_batch(request)

    print(prepared.summary())

    if not prepared.ok:
        raise SystemExit("Prepared batch is not OK.")

    if prepared.artifacts.generic_input_count != 3:
        raise SystemExit("Expected 3 generic input sets.")

    if prepared.artifacts.af3_input_count != 1:
        raise SystemExit("Expected 1 AF3 input set.")

    if prepared.artifacts.boltz2_input_count != 1:
        raise SystemExit("Expected 1 Boltz-2 input set.")

    if prepared.artifacts.foldcp_input_count != 1:
        raise SystemExit("Expected 1 FoldCP input set.")

    if len(prepared.slurm_scripts.scripts) != 3:
        raise SystemExit("Expected 3 Slurm script groups.")

    manifest = root / "job_manifest.tsv"
    if not manifest.exists():
        raise SystemExit(f"Missing manifest: {manifest}")

    af3_script = root / "slurm" / "submit_af3_small.slurm"
    if not af3_script.exists():
        raise SystemExit(f"Missing AF3 Slurm script: {af3_script}")

    af3_text = af3_script.read_text(encoding="utf-8")

    for fragment in [
        "module load alphafold3/1.0",
        'AF3_CMD="${AF3_SCRIPT:-af3}"',
        "AF3_ARGS=(",
        "--executor local",
        '"$AF3_CMD" "${AF3_ARGS[@]}"',
    ]:
        if fragment not in af3_text:
            raise SystemExit(f"Missing AF3 script fragment: {fragment}")

    print()
    print("OK: structural prediction service self-check passed.")
    print(f"root={root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
