#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

from ppigfinder.structure_prediction import (
    SequenceTarget,
    build_multibackend_jobs,
    plan_prediction_batch,
    validate_prediction_manifest,
    write_batch_slurm_scripts,
    write_prediction_batch_artifacts,
)


def require(path: Path) -> None:
    if not path.exists():
        raise SystemExit(f"Missing expected file: {path}")


def main() -> int:
    root = Path(tempfile.mkdtemp(prefix="ppigfinder_structure_backend_check_"))

    targets = [
        SequenceTarget("ORF00001", "M" * 90, chain_id="A"),
        SequenceTarget("ORF2601", "M" * 120, chain_id="B"),
        SequenceTarget("ORF00002", "M" * 70, chain_id="C"),
    ]

    jobs = build_multibackend_jobs(
        base_job_id="ORF00001_ORF2601_vs_ORF00002",
        targets=targets,
        backend_ids=("af3", "boltz2", "foldcp"),
    )

    batch = plan_prediction_batch(
        batch_id="structure_backend_check",
        root_dir=root,
        jobs=jobs,
        run_id="dry_run",
        create_dirs=False,
    )

    artifacts = write_prediction_batch_artifacts(batch)

    if not artifacts.valid:
        print(artifacts.summary())
        raise SystemExit("Batch artifact validation failed.")

    if artifacts.generic_input_count != 3:
        raise SystemExit(f"Expected 3 generic inputs, observed {artifacts.generic_input_count}.")

    if artifacts.af3_input_count != 1:
        raise SystemExit(f"Expected 1 AF3 input, observed {artifacts.af3_input_count}.")

    if artifacts.boltz2_input_count != 1:
        raise SystemExit(f"Expected 1 Boltz-2 input, observed {artifacts.boltz2_input_count}.")

    if artifacts.foldcp_input_count != 1:
        raise SystemExit(f"Expected 1 FoldCP input, observed {artifacts.foldcp_input_count}.")

    manifest = root / "job_manifest.tsv"
    require(manifest)

    validation = validate_prediction_manifest(manifest)
    if not validation.valid or validation.row_count != 3:
        print(validation.summary())
        raise SystemExit("Manifest validation failed.")

    af3_json = root / "af3/ORF00001_ORF2601_vs_ORF00002_af3_dry_run/input/af3_input.json"
    boltz_yaml = root / "boltz2/ORF00001_ORF2601_vs_ORF00002_boltz2_dry_run/input/boltz2_job_spec.yaml"
    foldcp_yaml = root / "foldcp/ORF00001_ORF2601_vs_ORF00002_foldcp_dry_run/input/foldcp_job_spec.yaml"

    require(af3_json)
    require(boltz_yaml)
    require(foldcp_yaml)

    payload = json.loads(af3_json.read_text(encoding="utf-8"))
    if not isinstance(payload, list):
        raise SystemExit("AF3 JSON should be list-wrapped.")

    if payload[0]["modelSeeds"] != [1]:
        raise SystemExit("Unexpected AF3 modelSeeds.")

    if len(payload[0]["sequences"]) != 3:
        raise SystemExit("AF3 JSON should contain 3 sequences.")

    for yaml_path in [boltz_yaml, foldcp_yaml]:
        text = yaml_path.read_text(encoding="utf-8")
        if "\\n" in text:
            raise SystemExit(f"Literal newline escape found in YAML: {yaml_path}")
        for token in ["job:", "targets:", "inputs:"]:
            if token not in text:
                raise SystemExit(f"Missing {token} in {yaml_path}")

    scripts = write_batch_slurm_scripts(batch, cluster="davinci")

    if len(scripts.scripts) != 3:
        raise SystemExit(f"Expected 3 Slurm scripts, observed {len(scripts.scripts)}.")

    for item in scripts.scripts:
        require(item.script_path)
        require(item.task_table_path)

        script = item.script_path.read_text(encoding="utf-8")
        subprocess.run(["bash", "-n", str(item.script_path)], check=True)

        if "TASK_TABLE=" not in script:
            raise SystemExit(f"Missing TASK_TABLE in {item.script_path}")

        if item.backend_id == "af3":
            for fragment in [
                "module load alphafold3/1.0",
                'AF3_CMD="${AF3_SCRIPT:-af3}"',
                "AF3_ARGS=(",
                "--executor local",
                '"$AF3_CMD" "${AF3_ARGS[@]}"',
            ]:
                if fragment not in script:
                    raise SystemExit(f"Missing AF3 fragment {fragment} in {item.script_path}")

    print("OK: structure_prediction backend smoke test passed.")
    print(f"root={root}")
    print(f"manifest={manifest}")
    print(f"slurm_dir={scripts.slurm_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
