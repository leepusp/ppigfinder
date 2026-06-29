#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

from ppigfinder.structure_prediction import (
    SequenceTarget,
    build_multibackend_jobs,
    plan_prediction_batch,
    write_batch_slurm_scripts,
    write_prediction_batch_artifacts,
)


def main() -> int:
    root = Path(tempfile.mkdtemp(prefix="ppigfinder_af3_slurm_check_"))

    targets = [
        SequenceTarget("ORF00001", "M" * 90, chain_id="A"),
        SequenceTarget("ORF2601", "M" * 120, chain_id="B"),
        SequenceTarget("ORF00002", "M" * 70, chain_id="C"),
    ]

    jobs = build_multibackend_jobs(
        base_job_id="ORF00001_ORF2601_vs_ORF00002",
        targets=targets,
        backend_ids=("af3",),
    )

    batch = plan_prediction_batch(
        batch_id="af3_slurm_check",
        root_dir=root,
        jobs=jobs,
        run_id="dry_run",
        create_dirs=False,
    )

    artifacts = write_prediction_batch_artifacts(batch)
    if not artifacts.valid:
        print(artifacts.summary())
        raise SystemExit("Artifact validation failed.")

    scripts = write_batch_slurm_scripts(batch, cluster="davinci")
    if len(scripts.scripts) != 1:
        raise SystemExit(f"Expected one AF3 Slurm script, observed {len(scripts.scripts)}.")

    script_path = scripts.scripts[0].script_path
    task_table_path = scripts.scripts[0].task_table_path

    script = script_path.read_text(encoding="utf-8")

    required_fragments = [
        "module load alphafold3/1.0",
        'AF3_CMD="${AF3_SCRIPT:-af3}"',
        "AF3_ARGS=(",
        '--json-path "$AF3_JSON"',
        "--stage all",
        "--executor local",
        '"$AF3_CMD" "${AF3_ARGS[@]}"',
        "#SBATCH --partition=max50",
        "#SBATCH --gres=shard:10",
        "TASK_TABLE=",
    ]

    for fragment in required_fragments:
        if fragment not in script:
            raise SystemExit(f"Missing expected fragment: {fragment}")

    forbidden_fragments = [
        "TODO: call the validated local AF3 runner here.",
        "--executor slurm",
    ]

    for fragment in forbidden_fragments:
        if fragment in script:
            raise SystemExit(f"Forbidden fragment found: {fragment}")

    subprocess.run(["bash", "-n", str(script_path)], check=True)

    if not task_table_path.exists():
        raise SystemExit(f"Missing task table: {task_table_path}")

    print("OK: AF3 Slurm renderer self-check passed.")
    print(f"root={root}")
    print(f"script={script_path}")
    print(f"tasks={task_table_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
