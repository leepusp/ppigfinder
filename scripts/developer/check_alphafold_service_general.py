#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from ppigfinder.hpc.profiles.generic import GENERIC_AF3
from ppigfinder.hpc.profiles.davinci import DAVINCI_AF3
from ppigfinder.services.alphafold_service import AlphaFoldService


def main() -> int:
    service = AlphaFoldService(command="/bin/echo")

    if service.profile is not GENERIC_AF3:
        raise SystemExit("AlphaFoldService default profile should be GENERIC_AF3.")

    generic = service.build_slurm_options(
        job_name="generic_job",
        json_path="/tmp/input.json",
        dry_run=True,
    )

    if generic.slurm_partition is not None:
        raise SystemExit("Generic profile should not impose a Slurm partition.")
    if generic.slurm_gres is not None:
        raise SystemExit("Generic profile should not impose a GRES value.")

    davinci = service.build_davinci_options(
        job_name="davinci_job",
        json_path="/tmp/input.json",
        dry_run=True,
    )

    if davinci.slurm_partition != DAVINCI_AF3.default_partition:
        raise SystemExit("DaVinci alias did not use DaVinci partition default.")
    if davinci.slurm_gres != DAVINCI_AF3.default_gres:
        raise SystemExit("DaVinci alias did not use DaVinci GRES default.")

    cmd = service.build_command(generic, adaptive=False)
    if "--job-name" not in cmd or "--json-path" not in cmd:
        raise SystemExit("Generic command did not include required AF3 arguments.")

    with tempfile.TemporaryDirectory(prefix="ppigfinder_af3_service_check_") as tmp:
        out = Path(tmp) / "server_jobs.json"
        service.export_server_json_from_sequence_pairs(
            out,
            [("A", "MAAA", "B", "MCCC")],
            model_seeds=["1"],
        )
        data = json.loads(out.read_text(encoding="utf-8"))
        if not isinstance(data, list) or not data:
            raise SystemExit("Server JSON export did not produce a non-empty list.")

    print("OK: AlphaFoldService generic profile self-check passed.")
    print("default_profile:", service.profile)
    print("davinci_partition:", davinci.slurm_partition)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
