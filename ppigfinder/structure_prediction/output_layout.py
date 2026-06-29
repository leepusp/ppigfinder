from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from ppigfinder.structure_prediction.models import PredictionJobSpec


@dataclass(frozen=True)
class PredictionOutputLayout:
    root_dir: Path
    backend_dir: Path
    job_dir: Path
    input_dir: Path
    result_dir: Path
    log_dir: Path
    manifest_dir: Path
    retry_dir: Path

    def create(self) -> "PredictionOutputLayout":
        for path in [
            self.backend_dir,
            self.job_dir,
            self.input_dir,
            self.result_dir,
            self.log_dir,
            self.manifest_dir,
            self.retry_dir,
        ]:
            path.mkdir(parents=True, exist_ok=True)
        return self

    def manifest_path(self, filename: str = "job_manifest.tsv") -> Path:
        return self.manifest_dir / filename

    def slurm_script_path(self, filename: str = "submit.slurm") -> Path:
        return self.job_dir / filename


def safe_job_name(name: str) -> str:
    cleaned = []
    for char in name:
        if char.isalnum() or char in {"_", "-", "."}:
            cleaned.append(char)
        else:
            cleaned.append("_")
    value = "".join(cleaned).strip("_")
    return value or "prediction_job"


def build_prediction_output_layout(
    root_dir: str | Path,
    job: PredictionJobSpec,
    run_id: Optional[str] = None,
) -> PredictionOutputLayout:
    root = Path(root_dir)
    backend = safe_job_name(job.backend_id)
    job_name = safe_job_name(job.job_id)

    if run_id:
        job_name = f"{job_name}_{safe_job_name(run_id)}"

    backend_dir = root / backend
    job_dir = backend_dir / job_name

    return PredictionOutputLayout(
        root_dir=root,
        backend_dir=backend_dir,
        job_dir=job_dir,
        input_dir=job_dir / "input",
        result_dir=job_dir / "results",
        log_dir=job_dir / "logs",
        manifest_dir=job_dir / "manifests",
        retry_dir=job_dir / "retries",
    )
