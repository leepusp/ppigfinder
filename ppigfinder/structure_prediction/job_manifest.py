from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List

from ppigfinder.structure_prediction.hpc_planner import HPCResourcePlan
from ppigfinder.structure_prediction.models import PredictionJobSpec
from ppigfinder.structure_prediction.output_layout import PredictionOutputLayout


@dataclass(frozen=True)
class PredictionManifestRow:
    job_id: str
    backend_id: str
    model_mode: str
    target_count: int
    estimated_tokens: int
    token_class: str
    partition_hint: str
    gpu_hint: str
    cpus_per_task: int
    memory_gb: int
    time_limit: str
    array_chunk_size: int
    job_dir: str
    input_dir: str
    result_dir: str
    log_dir: str
    retry_dir: str
    status: str = "planned"

    def as_tsv_row(self) -> List[str]:
        return [
            self.job_id,
            self.backend_id,
            self.model_mode,
            str(self.target_count),
            str(self.estimated_tokens),
            self.token_class,
            self.partition_hint,
            self.gpu_hint,
            str(self.cpus_per_task),
            str(self.memory_gb),
            self.time_limit,
            str(self.array_chunk_size),
            self.job_dir,
            self.input_dir,
            self.result_dir,
            self.log_dir,
            self.retry_dir,
            self.status,
        ]


MANIFEST_HEADER = [
    "job_id",
    "backend_id",
    "model_mode",
    "target_count",
    "estimated_tokens",
    "token_class",
    "partition_hint",
    "gpu_hint",
    "cpus_per_task",
    "memory_gb",
    "time_limit",
    "array_chunk_size",
    "job_dir",
    "input_dir",
    "result_dir",
    "log_dir",
    "retry_dir",
    "status",
]


def build_manifest_row(
    job: PredictionJobSpec,
    plan: HPCResourcePlan,
    layout: PredictionOutputLayout,
    status: str = "planned",
) -> PredictionManifestRow:
    return PredictionManifestRow(
        job_id=job.job_id,
        backend_id=job.backend_id,
        model_mode=job.model_mode,
        target_count=job.target_count(),
        estimated_tokens=plan.token_count,
        token_class=plan.token_class,
        partition_hint=plan.partition_hint,
        gpu_hint=plan.gpu_hint,
        cpus_per_task=plan.cpus_per_task,
        memory_gb=plan.memory_gb,
        time_limit=plan.time_limit,
        array_chunk_size=plan.array_chunk_size,
        job_dir=str(layout.job_dir),
        input_dir=str(layout.input_dir),
        result_dir=str(layout.result_dir),
        log_dir=str(layout.log_dir),
        retry_dir=str(layout.retry_dir),
        status=status,
    )


def write_prediction_manifest(
    rows: Iterable[PredictionManifestRow],
    output_path: str | Path,
) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(MANIFEST_HEADER) + "\n")
        for row in rows:
            handle.write("\t".join(row.as_tsv_row()) + "\n")

    return output
