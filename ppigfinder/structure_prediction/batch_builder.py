from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence

from ppigfinder.structure_prediction.hpc_planner import HPCResourcePlan, plan_hpc_resources
from ppigfinder.structure_prediction.job_manifest import (
    PredictionManifestRow,
    build_manifest_row,
    write_prediction_manifest,
)
from ppigfinder.structure_prediction.models import PredictionJobSpec, SequenceTarget
from ppigfinder.structure_prediction.output_layout import (
    PredictionOutputLayout,
    build_prediction_output_layout,
    safe_job_name,
)


@dataclass(frozen=True)
class PlannedPredictionJob:
    job: PredictionJobSpec
    plan: HPCResourcePlan
    layout: PredictionOutputLayout
    manifest_row: PredictionManifestRow


@dataclass(frozen=True)
class PredictionBatchPlan:
    batch_id: str
    root_dir: Path
    planned_jobs: List[PlannedPredictionJob]

    def manifest_rows(self) -> List[PredictionManifestRow]:
        return [item.manifest_row for item in self.planned_jobs]

    def create_layouts(self) -> "PredictionBatchPlan":
        for item in self.planned_jobs:
            item.layout.create()
        return self

    def write_manifest(self, filename: str = "job_manifest.tsv") -> Path:
        return write_prediction_manifest(
            self.manifest_rows(),
            self.root_dir / filename,
        )


def build_multibackend_jobs(
    base_job_id: str,
    targets: Sequence[SequenceTarget],
    backend_ids: Iterable[str] = ("af3", "boltz2", "foldcp"),
    model_mode: str = "multicomponent_complex",
    priority: str = "normal",
) -> List[PredictionJobSpec]:
    jobs: List[PredictionJobSpec] = []
    safe_base = safe_job_name(base_job_id)

    for backend_id in backend_ids:
        backend = backend_id.lower()
        job_id = f"{safe_base}_{backend}"
        jobs.append(
            PredictionJobSpec(
                job_id=job_id,
                backend_id=backend,
                targets=list(targets),
                model_mode=model_mode,
                priority=priority,
            )
        )

    return jobs


def plan_prediction_batch(
    batch_id: str,
    root_dir: str | Path,
    jobs: Sequence[PredictionJobSpec],
    run_id: str | None = None,
    create_dirs: bool = False,
) -> PredictionBatchPlan:
    root = Path(root_dir)
    planned: List[PlannedPredictionJob] = []

    for job in jobs:
        plan = plan_hpc_resources(job)
        layout = build_prediction_output_layout(root, job, run_id=run_id)

        if create_dirs:
            layout.create()

        row = build_manifest_row(job, plan, layout)

        planned.append(
            PlannedPredictionJob(
                job=job,
                plan=plan,
                layout=layout,
                manifest_row=row,
            )
        )

    return PredictionBatchPlan(
        batch_id=safe_job_name(batch_id),
        root_dir=root,
        planned_jobs=planned,
    )
