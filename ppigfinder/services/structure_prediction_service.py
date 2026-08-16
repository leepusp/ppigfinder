from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence

from ppigfinder.structure_prediction import (
    BatchSlurmScripts,
    PredictionBatchArtifacts,
    PredictionBatchPlan,
    SequenceTarget,
    build_multibackend_jobs,
    plan_prediction_batch,
    write_batch_slurm_scripts,
    write_prediction_batch_artifacts,
)


@dataclass(frozen=True)
class StructuralPredictionRequest:
    base_job_id: str
    root_dir: Path
    targets: Sequence[SequenceTarget]
    backend_ids: Sequence[str] = ("af3", "boltz2", "foldcp")
    model_mode: str = "multicomponent_complex"
    priority: str = "normal"
    run_id: str | None = None
    cluster: str = "generic_slurm"


@dataclass(frozen=True)
class StructuralPredictionPreparedBatch:
    request: StructuralPredictionRequest
    batch: PredictionBatchPlan
    artifacts: PredictionBatchArtifacts
    slurm_scripts: BatchSlurmScripts

    @property
    def ok(self) -> bool:
        return self.artifacts.valid

    def summary(self) -> str:
        lines = [
            f"base_job_id: {self.request.base_job_id}",
            f"root_dir: {self.request.root_dir}",
            f"backends: {','.join(self.request.backend_ids)}",
            f"cluster: {self.request.cluster}",
            f"ok: {self.ok}",
            "",
            "Artifacts:",
            self.artifacts.summary(),
            "",
            "Slurm:",
            self.slurm_scripts.summary(),
        ]

        return "\n".join(lines)


def prepare_structural_prediction_batch(
    request: StructuralPredictionRequest,
) -> StructuralPredictionPreparedBatch:
    jobs = build_multibackend_jobs(
        base_job_id=request.base_job_id,
        targets=request.targets,
        backend_ids=request.backend_ids,
        model_mode=request.model_mode,
        priority=request.priority,
    )

    batch = plan_prediction_batch(
        batch_id=request.base_job_id,
        root_dir=request.root_dir,
        jobs=jobs,
        run_id=request.run_id,
        create_dirs=False,
    )

    artifacts = write_prediction_batch_artifacts(batch)
    slurm_scripts = write_batch_slurm_scripts(
        batch,
        cluster=request.cluster,
    )

    return StructuralPredictionPreparedBatch(
        request=request,
        batch=batch,
        artifacts=artifacts,
        slurm_scripts=slurm_scripts,
    )


def make_sequence_target(
    target_id: str,
    sequence: str,
    chain_id: str | None = None,
    molecule_type: str = "protein",
    role: str = "target",
) -> SequenceTarget:
    return SequenceTarget(
        target_id=target_id,
        sequence=sequence,
        chain_id=chain_id,
        molecule_type=molecule_type,
        role=role,
    )


def make_sequence_targets_from_records(
    records: Iterable[dict],
) -> List[SequenceTarget]:
    targets: List[SequenceTarget] = []

    for index, record in enumerate(records, start=1):
        target_id = str(record.get("target_id") or record.get("id") or f"target_{index}")
        sequence = str(record.get("sequence") or "")
        chain_id = record.get("chain_id")

        targets.append(
            make_sequence_target(
                target_id=target_id,
                sequence=sequence,
                chain_id=None if chain_id is None else str(chain_id),
                molecule_type=str(record.get("molecule_type") or "protein"),
                role=str(record.get("role") or "target"),
            )
        )

    return targets
