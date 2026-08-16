from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List

from ppigfinder.structure_prediction.batch_builder import PredictionBatchPlan
from ppigfinder.structure_prediction.models import PredictionJobSpec, SequenceTarget
from ppigfinder.structure_prediction.output_layout import PredictionOutputLayout


@dataclass(frozen=True)
class GenericInputFiles:
    fasta_path: Path
    job_spec_path: Path


def _wrap_sequence(sequence: str, width: int = 80) -> str:
    clean = "".join((sequence or "").split())
    return "\n".join(clean[i : i + width] for i in range(0, len(clean), width))


def _safe_fasta_header(target: SequenceTarget, index: int) -> str:
    parts = [
        target.target_id or f"target_{index}",
        f"molecule_type={target.molecule_type}",
        f"role={target.role}",
    ]

    if target.chain_id:
        parts.append(f"chain_id={target.chain_id}")

    return " ".join(str(part).replace("\t", "_").replace("\n", "_") for part in parts)


def render_targets_fasta(targets: Iterable[SequenceTarget]) -> str:
    records: List[str] = []

    for index, target in enumerate(targets, start=1):
        header = _safe_fasta_header(target, index)
        sequence = _wrap_sequence(target.sequence)

        records.append(f">{header}")
        records.append(sequence)

    return "\n".join(records) + "\n"


def job_spec_to_dict(job: PredictionJobSpec) -> dict:
    return {
        "job_id": job.job_id,
        "backend_id": job.backend_id,
        "model_mode": job.model_mode,
        "priority": job.priority,
        "target_count": job.target_count(),
        "estimated_tokens": job.estimated_tokens(),
        "metadata": dict(job.metadata),
        "targets": [
            {
                "target_id": target.target_id,
                "molecule_type": target.molecule_type,
                "chain_id": target.chain_id,
                "role": target.role,
                "sequence_length": target.token_length(),
                "sequence": "".join((target.sequence or "").split()),
                "metadata": dict(target.metadata),
            }
            for target in job.targets
        ],
    }


def write_targets_fasta(job: PredictionJobSpec, output_path: str | Path) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(render_targets_fasta(job.targets), encoding="utf-8")
    return output


def write_job_spec_json(job: PredictionJobSpec, output_path: str | Path) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(job_spec_to_dict(job), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return output


def write_generic_backend_inputs(
    job: PredictionJobSpec,
    layout: PredictionOutputLayout,
) -> GenericInputFiles:
    layout.create()

    fasta_path = write_targets_fasta(
        job,
        layout.input_dir / "targets.fasta",
    )

    job_spec_path = write_job_spec_json(
        job,
        layout.input_dir / "job_spec.json",
    )

    return GenericInputFiles(
        fasta_path=fasta_path,
        job_spec_path=job_spec_path,
    )


def write_batch_generic_inputs(batch: PredictionBatchPlan) -> List[GenericInputFiles]:
    written: List[GenericInputFiles] = []

    for item in batch.planned_jobs:
        written.append(
            write_generic_backend_inputs(
                item.job,
                item.layout,
            )
        )

    return written
