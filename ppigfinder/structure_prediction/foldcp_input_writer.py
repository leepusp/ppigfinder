from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List

from ppigfinder.structure_prediction.batch_builder import PredictionBatchPlan
from ppigfinder.structure_prediction.input_writers import render_targets_fasta
from ppigfinder.structure_prediction.models import PredictionJobSpec
from ppigfinder.structure_prediction.output_layout import PredictionOutputLayout


@dataclass(frozen=True)
class FoldCPInputFiles:
    fasta_path: Path
    job_yaml_path: Path


def _yaml_scalar(value: object) -> str:
    text = "" if value is None else str(value)
    text = text.replace("\\", "\\\\").replace('"', '\\"')
    return f'"{text}"'


def render_foldcp_job_yaml(job: PredictionJobSpec) -> str:
    if job.backend_id.lower() != "foldcp":
        raise ValueError(
            f"FoldCP input writer received non-FoldCP backend: {job.backend_id}"
        )

    lines = [
        "job:",
        f"  job_id: {_yaml_scalar(job.job_id)}",
        f"  backend_id: {_yaml_scalar(job.backend_id)}",
        f"  model_mode: {_yaml_scalar(job.model_mode)}",
        f"  priority: {_yaml_scalar(job.priority)}",
        f"  target_count: {job.target_count()}",
        f"  estimated_tokens: {job.estimated_tokens()}",
        "targets:",
    ]

    for target in job.targets:
        lines.extend(
            [
                f"  - target_id: {_yaml_scalar(target.target_id)}",
                f"    molecule_type: {_yaml_scalar(target.molecule_type)}",
                f"    chain_id: {_yaml_scalar(target.chain_id or '')}",
                f"    role: {_yaml_scalar(target.role)}",
                f"    sequence_length: {target.token_length()}",
            ]
        )

    lines.extend(
        [
            "inputs:",
            "  fasta: foldcp_input.fasta",
            "expected_outputs:",
            "  - structural_comparison_or_prediction_support_results",
            "notes:",
            "  - This is a ppigFinder intermediate FoldCP input description.",
            "  - FoldCP execution should be adapted to the local workflow.",
            "  - If FoldCP requires structures instead of sequences, this file remains the traceability layer.",
        ]
    )

    return "\n".join(lines) + "\n"


def write_foldcp_input_fasta(
    job: PredictionJobSpec,
    output_path: str | Path,
) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(render_targets_fasta(job.targets), encoding="utf-8")
    return output


def write_foldcp_job_yaml(
    job: PredictionJobSpec,
    output_path: str | Path,
) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(render_foldcp_job_yaml(job), encoding="utf-8")
    return output


def write_foldcp_backend_inputs(
    job: PredictionJobSpec,
    layout: PredictionOutputLayout,
) -> FoldCPInputFiles:
    if job.backend_id.lower() != "foldcp":
        raise ValueError(
            f"FoldCP input writer received non-FoldCP backend: {job.backend_id}"
        )

    layout.create()

    fasta_path = write_foldcp_input_fasta(
        job,
        layout.input_dir / "foldcp_input.fasta",
    )

    job_yaml_path = write_foldcp_job_yaml(
        job,
        layout.input_dir / "foldcp_job_spec.yaml",
    )

    return FoldCPInputFiles(
        fasta_path=fasta_path,
        job_yaml_path=job_yaml_path,
    )


def write_batch_foldcp_inputs(batch: PredictionBatchPlan) -> List[FoldCPInputFiles]:
    written: List[FoldCPInputFiles] = []

    for item in batch.planned_jobs:
        if item.job.backend_id.lower() != "foldcp":
            continue

        written.append(
            write_foldcp_backend_inputs(
                job=item.job,
                layout=item.layout,
            )
        )

    return written
