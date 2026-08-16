from __future__ import annotations

import json
import string
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence

from ppigfinder.structure_prediction.batch_builder import PredictionBatchPlan
from ppigfinder.structure_prediction.models import PredictionJobSpec, SequenceTarget
from ppigfinder.structure_prediction.output_layout import PredictionOutputLayout


@dataclass(frozen=True)
class AF3InputFiles:
    af3_json_path: Path


def _clean_sequence(sequence: str) -> str:
    return "".join((sequence or "").split()).upper()


def _default_chain_id(index: int) -> str:
    letters = string.ascii_uppercase

    if index < len(letters):
        return letters[index]

    first = letters[(index // len(letters)) - 1]
    second = letters[index % len(letters)]
    return first + second


def _target_chain_id(target: SequenceTarget, index: int) -> str:
    chain_id = (target.chain_id or "").strip()
    if chain_id:
        return chain_id
    return _default_chain_id(index)


def _target_to_af3_sequence(target: SequenceTarget, index: int) -> dict:
    molecule_type = (target.molecule_type or "protein").lower()

    if molecule_type != "protein":
        raise ValueError(
            "Initial AF3 writer supports only protein targets. "
            f"Unsupported molecule_type for {target.target_id}: {target.molecule_type}"
        )

    sequence = _clean_sequence(target.sequence)
    if not sequence:
        raise ValueError(f"Empty sequence for target: {target.target_id}")

    return {
        "protein": {
            "id": _target_chain_id(target, index),
            "sequence": sequence,
        }
    }


def render_af3_job_dict(
    job: PredictionJobSpec,
    model_seeds: Sequence[int] = (1,),
    dialect: str = "alphafold3",
    version: int = 1,
) -> dict:
    if job.backend_id.lower() != "af3":
        raise ValueError(
            f"AF3 input writer received non-AF3 backend: {job.backend_id}"
        )

    return {
        "name": job.job_id,
        "modelSeeds": [int(seed) for seed in model_seeds],
        "sequences": [
            _target_to_af3_sequence(target, index)
            for index, target in enumerate(job.targets)
        ],
        "dialect": dialect,
        "version": int(version),
    }


def render_af3_input_json(
    job: PredictionJobSpec,
    model_seeds: Sequence[int] = (1,),
    wrap_as_list: bool = True,
    dialect: str = "alphafold3",
    version: int = 1,
) -> str:
    payload = render_af3_job_dict(
        job=job,
        model_seeds=model_seeds,
        dialect=dialect,
        version=version,
    )

    if wrap_as_list:
        payload = [payload]

    return json.dumps(payload, indent=2, ensure_ascii=False) + "\n"


def write_af3_input_json(
    job: PredictionJobSpec,
    output_path: str | Path,
    model_seeds: Sequence[int] = (1,),
    wrap_as_list: bool = True,
    dialect: str = "alphafold3",
    version: int = 1,
) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        render_af3_input_json(
            job=job,
            model_seeds=model_seeds,
            wrap_as_list=wrap_as_list,
            dialect=dialect,
            version=version,
        ),
        encoding="utf-8",
    )
    return output


def write_af3_backend_inputs(
    job: PredictionJobSpec,
    layout: PredictionOutputLayout,
    model_seeds: Sequence[int] = (1,),
    wrap_as_list: bool = True,
) -> AF3InputFiles:
    layout.create()

    af3_json_path = write_af3_input_json(
        job=job,
        output_path=layout.input_dir / "af3_input.json",
        model_seeds=model_seeds,
        wrap_as_list=wrap_as_list,
    )

    return AF3InputFiles(af3_json_path=af3_json_path)


def write_batch_af3_inputs(
    batch: PredictionBatchPlan,
    model_seeds: Sequence[int] = (1,),
    wrap_as_list: bool = True,
) -> List[AF3InputFiles]:
    written: List[AF3InputFiles] = []

    for item in batch.planned_jobs:
        if item.job.backend_id.lower() != "af3":
            continue

        written.append(
            write_af3_backend_inputs(
                job=item.job,
                layout=item.layout,
                model_seeds=model_seeds,
                wrap_as_list=wrap_as_list,
            )
        )

    return written
