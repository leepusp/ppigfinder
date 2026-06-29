from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from ppigfinder.structure_prediction.af3_input_writer import write_batch_af3_inputs
from ppigfinder.structure_prediction.batch_builder import PredictionBatchPlan
from ppigfinder.structure_prediction.boltz2_input_writer import write_batch_boltz2_inputs
from ppigfinder.structure_prediction.foldcp_input_writer import write_batch_foldcp_inputs
from ppigfinder.structure_prediction.input_writers import write_batch_generic_inputs
from ppigfinder.structure_prediction.manifest_validation import (
    ManifestValidationResult,
    validate_prediction_manifest,
)


@dataclass(frozen=True)
class PredictionBatchArtifacts:
    batch_id: str
    root_dir: Path
    manifest_path: Optional[Path]
    manifest_validation: Optional[ManifestValidationResult]
    generic_input_count: int
    af3_input_count: int
    boltz2_input_count: int
    foldcp_input_count: int

    @property
    def valid(self) -> bool:
        if self.manifest_validation is None:
            return True
        return self.manifest_validation.valid

    def summary(self) -> str:
        lines = [
            f"batch_id: {self.batch_id}",
            f"root_dir: {self.root_dir}",
            f"manifest_path: {self.manifest_path}",
            f"manifest_valid: {self.valid}",
            f"generic_input_count: {self.generic_input_count}",
            f"af3_input_count: {self.af3_input_count}",
            f"boltz2_input_count: {self.boltz2_input_count}",
            f"foldcp_input_count: {self.foldcp_input_count}",
        ]

        if self.manifest_validation is not None:
            lines.append("")
            lines.append(self.manifest_validation.summary())

        return "\n".join(lines)


def write_prediction_batch_artifacts(
    batch: PredictionBatchPlan,
    write_manifest: bool = True,
    write_inputs: bool = True,
    validate_manifest: bool = True,
) -> PredictionBatchArtifacts:
    batch.create_layouts()

    manifest_path: Optional[Path] = None
    validation: Optional[ManifestValidationResult] = None

    if write_manifest:
        manifest_path = batch.write_manifest()

        if validate_manifest:
            validation = validate_prediction_manifest(manifest_path)

    generic_count = 0
    af3_count = 0
    boltz2_count = 0
    foldcp_count = 0

    if write_inputs:
        generic_count = len(write_batch_generic_inputs(batch))
        af3_count = len(write_batch_af3_inputs(batch))
        boltz2_count = len(write_batch_boltz2_inputs(batch))
        foldcp_count = len(write_batch_foldcp_inputs(batch))

    return PredictionBatchArtifacts(
        batch_id=batch.batch_id,
        root_dir=batch.root_dir,
        manifest_path=manifest_path,
        manifest_validation=validation,
        generic_input_count=generic_count,
        af3_input_count=af3_count,
        boltz2_input_count=boltz2_count,
        foldcp_input_count=foldcp_count,
    )
