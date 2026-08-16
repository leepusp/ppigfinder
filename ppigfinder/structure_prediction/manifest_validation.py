from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List

from ppigfinder.structure_prediction.job_manifest import MANIFEST_HEADER


VALID_STATUSES = {
    "planned",
    "submitted",
    "running",
    "completed",
    "failed",
    "retry_planned",
    "skipped",
}


@dataclass(frozen=True)
class ManifestValidationResult:
    path: Path
    row_count: int
    valid: bool
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    def summary(self) -> str:
        status = "OK" if self.valid else "FAILED"
        lines = [
            f"manifest: {self.path}",
            f"status: {status}",
            f"rows: {self.row_count}",
            f"errors: {len(self.errors)}",
            f"warnings: {len(self.warnings)}",
        ]

        for error in self.errors:
            lines.append(f"ERROR: {error}")

        for warning in self.warnings:
            lines.append(f"WARNING: {warning}")

        return "\n".join(lines)


def read_manifest_rows(path: str | Path) -> List[Dict[str, str]]:
    manifest = Path(path)
    lines = manifest.read_text(encoding="utf-8").splitlines()

    if not lines:
        return []

    header = lines[0].split("\t")
    rows: List[Dict[str, str]] = []

    for line in lines[1:]:
        fields = line.split("\t")
        row = {key: value for key, value in zip(header, fields)}
        rows.append(row)

    return rows


def validate_prediction_manifest(path: str | Path) -> ManifestValidationResult:
    manifest = Path(path)
    errors: List[str] = []
    warnings: List[str] = []

    if not manifest.exists():
        return ManifestValidationResult(
            path=manifest,
            row_count=0,
            valid=False,
            errors=[f"Manifest does not exist: {manifest}"],
        )

    lines = manifest.read_text(encoding="utf-8").splitlines()

    if not lines:
        return ManifestValidationResult(
            path=manifest,
            row_count=0,
            valid=False,
            errors=["Manifest is empty."],
        )

    header = lines[0].split("\t")

    if header != MANIFEST_HEADER:
        errors.append(
            "Header does not match MANIFEST_HEADER. "
            f"Expected {len(MANIFEST_HEADER)} columns, observed {len(header)}."
        )

    expected_cols = len(MANIFEST_HEADER)
    seen_job_ids = set()

    for line_number, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")

        if len(fields) != expected_cols:
            errors.append(
                f"Line {line_number} has {len(fields)} columns; expected {expected_cols}."
            )
            continue

        row = dict(zip(MANIFEST_HEADER, fields))
        job_id = row["job_id"]

        if not job_id:
            errors.append(f"Line {line_number} has empty job_id.")

        if job_id in seen_job_ids:
            errors.append(f"Duplicate job_id at line {line_number}: {job_id}")

        seen_job_ids.add(job_id)

        status = row["status"]
        if status not in VALID_STATUSES:
            warnings.append(
                f"Line {line_number} has non-standard status '{status}'."
            )

        for col in ["job_dir", "input_dir", "result_dir", "log_dir", "retry_dir"]:
            if not row[col]:
                errors.append(f"Line {line_number} has empty {col}.")

    return ManifestValidationResult(
        path=manifest,
        row_count=max(0, len(lines) - 1),
        valid=not errors,
        errors=errors,
        warnings=warnings,
    )
