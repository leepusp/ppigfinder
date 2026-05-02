#!/usr/bin/env python3
"""
AlphaFold Server JSON builder.

This module generates JSON files compatible with the AlphaFold Server web
format. The top-level object is always a list of jobs, even for a single job.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import json
import re

from ppigfinder.alphafold.sequence_validation import clean_protein_sequence


def sanitize_job_name(name: str) -> str:
    """
    Return a readable but safe AlphaFold job name.
    """
    name = str(name or "").strip()
    name = re.sub(r"\s+", "_", name)
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name)
    name = name.strip("._-")
    return name[:120] or "af3_job"




@dataclass(slots=True)
class ServerProteinChain:
    """
    Protein entity for AlphaFold Server JSON.
    """

    sequence: str
    count: int = 1
    use_structure_template: bool | None = None
    max_template_date: str | None = None
    unpaired_msa: str | None = None

    def to_json(self) -> dict:
        payload: dict = {
            "sequence": clean_protein_sequence(self.sequence),
            "count": int(self.count),
        }

        if self.use_structure_template is not None:
            payload["useStructureTemplate"] = bool(self.use_structure_template)

        if self.max_template_date:
            payload["maxTemplateDate"] = self.max_template_date

        if self.unpaired_msa is not None:
            payload["unpairedMsa"] = self.unpaired_msa

        return {"proteinChain": payload}


@dataclass(slots=True)
class ServerJob:
    """
    One AlphaFold Server job.
    """

    name: str
    sequences: list[ServerProteinChain] = field(default_factory=list)
    model_seeds: list[str] = field(default_factory=list)
    version: int = 1

    def to_json(self) -> dict:
        if not self.sequences:
            raise ValueError(f"Job {self.name!r} has no sequences.")

        return {
            "name": sanitize_job_name(self.name),
            "modelSeeds": list(self.model_seeds),
            "sequences": [entity.to_json() for entity in self.sequences],
            "dialect": "alphafoldserver",
            "version": self.version,
        }


def build_server_job(
    name: str,
    protein_sequences: list[str],
    model_seeds: list[str] | None = None,
    use_structure_template: bool | None = None,
    max_template_date: str | None = None,
) -> ServerJob:
    """
    Build one AlphaFold Server job from protein sequences.
    """
    chains = [
        ServerProteinChain(
            sequence=sequence,
            count=1,
            use_structure_template=use_structure_template,
            max_template_date=max_template_date,
        )
        for sequence in protein_sequences
    ]

    return ServerJob(
        name=name,
        sequences=chains,
        model_seeds=model_seeds or [],
    )


def build_pair_job(
    name: str,
    protein_a: str,
    protein_b: str,
    model_seeds: list[str] | None = None,
    use_structure_template: bool | None = None,
    max_template_date: str | None = None,
) -> ServerJob:
    """
    Build a two-chain protein-protein interaction job.
    """
    return build_server_job(
        name=name,
        protein_sequences=[protein_a, protein_b],
        model_seeds=model_seeds,
        use_structure_template=use_structure_template,
        max_template_date=max_template_date,
    )


def build_pair_jobs_from_sequences(
    pairs: list[tuple[str, str, str, str]],
    model_seeds: list[str] | None = None,
    use_structure_template: bool | None = None,
    max_template_date: str | None = None,
) -> list[ServerJob]:
    """
    Build multiple PPI jobs.

    pairs format:
        [(name_a, sequence_a, name_b, sequence_b), ...]
    """
    jobs: list[ServerJob] = []

    for name_a, seq_a, name_b, seq_b in pairs:
        job_name = f"{name_a}_x_{name_b}"
        jobs.append(
            build_pair_job(
                name=job_name,
                protein_a=seq_a,
                protein_b=seq_b,
                model_seeds=model_seeds,
                use_structure_template=use_structure_template,
                max_template_date=max_template_date,
            )
        )

    return jobs


def build_pair_jobs_from_legacy_orfs(
    orfs: list[dict],
    pairs: list[tuple[int, int]],
    model_seeds: list[str] | None = None,
    use_structure_template: bool | None = None,
    max_template_date: str | None = None,
) -> list[ServerJob]:
    """
    Build AlphaFold Server jobs from legacy ORF dictionaries and ORF index pairs.

    pairs are zero-based ORF indexes:
        [(0, 1), (0, 2), ...]
    """
    jobs: list[ServerJob] = []

    for index_a, index_b in pairs:
        orf_a = orfs[index_a]
        orf_b = orfs[index_b]

        name_a = orf_a.get("id") or f"ORF{index_a + 1}"
        name_b = orf_b.get("id") or f"ORF{index_b + 1}"

        seq_a = str(orf_a.get("protein", "")).rstrip("*")
        seq_b = str(orf_b.get("protein", "")).rstrip("*")

        jobs.append(
            build_pair_job(
                name=f"{name_a}_x_{name_b}",
                protein_a=seq_a,
                protein_b=seq_b,
                model_seeds=model_seeds,
                use_structure_template=use_structure_template,
                max_template_date=max_template_date,
            )
        )

    return jobs


def jobs_to_json_payload(jobs: list[ServerJob]) -> list[dict]:
    """
    Convert jobs to AlphaFold Server JSON payload.

    The top-level payload is always a list.
    """
    return [job.to_json() for job in jobs]


def write_server_json(path: str | Path, jobs: list[ServerJob], indent: int = 2) -> None:
    """
    Write AlphaFold Server-compatible JSON.
    """
    path = Path(path)
    payload = jobs_to_json_payload(jobs)

    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=indent)
        handle.write("\n")


def validate_server_payload(payload: object) -> None:
    """
    Validate the minimal AlphaFold Server JSON structure.
    """
    if not isinstance(payload, list):
        raise ValueError("AlphaFold Server JSON must be a top-level list.")

    if not payload:
        raise ValueError("AlphaFold Server JSON contains no jobs.")

    for job in payload:
        if not isinstance(job, dict):
            raise ValueError("Each AlphaFold Server job must be a dictionary.")

        if job.get("dialect") != "alphafoldserver":
            raise ValueError("Each job must have dialect='alphafoldserver'.")

        if "name" not in job:
            raise ValueError("Each job must have a name.")

        if "modelSeeds" not in job:
            raise ValueError("Each job must have modelSeeds.")

        if not isinstance(job.get("sequences"), list) or not job["sequences"]:
            raise ValueError("Each job must have at least one sequence.")
