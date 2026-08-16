from __future__ import annotations

from dataclasses import dataclass, field
from typing import List

from ppigfinder.structure_prediction.models import PredictionJobSpec
from ppigfinder.structure_prediction.token_budget import classify_token_load


@dataclass
class HPCResourcePlan:
    job_id: str
    backend_id: str
    token_count: int
    token_class: str
    partition_hint: str
    cpus_per_task: int
    memory_gb: int
    gpu_hint: str
    time_limit: str
    array_chunk_size: int
    strategy: str
    notes: List[str] = field(default_factory=list)

    def as_sbatch_comment_block(self) -> str:
        lines = [
            "# ppigFinder structural prediction resource plan",
            f"# job_id: {self.job_id}",
            f"# backend: {self.backend_id}",
            f"# estimated_tokens: {self.token_count}",
            f"# token_class: {self.token_class}",
            f"# partition_hint: {self.partition_hint}",
            f"# cpus_per_task: {self.cpus_per_task}",
            f"# memory_gb: {self.memory_gb}",
            f"# gpu_hint: {self.gpu_hint}",
            f"# time_limit: {self.time_limit}",
            f"# array_chunk_size: {self.array_chunk_size}",
            f"# strategy: {self.strategy}",
        ]

        for note in self.notes:
            lines.append(f"# note: {note}")

        return "\n".join(lines)


def plan_hpc_resources(job: PredictionJobSpec) -> HPCResourcePlan:
    token_count = job.estimated_tokens()
    token_class = classify_token_load(token_count)
    backend = job.backend_id.lower()

    notes: List[str] = []

    if backend in {"af3", "boltz2"}:
        strategy = "cpu_preprocess_gpu_inference_cpu_postprocess"

        if token_class == "small":
            partition = "gpu_shared_or_standard"
            cpus = 6
            mem = 48
            gpu = "standard_gpu_or_slice"
            wall = "08:00:00"
            chunk = 8
        elif token_class == "medium":
            partition = "gpu_standard"
            cpus = 8
            mem = 80
            gpu = "one_standard_gpu"
            wall = "16:00:00"
            chunk = 4
        elif token_class == "large":
            partition = "gpu_high_memory_preferred"
            cpus = 12
            mem = 128
            gpu = "one_high_memory_gpu"
            wall = "24:00:00"
            chunk = 2
        else:
            partition = "gpu_high_memory_required"
            cpus = 16
            mem = 192
            gpu = "largest_available_gpu_memory"
            wall = "48:00:00"
            chunk = 1
            notes.append("Consider partitioning the complex or reducing token load.")

        notes.append("Use Slurm arrays for independent predictions.")
        notes.append("Retry failed jobs with higher memory or stronger GPU profile.")
        notes.append("Stage CPU preprocessing, GPU inference and CPU postprocessing.")

        if backend == "boltz2":
            notes.append("Use Boltz-2 as complementary backend for comparison with AF3.")

    elif backend == "foldcp":
        strategy = "cpu_or_backend_specific_structural_support"
        partition = "cpu_or_backend_specific"
        cpus = 12 if token_class in {"small", "medium"} else 20
        mem = 48 if token_class in {"small", "medium"} else 96
        gpu = "none_by_default"
        wall = "12:00:00" if token_class in {"small", "medium"} else "24:00:00"
        chunk = 8 if token_class in {"small", "medium"} else 2
        notes.append("FoldCP requirements should be adapted to the local deployment.")

    else:
        strategy = "unknown_backend_conservative_plan"
        partition = "review_required"
        cpus = 8
        mem = 64
        gpu = "review_required"
        wall = "12:00:00"
        chunk = 1
        notes.append("Unknown backend: manual review required.")

    return HPCResourcePlan(
        job_id=job.job_id,
        backend_id=job.backend_id,
        token_count=token_count,
        token_class=token_class,
        partition_hint=partition,
        cpus_per_task=cpus,
        memory_gb=mem,
        gpu_hint=gpu,
        time_limit=wall,
        array_chunk_size=chunk,
        strategy=strategy,
        notes=notes,
    )
