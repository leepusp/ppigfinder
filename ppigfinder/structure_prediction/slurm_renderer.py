from __future__ import annotations

from typing import Iterable, Optional

from ppigfinder.structure_prediction.hpc_planner import HPCResourcePlan


def render_sbatch_header(
    plan: HPCResourcePlan,
    job_name: Optional[str] = None,
    partition: Optional[str] = None,
    gres: Optional[str] = None,
    output_log: str = "logs/%x_%A_%a.out",
    error_log: str = "logs/%x_%A_%a.err",
) -> str:
    effective_job_name = job_name or plan.job_id
    effective_partition = partition or plan.partition_hint

    lines = [
        "#!/usr/bin/env bash",
        f"#SBATCH --job-name={effective_job_name}",
        f"#SBATCH --cpus-per-task={plan.cpus_per_task}",
        f"#SBATCH --mem={plan.memory_gb}G",
        f"#SBATCH --time={plan.time_limit}",
        f"#SBATCH --output={output_log}",
        f"#SBATCH --error={error_log}",
    ]

    if effective_partition and effective_partition != "review_required":
        lines.append(f"#SBATCH --partition={effective_partition}")

    if gres:
        lines.append(f"#SBATCH --gres={gres}")

    lines.extend(
        [
            "",
            "set -euo pipefail",
            "",
            plan.as_sbatch_comment_block(),
            "",
        ]
    )

    return "\n".join(lines)


def render_array_directive(job_count: int, chunk_size: int = 1) -> str:
    if job_count <= 0:
        raise ValueError("job_count must be positive")

    chunk_size = max(1, int(chunk_size))

    if job_count == 1:
        return "# Single job: no Slurm array directive required."

    return f"#SBATCH --array=1-{job_count}%{chunk_size}"


def render_submission_preview(
    plans: Iterable[HPCResourcePlan],
    backend_module_command: str,
    job_count: int,
    gres: Optional[str] = None,
) -> str:
    plans = list(plans)
    if not plans:
        raise ValueError("at least one plan is required")

    first = plans[0]
    header = render_sbatch_header(
        first,
        job_name=f"ppig_{first.backend_id}",
        gres=gres,
    )

    array_line = render_array_directive(
        job_count=job_count,
        chunk_size=first.array_chunk_size,
    )

    body = [
        array_line,
        "",
        "# Load backend environment here.",
        backend_module_command,
        "",
        "# Replace this placeholder with backend-specific execution.",
        'echo "Running ppigFinder structural prediction task ${SLURM_ARRAY_TASK_ID:-1}"',
        "",
    ]

    return header + "\n".join(body)
