from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence

from ppigfinder.structure_prediction.batch_builder import PlannedPredictionJob, PredictionBatchPlan
from ppigfinder.structure_prediction.cluster_profiles import resolve_slurm_overrides
from ppigfinder.structure_prediction.output_layout import safe_job_name
from ppigfinder.structure_prediction.slurm_renderer import render_sbatch_header


DEFAULT_BACKEND_MODULE_COMMANDS: Dict[str, str] = {
    "af3": "module load alphafold3/1.0",
    "boltz2": "echo 'TODO: load Boltz-2 environment/module'",
    "foldcp": "echo 'TODO: load FoldCP environment/module'",
}


@dataclass(frozen=True)
class BackendSlurmFiles:
    group_id: str
    backend_id: str
    token_class: str
    script_path: Path
    task_table_path: Path
    job_count: int


@dataclass(frozen=True)
class BatchSlurmScripts:
    root_dir: Path
    slurm_dir: Path
    scripts: List[BackendSlurmFiles]

    def summary(self) -> str:
        lines = [
            f"root_dir: {self.root_dir}",
            f"slurm_dir: {self.slurm_dir}",
            f"script_count: {len(self.scripts)}",
        ]

        for item in self.scripts:
            lines.append(
                f"- {item.group_id}: jobs={item.job_count} "
                f"script={item.script_path} tasks={item.task_table_path}"
            )

        return "\n".join(lines)


def _group_id(item: PlannedPredictionJob) -> str:
    backend = safe_job_name(item.job.backend_id.lower())
    token_class = safe_job_name(item.plan.token_class.lower())
    return f"{backend}_{token_class}"


def group_planned_jobs_for_slurm(
    planned_jobs: Iterable[PlannedPredictionJob],
) -> Dict[str, List[PlannedPredictionJob]]:
    grouped: Dict[str, List[PlannedPredictionJob]] = {}

    for item in planned_jobs:
        grouped.setdefault(_group_id(item), []).append(item)

    return grouped


def write_slurm_task_table(
    items: Sequence[PlannedPredictionJob],
    output_path: str | Path,
) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    header = [
        "array_index",
        "job_id",
        "backend_id",
        "token_class",
        "job_dir",
        "input_dir",
        "result_dir",
        "log_dir",
        "retry_dir",
    ]

    with output.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(header) + "\n")

        for index, item in enumerate(items, start=1):
            row = [
                str(index),
                item.job.job_id,
                item.job.backend_id,
                item.plan.token_class,
                str(item.layout.job_dir),
                str(item.layout.input_dir),
                str(item.layout.result_dir),
                str(item.layout.log_dir),
                str(item.layout.retry_dir),
            ]
            handle.write("\t".join(row) + "\n")

    return output


def _render_backend_body(
    backend_id: str,
    task_table_path: Path,
    module_command: str,
) -> str:
    return f'''TASK_TABLE="{task_table_path}"

TASK_LINE="$(awk -F '\\t' -v idx="${{SLURM_ARRAY_TASK_ID:-1}}" 'NR > 1 && $1 == idx {{print; exit}}' "$TASK_TABLE")"

if [ -z "$TASK_LINE" ]; then
  echo "ERROR: no task found for SLURM_ARRAY_TASK_ID=${{SLURM_ARRAY_TASK_ID:-1}}" >&2
  exit 2
fi

IFS=$'\\t' read -r ARRAY_INDEX JOB_ID BACKEND_ID TOKEN_CLASS JOB_DIR INPUT_DIR RESULT_DIR LOG_DIR RETRY_DIR <<< "$TASK_LINE"

mkdir -p "$RESULT_DIR" "$LOG_DIR" "$RETRY_DIR"

echo "=== ppigFinder structural prediction task ==="
echo "array_index=$ARRAY_INDEX"
echo "job_id=$JOB_ID"
echo "backend_id=$BACKEND_ID"
echo "token_class=$TOKEN_CLASS"
echo "job_dir=$JOB_DIR"
echo "input_dir=$INPUT_DIR"
echo "result_dir=$RESULT_DIR"

{module_command}

case "{backend_id}" in
  af3)
    AF3_JSON="$INPUT_DIR/af3_input.json"
    AF3_WORKDIR="$RESULT_DIR"
    AF3_CMD="${{AF3_SCRIPT:-af3}}"

    echo "AF3 input: $AF3_JSON"
    echo "AF3 workdir: $AF3_WORKDIR"
    echo "AF3 command: $AF3_CMD"

    "$AF3_CMD" \
      --json-path "$AF3_JSON" \
      --job-name "$JOB_ID" \
      --workdir "$AF3_WORKDIR" \
      --stage all \
      --executor local \
      --force \
      --image "${{AF3_IMAGE:-}}" \
      --model-dir "${{AF3_MODEL_DIR:-}}" \
      --db-dir "${{AF3_DB_DIR:-}}"
    ;;
  boltz2)
    echo "Boltz-2 FASTA: $INPUT_DIR/boltz2_input.fasta"
    echo "Boltz-2 YAML: $INPUT_DIR/boltz2_job_spec.yaml"
    echo "TODO: call the validated local Boltz-2 runner here."
    ;;
  foldcp)
    echo "FoldCP FASTA: $INPUT_DIR/foldcp_input.fasta"
    echo "FoldCP YAML: $INPUT_DIR/foldcp_job_spec.yaml"
    echo "TODO: call the validated local FoldCP workflow here."
    ;;
  *)
    echo "ERROR: unsupported backend: {backend_id}" >&2
    exit 3
    ;;
esac
'''


def write_batch_slurm_scripts(
    batch: PredictionBatchPlan,
    cluster: str = "davinci",
    module_commands: Mapping[str, str] | None = None,
    slurm_subdir: str = "slurm",
) -> BatchSlurmScripts:
    commands = dict(DEFAULT_BACKEND_MODULE_COMMANDS)
    if module_commands:
        commands.update({key.lower(): value for key, value in module_commands.items()})

    slurm_dir = batch.root_dir / slurm_subdir
    log_dir = slurm_dir / "logs"
    slurm_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    grouped = group_planned_jobs_for_slurm(batch.planned_jobs)
    written: List[BackendSlurmFiles] = []

    for group_id, items in sorted(grouped.items()):
        first = items[0]
        backend = first.job.backend_id.lower()
        token_class = first.plan.token_class

        task_table = write_slurm_task_table(
            items,
            slurm_dir / f"tasks_{group_id}.tsv",
        )

        overrides = resolve_slurm_overrides(first.plan, cluster=cluster)

        header = render_sbatch_header(
            first.plan,
            job_name=f"ppig_{group_id}",
            overrides=overrides,
            array_job_count=len(items),
            output_log=str(log_dir / "%x_%A_%a.out"),
            error_log=str(log_dir / "%x_%A_%a.err"),
        )

        body = _render_backend_body(
            backend_id=backend,
            task_table_path=task_table,
            module_command=commands.get(backend, "echo 'TODO: load backend environment'"),
        )

        script_path = slurm_dir / f"submit_{group_id}.slurm"
        script_path.write_text(header + body, encoding="utf-8")
        script_path.chmod(0o755)

        written.append(
            BackendSlurmFiles(
                group_id=group_id,
                backend_id=backend,
                token_class=token_class,
                script_path=script_path,
                task_table_path=task_table,
                job_count=len(items),
            )
        )

    return BatchSlurmScripts(
        root_dir=batch.root_dir,
        slurm_dir=slurm_dir,
        scripts=written,
    )
