#!/usr/bin/env python3
"""
AlphaFold 3 CLI integration for ppigFinder.

This module targets an `af3` command-line wrapper. It does not assume a
specific cluster; provide a command path or load the required environment first.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shlex
import shutil
import subprocess


@dataclass
class AF3CliOptions:
    """
    Options accepted by an AF3 command-line wrapper.
    """

    job_name: str
    fasta: str | None = None
    json_path: str | None = None
    input_dir: str | None = None
    workdir: str | None = None
    mode: str | None = None
    copies: int | None = None
    model_seeds: str | None = None
    resource_mode: str | None = None
    slurm_partition: str | None = None
    slurm_nodes: int | None = None
    slurm_ntasks: int | None = None
    slurm_mem: str | None = None
    slurm_gres: str | None = None
    slurm_time: str | None = None
    dry_run: bool = False
    force: bool = False


def find_af3_command(command: str = "af3") -> str | None:
    """
    Return the path to the af3 executable if available.
    """
    return shutil.which(command)


def get_af3_help(command: str = "af3", timeout: int = 10) -> str:
    """
    Return `af3 --help` output.
    """
    executable = find_af3_command(command)

    if executable is None:
        raise FileNotFoundError(
            f"Could not find '{command}' in PATH. "
            "Load the environment where af3 is available, or provide the full path."
        )

    result = subprocess.run(
        [executable, "--help"],
        capture_output=True,
        text=True,
        timeout=timeout,
    )

    return (result.stdout or "") + (result.stderr or "")


def detect_af3_flags(help_text: str) -> set[str]:
    """
    Detect supported --flags from the help output.
    """
    flags: set[str] = set()

    for token in help_text.replace(",", " ").split():
        if token.startswith("--"):
            flags.add(token.strip())

    return flags


def build_af3_command(
    options: AF3CliOptions,
    command: str = "af3",
    validate_executable: bool = True,
    supported_flags: set[str] | None = None,
) -> list[str]:
    """
    Build an af3 command safely as a list of arguments.
    """
    executable = find_af3_command(command) if validate_executable else command

    if executable is None:
        raise FileNotFoundError(f"Could not find '{command}' in PATH.")

    cmd = [executable, "--job-name", options.job_name]

    input_count = sum(
        value is not None
        for value in [options.fasta, options.json_path, options.input_dir]
    )

    if input_count != 1:
        raise ValueError("Exactly one of fasta, json_path or input_dir must be provided.")

    def add(flag: str, value):
        if value is None:
            return
        if supported_flags is not None and flag not in supported_flags:
            return
        cmd.extend([flag, str(value)])

    def add_bool(flag: str, enabled: bool):
        if not enabled:
            return
        if supported_flags is not None and flag not in supported_flags:
            return
        cmd.append(flag)

    add("--fasta", options.fasta)
    add("--json-path", options.json_path)
    add("--input-dir", options.input_dir)
    add("--workdir", options.workdir)
    add("--mode", options.mode)
    add("--copies", options.copies)
    add("--model-seeds", options.model_seeds)
    add("--resource-mode", options.resource_mode)
    add("--slurm-partition", options.slurm_partition)
    add("--slurm-nodes", options.slurm_nodes)
    add("--slurm-ntasks", options.slurm_ntasks)
    add("--slurm-mem", options.slurm_mem)
    add("--slurm-gres", options.slurm_gres)
    add("--slurm-time", options.slurm_time)
    add_bool("--dry-run", options.dry_run)
    add_bool("--force", options.force)

    return cmd


def shell_join(command: list[str]) -> str:
    """
    Convert command list into a shell-safe string for display or scripts.
    """
    return " ".join(shlex.quote(part) for part in command)


def write_af3_run_script(
    path: str | Path,
    command: list[str],
    title: str = "Run AlphaFold 3 job",
) -> None:
    """
    Write a standalone shell script that can be copied to a server and run.
    """
    path = Path(path)

    script = f"""#!/usr/bin/env bash
set -euo pipefail

echo "{title}"
echo "Host: $(hostname)"
echo "Date: $(date)"
echo

{shell_join(command)}
"""

    path.write_text(script, encoding="utf-8")
    path.chmod(0o755)
