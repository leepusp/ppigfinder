#!/usr/bin/env python3
"""
AlphaFold 3 service layer.

This service builds DaVinci-compatible AF3 commands and standalone scripts.
It is independent from the GUI and can be used by desktop, CLI or batch modes.
"""

from __future__ import annotations

from pathlib import Path
import subprocess

from ppigfinder.alphafold.af3_cli import (
    AF3CliOptions,
    build_af3_command,
    detect_af3_flags,
    get_af3_help,
    shell_join,
    write_af3_run_script,
)
from ppigfinder.hpc.profiles.davinci import DAVINCI_AF3
from ppigfinder.alphafold.server_json import (
    build_pair_jobs_from_sequences,
    build_pair_jobs_from_legacy_orfs,
    write_server_json,
)


class AlphaFoldService:
    """
    High-level AF3 orchestration service.
    """

    def __init__(self, command: str | None = None):
        self.command = command or DAVINCI_AF3.command
        self._supported_flags: set[str] | None = None

    def refresh_supported_flags(self) -> set[str]:
        """
        Inspect the local af3 CLI and cache supported flags.
        """
        help_text = get_af3_help(self.command)
        self._supported_flags = detect_af3_flags(help_text)
        return self._supported_flags

    def build_davinci_options(
        self,
        job_name: str,
        fasta: str | None = None,
        json_path: str | None = None,
        input_dir: str | None = None,
        workdir: str | None = None,
        mode: str | None = "auto",
        dry_run: bool = False,
        force: bool = False,
        partition: str | None = None,
        ntasks: int | None = None,
        mem: str | None = None,
        gres: str | None = None,
        time: str | None = None,
        resource_mode: str | None = None,
    ) -> AF3CliOptions:
        """
        Build options using DaVinci defaults.
        """
        return AF3CliOptions(
            job_name=job_name,
            fasta=fasta,
            json_path=json_path,
            input_dir=input_dir,
            workdir=workdir,
            mode=mode,
            resource_mode=resource_mode or DAVINCI_AF3.default_resource_mode,
            slurm_partition=partition or DAVINCI_AF3.default_partition,
            slurm_nodes=DAVINCI_AF3.default_nodes,
            slurm_ntasks=ntasks or DAVINCI_AF3.default_ntasks,
            slurm_mem=mem or DAVINCI_AF3.default_mem,
            slurm_gres=gres or DAVINCI_AF3.default_gres,
            slurm_time=time or DAVINCI_AF3.default_time,
            dry_run=dry_run,
            force=force,
        )

    def build_command(self, options: AF3CliOptions, adaptive: bool = True) -> list[str]:
        """
        Build an af3 command.

        adaptive=True means unsupported flags are skipped based on `af3 --help`.
        """
        supported = None

        if adaptive:
            if self._supported_flags is None:
                self.refresh_supported_flags()
            supported = self._supported_flags

        return build_af3_command(
            options,
            command=self.command,
            validate_executable=True,
            supported_flags=supported,
        )

    def command_preview(self, options: AF3CliOptions, adaptive: bool = True) -> str:
        """
        Return shell-safe command string.
        """
        return shell_join(self.build_command(options, adaptive=adaptive))

    def export_run_script(
        self,
        path: str | Path,
        options: AF3CliOptions,
        adaptive: bool = True,
    ) -> None:
        """
        Export a standalone script to run AF3 manually on the server.
        """
        command = self.build_command(options, adaptive=adaptive)
        write_af3_run_script(path, command, title=f"Run AF3 job: {options.job_name}")


    def export_server_json_from_sequence_pairs(
        self,
        path,
        pairs: list[tuple[str, str, str, str]],
        model_seeds: list[str] | None = None,
        use_structure_template: bool | None = None,
        max_template_date: str | None = None,
    ) -> None:
        """
        Export AlphaFold Server JSON from explicit protein sequence pairs.

        pairs format:
            [(name_a, sequence_a, name_b, sequence_b), ...]
        """
        jobs = build_pair_jobs_from_sequences(
            pairs,
            model_seeds=model_seeds,
            use_structure_template=use_structure_template,
            max_template_date=max_template_date,
        )
        write_server_json(path, jobs)

    def export_server_json_from_legacy_orf_pairs(
        self,
        path,
        orfs: list[dict],
        pairs: list[tuple[int, int]],
        model_seeds: list[str] | None = None,
        use_structure_template: bool | None = None,
        max_template_date: str | None = None,
    ) -> None:
        """
        Export AlphaFold Server JSON from ppigFinder ORFs and ORF index pairs.
        """
        jobs = build_pair_jobs_from_legacy_orfs(
            orfs,
            pairs,
            model_seeds=model_seeds,
            use_structure_template=use_structure_template,
            max_template_date=max_template_date,
        )
        write_server_json(path, jobs)

    def run(self, options: AF3CliOptions, adaptive: bool = True) -> subprocess.CompletedProcess:
        """
        Execute the AF3 command directly.
        """
        command = self.build_command(options, adaptive=adaptive)
        return subprocess.run(command, text=True)
