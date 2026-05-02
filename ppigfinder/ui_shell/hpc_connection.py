#!/usr/bin/env python3
"""
DaVinci/HPC connection helpers for the experimental guided UI shell.

This module intentionally starts simple:
- detects whether the user is already running on a DaVinci/GN node
- tests SSH connectivity using key/agent-based authentication
- avoids interactive password prompts in the GUI
"""

from __future__ import annotations

from dataclasses import dataclass
import getpass
import platform
import socket
import subprocess
from pathlib import Path


@dataclass
class HPCConnectionConfig:
    host: str = "davinci.icb.usp.br"
    user: str = ""
    port: int = 22
    profile: str = "DaVinci"
    key_path: str = ""


@dataclass
class HPCConnectionStatus:
    profile: str
    host: str
    user: str
    port: int
    local_hostname: str
    running_on_cluster: bool
    ssh_available: bool
    connection_ok: bool
    message: str
    command_preview: str


def default_config() -> HPCConnectionConfig:
    return HPCConnectionConfig(
        host="davinci.icb.usp.br",
        user=getpass.getuser(),
        port=22,
        profile="DaVinci",
        key_path=str(Path.home() / ".ssh" / "id_ed25519"),
    )


def local_hostname() -> str:
    return socket.gethostname()


def is_probably_davinci_node(hostname: str | None = None) -> bool:
    name = (hostname or local_hostname()).lower()
    return (
        "davinci" in name
        or name.startswith("gn")
        or name.startswith("n0")
        or name.startswith("n1")
    )


def ssh_command(config: HPCConnectionConfig, remote_command: str = "hostname") -> list[str]:
    destination = f"{config.user}@{config.host}" if config.user else config.host

    cmd = [
        "ssh",
        "-o",
        "BatchMode=yes",
        "-o",
        "ConnectTimeout=8",
        "-p",
        str(config.port),
    ]

    if config.key_path:
        key = Path(config.key_path).expanduser()
        if key.exists():
            cmd.extend(["-i", str(key)])

    cmd.extend([destination, remote_command])
    return cmd


def command_preview(config: HPCConnectionConfig) -> str:
    return " ".join(ssh_command(config, "hostname"))


def test_connection(config: HPCConnectionConfig) -> HPCConnectionStatus:
    hostname = local_hostname()
    running_here = is_probably_davinci_node(hostname)

    if running_here:
        return HPCConnectionStatus(
            profile=config.profile,
            host=config.host,
            user=config.user,
            port=config.port,
            local_hostname=hostname,
            running_on_cluster=True,
            ssh_available=True,
            connection_ok=True,
            message=(
                "This session appears to be running on a DaVinci/GN node. "
                "Guided workflows can use local cluster paths and Slurm commands directly."
            ),
            command_preview="local DaVinci/GN session detected",
        )

    cmd = ssh_command(config, "hostname")

    try:
        result = subprocess.run(
            cmd,
            check=False,
            capture_output=True,
            text=True,
            timeout=12,
        )
    except FileNotFoundError:
        return HPCConnectionStatus(
            profile=config.profile,
            host=config.host,
            user=config.user,
            port=config.port,
            local_hostname=hostname,
            running_on_cluster=False,
            ssh_available=False,
            connection_ok=False,
            message="ssh executable was not found in PATH.",
            command_preview=command_preview(config),
        )
    except subprocess.TimeoutExpired:
        return HPCConnectionStatus(
            profile=config.profile,
            host=config.host,
            user=config.user,
            port=config.port,
            local_hostname=hostname,
            running_on_cluster=False,
            ssh_available=True,
            connection_ok=False,
            message="SSH connection timed out.",
            command_preview=command_preview(config),
        )

    if result.returncode == 0:
        remote = result.stdout.strip() or "remote host"
        return HPCConnectionStatus(
            profile=config.profile,
            host=config.host,
            user=config.user,
            port=config.port,
            local_hostname=hostname,
            running_on_cluster=False,
            ssh_available=True,
            connection_ok=True,
            message=f"SSH connection OK. Remote hostname: {remote}",
            command_preview=command_preview(config),
        )

    stderr = result.stderr.strip() or result.stdout.strip() or "SSH command failed."

    return HPCConnectionStatus(
        profile=config.profile,
        host=config.host,
        user=config.user,
        port=config.port,
        local_hostname=hostname,
        running_on_cluster=False,
        ssh_available=True,
        connection_ok=False,
        message=(
            "SSH connection failed without interactive password prompt. "
            "Use SSH keys/agent or test manually in the terminal. Details: "
            + stderr
        ),
        command_preview=command_preview(config),
    )


def slurm_template_for_af3(job_name: str = "ppigfinder_af3_job") -> str:
    return f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition=max50
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output={job_name}.out
#SBATCH --error={job_name}.err

module purge

# Example placeholder.
# Replace this block with the configured AlphaFold 3 command or ppigFinder AF3 runner.

echo "Running {job_name} on $(hostname)"
date
"""
