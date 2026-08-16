#!/usr/bin/env python3
"""
External backend detection for ppigFinder.

Startup must be fast. Therefore, default backend detection only checks
whether executables/modules are available. Slow version probing should be
requested explicitly with detailed=True.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any


try:
    import pyrodigal

    PYRODIGAL_AVAILABLE = True
    PYRODIGAL_VERSION = getattr(pyrodigal, "__version__", "detected")
except Exception:
    pyrodigal = None
    PYRODIGAL_AVAILABLE = False
    PYRODIGAL_VERSION = None


try:
    import paramiko

    PARAMIKO_AVAILABLE = True
    PARAMIKO_VERSION = getattr(paramiko, "__version__", "detected")
except Exception:
    paramiko = None
    PARAMIKO_AVAILABLE = False
    PARAMIKO_VERSION = None


def _bundled_tools_dir() -> Path | None:
    try:
        base = Path(getattr(sys, "_MEIPASS", ""))
        if not base.is_dir():
            base = Path(sys.argv[0]).resolve().parent

        candidate = base / "tools"
        if candidate.is_dir():
            return candidate
    except Exception:
        return None

    return None


def _find_executable(command: str) -> str | None:
    bundled_tools = _bundled_tools_dir()

    if bundled_tools:
        executable_name = command + (".exe" if os.name == "nt" else "")
        candidate = bundled_tools / executable_name
        if candidate.is_file():
            return str(candidate)

    return shutil.which(command)


def _probe_version(path: str, version_args: list[str], timeout: int = 3) -> str:
    """
    Slow version check. Use only when explicitly requested.
    """
    try:
        result = subprocess.run(
            [path, *version_args],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        output = result.stdout or result.stderr or ""
        return output.splitlines()[0] if output else "detected"
    except Exception:
        return "detected"


def _detect_executable(
    tool_name: str,
    command: str,
    version_args: list[str],
    detailed: bool = False,
) -> dict[str, Any]:
    path = _find_executable(command)

    if path:
        version = _probe_version(path, version_args) if detailed else "available"
        return {
            "path": path,
            "version": version,
            "available": True,
            "wsl": False,
        }

    # Do not probe WSL during normal Linux/HPC startup.
    # WSL probing can be slow and should only be used in detailed mode on Windows.
    if detailed and os.name == "nt":
        try:
            result = subprocess.run(
                ["wsl", "bash", "-lc", f"command -v {command}"],
                capture_output=True,
                text=True,
                timeout=3,
            )
            if result.returncode == 0:
                return {
                    "path": f"wsl bash -lc {command}",
                    "version": "via WSL",
                    "available": True,
                    "wsl": True,
                }
        except Exception:
            pass

    return {
        "path": None,
        "version": None,
        "available": False,
        "wsl": False,
    }


def detect_backends(detailed: bool = False) -> dict[str, dict[str, Any]]:
    """
    Detect external and optional backends used by ppigFinder.

    detailed=False is intended for application startup.
    detailed=True may execute external commands to obtain versions.
    """
    backends: dict[str, dict[str, Any]] = {}

    backends["blast+"] = _detect_executable(
        tool_name="blast+",
        command="blastp",
        version_args=["-version"],
        detailed=detailed,
    )

    backends["hmmer3"] = _detect_executable(
        tool_name="hmmer3",
        command="hmmsearch",
        version_args=["-h"],
        detailed=detailed,
    )

    backends["pyrodigal"] = {
        "path": "pyrodigal Python module" if PYRODIGAL_AVAILABLE else None,
        "version": PYRODIGAL_VERSION if detailed else ("available" if PYRODIGAL_AVAILABLE else None),
        "available": PYRODIGAL_AVAILABLE,
        "wsl": False,
    }

    backends["paramiko"] = {
        "path": "paramiko Python module" if PARAMIKO_AVAILABLE else None,
        "version": PARAMIKO_VERSION if detailed else ("available" if PARAMIKO_AVAILABLE else None),
        "available": PARAMIKO_AVAILABLE,
        "wsl": False,
    }

    return backends


# Fast startup cache.
BACKENDS = detect_backends(detailed=False)


def refresh_backends(detailed: bool = True) -> dict[str, dict[str, Any]]:
    """
    Refresh backend information on demand.
    """
    global BACKENDS
    BACKENDS = detect_backends(detailed=detailed)
    return BACKENDS


def missing_optional_packages() -> list[str]:
    missing = []

    if not PYRODIGAL_AVAILABLE:
        missing.append("pyrodigal>=2.0")

    if not PARAMIKO_AVAILABLE:
        missing.append("paramiko>=2.9")

    return missing
