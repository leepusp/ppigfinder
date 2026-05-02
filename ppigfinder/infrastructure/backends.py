#!/usr/bin/env python3
"""
External backend detection for ppigFinder.

This module centralizes detection of command-line and optional Python
backends used by ppigFinder, such as BLAST+, HMMER3, Pyrodigal and Paramiko.
It is independent from the PyQt interface.
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
    """
    Return a bundled tools directory when running from a packaged build.

    Search order:
    1. PyInstaller runtime directory, if present.
    2. Directory next to the running script.
    """
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


def _detect_executable(
    tool_name: str,
    command: str,
    version_args: list[str],
    version_timeout: int = 5,
) -> dict[str, Any]:
    """
    Detect an executable in bundled tools, PATH, or WSL on Windows.
    """
    bundled_tools = _bundled_tools_dir()
    path = None
    version = None

    if bundled_tools:
        executable_name = command + (".exe" if os.name == "nt" else "")
        candidate = bundled_tools / executable_name
        if candidate.is_file():
            path = str(candidate)

    if path is None:
        path = shutil.which(command)

    if path:
        try:
            result = subprocess.run(
                [path, *version_args],
                capture_output=True,
                text=True,
                timeout=version_timeout,
            )
            output = result.stdout or result.stderr or ""
            version = output.splitlines()[0] if output else "detected"
        except Exception:
            version = "detected"

        return {
            "path": path,
            "version": version,
            "available": True,
            "wsl": False,
        }

    if os.name == "nt":
        try:
            result = subprocess.run(
                ["wsl", "bash", "-lc", f"{command} {' '.join(version_args)}"],
                capture_output=True,
                text=True,
                timeout=10,
            )
            output = result.stdout + result.stderr
            if result.returncode == 0 or "Usage" in output or "usage" in output:
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


def detect_backends() -> dict[str, dict[str, Any]]:
    """
    Detect external and optional backends used by ppigFinder.
    """
    backends: dict[str, dict[str, Any]] = {}

    backends["blast+"] = _detect_executable(
        tool_name="blast+",
        command="blastp",
        version_args=["-version"],
    )

    backends["hmmer3"] = _detect_executable(
        tool_name="hmmer3",
        command="hmmsearch",
        version_args=["-h"],
    )

    backends["pyrodigal"] = {
        "path": "pyrodigal Python module" if PYRODIGAL_AVAILABLE else None,
        "version": PYRODIGAL_VERSION,
        "available": PYRODIGAL_AVAILABLE,
        "wsl": False,
    }

    backends["paramiko"] = {
        "path": "paramiko Python module" if PARAMIKO_AVAILABLE else None,
        "version": PARAMIKO_VERSION,
        "available": PARAMIKO_AVAILABLE,
        "wsl": False,
    }

    return backends


BACKENDS = detect_backends()


def missing_optional_packages() -> list[str]:
    """
    Return optional Python packages that are currently unavailable.
    """
    missing = []

    if not PYRODIGAL_AVAILABLE:
        missing.append("pyrodigal>=2.0")

    if not PARAMIKO_AVAILABLE:
        missing.append("paramiko>=2.9")

    return missing
