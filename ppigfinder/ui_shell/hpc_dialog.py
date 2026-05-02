#!/usr/bin/env python3
"""
Flexible DaVinci/HPC dialog for the guided ppigFinder shell.

This dialog is intentionally optional. Users without a server can still export
AF3 JSON, reports, tables and figures. Users with DaVinci/HPC access can prepare
or adapt Slurm templates for several tools.
"""

from __future__ import annotations

from dataclasses import dataclass
import getpass
import os
import socket
import subprocess
from pathlib import Path

try:
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import (
        QDialog,
        QVBoxLayout,
        QHBoxLayout,
        QGridLayout,
        QLabel,
        QLineEdit,
        QSpinBox,
        QPushButton,
        QTextEdit,
        QComboBox,
        QFileDialog,
        QMessageBox,
        QTabWidget,
        QWidget,
        QGroupBox,
    )
    QT6 = True
except Exception:
    from PyQt5.QtCore import Qt
    from PyQt5.QtWidgets import (
        QDialog,
        QVBoxLayout,
        QHBoxLayout,
        QGridLayout,
        QLabel,
        QLineEdit,
        QSpinBox,
        QPushButton,
        QTextEdit,
        QComboBox,
        QFileDialog,
        QMessageBox,
        QTabWidget,
        QWidget,
        QGroupBox,
    )
    QT6 = False

from ppigfinder.ui_shell.branding import apply_ppigfinder_branding


@dataclass
class HPCConnectionConfig:
    profile: str = "DaVinci"
    host: str = "davinci.icb.usp.br"
    user: str = ""
    port: int = 22
    ssh_key_path: str = ""
    scheduler: str = "Slurm"
    partition: str = "max50"
    cpus: int = 8
    memory: str = "64G"
    time_limit: str = "24:00:00"
    gres: str = ""
    workdir: str = "$PWD"


@dataclass
class HPCConnectionStatus:
    connection_ok: bool = False
    running_on_cluster: bool = False
    scheduler_available: bool = False
    hostname: str = ""
    message: str = ""


def _window_flags():
    flags = Qt.WindowType.Window if QT6 else Qt.Window
    flags |= Qt.WindowType.WindowMinimizeButtonHint if QT6 else Qt.WindowMinimizeButtonHint
    flags |= Qt.WindowType.WindowMaximizeButtonHint if QT6 else Qt.WindowMaximizeButtonHint
    flags |= Qt.WindowType.WindowCloseButtonHint if QT6 else Qt.WindowCloseButtonHint
    return flags


def _default_user() -> str:
    return os.environ.get("USER") or getpass.getuser() or ""


def _default_key() -> str:
    home = Path.home()
    for name in ("id_ed25519", "id_rsa"):
        candidate = home / ".ssh" / name
        if candidate.exists():
            return str(candidate)
    return str(home / ".ssh" / "id_ed25519")


def _detect_local_cluster() -> tuple[bool, bool, str]:
    hostname = socket.gethostname()
    looks_like_cluster = any(token in hostname.lower() for token in ("davinci", "gn", "n01"))

    try:
        result = subprocess.run(
            ["bash", "-lc", "command -v sbatch >/dev/null 2>&1 && command -v squeue >/dev/null 2>&1"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=3,
        )
        slurm_available = result.returncode == 0
    except Exception:
        slurm_available = False

    return looks_like_cluster, slurm_available, hostname


def _template_for(tool: str, cfg: HPCConnectionConfig) -> str:
    gres_line = f"#SBATCH --gres={cfg.gres}\n" if cfg.gres.strip() else ""

    header = f"""#!/bin/bash
#SBATCH --job-name=ppigfinder_{tool.lower().replace(' ', '_').replace('/', '_')}
#SBATCH --partition={cfg.partition}
#SBATCH --nodes=1
#SBATCH --cpus-per-task={cfg.cpus}
#SBATCH --mem={cfg.memory}
#SBATCH --time={cfg.time_limit}
{gres_line}#SBATCH --output=ppigfinder_%j.out
#SBATCH --error=ppigfinder_%j.err

set -euo pipefail

echo "Running on: $(hostname)"
echo "Started at: $(date)"
echo "Workdir: {cfg.workdir}"

cd {cfg.workdir}

module purge

"""

    if tool == "AlphaFold 3":
        body = """# AlphaFold 3 example for DaVinci.
# Adjust module/image/database paths according to the local installation.

module use /home/public/davinci/etc/lmod/modules || true
module load alphafold3 || true

# Example using the ppigFinder/DaVinci AF3 helper:
# af3 --json-path input.json --job-name ppigfinder_af3 --resource-mode shared

echo "Replace this block with the configured AlphaFold 3 command."
"""
    elif tool == "GROMACS":
        body = """# GROMACS example.
# Adjust module name and input files according to your environment.

module load gromacs || true

# Example:
# gmx mdrun -s topol.tpr -deffnm md_run -ntmpi 1 -ntomp ${SLURM_CPUS_PER_TASK}

echo "Replace this block with the configured GROMACS command."
"""
    elif tool == "BLAST+":
        body = """# BLAST+ example.
# Use this for local protein/nucleotide search workflows.

module load blast || true

# Example:
# blastp -query query.faa -db database -out results.tsv -outfmt 6 -num_threads ${SLURM_CPUS_PER_TASK}

echo "Replace this block with the configured BLAST+ command."
"""
    elif tool == "HMMER3":
        body = """# HMMER3 example.
# Use this for domain annotation workflows.

module load hmmer || true

# Example:
# hmmsearch --cpu ${SLURM_CPUS_PER_TASK} --tblout hmm_hits.tbl profiles.hmm proteins.faa

echo "Replace this block with the configured HMMER3 command."
"""
    elif tool == "CryoSPARC":
        body = """# CryoSPARC note.
# CryoSPARC usually submits jobs through its own scheduler integration.
# This template is only a placeholder for cluster-side helper tasks.

echo "Use CryoSPARC's configured lane/cluster submission when available."
"""
    elif tool == "Python / Custom":
        body = """# Python/custom workflow example.

module load conda || true
# source activate your_environment

# Example:
# python your_script.py --input input.dat --output results/

echo "Replace this block with your custom command."
"""
    else:
        body = """# Custom command block.

echo "Add your software command here."
"""

    footer = """
echo "Finished at: $(date)"
"""

    return header + body + footer


class HPCDialog(QDialog):
    """
    Optional server/HPC execution center.
    """

    def __init__(self, parent=None):
        super().__init__(parent)

        self.status: HPCConnectionStatus | None = None
        self.config = HPCConnectionConfig(
            user=_default_user(),
            ssh_key_path=_default_key(),
        )

        self.setWindowTitle("DaVinci / HPC Execution Center")
        self.setWindowFlags(_window_flags())
        self.resize(1120, 820)
        self.setMinimumSize(920, 680)
        self.setSizeGripEnabled(True)
        apply_ppigfinder_branding(self)

        root = QVBoxLayout(self)
        root.setContentsMargins(18, 18, 18, 18)
        root.setSpacing(12)

        title = QLabel("DaVinci / HPC Execution Center")
        title.setObjectName("SectionTitle")
        root.addWidget(title)

        subtitle = QLabel(
            "Optional execution layer for users with server access. "
            "If no server is available, ppigFinder can still export AF3 JSON files, HTML reports, tables and vector figures for external use."
        )
        subtitle.setWordWrap(True)
        root.addWidget(subtitle)

        self.tabs = QTabWidget()
        root.addWidget(self.tabs, 1)

        self._build_connection_tab()
        self._build_software_tab()
        self._build_template_tab()
        self._build_outputs_tab()

        buttons = QHBoxLayout()

        btn_test = QPushButton("Test connection / local mode")
        btn_test.clicked.connect(self._test_connection)
        buttons.addWidget(btn_test)

        btn_update = QPushButton("Update template")
        btn_update.clicked.connect(self._update_template)
        buttons.addWidget(btn_update)

        btn_export = QPushButton("Export Slurm script")
        btn_export.clicked.connect(self._export_template)
        buttons.addWidget(btn_export)

        buttons.addStretch(1)

        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.accept)
        buttons.addWidget(btn_close)

        root.addLayout(buttons)

        self._update_template()

    # --------------------------------------------------------
    # Tabs
    # --------------------------------------------------------

    def _build_connection_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)

        box = QGroupBox("Connection profile")
        grid = QGridLayout(box)

        self.profile_edit = QLineEdit(self.config.profile)
        self.host_edit = QLineEdit(self.config.host)
        self.user_edit = QLineEdit(self.config.user)
        self.port_spin = QSpinBox()
        self.port_spin.setRange(1, 65535)
        self.port_spin.setValue(self.config.port)
        self.key_edit = QLineEdit(self.config.ssh_key_path)

        btn_key = QPushButton("Browse")
        btn_key.clicked.connect(self._browse_key)

        grid.addWidget(QLabel("Profile"), 0, 0)
        grid.addWidget(self.profile_edit, 0, 1, 1, 2)

        grid.addWidget(QLabel("Host"), 1, 0)
        grid.addWidget(self.host_edit, 1, 1, 1, 2)

        grid.addWidget(QLabel("User"), 2, 0)
        grid.addWidget(self.user_edit, 2, 1, 1, 2)

        grid.addWidget(QLabel("Port"), 3, 0)
        grid.addWidget(self.port_spin, 3, 1, 1, 2)

        grid.addWidget(QLabel("SSH key path"), 4, 0)
        grid.addWidget(self.key_edit, 4, 1)
        grid.addWidget(btn_key, 4, 2)

        layout.addWidget(box)

        self.status_text = QTextEdit()
        self.status_text.setReadOnly(True)
        self.status_text.setPlainText(
            "Connection not tested yet.\n\n"
            "Local DaVinci/GN mode is detected when the GUI is running directly on a cluster node. "
            "Remote SSH mode uses BatchMode=yes to avoid password prompts freezing the GUI."
        )
        layout.addWidget(self.status_text, 1)

        self.tabs.addTab(tab, "Connection")

    def _build_software_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)

        box = QGroupBox("Execution target")
        grid = QGridLayout(box)

        self.scheduler_combo = QComboBox()
        self.scheduler_combo.addItems(["Slurm", "PBS/Torque", "LSF", "External/manual"])

        self.software_combo = QComboBox()
        self.software_combo.addItems(
            [
                "AlphaFold 3",
                "GROMACS",
                "BLAST+",
                "HMMER3",
                "CryoSPARC",
                "Python / Custom",
            ]
        )
        self.software_combo.currentIndexChanged.connect(self._update_template)

        self.partition_edit = QLineEdit(self.config.partition)

        self.cpus_spin = QSpinBox()
        self.cpus_spin.setRange(1, 512)
        self.cpus_spin.setValue(self.config.cpus)

        self.memory_edit = QLineEdit(self.config.memory)
        self.time_edit = QLineEdit(self.config.time_limit)
        self.gres_edit = QLineEdit(self.config.gres)
        self.workdir_edit = QLineEdit(self.config.workdir)

        for widget in (
            self.scheduler_combo,
            self.partition_edit,
            self.cpus_spin,
            self.memory_edit,
            self.time_edit,
            self.gres_edit,
            self.workdir_edit,
        ):
            try:
                if hasattr(widget, "textChanged"):
                    widget.textChanged.connect(self._update_template)
                elif hasattr(widget, "valueChanged"):
                    widget.valueChanged.connect(self._update_template)
            except Exception:
                pass

        self.scheduler_combo.currentIndexChanged.connect(self._update_template)

        grid.addWidget(QLabel("Scheduler"), 0, 0)
        grid.addWidget(self.scheduler_combo, 0, 1)

        grid.addWidget(QLabel("Software/workflow"), 1, 0)
        grid.addWidget(self.software_combo, 1, 1)

        grid.addWidget(QLabel("Partition/queue"), 2, 0)
        grid.addWidget(self.partition_edit, 2, 1)

        grid.addWidget(QLabel("CPUs"), 3, 0)
        grid.addWidget(self.cpus_spin, 3, 1)

        grid.addWidget(QLabel("Memory"), 4, 0)
        grid.addWidget(self.memory_edit, 4, 1)

        grid.addWidget(QLabel("Time"), 5, 0)
        grid.addWidget(self.time_edit, 5, 1)

        grid.addWidget(QLabel("GRES/GPU"), 6, 0)
        grid.addWidget(self.gres_edit, 6, 1)

        grid.addWidget(QLabel("Workdir"), 7, 0)
        grid.addWidget(self.workdir_edit, 7, 1)

        layout.addWidget(box)

        note = QLabel(
            "This module is optional. It prepares editable templates for external execution. "
            "For DaVinci, examples include AF3, GROMACS, BLAST+, HMMER3 and custom scripts. "
            "CryoSPARC is shown as a supported workflow context, but usually submits through its own scheduler integration."
        )
        note.setWordWrap(True)
        layout.addWidget(note)

        layout.addStretch(1)

        self.tabs.addTab(tab, "Software / Slurm")

    def _build_template_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)

        label = QLabel("Editable scheduler template")
        label.setObjectName("SectionSubTitle")
        layout.addWidget(label)

        self.template_text = QTextEdit()
        self.template_text.setAcceptRichText(False)
        layout.addWidget(self.template_text, 1)

        self.tabs.addTab(tab, "Template")

    def _build_outputs_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)

        text = QTextEdit()
        text.setReadOnly(True)
        text.setPlainText(
            "Server execution is optional.\n\n"
            "Recommended local/export outputs when no server is available:\n\n"
            "1. AlphaFold 3 Server JSON\n"
            "   - Can be uploaded to AlphaFold Server or used by a local AF3 installation.\n\n"
            "2. HTML report\n"
            "   - Should become an interactive report with filters, sortable tables and embedded figures.\n\n"
            "3. Vector figures\n"
            "   - SVG/PDF should be preferred for genome maps, neighbourhood diagrams and workflow figures.\n"
            "   - Users can open these files in Inkscape, Illustrator, Affinity Designer or similar tools.\n\n"
            "4. Tables\n"
            "   - TSV/CSV for ORFs, BLAST hits, HMM hits, neighbourhood candidates and AF3 metrics.\n\n"
            "5. Project Snapshot\n"
            "   - Portable JSON storing the analysis state for reproducibility."
        )
        layout.addWidget(text, 1)

        self.tabs.addTab(tab, "Outputs without server")

    # --------------------------------------------------------
    # Helpers
    # --------------------------------------------------------

    def _browse_key(self):
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select SSH key",
            str(Path.home() / ".ssh"),
            "SSH keys (*)",
        )
        if path:
            self.key_edit.setText(path)

    def _read_config(self) -> HPCConnectionConfig:
        return HPCConnectionConfig(
            profile=self.profile_edit.text().strip() or "HPC",
            host=self.host_edit.text().strip(),
            user=self.user_edit.text().strip(),
            port=int(self.port_spin.value()),
            ssh_key_path=self.key_edit.text().strip(),
            scheduler=self.scheduler_combo.currentText(),
            partition=self.partition_edit.text().strip() or "max50",
            cpus=int(self.cpus_spin.value()),
            memory=self.memory_edit.text().strip() or "64G",
            time_limit=self.time_edit.text().strip() or "24:00:00",
            gres=self.gres_edit.text().strip(),
            workdir=self.workdir_edit.text().strip() or "$PWD",
        )

    def _update_template(self):
        self.config = self._read_config() if hasattr(self, "profile_edit") else self.config

        tool = self.software_combo.currentText() if hasattr(self, "software_combo") else "AlphaFold 3"

        if self.config.scheduler != "Slurm":
            self.template_text.setPlainText(
                f"# {self.config.scheduler} selected.\n"
                "# Slurm template generation is currently the primary supported mode.\n"
                "# Use this tab as an editable manual script area.\n\n"
                "# Add your external execution commands here.\n"
            )
            return

        self.template_text.setPlainText(_template_for(tool, self.config))

    def _test_connection(self):
        self.config = self._read_config()

        local_cluster, slurm_available, hostname = _detect_local_cluster()

        if local_cluster or slurm_available:
            self.status = HPCConnectionStatus(
                connection_ok=True,
                running_on_cluster=True,
                scheduler_available=slurm_available,
                hostname=hostname,
                message="Local cluster mode detected.",
            )
            self.status_text.setPlainText(
                f"Local cluster mode detected.\n\n"
                f"Hostname: {hostname}\n"
                f"Slurm available: {'yes' if slurm_available else 'no'}\n\n"
                f"You can export or adapt Slurm scripts directly on this node."
            )
            return

        cmd = [
            "ssh",
            "-o",
            "BatchMode=yes",
            "-o",
            "ConnectTimeout=5",
            "-p",
            str(self.config.port),
        ]

        if self.config.ssh_key_path:
            cmd.extend(["-i", self.config.ssh_key_path])

        target = f"{self.config.user}@{self.config.host}" if self.config.user else self.config.host
        cmd.extend([target, "hostname && command -v sbatch || true"])

        try:
            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=8,
            )

            ok = result.returncode == 0
            scheduler_available = "sbatch" in result.stdout or "/sbatch" in result.stdout

            self.status = HPCConnectionStatus(
                connection_ok=ok,
                running_on_cluster=False,
                scheduler_available=scheduler_available,
                hostname=result.stdout.strip().splitlines()[0] if result.stdout.strip() else "",
                message="SSH connection OK." if ok else "SSH connection failed.",
            )

            self.status_text.setPlainText(
                f"SSH command:\n{' '.join(cmd)}\n\n"
                f"Return code: {result.returncode}\n\n"
                f"STDOUT:\n{result.stdout}\n\n"
                f"STDERR:\n{result.stderr}\n"
            )

        except Exception as exc:
            self.status = HPCConnectionStatus(
                connection_ok=False,
                running_on_cluster=False,
                scheduler_available=False,
                hostname="",
                message=str(exc),
            )
            self.status_text.setPlainText(f"Connection test failed:\n\n{exc}")

    def _export_template(self):
        self.config = self._read_config()

        default_name = f"ppigfinder_{self.software_combo.currentText().lower().replace(' ', '_').replace('/', '_')}.slurm"

        path, _ = QFileDialog.getSaveFileName(
            self,
            "Export scheduler script",
            default_name,
            "Slurm scripts (*.slurm *.sh);;Shell scripts (*.sh);;All files (*)",
        )

        if not path:
            return

        Path(path).write_text(self.template_text.toPlainText(), encoding="utf-8")
        QMessageBox.information(self, "Exported", f"Scheduler script exported:\n\n{path}")


def open_hpc_connection_dialog(parent=None):
    dialog = HPCDialog(parent=parent)
    dialog.exec() if hasattr(dialog, "exec") else dialog.exec_()

    return dialog.status, dialog.config
