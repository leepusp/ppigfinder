#!/usr/bin/env python3
"""
DaVinci/HPC setup dialog for the experimental guided UI shell.
"""

from __future__ import annotations

try:
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import (
        QDialog,
        QVBoxLayout,
        QHBoxLayout,
        QFormLayout,
        QLabel,
        QLineEdit,
        QSpinBox,
        QPushButton,
        QTextEdit,
    )
    QT6 = True
except Exception:
    from PyQt5.QtCore import Qt
    from PyQt5.QtWidgets import (
        QDialog,
        QVBoxLayout,
        QHBoxLayout,
        QFormLayout,
        QLabel,
        QLineEdit,
        QSpinBox,
        QPushButton,
        QTextEdit,
    )
    QT6 = False

from ppigfinder.ui_shell.branding import apply_ppigfinder_branding
from ppigfinder.ui_shell.hpc_connection import (
    HPCConnectionConfig,
    default_config,
    test_connection,
    slurm_template_for_af3,
)


def _window_flags():
    flags = Qt.WindowType.Window if QT6 else Qt.Window
    flags |= Qt.WindowType.WindowMinimizeButtonHint if QT6 else Qt.WindowMinimizeButtonHint
    flags |= Qt.WindowType.WindowMaximizeButtonHint if QT6 else Qt.WindowMaximizeButtonHint
    flags |= Qt.WindowType.WindowCloseButtonHint if QT6 else Qt.WindowCloseButtonHint
    return flags


class HPCConnectionDialog(QDialog):
    """
    Dialog for configuring and testing DaVinci/HPC connectivity.
    """

    def __init__(self, initial_config: HPCConnectionConfig | None = None, parent=None):
        super().__init__(parent)

        self.config = initial_config or default_config()
        self.status = None

        self.setWindowTitle("DaVinci / HPC Connection")
        self.setWindowFlags(_window_flags())
        self.resize(900, 680)
        self.setMinimumSize(760, 520)
        self.setSizeGripEnabled(True)
        apply_ppigfinder_branding(self)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(18, 18, 18, 18)
        layout.setSpacing(12)

        title = QLabel("DaVinci / HPC connection")
        title.setObjectName("SectionTitle")
        layout.addWidget(title)

        description = QLabel(
            "Configure optional access to DaVinci/HPC workflows. "
            "The test uses SSH BatchMode to avoid freezing the GUI with password prompts. "
            "If you are already running on a DaVinci/GN node, the dialog detects local cluster mode."
        )
        description.setWordWrap(True)
        layout.addWidget(description)

        form = QFormLayout()

        self.profile_input = QLineEdit(self.config.profile)
        self.host_input = QLineEdit(self.config.host)
        self.user_input = QLineEdit(self.config.user)
        self.port_input = QSpinBox()
        self.port_input.setRange(1, 65535)
        self.port_input.setValue(int(self.config.port))
        self.key_input = QLineEdit(self.config.key_path)

        form.addRow("Profile", self.profile_input)
        form.addRow("Host", self.host_input)
        form.addRow("User", self.user_input)
        form.addRow("Port", self.port_input)
        form.addRow("SSH key path", self.key_input)

        layout.addLayout(form)

        buttons = QHBoxLayout()

        test_button = QPushButton("Test connection")
        test_button.clicked.connect(self._test_connection)
        buttons.addWidget(test_button)

        template_button = QPushButton("Show AF3 Slurm template")
        template_button.clicked.connect(self._show_slurm_template)
        buttons.addWidget(template_button)

        buttons.addStretch(1)
        layout.addLayout(buttons)

        self.output = QTextEdit()
        self.output.setReadOnly(True)
        self.output.setPlaceholderText("Connection status and generated templates will appear here.")
        layout.addWidget(self.output, 1)

    def current_config(self) -> HPCConnectionConfig:
        return HPCConnectionConfig(
            profile=self.profile_input.text().strip() or "DaVinci",
            host=self.host_input.text().strip() or "davinci.icb.usp.br",
            user=self.user_input.text().strip(),
            port=int(self.port_input.value()),
            key_path=self.key_input.text().strip(),
        )

    def _test_connection(self) -> None:
        self.config = self.current_config()
        self.status = test_connection(self.config)

        text = [
            "Connection test",
            "",
            f"Profile: {self.status.profile}",
            f"Host: {self.status.host}",
            f"User: {self.status.user}",
            f"Port: {self.status.port}",
            f"Local hostname: {self.status.local_hostname}",
            f"Running on cluster: {self.status.running_on_cluster}",
            f"SSH available: {self.status.ssh_available}",
            f"Connection OK: {self.status.connection_ok}",
            "",
            "Command preview:",
            self.status.command_preview,
            "",
            "Message:",
            self.status.message,
        ]

        self.output.setPlainText("\n".join(text))

    def _show_slurm_template(self) -> None:
        self.output.setPlainText(slurm_template_for_af3())


def open_hpc_connection_dialog(parent=None):
    dialog = HPCConnectionDialog(parent=parent)
    dialog.exec() if hasattr(dialog, "exec") else dialog.exec_()
    return dialog.status, dialog.current_config()
