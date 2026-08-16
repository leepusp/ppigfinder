#!/usr/bin/env python3
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

try:
    from PyQt6.QtCore import Qt, QTimer
    from PyQt6.QtWidgets import QApplication, QLabel, QProgressBar, QVBoxLayout, QWidget
    QT6 = True
except Exception:
    from PyQt5.QtCore import Qt, QTimer
    from PyQt5.QtWidgets import QApplication, QLabel, QProgressBar, QVBoxLayout, QWidget
    QT6 = False


REPO_ROOT = Path(__file__).resolve().parents[2]


def _align_center():
    return Qt.AlignmentFlag.AlignCenter if QT6 else Qt.AlignCenter


class PpigFinderStarter(QWidget):
    """
    Lightweight splash launcher for the original ppigFinder interface.

    The splash stays visible until the real interface signals that its main
    window has been shown. This avoids a visual gap between "loading finished"
    and the actual GUI opening.
    """

    def __init__(self):
        super().__init__()

        self.child_process: subprocess.Popen | None = None
        self.started_at = time.monotonic()
        self.ready_file = Path(tempfile.gettempdir()) / f"ppigfinder_ready_{os.getpid()}.txt"

        try:
            self.ready_file.unlink()
        except FileNotFoundError:
            pass
        except Exception:
            pass

        self.setWindowTitle("ppigFinder — Loading")
        self.resize(560, 360)

        try:
            from ppigfinder.ui.icon_provider import set_window_icon
            set_window_icon(self)
        except Exception:
            pass

        layout = QVBoxLayout(self)
        layout.setContentsMargins(34, 34, 34, 34)
        layout.setSpacing(14)

        icon = QLabel("P")
        icon.setAlignment(_align_center())
        icon.setStyleSheet(
            """
            QLabel {
                background: #17384d;
                color: white;
                border-radius: 18px;
                font-size: 48px;
                font-weight: 800;
                min-width: 96px;
                min-height: 96px;
                max-width: 96px;
                max-height: 96px;
            }
            """
        )

        title = QLabel("ppigFinder")
        title.setAlignment(_align_center())
        title.setStyleSheet("font-size: 30px; font-weight: 800; color: #17384d;")

        subtitle = QLabel("Protein-Protein Interaction Genomic Finder")
        subtitle.setAlignment(_align_center())
        subtitle.setStyleSheet("font-size: 12px; color: #60717f;")

        self.status_label = QLabel("Loading original ppigFinder interface...")
        self.status_label.setAlignment(_align_center())
        self.status_label.setWordWrap(True)
        self.status_label.setStyleSheet("font-size: 12px; color: #263238;")

        self.progress = QProgressBar()
        self.progress.setRange(0, 0)
        self.progress.setTextVisible(False)
        self.progress.setMaximumHeight(8)

        layout.addStretch(1)
        layout.addWidget(icon, alignment=_align_center())
        layout.addWidget(title)
        layout.addWidget(subtitle)
        layout.addSpacing(8)
        layout.addWidget(self.status_label)
        layout.addWidget(self.progress)
        layout.addStretch(1)

        QTimer.singleShot(250, self.open_original_interface)

        self.poll_timer = QTimer(self)
        self.poll_timer.timeout.connect(self._poll_ready)
        self.poll_timer.start(150)

    def open_original_interface(self):
        env = os.environ.copy()
        env["PPIG_READY_FILE"] = str(self.ready_file)

        self.status_label.setText("Starting full interface...")

        self.child_process = subprocess.Popen(
            [sys.executable, "main.py"],
            cwd=str(REPO_ROOT),
            env=env,
        )

    def _poll_ready(self):
        if self.ready_file.exists():
            self.status_label.setText("Interface ready.")
            self.poll_timer.stop()
            QTimer.singleShot(250, self.close)
            return

        if self.child_process is not None and self.child_process.poll() is not None:
            self.status_label.setText(
                "The original interface closed before reporting readiness. "
                "Check the terminal for errors."
            )
            self.poll_timer.stop()
            QTimer.singleShot(2500, self.close)
            return

        elapsed = time.monotonic() - self.started_at
        if elapsed > 45:
            self.status_label.setText("Still loading the original interface...")
            # Keep waiting; large HPC filesystems can be slow.


def main() -> int:
    app = QApplication(sys.argv)
    app.setApplicationName("ppigFinder")
    app.setApplicationDisplayName("ppigFinder — Loading")

    window = PpigFinderStarter()
    window.show()

    return app.exec() if QT6 else app.exec_()


if __name__ == "__main__":
    raise SystemExit(main())
