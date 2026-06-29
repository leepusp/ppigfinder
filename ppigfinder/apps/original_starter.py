#!/usr/bin/env python3
from __future__ import annotations

import os
import subprocess
import sys
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

    It keeps startup visually responsive and starts the full original interface
    automatically, without requiring an intermediate button click.
    """

    def __init__(self):
        super().__init__()

        self.child_process: subprocess.Popen | None = None

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

        QTimer.singleShot(350, self.open_original_interface)

    def open_original_interface(self):
        env = os.environ.copy()

        self.status_label.setText("Starting full interface...")

        self.child_process = subprocess.Popen(
            [sys.executable, "main.py"],
            cwd=str(REPO_ROOT),
            env=env,
        )

        QTimer.singleShot(1400, self.close)


def main() -> int:
    app = QApplication(sys.argv)
    app.setApplicationName("ppigFinder")
    app.setApplicationDisplayName("ppigFinder — Loading")

    window = PpigFinderStarter()
    window.show()

    return app.exec() if QT6 else app.exec_()


if __name__ == "__main__":
    raise SystemExit(main())
