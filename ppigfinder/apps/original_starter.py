#!/usr/bin/env python3
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

try:
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import QApplication, QLabel, QPushButton, QVBoxLayout, QWidget
    QT6 = True
except Exception:
    from PyQt5.QtCore import Qt
    from PyQt5.QtWidgets import QApplication, QLabel, QPushButton, QVBoxLayout, QWidget
    QT6 = False


REPO_ROOT = Path(__file__).resolve().parents[2]


def _align_center():
    return Qt.AlignmentFlag.AlignCenter if QT6 else Qt.AlignCenter


class PpigFinderStarter(QWidget):
    """
    Lightweight ppigFinder starter.

    This keeps startup fast because it does not instantiate the full legacy
    interface until the user explicitly opens it.
    """

    def __init__(self):
        super().__init__()

        self.setWindowTitle("ppigFinder — Start")
        self.resize(520, 340)

        try:
            from ppigfinder.ui.icon_provider import set_window_icon
            set_window_icon(self)
        except Exception:
            pass

        layout = QVBoxLayout(self)
        layout.setContentsMargins(28, 28, 28, 28)
        layout.setSpacing(16)

        icon = QLabel("P")
        icon.setAlignment(_align_center())
        icon.setStyleSheet(
            """
            QLabel {
                background: #17384d;
                color: white;
                border-radius: 18px;
                font-size: 44px;
                font-weight: 800;
                min-width: 86px;
                min-height: 86px;
                max-width: 86px;
                max-height: 86px;
            }
            """
        )

        title = QLabel("ppigFinder")
        title.setAlignment(_align_center())
        title.setStyleSheet("font-size: 28px; font-weight: 800; color: #17384d;")

        subtitle = QLabel("Protein-Protein Interaction Genomic Finder")
        subtitle.setAlignment(_align_center())
        subtitle.setStyleSheet("font-size: 12px; color: #60717f;")

        description = QLabel(
            "Fast launcher for the original ppigFinder interface. "
            "The full application is loaded only after clicking Open."
        )
        description.setWordWrap(True)
        description.setAlignment(_align_center())

        open_button = QPushButton("Open original ppigFinder interface")
        open_button.clicked.connect(self.open_original_interface)

        layout.addWidget(icon, alignment=_align_center())
        layout.addWidget(title)
        layout.addWidget(subtitle)
        layout.addWidget(description)
        layout.addStretch(1)
        layout.addWidget(open_button)

    def open_original_interface(self):
        env = os.environ.copy()
        subprocess.Popen(
            [sys.executable, "main.py"],
            cwd=str(REPO_ROOT),
            env=env,
        )
        self.close()


def main() -> int:
    app = QApplication(sys.argv)
    app.setApplicationName("ppigFinder")
    app.setApplicationDisplayName("ppigFinder — Start")

    window = PpigFinderStarter()
    window.show()

    return app.exec() if QT6 else app.exec_()


if __name__ == "__main__":
    raise SystemExit(main())
