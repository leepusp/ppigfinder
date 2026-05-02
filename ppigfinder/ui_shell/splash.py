#!/usr/bin/env python3
"""
Splash/loading screen for ppigFinder.
"""

from __future__ import annotations

from ppigfinder.ui_shell.qt import (
    QWidget,
    QVBoxLayout,
    QLabel,
    QProgressBar,
    QTimer,
    align_center,
)
from ppigfinder.ui_shell.theme import APP_TITLE, APP_SUBTITLE
from ppigfinder.ui_shell.branding import create_ppigfinder_pixmap, apply_ppigfinder_branding


class SplashWindow(QWidget):
    """
    Lightweight startup splash window.
    """

    def __init__(self):
        super().__init__()

        self.setWindowTitle(APP_TITLE)
        apply_ppigfinder_branding(self)
        self.resize(520, 280)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(36, 36, 36, 36)
        layout.setSpacing(14)

        logo = QLabel()
        logo.setPixmap(create_ppigfinder_pixmap(72))
        logo.setAlignment(align_center())
        layout.addWidget(logo)

        title = QLabel(APP_TITLE)
        title.setObjectName("HeroTitle")
        title.setAlignment(align_center())
        layout.addWidget(title)

        subtitle = QLabel(APP_SUBTITLE)
        subtitle.setObjectName("HeroSubtitle")
        subtitle.setAlignment(align_center())
        layout.addWidget(subtitle)

        self.message = QLabel("Preparing analysis environment...")
        self.message.setAlignment(align_center())
        layout.addWidget(self.message)

        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        self.progress.setValue(0)
        layout.addWidget(self.progress)

        self._value = 0
        self._timer = QTimer(self)
        self._timer.timeout.connect(self._advance)

    def start(self) -> None:
        self._timer.start(60)

    def _advance(self) -> None:
        self._value = min(100, self._value + 4)
        self.progress.setValue(self._value)

        if self._value < 35:
            self.message.setText("Checking backend tools...")
        elif self._value < 70:
            self.message.setText("Preparing genome and protein modules...")
        elif self._value < 100:
            self.message.setText("Loading interface...")
        else:
            self.message.setText("Ready")
            self._timer.stop()
