#!/usr/bin/env python3
"""
Qt compatibility helpers for the future ppigFinder UI shell.
"""

from __future__ import annotations

try:
    from PyQt6.QtCore import Qt, QSize, QTimer
    from PyQt6.QtGui import QFont, QPixmap, QColor
    from PyQt6.QtWidgets import (
        QApplication,
        QMainWindow,
        QWidget,
        QVBoxLayout,
        QHBoxLayout,
        QGridLayout,
        QLabel,
        QPushButton,
        QFrame,
        QStackedWidget,
        QListWidget,
        QListWidgetItem,
        QMessageBox,
        QProgressBar,
        QSizePolicy,
        QFileDialog,
    )

    QT6 = True

except Exception:
    from PyQt5.QtCore import Qt, QSize, QTimer
    from PyQt5.QtGui import QFont, QPixmap, QColor
    from PyQt5.QtWidgets import (
        QApplication,
        QMainWindow,
        QWidget,
        QVBoxLayout,
        QHBoxLayout,
        QGridLayout,
        QLabel,
        QPushButton,
        QFrame,
        QStackedWidget,
        QListWidget,
        QListWidgetItem,
        QMessageBox,
        QProgressBar,
        QSizePolicy,
        QFileDialog,
    )

    QT6 = False


def align_center():
    return Qt.AlignmentFlag.AlignCenter if QT6 else Qt.AlignCenter


def align_left():
    return Qt.AlignmentFlag.AlignLeft if QT6 else Qt.AlignLeft


def pointing_hand_cursor():
    return Qt.CursorShape.PointingHandCursor if QT6 else Qt.PointingHandCursor


def expanding_policy():
    return QSizePolicy.Policy.Expanding if QT6 else QSizePolicy.Expanding


def preferred_policy():
    return QSizePolicy.Policy.Preferred if QT6 else QSizePolicy.Preferred


def exec_app(app: QApplication) -> int:
    return app.exec() if QT6 else app.exec_()
