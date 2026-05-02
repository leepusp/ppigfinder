#!/usr/bin/env python3
"""
Visual theme for the future ppigFinder UI shell.
"""

from __future__ import annotations


APP_TITLE = "ppigFinder"
APP_SUBTITLE = "Protein-Protein Interaction Genomic Finder"


def shell_stylesheet() -> str:
    """
    Central stylesheet for the new UI shell.
    """
    return """
QMainWindow {
    background: #f5f7f8;
}

QWidget {
    font-family: Arial, Helvetica, sans-serif;
    font-size: 12px;
    color: #263238;
}

QLabel#HeroTitle {
    font-size: 28px;
    font-weight: 700;
    color: #1b3a4b;
}

QLabel#HeroSubtitle {
    font-size: 14px;
    color: #607d8b;
}

QLabel#SectionTitle {
    font-size: 16px;
    font-weight: 700;
    color: #1b3a4b;
}

QFrame#Card {
    background: #ffffff;
    border: 1px solid #d7dee2;
    border-radius: 12px;
}

QFrame#Card:hover {
    border: 1px solid #90a4ae;
}

QFrame#VisualizationPlaceholder {
    background: #f8fbfc;
    border: 1px dashed #90a4ae;
    border-radius: 12px;
}

QLabel#CardTitle {
    font-size: 14px;
    font-weight: 700;
    color: #263238;
}

QLabel#CardDescription {
    color: #607d8b;
}

QLabel#MetricValue {
    font-size: 18px;
    font-weight: 700;
    color: #1b3a4b;
}

QLabel#FlowStep {
    background: #e8f0f4;
    color: #455a64;
    border-radius: 10px;
    padding: 8px 10px;
}

QLabel#FlowStepActive {
    background: #1b3a4b;
    color: #ffffff;
    border-radius: 10px;
    padding: 8px 10px;
    font-weight: 700;
}

QLabel#FlowArrow {
    color: #607d8b;
    font-weight: 700;
    padding: 8px 0;
}

QPushButton {
    background: #1b3a4b;
    color: white;
    border-radius: 8px;
    padding: 8px 12px;
}

QPushButton:hover {
    background: #24536a;
}

QPushButton#SecondaryButton {
    background: #e8f0f4;
    color: #1b3a4b;
}

QPushButton#SecondaryButton:hover {
    background: #d7e6ee;
}

QPushButton:disabled {
    background: #b0bec5;
    color: #eceff1;
}

QListWidget {
    background: #ffffff;
    border: 1px solid #d7dee2;
    border-radius: 10px;
    padding: 6px;
}

QListWidget::item {
    padding: 10px;
    border-radius: 8px;
}

QListWidget::item:selected {
    background: #d7e6ee;
    color: #1b3a4b;
}
"""
