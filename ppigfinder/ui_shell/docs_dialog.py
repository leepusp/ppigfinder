#!/usr/bin/env python3
"""
Documentation dialog for the experimental UI shell.
"""

from __future__ import annotations

from html import escape
import re

try:
    from PyQt6.QtWidgets import QDialog, QVBoxLayout, QTextBrowser, QPushButton
except Exception:
    from PyQt5.QtWidgets import QDialog, QVBoxLayout, QTextBrowser, QPushButton

from ppigfinder.ui_shell.doc_loader import load_module_markdown


def _basic_markdown_to_html(text: str) -> str:
    """
    Minimal Markdown-to-HTML converter for local documentation.
    """
    html_lines = []

    for raw in text.splitlines():
        line = raw.rstrip()

        if line.startswith("# "):
            html_lines.append(f"<h1>{escape(line[2:].strip())}</h1>")
        elif line.startswith("## "):
            html_lines.append(f"<h2>{escape(line[3:].strip())}</h2>")
        elif line.startswith("### "):
            html_lines.append(f"<h3>{escape(line[4:].strip())}</h3>")
        elif line.startswith("- "):
            html_lines.append(f"<li>{escape(line[2:].strip())}</li>")
        elif not line.strip():
            html_lines.append("<br>")
        else:
            escaped = escape(line)
            escaped = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", escaped)
            html_lines.append(f"<p>{escaped}</p>")

    return """
<html>
<head>
<style>
body {
    font-family: Arial, Helvetica, sans-serif;
    color: #263238;
    background: #ffffff;
}
h1, h2, h3 {
    color: #1b3a4b;
}
li {
    margin-bottom: 4px;
}
p {
    line-height: 1.35;
}
</style>
</head>
<body>
""" + "\n".join(html_lines) + """
</body>
</html>
"""


class ModuleDocumentationDialog(QDialog):
    """
    Dialog showing module documentation.
    """

    def __init__(self, module_id: str, title: str, parent=None):
        super().__init__(parent)

        self.setWindowTitle(f"{title} — Guide")
        self.resize(760, 620)

        layout = QVBoxLayout(self)

        browser = QTextBrowser(self)
        markdown = load_module_markdown(module_id)
        browser.setHtml(_basic_markdown_to_html(markdown))
        layout.addWidget(browser)

        close_button = QPushButton("Close")
        close_button.clicked.connect(self.accept)
        layout.addWidget(close_button)


def show_module_documentation(module_id: str, title: str, parent=None) -> None:
    dialog = ModuleDocumentationDialog(module_id, title, parent=parent)
    dialog.exec() if hasattr(dialog, "exec") else dialog.exec_()
