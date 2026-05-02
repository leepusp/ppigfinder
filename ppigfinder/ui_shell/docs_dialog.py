#!/usr/bin/env python3
"""
Documentation dialog for the experimental UI shell.
"""

from __future__ import annotations

from html import escape
import re

try:
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import QDialog, QVBoxLayout, QTextBrowser
    QT6 = True
except Exception:
    from PyQt5.QtCore import Qt
    from PyQt5.QtWidgets import QDialog, QVBoxLayout, QTextBrowser
    QT6 = False

from ppigfinder.ui_shell.doc_loader import load_module_markdown
from ppigfinder.ui_shell.branding import apply_ppigfinder_branding


def _window_flags():
    flags = Qt.WindowType.Window if QT6 else Qt.Window
    flags |= Qt.WindowType.WindowMinimizeButtonHint if QT6 else Qt.WindowMinimizeButtonHint
    flags |= Qt.WindowType.WindowMaximizeButtonHint if QT6 else Qt.WindowMaximizeButtonHint
    flags |= Qt.WindowType.WindowCloseButtonHint if QT6 else Qt.WindowCloseButtonHint
    return flags


def _inline_markup(text: str) -> str:
    text = escape(text)
    text = re.sub(r"`([^`]+)`", r"<code>\1</code>", text)
    text = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", text)
    return text


def _basic_markdown_to_html(text: str) -> str:
    """
    Minimal but structured Markdown-to-HTML converter for local documentation.
    Supports headings, bullet lists and fenced code blocks.
    """
    html_lines = []
    in_list = False
    in_code = False
    code_lines = []

    def close_list():
        nonlocal in_list
        if in_list:
            html_lines.append("</ul>")
            in_list = False

    def close_code():
        nonlocal in_code, code_lines
        if in_code:
            html_lines.append("<pre><code>")
            html_lines.append(escape("\n".join(code_lines)))
            html_lines.append("</code></pre>")
            code_lines = []
            in_code = False

    for raw in text.splitlines():
        line = raw.rstrip()

        if line.startswith("```"):
            if in_code:
                close_code()
            else:
                close_list()
                in_code = True
                code_lines = []
            continue

        if in_code:
            code_lines.append(line)
            continue

        if line.startswith("# "):
            close_list()
            html_lines.append(f"<h1>{_inline_markup(line[2:].strip())}</h1>")
        elif line.startswith("## "):
            close_list()
            html_lines.append(f"<h2>{_inline_markup(line[3:].strip())}</h2>")
        elif line.startswith("### "):
            close_list()
            html_lines.append(f"<h3>{_inline_markup(line[4:].strip())}</h3>")
        elif line.startswith("- "):
            if not in_list:
                html_lines.append("<ul>")
                in_list = True
            html_lines.append(f"<li>{_inline_markup(line[2:].strip())}</li>")
        elif not line.strip():
            close_list()
        else:
            close_list()
            html_lines.append(f"<p>{_inline_markup(line)}</p>")

    close_code()
    close_list()

    return """
<html>
<head>
<style>
body {
    font-family: Arial, Helvetica, sans-serif;
    color: #263238;
    background: #ffffff;
    font-size: 13px;
    margin: 22px;
}
h1 {
    color: #1b3a4b;
    font-size: 26px;
    margin-bottom: 18px;
}
h2 {
    color: #1b3a4b;
    font-size: 18px;
    margin-top: 22px;
    margin-bottom: 8px;
    border-bottom: 1px solid #d7dee2;
    padding-bottom: 4px;
}
h3 {
    color: #263238;
    font-size: 15px;
    margin-top: 14px;
    margin-bottom: 5px;
}
p {
    line-height: 1.45;
    margin-top: 6px;
    margin-bottom: 8px;
}
ul {
    margin-top: 6px;
    margin-bottom: 10px;
}
li {
    margin-bottom: 5px;
}
code {
    background: #eef3f5;
    color: #1b3a4b;
    padding: 2px 4px;
    border-radius: 4px;
}
pre {
    background: #eef3f5;
    border: 1px solid #d7dee2;
    border-radius: 8px;
    padding: 10px;
    overflow-x: auto;
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
        self.setWindowFlags(_window_flags())
        self.resize(980, 760)
        self.setMinimumSize(760, 560)
        self.setSizeGripEnabled(True)
        apply_ppigfinder_branding(self)

        layout = QVBoxLayout(self)

        browser = QTextBrowser(self)
        browser.setOpenExternalLinks(True)
        markdown = load_module_markdown(module_id)
        browser.setHtml(_basic_markdown_to_html(markdown))
        layout.addWidget(browser)


def show_module_documentation(module_id: str, title: str, parent=None) -> None:
    dialog = ModuleDocumentationDialog(module_id, title, parent=parent)
    dialog.showMaximized()
    dialog.exec() if hasattr(dialog, "exec") else dialog.exec_()
