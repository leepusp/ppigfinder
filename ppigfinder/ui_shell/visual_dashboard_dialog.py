#!/usr/bin/env python3
"""
Visual workflow dashboard for the experimental guided UI shell.

This is the first explicit illustration layer:
- workflow graph
- input/operation/output process diagram
- ORF map preview
- ORF length distribution

It uses Qt/QPainter only, so it does not depend on QWebEngine or internet access.
"""

from __future__ import annotations

from collections import Counter

try:
    from PyQt6.QtCore import Qt, QPointF, QRectF
    from PyQt6.QtGui import QPainter, QColor, QPen, QBrush, QFont, QPolygonF
    from PyQt6.QtWidgets import (
        QDialog,
        QVBoxLayout,
        QLabel,
        QWidget,
        QScrollArea,
    )
    QT6 = True
except Exception:
    from PyQt5.QtCore import Qt, QPointF, QRectF
    from PyQt5.QtGui import QPainter, QColor, QPen, QBrush, QFont, QPolygonF
    from PyQt5.QtWidgets import (
        QDialog,
        QVBoxLayout,
        QLabel,
        QWidget,
        QScrollArea,
    )
    QT6 = False

from ppigfinder.ui_shell.branding import apply_ppigfinder_branding


WORKFLOW_NODES = [
    ("data", "Input data"),
    ("genome", "Genome"),
    ("orfs", "ORFs"),
    ("annotation", "Annotation"),
    ("alphafold", "AF3 / PPI"),
    ("hpc", "DaVinci / HPC"),
    ("reports", "Reports"),
]


def _window_flags():
    flags = Qt.WindowType.Window if QT6 else Qt.Window
    flags |= Qt.WindowType.WindowMinimizeButtonHint if QT6 else Qt.WindowMinimizeButtonHint
    flags |= Qt.WindowType.WindowMaximizeButtonHint if QT6 else Qt.WindowMaximizeButtonHint
    flags |= Qt.WindowType.WindowCloseButtonHint if QT6 else Qt.WindowCloseButtonHint
    return flags


def _align_center():
    return Qt.AlignmentFlag.AlignCenter if QT6 else Qt.AlignCenter


def _safe_int(value, default=0):
    try:
        return int(value)
    except Exception:
        return default


def _state_get(state, key, default=None):
    if state is None:
        return default
    if hasattr(state, "get"):
        return state.get(key, default)
    if isinstance(state, dict):
        return state.get(key, default)
    return default


def _completed_steps(state) -> set[str]:
    if state is None:
        return set()
    if hasattr(state, "completed_steps"):
        try:
            return set(state.completed_steps())
        except Exception:
            return set()
    return set()


def _current_route(state) -> str:
    return getattr(state, "current_route", "data")


def _orf_length(orf) -> int:
    if hasattr(orf, "aa_length"):
        return _safe_int(getattr(orf, "aa_length"))
    return len(getattr(orf, "protein_sequence", "") or "")


class VisualDashboardCanvas(QWidget):
    """
    Large illustrated dashboard canvas.
    """

    def __init__(self, state=None, orfs=None, parent=None):
        super().__init__(parent)
        self.state = state
        self.orfs = list(orfs or [])
        self.setMinimumSize(1280, 1050)

    def set_data(self, state=None, orfs=None) -> None:
        self.state = state
        self.orfs = list(orfs or [])
        self.update()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing if QT6 else QPainter.Antialiasing)
        painter.fillRect(self.rect(), QColor("#f7fafc"))

        self._draw_title(painter)
        self._draw_workflow_graph(painter, 40, 80, self.width() - 80, 150)
        self._draw_data_process(painter, 40, 260, self.width() - 80, 220)
        self._draw_orf_summary(painter, 40, 520, self.width() - 80, 260)
        self._draw_downstream_panel(painter, 40, 820, self.width() - 80, 190)

        painter.end()

    # --------------------------------------------------------
    # Helpers
    # --------------------------------------------------------

    def _font(self, size=10, bold=False):
        font = QFont("Arial")
        font.setPointSize(size)
        font.setBold(bold)
        return font

    def _draw_rounded_rect(self, painter, rect, color, border="#cfd8dc", radius=14):
        painter.setPen(QPen(QColor(border), 1))
        painter.setBrush(QBrush(QColor(color)))
        painter.drawRoundedRect(QRectF(*rect), radius, radius)

    def _draw_text(self, painter, x, y, w, h, text, size=10, color="#263238", bold=False, align=None):
        painter.setFont(self._font(size=size, bold=bold))
        painter.setPen(QColor(color))
        painter.drawText(
            QRectF(x, y, w, h),
            align if align is not None else (_align_center()),
            text,
        )

    def _draw_section_title(self, painter, x, y, title, subtitle=""):
        painter.setFont(self._font(16, True))
        painter.setPen(QColor("#17384d"))
        painter.drawText(QRectF(x, y, 720, 30), title)

        if subtitle:
            painter.setFont(self._font(9, False))
            painter.setPen(QColor("#60717f"))
            painter.drawText(QRectF(x, y + 28, 900, 40), subtitle)

    # --------------------------------------------------------
    # Sections
    # --------------------------------------------------------

    def _draw_title(self, painter):
        painter.setFont(self._font(20, True))
        painter.setPen(QColor("#17384d"))
        painter.drawText(QRectF(40, 22, 600, 40), "ppigFinder visual workflow dashboard")

        painter.setFont(self._font(10, False))
        painter.setPen(QColor("#60717f"))
        painter.drawText(
            QRectF(40, 54, 900, 24),
            "Illustrated view of data, operations, generated outputs and next decisions.",
        )

    def _draw_workflow_graph(self, painter, x, y, w, h):
        self._draw_section_title(
            painter,
            x,
            y,
            "Workflow progression",
            "Each node becomes active as its required data or output is generated.",
        )

        completed = _completed_steps(self.state)
        current = _current_route(self.state)

        node_w = 145
        node_h = 58
        gap = max(20, (w - len(WORKFLOW_NODES) * node_w) / (len(WORKFLOW_NODES) - 1))
        start_y = y + 75

        for i, (node_id, label) in enumerate(WORKFLOW_NODES):
            nx = x + i * (node_w + gap)
            active = node_id == current
            done = node_id in completed

            if active:
                fill = "#17384d"
                text = "#ffffff"
                border = "#17384d"
            elif done:
                fill = "#d7ecf5"
                text = "#17384d"
                border = "#8bb8c9"
            else:
                fill = "#ffffff"
                text = "#60717f"
                border = "#cfd8dc"

            self._draw_rounded_rect(painter, (nx, start_y, node_w, node_h), fill, border, 16)
            self._draw_text(painter, nx + 8, start_y + 6, node_w - 16, node_h - 12, label, 10, text, True)

            if i < len(WORKFLOW_NODES) - 1:
                ax1 = nx + node_w + 4
                ax2 = nx + node_w + gap - 8
                ay = start_y + node_h / 2
                painter.setPen(QPen(QColor("#90a4ae"), 2))
                painter.drawLine(int(ax1), int(ay), int(ax2), int(ay))

                arrow = QPolygonF(
                    [
                        QPointF(ax2, ay),
                        QPointF(ax2 - 8, ay - 5),
                        QPointF(ax2 - 8, ay + 5),
                    ]
                )
                painter.setBrush(QBrush(QColor("#90a4ae")))
                painter.drawPolygon(arrow)

    def _draw_data_process(self, painter, x, y, w, h):
        self._draw_section_title(
            painter,
            x,
            y,
            "Data → operation → output",
            "The interface should show what data entered, what operation is possible, and what artifact is produced.",
        )

        box_y = y + 70
        box_h = 115
        box_w = (w - 60) / 3

        genome = _state_get(self.state, "genome_file", "No genome loaded")
        protein = _state_get(self.state, "protein_query_file", "No protein query")
        hmm = _state_get(self.state, "hmm_profile_file", "No HMM profiles")

        input_text = (
            "Genome: " + str(genome).split("/")[-1] + "\n"
            "Protein query: " + str(protein).split("/")[-1] + "\n"
            "HMM profiles: " + str(hmm).split("/")[-1]
        )

        orf_count = _state_get(self.state, "guided_orf_count", 0)
        operation_text = (
            "Available operations\n"
            "• ORF prediction\n"
            "• BLAST query vs ORFs\n"
            "• HMM/domain annotation\n"
            "• Neighbourhood selection"
        )

        output_text = (
            f"Generated outputs\n"
            f"• ORFs: {orf_count}\n"
            f"• AF3 pairs: {_state_get(self.state, 'af3_pair_count', 0)}\n"
            f"• HPC: {_state_get(self.state, 'hpc_status', 'not configured')}"
        )

        labels = [
            ("Input data", input_text, "#ffffff"),
            ("Operations", operation_text, "#eef7fb"),
            ("Outputs", output_text, "#ffffff"),
        ]

        for i, (title, body, color) in enumerate(labels):
            bx = x + i * (box_w + 30)
            self._draw_rounded_rect(painter, (bx, box_y, box_w, box_h), color, "#cfd8dc", 14)

            painter.setFont(self._font(11, True))
            painter.setPen(QColor("#17384d"))
            painter.drawText(QRectF(bx + 16, box_y + 12, box_w - 32, 24), title)

            painter.setFont(self._font(9, False))
            painter.setPen(QColor("#37474f"))
            painter.drawText(QRectF(bx + 16, box_y + 38, box_w - 32, box_h - 45), body)

            if i < 2:
                ax1 = bx + box_w + 6
                ax2 = bx + box_w + 24
                ay = box_y + box_h / 2
                painter.setPen(QPen(QColor("#90a4ae"), 2))
                painter.drawLine(int(ax1), int(ay), int(ax2), int(ay))
                painter.setBrush(QBrush(QColor("#90a4ae")))
                painter.drawPolygon(
                    QPolygonF(
                        [
                            QPointF(ax2, ay),
                            QPointF(ax2 - 7, ay - 5),
                            QPointF(ax2 - 7, ay + 5),
                        ]
                    )
                )

    def _draw_orf_summary(self, painter, x, y, w, h):
        self._draw_section_title(
            painter,
            x,
            y,
            "ORF discovery illustration",
            "After ORF prediction, this panel becomes a graphical summary of the generated protein set.",
        )

        panel_y = y + 70
        map_x = x
        map_w = w * 0.58
        hist_x = x + map_w + 30
        hist_w = w - map_w - 30
        panel_h = 160

        self._draw_rounded_rect(painter, (map_x, panel_y, map_w, panel_h), "#ffffff", "#cfd8dc", 14)
        self._draw_rounded_rect(painter, (hist_x, panel_y, hist_w, panel_h), "#ffffff", "#cfd8dc", 14)

        painter.setFont(self._font(11, True))
        painter.setPen(QColor("#17384d"))
        painter.drawText(QRectF(map_x + 16, panel_y + 12, map_w - 32, 24), "Genome / ORF map preview")
        painter.drawText(QRectF(hist_x + 16, panel_y + 12, hist_w - 32, 24), "ORF length distribution")

        if not self.orfs:
            painter.setFont(self._font(10, False))
            painter.setPen(QColor("#60717f"))
            painter.drawText(
                QRectF(map_x + 20, panel_y + 55, map_w - 40, 70),
                _align_center(),
                "Run Predict ORFs to display gene arrows along the genome.",
            )
            painter.drawText(
                QRectF(hist_x + 20, panel_y + 55, hist_w - 40, 70),
                _align_center(),
                "No ORF length distribution yet.",
            )
            return

        self._draw_orf_map(painter, map_x + 18, panel_y + 50, map_w - 36, 90)
        self._draw_orf_histogram(painter, hist_x + 18, panel_y + 45, hist_w - 36, 95)

    def _draw_orf_map(self, painter, x, y, w, h):
        genome_length = max(_safe_int(getattr(orf, "end", 0)) for orf in self.orfs)
        axis_y = y + h / 2

        painter.setPen(QPen(QColor("#607d8b"), 2))
        painter.drawLine(int(x), int(axis_y), int(x + w), int(axis_y))

        visible = self.orfs[:100]

        for orf in visible:
            start = _safe_int(getattr(orf, "start", 0))
            end = _safe_int(getattr(orf, "end", 0))
            strand = getattr(orf, "strand", "+")

            if end < start:
                start, end = end, start

            x1 = x + (start / genome_length) * w
            x2 = x + (end / genome_length) * w

            if x2 - x1 < 4:
                x2 = x1 + 4

            yy = axis_y - 28 if strand == "+" else axis_y + 12
            hh = 12
            arrow = min(10, max(4, (x2 - x1) * 0.35))

            if strand == "+":
                points = [
                    QPointF(x1, yy),
                    QPointF(x2 - arrow, yy),
                    QPointF(x2, yy + hh / 2),
                    QPointF(x2 - arrow, yy + hh),
                    QPointF(x1, yy + hh),
                ]
                color = "#43a047"
            else:
                points = [
                    QPointF(x2, yy),
                    QPointF(x1 + arrow, yy),
                    QPointF(x1, yy + hh / 2),
                    QPointF(x1 + arrow, yy + hh),
                    QPointF(x2, yy + hh),
                ]
                color = "#1e88e5"

            painter.setBrush(QBrush(QColor(color)))
            painter.setPen(QPen(QColor("#263238"), 0.4))
            painter.drawPolygon(QPolygonF(points))

        painter.setFont(self._font(8, False))
        painter.setPen(QColor("#60717f"))
        painter.drawText(QRectF(x, y + h - 16, w, 16), f"Showing {len(visible)} of {len(self.orfs)} ORFs")

    def _draw_orf_histogram(self, painter, x, y, w, h):
        lengths = [_orf_length(orf) for orf in self.orfs]
        if not lengths:
            return

        bins = [0, 50, 100, 200, 400, 800, 1600, 3200]
        labels = ["<50", "50", "100", "200", "400", "800", "1600+"]

        counts = [0] * (len(bins) - 1)

        for length in lengths:
            placed = False
            for i in range(len(bins) - 1):
                if bins[i] <= length < bins[i + 1]:
                    counts[i] += 1
                    placed = True
                    break
            if not placed:
                counts[-1] += 1

        max_count = max(counts) if counts else 1
        bar_gap = 6
        bar_w = (w - bar_gap * (len(counts) - 1)) / len(counts)
        base_y = y + h - 22

        painter.setPen(QPen(QColor("#cfd8dc"), 1))
        painter.drawLine(int(x), int(base_y), int(x + w), int(base_y))

        for i, count in enumerate(counts):
            bh = 0 if max_count == 0 else (count / max_count) * (h - 45)
            bx = x + i * (bar_w + bar_gap)
            by = base_y - bh

            painter.setBrush(QBrush(QColor("#1e88e5")))
            painter.setPen(QPen(QColor("#1565c0"), 0.5))
            painter.drawRoundedRect(QRectF(bx, by, bar_w, bh), 3, 3)

            painter.setFont(self._font(7, False))
            painter.setPen(QColor("#37474f"))
            painter.drawText(QRectF(bx, base_y + 3, bar_w, 16), _align_center(), labels[i])

    def _draw_downstream_panel(self, painter, x, y, w, h):
        self._draw_section_title(
            painter,
            x,
            y,
            "Downstream decision points",
            "The next visual layers should resemble genome-browser and dashboard views, not only forms.",
        )

        box_y = y + 70
        box_w = (w - 60) / 4
        box_h = 85

        items = [
            ("Annotation", "BLAST hits\nHMM domains\nCandidate table", "#eef7fb"),
            ("Neighbourhood", "Gene arrows\nDomain colors\nLocal context", "#ffffff"),
            ("AlphaFold / PPI", "Pair network\nAF3 JSON\nipTM/PAE metrics", "#eef7fb"),
            ("Reports", "Tables\nFigures\nWorkflow summary", "#ffffff"),
        ]

        for i, (title, body, fill) in enumerate(items):
            bx = x + i * (box_w + 20)
            self._draw_rounded_rect(painter, (bx, box_y, box_w, box_h), fill, "#cfd8dc", 14)

            painter.setFont(self._font(11, True))
            painter.setPen(QColor("#17384d"))
            painter.drawText(QRectF(bx + 14, box_y + 10, box_w - 28, 22), title)

            painter.setFont(self._font(8, False))
            painter.setPen(QColor("#37474f"))
            painter.drawText(QRectF(bx + 14, box_y + 34, box_w - 28, 44), body)


class VisualDashboardDialog(QDialog):
    """
    Dialog wrapping the visual dashboard canvas.
    """

    def __init__(self, state=None, orfs=None, parent=None):
        super().__init__(parent)

        self.setWindowTitle("ppigFinder Visual Dashboard")
        self.setWindowFlags(_window_flags())
        self.resize(1360, 920)
        self.setMinimumSize(1000, 720)
        self.setSizeGripEnabled(True)
        apply_ppigfinder_branding(self)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 16, 16, 16)
        layout.setSpacing(10)

        title = QLabel("Visual workflow dashboard")
        title.setObjectName("SectionTitle")
        layout.addWidget(title)

        subtitle = QLabel(
            "This dashboard illustrates how inputs, operations, generated outputs and next steps connect across the ppigFinder workflow."
        )
        subtitle.setWordWrap(True)
        layout.addWidget(subtitle)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)

        self.canvas = VisualDashboardCanvas(state=state, orfs=orfs)
        scroll.setWidget(self.canvas)

        layout.addWidget(scroll, 1)


def open_visual_dashboard(state=None, orfs=None, parent=None):
    dialog = VisualDashboardDialog(state=state, orfs=orfs, parent=parent)
    dialog.showMaximized()
    dialog.exec() if hasattr(dialog, "exec") else dialog.exec_()
