#!/usr/bin/env python3
"""
Interactive ORF Discovery browser for the experimental guided UI shell.

This view is inspired by locus/gene-arrow visualizers such as LoVis4u:
- searchable and filterable ORF table
- selectable ORFs
- zoomable/pannable ORF map
- focus selected ORF
- sequence/details panel
"""

from __future__ import annotations

from collections import Counter
from statistics import mean
import re

try:
    from PyQt6.QtCore import Qt, QRectF, QPointF
    from PyQt6.QtGui import QPainter, QColor, QPen, QBrush, QPolygonF, QFont
    from PyQt6.QtWidgets import (
        QDialog,
        QVBoxLayout,
        QHBoxLayout,
        QGridLayout,
        QLabel,
        QPushButton,
        QFrame,
        QSplitter,
        QWidget,
        QTableWidget,
        QTableWidgetItem,
        QHeaderView,
        QLineEdit,
        QComboBox,
        QSpinBox,
        QCheckBox,
        QPlainTextEdit,
        QTabWidget,
    )
    QT6 = True
except Exception:
    from PyQt5.QtCore import Qt, QRectF, QPointF
    from PyQt5.QtGui import QPainter, QColor, QPen, QBrush, QPolygonF, QFont
    from PyQt5.QtWidgets import (
        QDialog,
        QVBoxLayout,
        QHBoxLayout,
        QGridLayout,
        QLabel,
        QPushButton,
        QFrame,
        QSplitter,
        QWidget,
        QTableWidget,
        QTableWidgetItem,
        QHeaderView,
        QLineEdit,
        QComboBox,
        QSpinBox,
        QCheckBox,
        QPlainTextEdit,
        QTabWidget,
    )
    QT6 = False

from ppigfinder.ui_shell.branding import apply_ppigfinder_branding
from ppigfinder.ui_shell.genomic_tracks import read_genome_sequence, compute_gc_and_skew_bins


def _window_flags():
    flags = Qt.WindowType.Window if QT6 else Qt.Window
    flags |= Qt.WindowType.WindowMinimizeButtonHint if QT6 else Qt.WindowMinimizeButtonHint
    flags |= Qt.WindowType.WindowMaximizeButtonHint if QT6 else Qt.WindowMaximizeButtonHint
    flags |= Qt.WindowType.WindowCloseButtonHint if QT6 else Qt.WindowCloseButtonHint
    return flags


def _align_center():
    return Qt.AlignmentFlag.AlignCenter if QT6 else Qt.AlignCenter


def _align_left_vcenter():
    if QT6:
        return Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter
    return Qt.AlignLeft | Qt.AlignVCenter


def _safe_int(value, default=0):
    try:
        return int(value)
    except Exception:
        return default


def _orf_id(orf, fallback=""):
    return str(getattr(orf, "id", fallback) or fallback)


def _orf_start(orf):
    return _safe_int(getattr(orf, "start", 0))


def _orf_end(orf):
    return _safe_int(getattr(orf, "end", 0))


def _orf_strand(orf):
    return str(getattr(orf, "strand", "?") or "?")


def _orf_frame(orf):
    return str(getattr(orf, "frame", "?") or "?")


def _orf_length(orf) -> int:
    if hasattr(orf, "aa_length"):
        return _safe_int(getattr(orf, "aa_length"))
    return len(getattr(orf, "protein_sequence", "") or "")


def _nt_length(orf) -> int:
    seq = getattr(orf, "nt_sequence", "") or getattr(orf, "dna_sequence", "") or ""
    if seq:
        return len(seq)
    start = _orf_start(orf)
    end = _orf_end(orf)
    return abs(end - start) + 1 if start and end else 0


def _protein_sequence(orf):
    return str(getattr(orf, "protein_sequence", "") or getattr(orf, "sequence", "") or "")


def _nt_sequence(orf):
    return str(getattr(orf, "nt_sequence", "") or getattr(orf, "dna_sequence", "") or "")


def _protein_preview(orf, length: int = 58) -> str:
    seq = _protein_sequence(orf)
    if len(seq) <= length:
        return seq
    return seq[:length] + "..."



def _load_genome_sequence_from_parent(parent) -> str:
    """
    Recover loaded genome sequence from the workflow parent when available.
    """
    try:
        if parent is not None and hasattr(parent, "loaded_genome_sequence"):
            seq = getattr(parent, "loaded_genome_sequence") or ""
            if seq:
                return seq

        if parent is not None and hasattr(parent, "workflow_state"):
            genome_file = parent.workflow_state.get("genome_file")
            if genome_file:
                return read_genome_sequence(genome_file)
    except Exception:
        return ""

    return ""


def _format_sequence(seq: str, width: int = 70, block: int = 10) -> str:
    seq = re.sub(r"\s+", "", seq or "").upper()
    if not seq:
        return "No sequence available."

    lines = []
    for i in range(0, len(seq), width):
        chunk = seq[i:i + width]
        grouped = " ".join(chunk[j:j + block] for j in range(0, len(chunk), block))
        lines.append(f"{i + 1:>8}  {grouped:<90}  {i + len(chunk):>8}")
    return "\n".join(lines)


class SummaryCard(QFrame):
    def __init__(self, title: str, value: str, subtitle: str = "", parent=None):
        super().__init__(parent)
        self.setObjectName("InfoCard")

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 10, 12, 10)
        layout.setSpacing(5)

        title_label = QLabel(title)
        title_label.setObjectName("InfoCardTitle")

        value_label = QLabel(value)
        value_label.setObjectName("InfoCardValue")
        value_label.setWordWrap(True)

        subtitle_label = QLabel(subtitle)
        subtitle_label.setObjectName("InfoCardSubtitle")
        subtitle_label.setWordWrap(True)

        layout.addWidget(title_label)
        layout.addWidget(value_label)
        layout.addWidget(subtitle_label)


class ORFMapWidget(QWidget):
    """
    Zoomable/pannable ORF map preview.

    It is intentionally lightweight and local, but follows the same visual idea
    as locus browsers: directional arrows along a genomic coordinate axis.
    """

    def __init__(self, orfs, genome_sequence: str = "", parent=None):
        super().__init__(parent)

        self.orfs = sorted(list(orfs or []), key=lambda o: _orf_start(o))
        self.genome_sequence = genome_sequence or ""
        self.genome_length = max(self._infer_genome_length(), len(self.genome_sequence))
        self.view_start = 1
        self.view_end = max(1, self.genome_length)
        self.selected_orf_id = ""
        self.setMinimumHeight(360)

    def _infer_genome_length(self) -> int:
        if not self.orfs:
            return 0
        return max(_orf_end(orf) for orf in self.orfs)

    def set_selected_orf(self, orf_id: str, focus: bool = False) -> None:
        self.selected_orf_id = str(orf_id or "")
        if focus:
            orf = self._find_orf(orf_id)
            if orf is not None:
                self.focus_orf(orf)
        self.update()

    def _find_orf(self, orf_id: str):
        for orf in self.orfs:
            if _orf_id(orf) == orf_id:
                return orf
        return None

    def reset_view(self) -> None:
        self.view_start = 1
        self.view_end = max(1, self.genome_length)
        self.update()

    def focus_orf(self, orf, flank: int = 8000) -> None:
        start = _orf_start(orf)
        end = _orf_end(orf)

        if not start or not end:
            return

        if end < start:
            start, end = end, start

        self.view_start = max(1, start - flank)
        self.view_end = min(max(1, self.genome_length), end + flank)

        if self.view_end <= self.view_start:
            self.view_end = min(max(1, self.genome_length), self.view_start + 1000)

        self.update()

    def zoom(self, factor: float) -> None:
        if self.genome_length <= 0:
            return

        center = (self.view_start + self.view_end) / 2
        span = max(100, (self.view_end - self.view_start) * factor)
        self.view_start = max(1, int(center - span / 2))
        self.view_end = min(self.genome_length, int(center + span / 2))

        if self.view_end <= self.view_start:
            self.view_end = min(self.genome_length, self.view_start + 100)

        self.update()

    def pan(self, fraction: float) -> None:
        if self.genome_length <= 0:
            return

        span = self.view_end - self.view_start
        shift = int(span * fraction)

        new_start = self.view_start + shift
        new_end = self.view_end + shift

        if new_start < 1:
            new_start = 1
            new_end = new_start + span

        if new_end > self.genome_length:
            new_end = self.genome_length
            new_start = max(1, new_end - span)

        self.view_start = int(new_start)
        self.view_end = int(new_end)
        self.update()

    def _visible_orfs(self):
        visible = []
        for orf in self.orfs:
            start = _orf_start(orf)
            end = _orf_end(orf)
            if end < start:
                start, end = end, start
            if end >= self.view_start and start <= self.view_end:
                visible.append(orf)

        if len(visible) <= 700:
            return visible

        selected = [orf for orf in visible if _orf_id(orf) == self.selected_orf_id]
        step = max(1, len(visible) // 650)
        sampled = visible[::step]
        for orf in selected:
            if orf not in sampled:
                sampled.append(orf)

        return sorted(sampled, key=lambda o: _orf_start(o))

    def _x_for_pos(self, pos: int, left: float, width: float) -> float:
        span = max(1, self.view_end - self.view_start)
        return left + ((pos - self.view_start) / span) * width

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing if QT6 else QPainter.Antialiasing)

        painter.fillRect(self.rect(), QColor("#f8fbfc"))

        width = self.width()
        height = self.height()
        left = 44
        right = width - 26
        usable = max(1, right - left)
        axis_y = 112 if self.genome_sequence else height // 2

        painter.setPen(QPen(QColor("#17384d"), 2))
        painter.drawLine(left, axis_y, right, axis_y)

        if not self.orfs or self.genome_length <= 0:
            painter.setPen(QColor("#60717f"))
            painter.drawText(self.rect(), _align_center(), "No ORF map available.")
            painter.end()
            return

        # Ticks
        tick_count = 6
        painter.setFont(QFont("Arial", 8))
        painter.setPen(QColor("#60717f"))
        for i in range(tick_count + 1):
            frac = i / tick_count
            x = left + usable * frac
            pos = int(self.view_start + (self.view_end - self.view_start) * frac)
            painter.drawLine(int(x), int(axis_y - 5), int(x), int(axis_y + 5))
            painter.drawText(QRectF(x - 50, axis_y + 9, 100, 18), _align_center(), f"{pos:,}")

        visible = self._visible_orfs()

        plus_color = QColor("#43a047")
        minus_color = QColor("#1e88e5")
        selected_color = QColor("#ffb300")
        border = QColor("#263238")

        label_count = 0

        for orf in visible:
            start = _orf_start(orf)
            end = _orf_end(orf)
            strand = _orf_strand(orf)
            oid = _orf_id(orf)

            if end < start:
                start, end = end, start

            x1 = self._x_for_pos(start, left, usable)
            x2 = self._x_for_pos(end, left, usable)

            if x2 < left or x1 > right:
                continue

            x1 = max(left, x1)
            x2 = min(right, x2)

            if x2 - x1 < 5:
                x2 = x1 + 5

            is_selected = oid == self.selected_orf_id
            y = axis_y - 48 if strand == "+" else axis_y + 28
            h = 18
            arrow = min(13, max(5, (x2 - x1) * 0.35))

            if strand == "+":
                points = [
                    QPointF(x1, y),
                    QPointF(x2 - arrow, y),
                    QPointF(x2, y + h / 2),
                    QPointF(x2 - arrow, y + h),
                    QPointF(x1, y + h),
                ]
                color = selected_color if is_selected else plus_color
            else:
                points = [
                    QPointF(x2, y),
                    QPointF(x1 + arrow, y),
                    QPointF(x1, y + h / 2),
                    QPointF(x1 + arrow, y + h),
                    QPointF(x2, y + h),
                ]
                color = selected_color if is_selected else minus_color

            painter.setBrush(QBrush(color))
            painter.setPen(QPen(border, 1.2 if is_selected else 0.5))
            painter.drawPolygon(QPolygonF(points))

            if is_selected or (x2 - x1 > 80 and label_count < 16):
                painter.setFont(QFont("Arial", 8, QFont.Weight.Bold if QT6 else QFont.Bold))
                painter.setPen(QColor("#263238"))
                label_y = y - 18 if strand == "+" else y + h + 2
                painter.drawText(QRectF(x1 - 20, label_y, max(60, x2 - x1 + 40), 16), _align_center(), oid)
                label_count += 1

        # Header/legend
        painter.setFont(QFont("Arial", 9))
        painter.setPen(QColor("#17384d"))
        painter.drawText(44, 22, f"View: {self.view_start:,} - {self.view_end:,} nt")
        painter.drawText(44, 42, f"Visible ORFs: {len(visible):,} / total ORFs: {len(self.orfs):,}")

        painter.setPen(QColor("#43a047"))
        painter.drawText(width - 210, 22, "+ strand")
        painter.setPen(QColor("#1e88e5"))
        painter.drawText(width - 130, 22, "- strand")
        painter.setPen(QColor("#ff8f00"))
        painter.drawText(width - 70, 22, "selected")

        if self.genome_sequence:
            self._draw_computed_signal_tracks(
                painter,
                left,
                right,
                axis_y + 70,
                max(80, height - axis_y - 95),
            )

        painter.end()

    def _draw_track_axis(self, painter, x: float, y: float, w: float, label: str, value_text: str = ""):
        painter.setPen(QPen(QColor("#cfd8dc"), 1))
        painter.drawLine(int(x), int(y), int(x + w), int(y))

        painter.setFont(QFont("Arial", 8))
        painter.setPen(QColor("#37474f"))
        painter.drawText(QRectF(x, y - 18, 260, 16), label)

        if value_text:
            painter.drawText(QRectF(x + w - 280, y - 18, 280, 16), _align_center(), value_text)

    def _draw_computed_signal_tracks(self, painter, left: float, right: float, top: float, height: float):
        """
        Draw computed GC content and GC skew tracks below the ORF map.

        These are available from genome sequence alone. DNA-Seq/RNA-Seq coverage
        requires external signal files and will be integrated separately.
        """
        width = max(1, right - left)
        bins = compute_gc_and_skew_bins(
            self.genome_sequence,
            self.view_start,
            self.view_end,
            bins=max(60, min(260, int(width / 5))),
        )

        if not bins:
            painter.setFont(QFont("Arial", 9))
            painter.setPen(QColor("#60717f"))
            painter.drawText(
                QRectF(left, top, width, 40),
                _align_center(),
                "Genome sequence not available for GC tracks.",
            )
            return

        gc_values = [row["gc"] for row in bins]
        skew_values = [row["skew"] for row in bins]

        avg_gc = sum(gc_values) / len(gc_values)
        avg_skew = sum(skew_values) / len(skew_values)

        gc_y = top + 30
        skew_y = top + 95
        track_h = 42

        self._draw_track_axis(
            painter,
            left,
            gc_y + track_h / 2,
            width,
            "GC content (%)",
            f"Average: {avg_gc * 100:.2f}%",
        )

        self._draw_track_axis(
            painter,
            left,
            skew_y + track_h / 2,
            width,
            "GC skew (G-C)/(G+C)",
            f"Average: {avg_skew:.3f}",
        )

        self._draw_gc_content_area(painter, bins, left, gc_y, width, track_h)
        self._draw_gc_skew_area(painter, bins, left, skew_y, width, track_h)

        painter.setFont(QFont("Arial", 8))
        painter.setPen(QColor("#60717f"))
        painter.drawText(
            QRectF(left, skew_y + track_h + 10, width, 18),
            "Computed tracks from genome sequence. Experimental DNA/RNA/proteomics tracks require external signal files.",
        )

    def _x_for_track_bin(self, pos: int, left: float, width: float) -> float:
        span = max(1, self.view_end - self.view_start)
        return left + ((pos - self.view_start) / span) * width

    def _draw_gc_content_area(self, painter, bins, left: float, y: float, width: float, height: float):
        baseline = y + height
        max_h = height

        painter.setPen(QPen(QColor("#fb8c00"), 1))
        painter.setBrush(QBrush(QColor(255, 183, 77, 120)))

        points = [QPointF(left, baseline)]

        for row in bins:
            x = self._x_for_track_bin(row["start"], left, width)
            value = row["gc"]
            py = baseline - value * max_h
            points.append(QPointF(x, py))

        points.append(QPointF(left + width, baseline))
        painter.drawPolygon(QPolygonF(points))

    def _draw_gc_skew_area(self, painter, bins, left: float, y: float, width: float, height: float):
        mid = y + height / 2
        scale = height / 2

        painter.setPen(QPen(QColor("#9e9e9e"), 1))
        painter.drawLine(int(left), int(mid), int(left + width), int(mid))

        for row in bins:
            x1 = self._x_for_track_bin(row["start"], left, width)
            x2 = self._x_for_track_bin(row["end"], left, width)
            skew = max(-1.0, min(1.0, row["skew"]))
            y2 = mid - skew * scale

            if skew >= 0:
                color = QColor(76, 175, 80, 130)
                border = QColor("#43a047")
            else:
                color = QColor(156, 39, 176, 120)
                border = QColor("#8e24aa")

            painter.setBrush(QBrush(color))
            painter.setPen(QPen(border, 0.4))
            top = min(mid, y2)
            bottom = max(mid, y2)
            painter.drawRect(QRectF(x1, top, max(1, x2 - x1), max(1, bottom - top)))




class LengthHistogramWidget(QWidget):
    """
    Simple ORF length distribution plot.
    """

    def __init__(self, orfs, parent=None):
        super().__init__(parent)
        self.orfs = list(orfs or [])
        self.setMinimumHeight(160)

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing if QT6 else QPainter.Antialiasing)
        painter.fillRect(self.rect(), QColor("#ffffff"))

        lengths = [_orf_length(orf) for orf in self.orfs]
        if not lengths:
            painter.setPen(QColor("#60717f"))
            painter.drawText(self.rect(), _align_center(), "No length distribution available.")
            painter.end()
            return

        bins = [0, 50, 100, 200, 400, 800, 1600, 3200, 6400]
        labels = ["<50", "50", "100", "200", "400", "800", "1600", "3200+"]

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

        left = 38
        top = 20
        bottom = self.height() - 34
        right = self.width() - 20
        usable_w = max(1, right - left)
        usable_h = max(1, bottom - top)

        painter.setPen(QPen(QColor("#cfd8dc"), 1))
        painter.drawLine(left, bottom, right, bottom)
        painter.drawLine(left, top, left, bottom)

        gap = 7
        bar_w = (usable_w - gap * (len(counts) - 1)) / len(counts)

        for i, count in enumerate(counts):
            h = 0 if max_count == 0 else (count / max_count) * usable_h
            x = left + i * (bar_w + gap)
            y = bottom - h

            painter.setBrush(QBrush(QColor("#1e88e5")))
            painter.setPen(QPen(QColor("#1565c0"), 0.5))
            painter.drawRoundedRect(QRectF(x, y, bar_w, h), 4, 4)

            painter.setFont(QFont("Arial", 8))
            painter.setPen(QColor("#37474f"))
            painter.drawText(QRectF(x, bottom + 4, bar_w, 18), _align_center(), labels[i])

        painter.setFont(QFont("Arial", 8))
        painter.setPen(QColor("#60717f"))
        painter.drawText(QRectF(4, top - 2, 30, 14), _align_center(), str(max_count))
        painter.drawText(QRectF(4, bottom - 8, 30, 14), _align_center(), "0")

        painter.end()



class FigureViewerDialog(QDialog):
    """
    Generic full-screen figure viewer for guided visual panels.

    For ORF maps, this viewer adds zoom, pan and region controls.
    """

    def __init__(self, title: str, content_widget: QWidget, parent=None):
        super().__init__(parent)

        self.content_widget = content_widget

        self.setWindowTitle(title)
        self.setWindowFlags(_window_flags())
        self.resize(1280, 850)
        self.setMinimumSize(980, 680)
        self.setSizeGripEnabled(True)
        apply_ppigfinder_branding(self)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 16, 16, 16)
        layout.setSpacing(10)

        title_label = QLabel(title)
        title_label.setObjectName("SectionTitle")
        layout.addWidget(title_label)

        if isinstance(content_widget, ORFMapWidget):
            controls = QHBoxLayout()

            btn_full = QPushButton("Full genome")
            btn_focus = QPushButton("Focus selected")
            btn_zoom_in = QPushButton("Zoom +")
            btn_zoom_out = QPushButton("Zoom -")
            btn_left = QPushButton("←")
            btn_right = QPushButton("→")

            self.region_start = QSpinBox()
            self.region_end = QSpinBox()
            self.region_start.setRange(1, max(1, content_widget.genome_length))
            self.region_end.setRange(1, max(1, content_widget.genome_length))
            self.region_start.setPrefix("start ")
            self.region_end.setPrefix("end ")
            self.region_start.setValue(max(1, content_widget.view_start))
            self.region_end.setValue(max(1, content_widget.view_end))

            btn_go = QPushButton("Go region")

            btn_full.clicked.connect(self._map_full_genome)
            btn_focus.clicked.connect(self._map_focus_selected)
            btn_zoom_in.clicked.connect(lambda: self._map_zoom(0.5))
            btn_zoom_out.clicked.connect(lambda: self._map_zoom(2.0))
            btn_left.clicked.connect(lambda: self._map_pan(-0.35))
            btn_right.clicked.connect(lambda: self._map_pan(0.35))
            btn_go.clicked.connect(self._map_go_region)

            controls.addWidget(btn_full)
            controls.addWidget(btn_focus)
            controls.addWidget(btn_zoom_in)
            controls.addWidget(btn_zoom_out)
            controls.addWidget(btn_left)
            controls.addWidget(btn_right)
            controls.addWidget(QLabel("Region"))
            controls.addWidget(self.region_start)
            controls.addWidget(self.region_end)
            controls.addWidget(btn_go)

            layout.addLayout(controls)

        layout.addWidget(content_widget, 1)

    def _sync_region_inputs(self):
        if not isinstance(self.content_widget, ORFMapWidget):
            return

        self.region_start.setValue(max(1, int(self.content_widget.view_start)))
        self.region_end.setValue(max(1, int(self.content_widget.view_end)))

    def _map_full_genome(self):
        if not isinstance(self.content_widget, ORFMapWidget):
            return

        self.content_widget.reset_view()
        self._sync_region_inputs()

    def _map_focus_selected(self):
        if not isinstance(self.content_widget, ORFMapWidget):
            return

        orf = self.content_widget._find_orf(self.content_widget.selected_orf_id)

        if orf is not None:
            self.content_widget.focus_orf(orf, flank=12000)
            self._sync_region_inputs()

    def _map_zoom(self, factor: float):
        if not isinstance(self.content_widget, ORFMapWidget):
            return

        self.content_widget.zoom(factor)
        self._sync_region_inputs()

    def _map_pan(self, fraction: float):
        if not isinstance(self.content_widget, ORFMapWidget):
            return

        self.content_widget.pan(fraction)
        self._sync_region_inputs()

    def _map_go_region(self):
        if not isinstance(self.content_widget, ORFMapWidget):
            return

        start = int(self.region_start.value())
        end = int(self.region_end.value())

        if end < start:
            start, end = end, start

        if end == start:
            end = min(self.content_widget.genome_length, start + 200)

        self.content_widget.view_start = max(1, start)
        self.content_widget.view_end = min(self.content_widget.genome_length, end)
        self.content_widget.update()
        self._sync_region_inputs()


class GuidedORFResultsDialog(QDialog):
    """
    Interactive ORF browser.
    """

    MAX_ROWS = 2500

    def __init__(self, orfs, genome_sequence: str = "", parent=None):
        super().__init__(parent)

        self.orfs = sorted(list(orfs or []), key=lambda o: _orf_start(o))
        self.genome_sequence = _load_genome_sequence_from_parent(parent)
        self.filtered_orfs = list(self.orfs)
        self.workflow_parent = parent
        self.selected_orf = None

        lengths = [_orf_length(orf) for orf in self.orfs]
        plus = sum(1 for orf in self.orfs if _orf_strand(orf) == "+")
        minus = sum(1 for orf in self.orfs if _orf_strand(orf) == "-")

        self.setWindowTitle("ORF Discovery Results")
        self.setWindowFlags(_window_flags())
        self.resize(1450, 900)
        self.setMinimumSize(1100, 720)
        self.setSizeGripEnabled(True)
        apply_ppigfinder_branding(self)

        root = QVBoxLayout(self)
        root.setContentsMargins(18, 18, 18, 18)
        root.setSpacing(10)

        title = QLabel("ORF Discovery Results")
        title.setObjectName("SectionTitle")
        root.addWidget(title)

        subtitle = QLabel(
            "Predicted ORFs are available for BLAST, HMM/domain annotation, neighbourhood analysis "
            "and AlphaFold/PPI candidate selection. Use filters, select an ORF, and zoom the map to inspect local context."
        )
        subtitle.setWordWrap(True)
        root.addWidget(subtitle)

        cards = QGridLayout()
        cards.setHorizontalSpacing(10)
        cards.setVerticalSpacing(10)
        cards.addWidget(SummaryCard("ORFs predicted", f"{len(self.orfs):,}", "Protein-coding candidates"), 0, 0)
        cards.addWidget(SummaryCard("Mean length", f"{mean(lengths):.1f} aa" if lengths else "N/A", "Average protein length"), 0, 1)
        cards.addWidget(SummaryCard("Longest ORF", f"{max(lengths)} aa" if lengths else "N/A", "Maximum protein length"), 0, 2)
        cards.addWidget(SummaryCard("Strands", f"+ {plus:,} / - {minus:,}", "Predicted orientation"), 0, 3)
        root.addLayout(cards)

        main_splitter = QSplitter()
        main_splitter.setChildrenCollapsible(False)
        main_splitter.setOrientation(Qt.Orientation.Horizontal if QT6 else Qt.Horizontal)

        # ----------------------------------------------------
        # Left panel: controls + table
        # ----------------------------------------------------
        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 8, 0)
        left_layout.setSpacing(8)

        controls = QFrame()
        controls.setObjectName("InfoCard")
        controls_layout = QGridLayout(controls)
        controls_layout.setContentsMargins(10, 8, 10, 8)
        controls_layout.setHorizontalSpacing(8)
        controls_layout.setVerticalSpacing(6)

        controls_layout.addWidget(QLabel("Search"), 0, 0)
        self.search_input = QLineEdit()
        self.search_input.setPlaceholderText("ORF ID or sequence fragment")
        self.search_input.textChanged.connect(self._apply_filters)
        controls_layout.addWidget(self.search_input, 0, 1, 1, 3)

        controls_layout.addWidget(QLabel("Strand"), 1, 0)
        self.strand_filter = QComboBox()
        self.strand_filter.addItems(["All", "+", "-"])
        self.strand_filter.currentIndexChanged.connect(self._apply_filters)
        controls_layout.addWidget(self.strand_filter, 1, 1)

        controls_layout.addWidget(QLabel("Min aa"), 1, 2)
        self.min_aa_filter = QSpinBox()
        self.min_aa_filter.setRange(0, 100000)
        self.min_aa_filter.setValue(0)
        self.min_aa_filter.valueChanged.connect(self._apply_filters)
        controls_layout.addWidget(self.min_aa_filter, 1, 3)

        self.auto_focus_checkbox = QCheckBox("Auto-focus selected ORF on map")
        self.auto_focus_checkbox.setChecked(True)
        controls_layout.addWidget(self.auto_focus_checkbox, 2, 0, 1, 4)

        left_layout.addWidget(controls)

        self.status_label = QLabel("")
        self.status_label.setObjectName("InfoFooter")
        left_layout.addWidget(self.status_label)

        self.table = QTableWidget(0, 8, self)
        self.table.setHorizontalHeaderLabels(
            [
                "ORF ID",
                "Start",
                "End",
                "Strand",
                "Frame",
                "AA length",
                "NT length",
                "Protein preview",
            ]
        )
        self.table.verticalHeader().setVisible(False)
        self.table.setAlternatingRowColors(True)
        self.table.setSortingEnabled(False)
        self.table.itemSelectionChanged.connect(self._on_table_selection_changed)

        header = self.table.horizontalHeader()
        try:
            for col in range(7):
                header.setSectionResizeMode(col, QHeaderView.ResizeMode.ResizeToContents if QT6 else QHeaderView.ResizeToContents)
            header.setSectionResizeMode(7, QHeaderView.ResizeMode.Stretch if QT6 else QHeaderView.Stretch)
        except Exception:
            pass

        font = QFont("Courier New")
        font.setPointSize(9)
        self.table.setFont(font)

        left_layout.addWidget(self.table, 1)

        # ----------------------------------------------------
        # Right panel: map, histogram, details, actions
        # ----------------------------------------------------
        right = QWidget()
        right_layout = QVBoxLayout(right)
        right_layout.setContentsMargins(8, 0, 0, 0)
        right_layout.setSpacing(8)

        map_title = QLabel("Interactive ORF map")
        map_title.setObjectName("SectionSubTitle")
        right_layout.addWidget(map_title)

        self.map_widget = ORFMapWidget(self.orfs, genome_sequence=self.genome_sequence)
        right_layout.addWidget(self.map_widget, 2)

        map_controls = QHBoxLayout()

        btn_all = QPushButton("Full genome")
        btn_all.clicked.connect(self.map_widget.reset_view)
        map_controls.addWidget(btn_all)

        btn_focus = QPushButton("Focus selected")
        btn_focus.clicked.connect(self._focus_selected_on_map)
        map_controls.addWidget(btn_focus)

        btn_zoom_in = QPushButton("Zoom +")
        btn_zoom_in.clicked.connect(lambda: self.map_widget.zoom(0.5))
        map_controls.addWidget(btn_zoom_in)

        btn_zoom_out = QPushButton("Zoom -")
        btn_zoom_out.clicked.connect(lambda: self.map_widget.zoom(2.0))
        map_controls.addWidget(btn_zoom_out)

        btn_left = QPushButton("←")
        btn_left.clicked.connect(lambda: self.map_widget.pan(-0.35))
        map_controls.addWidget(btn_left)

        btn_right = QPushButton("→")
        btn_right.clicked.connect(lambda: self.map_widget.pan(0.35))
        map_controls.addWidget(btn_right)

        btn_map_full = QPushButton("Open map full screen")
        btn_map_full.clicked.connect(self._open_map_fullscreen)
        map_controls.addWidget(btn_map_full)

        right_layout.addLayout(map_controls)

        hist_title = QLabel("ORF length distribution")
        hist_title.setObjectName("SectionSubTitle")
        right_layout.addWidget(hist_title)

        self.histogram = LengthHistogramWidget(self.orfs)
        right_layout.addWidget(self.histogram, 1)

        hist_controls = QHBoxLayout()
        hist_controls.addStretch(1)

        btn_hist_full = QPushButton("Open histogram full screen")
        btn_hist_full.clicked.connect(self._open_histogram_fullscreen)
        hist_controls.addWidget(btn_hist_full)

        right_layout.addLayout(hist_controls)

        self.details_tabs = QTabWidget()
        self.details_meta = QPlainTextEdit()
        self.details_protein = QPlainTextEdit()
        self.details_dna = QPlainTextEdit()

        for widget in (self.details_meta, self.details_protein, self.details_dna):
            widget.setReadOnly(True)
            widget.setFont(QFont("Courier New", 9))

        self.details_tabs.addTab(self.details_meta, "Selected ORF")
        self.details_tabs.addTab(self.details_protein, "Protein")
        self.details_tabs.addTab(self.details_dna, "DNA")
        right_layout.addWidget(self.details_tabs, 2)

        next_title = QLabel("Next workflow options")
        next_title.setObjectName("SectionSubTitle")
        right_layout.addWidget(next_title)

        action_row = QHBoxLayout()

        btn_annotation = QPushButton("Annotation")
        btn_annotation.clicked.connect(self._go_to_annotation)
        action_row.addWidget(btn_annotation)

        btn_candidates = QPushButton("Candidates")
        btn_candidates.clicked.connect(self._review_candidates)
        action_row.addWidget(btn_candidates)

        btn_export = QPushButton("Export FASTA")
        btn_export.clicked.connect(self._export_fasta)
        action_row.addWidget(btn_export)

        btn_af3 = QPushButton("AlphaFold / PPI")
        btn_af3.clicked.connect(self._go_to_alphafold)
        action_row.addWidget(btn_af3)

        right_layout.addLayout(action_row)

        main_splitter.addWidget(left)
        main_splitter.addWidget(right)
        main_splitter.setStretchFactor(0, 3)
        main_splitter.setStretchFactor(1, 2)

        root.addWidget(main_splitter, 1)

        self._apply_filters()
        if self.orfs:
            self._select_first_row()

    def _open_map_fullscreen(self):
        """
        Open the current ORF map view in a maximized figure window.
        """
        clone = ORFMapWidget(self.orfs, genome_sequence=self.genome_sequence)
        clone.view_start = self.map_widget.view_start
        clone.view_end = self.map_widget.view_end
        clone.selected_orf_id = self.map_widget.selected_orf_id
        clone.setMinimumHeight(650)

        dialog = FigureViewerDialog(
            "Interactive ORF map — full screen",
            clone,
            parent=self,
        )
        dialog.showMaximized()
        dialog.exec() if hasattr(dialog, "exec") else dialog.exec_()

    def _open_histogram_fullscreen(self):
        """
        Open ORF length distribution in a maximized figure window.
        """
        clone = LengthHistogramWidget(self.orfs)
        clone.setMinimumHeight(650)

        dialog = FigureViewerDialog(
            "ORF length distribution — full screen",
            clone,
            parent=self,
        )
        dialog.showMaximized()
        dialog.exec() if hasattr(dialog, "exec") else dialog.exec_()

    # --------------------------------------------------------
    # Filtering/table
    # --------------------------------------------------------

    def _apply_filters(self):
        query = self.search_input.text().strip().upper()
        strand = self.strand_filter.currentText()
        min_aa = self.min_aa_filter.value()

        filtered = []

        for orf in self.orfs:
            oid = _orf_id(orf).upper()
            protein = _protein_sequence(orf).upper()
            s = _orf_strand(orf)
            aa = _orf_length(orf)

            if strand != "All" and s != strand:
                continue

            if aa < min_aa:
                continue

            if query and query not in oid and query not in protein:
                continue

            filtered.append(orf)

        self.filtered_orfs = filtered
        self._populate_table()

    def _populate_table(self):
        self.table.setSortingEnabled(False)
        visible = self.filtered_orfs[: self.MAX_ROWS]

        self.table.setRowCount(len(visible))

        for row, orf in enumerate(visible):
            values = [
                _orf_id(orf, f"orf_{row + 1}"),
                str(_orf_start(orf)),
                str(_orf_end(orf)),
                _orf_strand(orf),
                _orf_frame(orf),
                str(_orf_length(orf)),
                str(_nt_length(orf)),
                _protein_preview(orf),
            ]

            for col, value in enumerate(values):
                item = QTableWidgetItem(value)
                if col in {1, 2, 3, 4, 5, 6}:
                    item.setTextAlignment(_align_center())
                else:
                    item.setTextAlignment(_align_left_vcenter())
                self.table.setItem(row, col, item)

        self.status_label.setText(
            f"Showing {len(visible):,} of {len(self.filtered_orfs):,} filtered ORFs "
            f"(total {len(self.orfs):,}). Use filters to narrow the table."
        )

        self.table.resizeRowsToContents()
        self.table.setSortingEnabled(True)

    def _select_first_row(self):
        if self.table.rowCount() > 0:
            self.table.selectRow(0)
            self._on_table_selection_changed()

    def _orf_from_table_row(self, row: int):
        item = self.table.item(row, 0)
        if item is None:
            return None

        oid = item.text()
        for orf in self.orfs:
            if _orf_id(orf) == oid:
                return orf
        return None

    def _on_table_selection_changed(self):
        indexes = self.table.selectedIndexes()
        if not indexes:
            return

        row = indexes[0].row()
        orf = self._orf_from_table_row(row)

        if orf is None:
            return

        self.selected_orf = orf
        oid = _orf_id(orf)

        self.map_widget.set_selected_orf(
            oid,
            focus=self.auto_focus_checkbox.isChecked(),
        )

        self._update_selected_details(orf)

    def _update_selected_details(self, orf):
        oid = _orf_id(orf)
        meta = [
            f"ORF ID:      {oid}",
            f"Start:       {_orf_start(orf)}",
            f"End:         {_orf_end(orf)}",
            f"Strand:      {_orf_strand(orf)}",
            f"Frame:       {_orf_frame(orf)}",
            f"AA length:   {_orf_length(orf)}",
            f"NT length:   {_nt_length(orf)}",
            "",
            "This selected ORF can be used as a protein of interest for",
            "annotation, neighbourhood analysis and AlphaFold/PPI candidate building.",
        ]

        self.details_meta.setPlainText("\n".join(meta))
        self.details_protein.setPlainText(_format_sequence(_protein_sequence(orf), width=70, block=10))
        self.details_dna.setPlainText(_format_sequence(_nt_sequence(orf), width=70, block=10))

    # --------------------------------------------------------
    # Map/actions
    # --------------------------------------------------------

    def _focus_selected_on_map(self):
        if self.selected_orf is not None:
            self.map_widget.focus_orf(self.selected_orf)

    def _go_to_annotation(self):
        if hasattr(self.workflow_parent, "show_route"):
            self.workflow_parent.show_route("annotation")

    def _go_to_alphafold(self):
        if hasattr(self.workflow_parent, "show_route"):
            self.workflow_parent.show_route("alphafold")

    def _review_candidates(self):
        if hasattr(self.workflow_parent, "_show_annotation_candidates"):
            self.workflow_parent._show_annotation_candidates()

    def _export_fasta(self):
        if hasattr(self.workflow_parent, "_export_guided_orfs_fasta"):
            self.workflow_parent._export_guided_orfs_fasta()


def show_guided_orf_results(orfs, parent=None) -> None:
    dialog = GuidedORFResultsDialog(orfs, parent=parent)
    dialog.showMaximized()
    dialog.exec() if hasattr(dialog, "exec") else dialog.exec_()
