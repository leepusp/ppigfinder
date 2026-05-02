#!/usr/bin/env python3
"""
Home/start window for the future ppigFinder interface.
"""

from __future__ import annotations

from ppigfinder.ui_shell.branding import (
    apply_ppigfinder_branding,
    create_ppigfinder_pixmap,
)
from ppigfinder.ui_shell.bridge import PreviewActionBridge
from ppigfinder.ui_shell.help_content import help_for
from ppigfinder.ui_shell.models import HomeAction
from ppigfinder.ui_shell.qt import (
    QMainWindow,
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QGridLayout,
    QLabel,
    QPushButton,
    QFrame,
    align_center,
    pointing_hand_cursor,
)
from ppigfinder.ui_shell.theme import APP_TITLE, APP_SUBTITLE, shell_stylesheet


DEFAULT_HOME_ACTIONS = [
    HomeAction(
        id="open_genome",
        title="Open genome",
        description="Start from a nucleotide genome file.",
        input_data="DNA / genome sequence",
        output_data="Genome workspace",
        action_name="load_fasta",
    ),
    HomeAction(
        id="open_project",
        title="Open project",
        description="Resume a previous ppigFinder session.",
        input_data="Project file",
        output_data="Restored session",
        action_name="load_project",
    ),
    HomeAction(
        id="predict_orfs",
        title="Predict ORFs",
        description="Identify protein-coding regions.",
        input_data="Loaded DNA / genome",
        output_data="ORF and protein table",
        action_name="analyze_orfs",
    ),
    HomeAction(
        id="annotation",
        title="Annotation",
        description="BLAST, HMM/domain and neighborhood analysis.",
        input_data="Predicted proteins / ORFs",
        output_data="Functional annotation",
        action_name=None,
    ),
    HomeAction(
        id="alphafold",
        title="AlphaFold / PPI",
        description="Prepare AF3 jobs and interpret interaction metrics.",
        input_data="Protein pairs / AF3 output",
        output_data="Interaction metrics",
        action_name=None,
    ),
    HomeAction(
        id="reports",
        title="Reports",
        description="Export reports, snapshots and tables.",
        input_data="Current project state",
        output_data="HTML / JSON / TSV",
        action_name="export_html_report",
    ),
]


class HomeWindow(QMainWindow):
    """
    Future ppigFinder home screen.
    """

    def __init__(self, bridge=None, actions: list[HomeAction] | None = None):
        super().__init__()

        self.bridge = bridge or PreviewActionBridge()
        self.actions = actions or DEFAULT_HOME_ACTIONS

        self.setWindowTitle(f"{APP_TITLE} — Start")
        self.resize(1180, 760)
        self.setStyleSheet(shell_stylesheet())
        apply_ppigfinder_branding(self)

        central = QWidget()
        self.setCentralWidget(central)

        root = QVBoxLayout(central)
        root.setContentsMargins(32, 26, 32, 26)
        root.setSpacing(18)

        root.addLayout(self._build_header())

        section = QLabel("Start or continue an analysis")
        section.setObjectName("SectionTitle")
        root.addWidget(section)

        body = QHBoxLayout()
        body.setSpacing(18)
        root.addLayout(body)

        grid = QGridLayout()
        grid.setSpacing(14)
        body.addLayout(grid, 3)

        for index, action in enumerate(self.actions):
            card = self._create_card(action)
            row = index // 2
            col = index % 2
            grid.addWidget(card, row, col)

        self.detail_panel = self._create_detail_panel()
        body.addWidget(self.detail_panel, 2)

        self._show_details(self.actions[0])

    def _build_header(self) -> QHBoxLayout:
        header = QHBoxLayout()
        header.setSpacing(16)

        logo = QLabel()
        logo.setPixmap(create_ppigfinder_pixmap(72))
        logo.setFixedSize(80, 80)
        logo.setAlignment(align_center())
        header.addWidget(logo)

        title_box = QVBoxLayout()

        title = QLabel(APP_TITLE)
        title.setObjectName("HeroTitle")
        title_box.addWidget(title)

        subtitle = QLabel(APP_SUBTITLE)
        subtitle.setObjectName("HeroSubtitle")
        title_box.addWidget(subtitle)

        caption = QLabel(
            "Guided workflow for genome loading, ORF prediction, annotation, "
            "AlphaFold/PPI analysis and final reporting."
        )
        caption.setWordWrap(True)
        title_box.addWidget(caption)

        header.addLayout(title_box)
        header.addStretch(1)

        return header

    def _create_card(self, action: HomeAction) -> QFrame:
        card = QFrame()
        card.setObjectName("Card")
        card.setCursor(pointing_hand_cursor())
        card.setMinimumHeight(160)
        card.setMaximumHeight(190)

        layout = QVBoxLayout(card)
        layout.setContentsMargins(16, 14, 16, 14)
        layout.setSpacing(8)

        title = QLabel(action.title)
        title.setObjectName("CardTitle")
        layout.addWidget(title)

        description = QLabel(action.description)
        description.setObjectName("CardDescription")
        description.setWordWrap(True)
        layout.addWidget(description)

        io_text = QLabel(
            f"<b>Input:</b> {action.input_data}<br>"
            f"<b>Output:</b> {action.output_data}"
        )
        io_text.setWordWrap(True)
        layout.addWidget(io_text)

        layout.addStretch(1)

        button = QPushButton("Details / Run")
        button.clicked.connect(lambda checked=False, a=action: self._activate(a))
        layout.addWidget(button)

        return card

    def _create_detail_panel(self) -> QFrame:
        panel = QFrame()
        panel.setObjectName("Card")
        panel.setMinimumWidth(360)

        layout = QVBoxLayout(panel)
        layout.setContentsMargins(18, 18, 18, 18)
        layout.setSpacing(12)

        self.detail_title = QLabel("Analysis guidance")
        self.detail_title.setObjectName("SectionTitle")
        layout.addWidget(self.detail_title)

        self.detail_body = QLabel("")
        self.detail_body.setWordWrap(True)
        layout.addWidget(self.detail_body)

        self.detail_next = QLabel("")
        self.detail_next.setWordWrap(True)
        layout.addWidget(self.detail_next)

        layout.addStretch(1)

        self.detail_button = QPushButton("Run selected step")
        self.detail_button.clicked.connect(self._run_selected_action)
        layout.addWidget(self.detail_button)

        self.selected_action: HomeAction | None = None

        return panel

    def _show_details(self, action: HomeAction) -> None:
        self.selected_action = action
        help_data = help_for(action.id)

        self.detail_title.setText(action.title)
        self.detail_body.setText(
            f"<b>Purpose</b><br>{help_data['purpose']}<br><br>"
            f"<b>When to use</b><br>{help_data['when']}<br><br>"
            f"<b>Input data</b><br>{help_data['data']}<br><br>"
            f"<b>Expected output</b><br>{action.output_data}"
        )
        self.detail_next.setText(f"<b>Next step</b><br>{help_data['next']}")

        if action.action_name and self.bridge.available(action.action_name):
            self.detail_button.setText("Run selected step")
            self.detail_button.setEnabled(True)
        else:
            self.detail_button.setText("Preview only")
            self.detail_button.setEnabled(False)

    def _activate(self, action: HomeAction) -> None:
        self._show_details(action)

        if action.action_name and self.bridge.available(action.action_name):
            self.bridge.call(action.action_name)

    def _run_selected_action(self) -> None:
        if not self.selected_action:
            return

        if self.selected_action.action_name:
            self.bridge.call(self.selected_action.action_name)
