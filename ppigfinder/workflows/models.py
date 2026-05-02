#!/usr/bin/env python3
"""
Workflow models for ppigFinder.

These models describe the scientific analysis flow independently from the GUI.
They can later be used by a wizard, start screen, CLI, batch runner or the
current legacy interface.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal


DataType = Literal[
    "project",
    "dna",
    "genome",
    "protein",
    "orf",
    "annotation",
    "alphafold",
    "report",
    "system",
]


@dataclass(slots=True)
class WorkflowStep:
    """
    One actionable step in the ppigFinder analysis workflow.
    """

    id: str
    title: str
    stage: str
    data_type: DataType
    description: str
    input_data: str
    output_data: str
    action_name: str | None = None
    tab_hint: str | None = None
    requires: list[str] = field(default_factory=list)
    optional: bool = False


@dataclass(slots=True)
class WorkflowStage:
    """
    A group of related workflow steps.
    """

    id: str
    title: str
    description: str
    steps: list[WorkflowStep] = field(default_factory=list)


@dataclass(slots=True)
class WorkflowDefinition:
    """
    Full ppigFinder workflow definition.
    """

    id: str
    title: str
    description: str
    stages: list[WorkflowStage] = field(default_factory=list)

    def all_steps(self) -> list[WorkflowStep]:
        steps: list[WorkflowStep] = []
        for stage in self.stages:
            steps.extend(stage.steps)
        return steps

    def step_by_id(self, step_id: str) -> WorkflowStep | None:
        for step in self.all_steps():
            if step.id == step_id:
                return step
        return None


@dataclass(slots=True)
class WorkflowContext:
    """
    Current analysis state used to decide which steps are ready.
    """

    has_project: bool = False
    has_genome: bool = False
    has_orfs: bool = False
    has_blast_results: bool = False
    has_hmm_results: bool = False
    has_af3_jobs: bool = False
    has_af3_results: bool = False
    has_report_data: bool = False

    blast_available: bool = False
    hmmer_available: bool = False
    pyrodigal_available: bool = False


@dataclass(slots=True)
class WorkflowStepStatus:
    """
    Runtime status of a workflow step.
    """

    step_id: str
    enabled: bool
    completed: bool
    reason: str = ""
