#!/usr/bin/env python3
"""
Workflow service for ppigFinder.

This service evaluates the current project/application state and determines
which workflow steps are available or completed.
"""

from __future__ import annotations

from ppigfinder.infrastructure.backends import BACKENDS
from ppigfinder.workflows.default_workflow import build_default_workflow
from ppigfinder.workflows.models import (
    WorkflowContext,
    WorkflowDefinition,
    WorkflowStep,
    WorkflowStepStatus,
)


class WorkflowService:
    """
    High-level service for workflow definitions and runtime status.
    """

    def __init__(self, workflow: WorkflowDefinition | None = None):
        self.workflow = workflow or build_default_workflow()

    def context_from_legacy_window(self, window) -> WorkflowContext:
        dna_sequence = getattr(window, "dna_sequence", "") or ""
        orfs = getattr(window, "orfs", []) or []
        blast_results = getattr(window, "blast_results", []) or []
        hmm_hits = getattr(window, "hmm_hits_all", []) or []
        af3_jobs = getattr(window, "af3_jobs", []) or []
        af3_results = getattr(window, "af3_results", []) or []

        return WorkflowContext(
            has_project=bool(dna_sequence or orfs or af3_results),
            has_genome=bool(dna_sequence),
            has_orfs=bool(orfs),
            has_blast_results=bool(blast_results),
            has_hmm_results=bool(hmm_hits),
            has_af3_jobs=bool(af3_jobs),
            has_af3_results=bool(af3_results),
            has_report_data=bool(dna_sequence or orfs or blast_results or hmm_hits or af3_results),
            blast_available=bool(BACKENDS.get("blast+", {}).get("available")),
            hmmer_available=bool(BACKENDS.get("hmmer3", {}).get("available")),
            pyrodigal_available=bool(BACKENDS.get("pyrodigal", {}).get("available")),
        )

    def completed_step_ids(self, context: WorkflowContext) -> set[str]:
        completed: set[str] = set()

        if context.has_genome:
            completed.update({"open_genome", "genome_overview"})

        if context.has_orfs:
            completed.update({"predict_orfs", "review_orfs"})

        if context.has_blast_results:
            completed.add("blast_query")

        if context.has_hmm_results:
            completed.add("hmm_domains")

        if context.has_af3_jobs:
            completed.add("export_af3_server_json")

        if context.has_af3_results:
            completed.update({"import_af3_results", "export_af3_table"})

        if context.has_report_data:
            completed.add("html_report")

        return completed

    def status_for_step(
        self,
        step: WorkflowStep,
        context: WorkflowContext,
    ) -> WorkflowStepStatus:
        completed = self.completed_step_ids(context)
        missing = [required for required in step.requires if required not in completed]

        if missing:
            return WorkflowStepStatus(
                step_id=step.id,
                enabled=False,
                completed=step.id in completed,
                reason="Requires: " + ", ".join(missing),
            )

        return WorkflowStepStatus(
            step_id=step.id,
            enabled=True,
            completed=step.id in completed,
            reason="Ready",
        )

    def status_by_step_id(self, context: WorkflowContext) -> dict[str, WorkflowStepStatus]:
        return {
            step.id: self.status_for_step(step, context)
            for step in self.workflow.all_steps()
        }
