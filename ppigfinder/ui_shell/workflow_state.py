#!/usr/bin/env python3
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


WORKFLOW_ORDER = [
    "overview",
    "data",
    "genome",
    "orfs",
    "annotation",
    "alphafold",
    "hpc",
    "reports",
]


@dataclass
class WorkflowEvent:
    step: str
    action: str
    message: str = ""


@dataclass
class WorkflowState:
    current_route: str = "overview"
    loaded_inputs: dict[str, Any] = field(default_factory=dict)
    metrics: dict[str, Any] = field(default_factory=dict)
    flags: dict[str, Any] = field(default_factory=dict)
    events: list[WorkflowEvent] = field(default_factory=list)

    def set_current_route(self, route_id: str) -> None:
        self.current_route = route_id

    def set_input(self, key: str, value: Any) -> None:
        self.loaded_inputs[key] = value

    def set_metric(self, key: str, value: Any) -> None:
        self.metrics[key] = value

    def set_flag(self, key: str, value: Any = True) -> None:
        self.flags[key] = value

    def get(self, key: str, default=None):
        if key in self.loaded_inputs:
            return self.loaded_inputs.get(key, default)
        if key in self.metrics:
            return self.metrics.get(key, default)
        if key in self.flags:
            return self.flags.get(key, default)
        return default

    def add_event(self, step: str, action: str, message: str = "") -> None:
        self.events.append(WorkflowEvent(step=step, action=action, message=message))

    def completed_steps(self) -> set[str]:
        completed = set()

        if self.loaded_inputs.get("genome_file"):
            completed.update({"data", "genome"})

        if self.metrics.get("guided_orf_count"):
            completed.add("orfs")

        if (
            self.flags.get("guided_blast_planned")
            or self.flags.get("guided_hmm_planned")
            or self.flags.get("guided_neighborhood_planned")
        ):
            completed.add("annotation")

        if self.metrics.get("af3_pair_count") or self.loaded_inputs.get("af3_json_path"):
            completed.add("alphafold")

        if self.flags.get("hpc_connected") or self.loaded_inputs.get("hpc_profile"):
            completed.add("hpc")

        if self.flags.get("guided_summary_exported") or self.flags.get("html_report_exported"):
            completed.add("reports")

        return completed

    def next_recommended_step(self) -> str:
        completed = self.completed_steps()

        for step in WORKFLOW_ORDER:
            if step == "overview":
                continue
            if step not in completed:
                return step

        return "reports"
