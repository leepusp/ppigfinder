#!/usr/bin/env python3
"""
Explicit guided workflow model for ppigFinder.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class WorkflowStep:
    id: str
    title: str
    objective: str
    required_inputs: tuple[str, ...] = field(default_factory=tuple)
    operations: tuple[str, ...] = field(default_factory=tuple)
    outputs: tuple[str, ...] = field(default_factory=tuple)
    visualizations: tuple[str, ...] = field(default_factory=tuple)
    next_steps: tuple[str, ...] = field(default_factory=tuple)
    optional: bool = False


WORKFLOW_STEPS: tuple[WorkflowStep, ...] = (
    WorkflowStep(
        id="data",
        title="Data / Project",
        objective="Start a new analysis, validate genome input, or restore a previous project.",
        required_inputs=("FASTA / GenBank / SnapGene / Project / Snapshot",),
        operations=("Input validation", "Metadata extraction", "Project restoration"),
        outputs=("Loaded genome/project state", "Input metadata", "Validation status"),
        visualizations=("Input summary cards", "Validation status", "File metadata"),
        next_steps=("genome", "orfs"),
    ),
    WorkflowStep(
        id="genome",
        title="DNA / Genome",
        objective="Inspect genome-level information before protein-level analysis.",
        required_inputs=("Loaded genome sequence",),
        operations=("Genome inspection", "Coordinate preparation", "Genome map preparation"),
        outputs=("Genome metadata", "Coordinate system", "Genome map data"),
        visualizations=("Genome summary", "Coordinate preview", "Future lovis4u-like map"),
        next_steps=("orfs",),
    ),
    WorkflowStep(
        id="orfs",
        title="Protein / ORFs",
        objective="Predict protein-coding ORFs and generate protein sequences.",
        required_inputs=("Loaded genome sequence",),
        operations=("Pyrodigal", "Six-frame ORF scan", "Hybrid ORF prediction", "Protein FASTA export"),
        outputs=("ORF table", "Protein sequences", "ORF FASTA", "ORF map"),
        visualizations=("ORF table", "ORF map preview", "Length/frame/strand summaries"),
        next_steps=("annotation",),
    ),
    WorkflowStep(
        id="annotation",
        title="Annotation",
        objective="Annotate ORFs and prioritize biologically plausible candidates.",
        required_inputs=("Predicted ORFs", "Protein sequences"),
        operations=("BLAST", "HMM/domain search", "Candidate review", "Neighbourhood preparation"),
        outputs=("BLAST hits", "HMM/domain hits", "Candidate ORFs", "Annotation tables"),
        visualizations=("Candidate table", "Domain status", "BLAST/HMM summaries"),
        next_steps=("alphafold",),
    ),
    WorkflowStep(
        id="alphafold",
        title="AlphaFold / PPI",
        objective="Prepare and review candidate protein-protein interaction predictions.",
        required_inputs=("Candidate ORFs", "Protein sequences"),
        operations=("Pair generation", "AF3 Server JSON export", "AF3 result import"),
        outputs=("AF3 candidate pairs", "AF3 JSON", "AF3 result folder", "Interaction metrics"),
        visualizations=("Pair table", "AF3 job summary", "Future ipTM × PAE dashboard"),
        next_steps=("hpc", "reports"),
    ),
    WorkflowStep(
        id="hpc",
        title="DaVinci / HPC",
        objective="Optionally prepare or execute computational workflows on DaVinci/HPC.",
        required_inputs=("AF3 jobs", "SSH/HPC profile"),
        operations=("Local cluster detection", "SSH test", "Slurm template preparation"),
        outputs=("Connection status", "Execution mode", "Slurm template"),
        visualizations=("HPC status cards", "Connection report", "Submission template"),
        next_steps=("reports",),
        optional=True,
    ),
    WorkflowStep(
        id="reports",
        title="Reports",
        objective="Export reproducible summaries, tables, figures and project state.",
        required_inputs=("Project state",),
        operations=("HTML report export", "Snapshot export", "Table export", "Guided summary export"),
        outputs=("HTML report", "Project Snapshot JSON", "TSV/CSV", "Markdown summary"),
        visualizations=("Report status", "Export checklist", "Project provenance"),
        next_steps=(),
    ),
)


STEP_BY_ID = {step.id: step for step in WORKFLOW_STEPS}


def get_step(step_id: str) -> WorkflowStep | None:
    return STEP_BY_ID.get(step_id)


def ordered_steps() -> tuple[WorkflowStep, ...]:
    return WORKFLOW_STEPS
