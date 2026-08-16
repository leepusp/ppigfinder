#!/usr/bin/env python3
"""
HTML report generation for ppigFinder.

Reports are generated from serializable data structures, especially the
versioned ProjectState snapshot.
"""

from __future__ import annotations

from html import escape
from pathlib import Path
from datetime import datetime

from ppigfinder.domain.project import ProjectState
from ppigfinder.io.project_json import read_project_json


def _format_value(value) -> str:
    if value is None:
        return ""

    if isinstance(value, float):
        return f"{value:.4g}"

    return str(value)


def _table(headers: list[str], rows: list[list[object]]) -> str:
    if not rows:
        return "<p class=\"empty\">No records available.</p>"

    th = "".join(f"<th>{escape(str(header))}</th>" for header in headers)

    body = []
    for row in rows:
        tds = "".join(f"<td>{escape(_format_value(value))}</td>" for value in row)
        body.append(f"<tr>{tds}</tr>")

    return f"<table><thead><tr>{th}</tr></thead><tbody>{''.join(body)}</tbody></table>"


def _section(title: str, body: str) -> str:
    return f"""
<section class="section">
<h2>{escape(title)}</h2>
{body}
</section>
"""


def _css() -> str:
    return """
:root {
    --bg: #f5f7f8;
    --panel: #ffffff;
    --text: #263238;
    --muted: #607d8b;
    --line: #d7dee2;
    --head: #e8f0f4;
    --accent: #1b3a4b;
}
body {
    font-family: Arial, Helvetica, sans-serif;
    margin: 32px;
    color: var(--text);
    background: var(--bg);
}
h1, h2, h3 {
    color: var(--accent);
}
.header {
    margin-bottom: 24px;
}
.meta {
    color: var(--muted);
    font-size: 13px;
}
.grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
    gap: 12px;
}
.card {
    background: var(--panel);
    border: 1px solid var(--line);
    border-radius: 10px;
    padding: 14px;
}
.card .label {
    color: var(--muted);
    font-size: 12px;
}
.card .value {
    font-size: 22px;
    font-weight: 700;
    margin-top: 4px;
}
.section {
    background: var(--panel);
    border: 1px solid var(--line);
    border-radius: 10px;
    padding: 18px;
    margin-bottom: 20px;
}
table {
    border-collapse: collapse;
    width: 100%;
    font-size: 13px;
}
th {
    background: var(--head);
    text-align: left;
}
th, td {
    border: 1px solid var(--line);
    padding: 6px 8px;
}
.empty {
    color: var(--muted);
    font-style: italic;
}
"""


def render_project_report(project: ProjectState, title: str = "ppigFinder Report") -> str:
    """
    Render a standalone HTML report from ProjectState.
    """
    genome = project.genome or {}
    orfs = project.orfs or []
    blast = project.blast or {}
    hmm = project.hmm or {}
    alphafold = project.alphafold or {}

    genome_name = genome.get("name", "")
    genome_sequence = genome.get("sequence", "") or ""
    genome_length = len(genome_sequence)

    blast_results = blast.get("results", []) or []
    hmm_hits = hmm.get("hits", []) or []
    af3_jobs = alphafold.get("jobs", []) or []
    af3_results = alphafold.get("results", []) or []

    summary_cards = f"""
<div class="grid">
  <div class="card"><div class="label">Genome</div><div class="value">{escape(str(genome_name or "-"))}</div></div>
  <div class="card"><div class="label">Genome length</div><div class="value">{genome_length}</div></div>
  <div class="card"><div class="label">ORFs</div><div class="value">{len(orfs)}</div></div>
  <div class="card"><div class="label">BLAST hits</div><div class="value">{len(blast_results)}</div></div>
  <div class="card"><div class="label">HMM/domain hits</div><div class="value">{len(hmm_hits)}</div></div>
  <div class="card"><div class="label">AF3 jobs/results</div><div class="value">{len(af3_jobs)} / {len(af3_results)}</div></div>
</div>
"""

    orf_rows = []
    for i, orf in enumerate(orfs, start=1):
        orf_rows.append([
            orf.get("id") or f"ORF{i}",
            orf.get("start", ""),
            orf.get("end", ""),
            orf.get("strand", ""),
            orf.get("frame", ""),
            orf.get("length", ""),
            orf.get("gc", ""),
            orf.get("source", ""),
            orf.get("candidate_score", ""),
        ])

    blast_rows = []
    for hit in blast_results:
        if isinstance(hit, dict):
            blast_rows.append([
                hit.get("orf_index", hit.get("index", "")),
                hit.get("identity", hit.get("pident", "")),
                hit.get("score", hit.get("bitscore", "")),
                hit.get("evalue", ""),
                hit.get("method", ""),
            ])

    hmm_rows = []
    for hit in hmm_hits:
        if isinstance(hit, dict):
            hmm_rows.append([
                hit.get("orf_index", ""),
                hit.get("domain", hit.get("name", "")),
                hit.get("start", ""),
                hit.get("end", ""),
                hit.get("score", ""),
                hit.get("evalue", ""),
                hit.get("source", ""),
            ])

    af3_rows = []
    for item in af3_results:
        if isinstance(item, dict):
            af3_rows.append([
                item.get("name", ""),
                item.get("iptm", item.get("ipTM", "")),
                item.get("ptm", item.get("pTM", "")),
                item.get("ranking_score", ""),
                item.get("pae_inter", ""),
                item.get("pae_min", ""),
                item.get("classification", ""),
            ])

    sections = [
        _section("Summary", summary_cards),
        _section(
            "ORFs",
            _table(
                ["ID", "Start", "End", "Strand", "Frame", "Length", "GC %", "Source", "Score"],
                orf_rows,
            ),
        ),
        _section(
            "BLAST results",
            _table(["ORF", "Identity", "Score", "E-value", "Method"], blast_rows),
        ),
        _section(
            "HMM/domain results",
            _table(["ORF", "Domain", "Start", "End", "Score", "E-value", "Source"], hmm_rows),
        ),
        _section(
            "AlphaFold results",
            _table(["Name", "ipTM", "pTM", "Ranking", "PAE inter", "PAE min", "Class"], af3_rows),
        ),
    ]

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{escape(title)}</title>
<style>{_css()}</style>
</head>
<body>
<div class="header">
<h1>{escape(title)}</h1>
<p class="meta">Generated by ppigFinder on {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
<p class="meta">Project schema version: {project.metadata.schema_version}</p>
</div>
{''.join(sections)}
</body>
</html>
"""


def render_basic_report(
    title: str,
    genome_name: str | None = None,
    genome_length: int | None = None,
    orfs: list[dict] | None = None,
    interaction_results: list[dict] | None = None,
) -> str:
    """
    Backward-compatible helper that builds a small ProjectState internally.
    """
    project = ProjectState(
        genome={
            "name": genome_name or "",
            "sequence": "N" * int(genome_length or 0),
        },
        orfs=orfs or [],
        alphafold={"results": interaction_results or []},
    )
    return render_project_report(project, title=title)


def write_project_report(path: str | Path, project: ProjectState, title: str = "ppigFinder Report") -> None:
    path = Path(path)
    html = render_project_report(project, title=title)

    with path.open("w", encoding="utf-8") as handle:
        handle.write(html)


def write_report_from_project_json(
    project_json_path: str | Path,
    html_path: str | Path,
    title: str = "ppigFinder Report",
) -> None:
    project = read_project_json(project_json_path)
    write_project_report(html_path, project, title=title)


def write_basic_report(path: str | Path, **kwargs) -> None:
    path = Path(path)
    html = render_basic_report(**kwargs)

    with path.open("w", encoding="utf-8") as handle:
        handle.write(html)
