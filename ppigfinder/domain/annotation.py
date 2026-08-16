#!/usr/bin/env python3
"""
Domain models and adapters for annotation objects.
"""

from __future__ import annotations

from .orf import DomainHit, ORF


def domain_hit_from_legacy_dict(data: dict) -> DomainHit:
    """
    Convert a legacy domain dictionary into a DomainHit object.
    """
    return DomainHit(
        name=data.get("domain") or data.get("name") or "unknown",
        start=data.get("start"),
        end=data.get("end"),
        score=data.get("score"),
        evalue=data.get("evalue"),
        source=data.get("source", "unknown"),
    )


def orf_from_legacy_dict(data: dict, index: int | None = None) -> ORF:
    """
    Convert a legacy ORF dictionary into an ORF dataclass.
    """
    orf_id = data.get("id") or data.get("name")

    if not orf_id:
        orf_id = f"ORF{index}" if index is not None else "ORF"

    domains = [
        domain_hit_from_legacy_dict(item)
        for item in data.get("domains", [])
        if isinstance(item, dict)
    ]

    return ORF(
        id=orf_id,
        start=int(data.get("start", 0)),
        end=int(data.get("end", 0)),
        strand=data.get("strand", "+"),
        frame=int(data.get("frame", 0)),
        dna=data.get("dna", ""),
        protein=data.get("protein", ""),
        length=int(data.get("length", 0)),
        gc=float(data.get("gc", 0.0)),
        source=data.get("source", "unknown"),
        candidate_score=float(data.get("candidate_score", 0.0)),
        domains=domains,
        neighborhood=list(data.get("neighborhood", [])),
        gene_name=data.get("gene_name"),
        putative_function=data.get("putative_function"),
    )


def orf_to_legacy_dict(orf: ORF) -> dict:
    """
    Convert an ORF dataclass back to the legacy dictionary format.
    """
    return {
        "id": orf.id,
        "start": orf.start,
        "end": orf.end,
        "strand": orf.strand,
        "frame": orf.frame,
        "dna": orf.dna,
        "protein": orf.protein,
        "length": orf.length,
        "gc": orf.gc,
        "source": orf.source,
        "candidate_score": orf.candidate_score,
        "domains": [
            {
                "domain": domain.name,
                "start": domain.start,
                "end": domain.end,
                "score": domain.score,
                "evalue": domain.evalue,
                "source": domain.source,
            }
            for domain in orf.domains
        ],
        "neighborhood": list(orf.neighborhood),
        "gene_name": orf.gene_name,
        "putative_function": orf.putative_function,
    }
