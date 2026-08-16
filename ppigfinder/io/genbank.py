#!/usr/bin/env python3
"""
GenBank reader/writer for ppigFinder.

This module is intentionally independent from the PyQt GUI.
"""

from pathlib import Path
from datetime import datetime
import re


def parse_genbank(filepath: str) -> dict:
    """
    Pure-stdlib parser for GenBank flat files (INSDC format).

    Returns the same dict as parse_snapgene_dna:
      sequence, topology, name, features, primers (empty), notes
    """
    result = {
        'sequence': '', 'topology': 'linear',
        'name': Path(filepath).stem,
        'features': [], 'primers': [], 'notes': {},
    }

    # Color palette by feature type (matching SnapGene defaults)
    _TYPE_COLORS = {
        'CDS':          '#99ccff',
        'gene':         '#ffcc99',
        'mRNA':         '#ff9999',
        'rRNA':         '#ccff99',
        'tRNA':         '#ffff99',
        'ncRNA':        '#ff99ff',
        'regulatory':   '#99ffff',
        'rep_origin':   '#ff9966',
        'misc_feature': '#cccccc',
        'promoter':     '#ffcc00',
        'terminator':   '#cc99ff',
        'primer_bind':  '#ff66cc',
    }

    with open(filepath, 'r', encoding='utf-8', errors='replace') as fh:
        content = fh.read()

    # ── Split into records (multi-record support) ───────────────
    # Use only the first record
    record_text = content.split('//')[0]

    # ── LOCUS ───────────────────────────────────────────────────
    locus_m = re.search(r'^LOCUS\s+(\S+)\s+.*?(circular|linear)', record_text,
                        re.MULTILINE | re.IGNORECASE)
    if locus_m:
        result['name']     = locus_m.group(1)
        result['topology'] = locus_m.group(2).lower()

    # ── DEFINITION / ACCESSION for notes ───────────────────────
    for tag in ('DEFINITION', 'ACCESSION', 'VERSION', 'ORGANISM'):
        m = re.search(rf'^{tag}\s+(.+?)(?=\n[A-Z])', record_text,
                      re.MULTILINE | re.DOTALL)
        if m:
            result['notes'][tag] = ' '.join(m.group(1).split())

    # ── FEATURES ────────────────────────────────────────────────
    feat_block_m = re.search(r'^FEATURES\s+.*?\n(.*?)^(?:ORIGIN|CONTIG)',
                              record_text, re.MULTILINE | re.DOTALL)
    if feat_block_m:
        feat_text = feat_block_m.group(1)

        # Split into individual feature entries
        # Each starts at col 5 with a keyword, qualifiers at col 21
        feat_entries = re.split(r'\n(?=     \S)', feat_text)

        for entry in feat_entries:
            lines = entry.split('\n')
            if not lines:
                continue
            first = lines[0].strip()
            if not first:
                continue

            # Feature type is the first token
            parts = first.split()
            if not parts:
                continue
            ftype = parts[0]
            loc_str = parts[1] if len(parts) > 1 else ''

            # Collect remaining qualifier lines
            for ln in lines[1:]:
                s = ln.strip()
                if s and not s.startswith('/'):
                    loc_str += s   # continuation of location
                else:
                    break

            # Parse location string → start, end, strand
            strand = '-' if loc_str.startswith('complement') else '+'
            # Strip complement(...) and join(...)
            loc_clean = re.sub(r'(complement|join|order)\(', '', loc_str).rstrip(')')
            # Grab all numeric pairs
            ranges = re.findall(r'(\d+)\.\.(\d+)', loc_clean)
            if not ranges:
                # single position
                m_single = re.search(r'(\d+)', loc_clean)
                if m_single:
                    pos = int(m_single.group(1))
                    ranges = [(str(pos), str(pos))]

            if not ranges:
                continue

            # For joined features keep first + last segment
            start = int(ranges[0][0]) - 1   # 0-based
            end   = int(ranges[-1][1])

            # Collect qualifiers
            qual_text = '\n'.join(lines[1:])
            name_m = re.search(r'/(?:gene|locus_tag|product|label)="([^"]+)"', qual_text)
            feat_name = name_m.group(1) if name_m else ftype

            color = _TYPE_COLORS.get(ftype, '#cccccc')

            if ftype == 'primer_bind':
                result['primers'].append({
                    'name': feat_name,
                    'start': start, 'end': end,
                    'strand': strand,
                })
            else:
                result['features'].append({
                    'name': feat_name,
                    'type': ftype,
                    'color': color,
                    'start': start,
                    'end':   end,
                    'strand': strand,
                })

    # ── ORIGIN (sequence) ───────────────────────────────────────
    origin_m = re.search(r'^ORIGIN\s*\n(.*)', record_text,
                         re.MULTILINE | re.DOTALL)
    if origin_m:
        raw_seq = origin_m.group(1)
        result['sequence'] = re.sub(r'[^a-zA-Z]', '', raw_seq).upper()

    return result


def write_genbank(filepath: str, sequence: str, orfs: list,
                  sg_features: list = None, name: str = 'sequence',
                  topology: str = 'linear'):
    """
    Write a GenBank flat-file (.gb) containing:
      - CDS features for each analysed ORF
      - misc_feature for each imported SnapGene feature (if any)
    """
    sl = len(sequence)
    now = datetime.now().strftime('%d-%b-%Y').upper()
    topo_str = 'circular' if topology == 'circular' else 'linear  '
    mol = 'DNA'

    lines = []

    # ── LOCUS ──────────────────────────────────────────────────
    lines.append(f"LOCUS       {name[:16]:<16} {sl:>9} bp    {mol}     {topo_str} {now}")
    lines.append(f"DEFINITION  {name} — exported by ORF Secretion Pipeline v20.")
    lines.append( "ACCESSION   .")
    lines.append( "VERSION     .")
    lines.append( "KEYWORDS    .")
    lines.append( "SOURCE      .")
    lines.append( "  ORGANISM  .")
    lines.append( "            .")
    lines.append( "FEATURES             Location/Qualifiers")

    def _loc(start, end, strand):
        """Format a GenBank location string (1-based, closed)."""
        s = start + 1
        e = end
        loc = f"{s}..{e}"
        return f"complement({loc})" if strand == '-' else loc

    def _wrap_qual(key, val, indent=21):
        """Wrap a qualifier value at col 80."""
        prefix = ' ' * indent
        full = f'{prefix}/{key}="{val}"'
        out = []
        while len(full) > 79:
            out.append(full[:79])
            full = prefix + full[79:]
        out.append(full)
        return '\n'.join(out)

    feat_indent = '     '

    # ── SnapGene features (preserved from original file) ────────
    for feat in (sg_features or []):
        ftype = feat.get('type', 'misc_feature')
        lines.append(f"{feat_indent}{ftype:<16}{_loc(feat['start'], feat['end'], feat.get('strand', '+'))}")
        lines.append(_wrap_qual('label', feat.get('name', ftype)))
        lines.append(_wrap_qual('note', f"Imported from SnapGene; color: {feat.get('color','#aaaaaa')}"))

    # ── ORFs as CDS ─────────────────────────────────────────────
    for i, orf in enumerate(orfs):
        loc = _loc(orf['start'], orf['end'], orf['strand'])
        lines.append(f"{feat_indent}{'CDS':<16}{loc}")
        gene_name = orf.get('gene_name') or f"orf{i+1}"
        lines.append(_wrap_qual('gene', gene_name))
        prot = orf['protein'].rstrip('*')
        lines.append(_wrap_qual('product',
            orf.get('putative_function') or orf.get('observation') or 'hypothetical protein'))
        lines.append(_wrap_qual('protein_id', f"ORF{i+1}"))
        domains = '; '.join(d['domain'] for d in orf.get('domains', []))
        if domains:
            lines.append(_wrap_qual('note', f"HMM domains: {domains}"))
        # Translate qualifier — wrap at 60 aa per line inside the qualifier
        prot_chunks = [prot[j:j+60] for j in range(0, len(prot), 60)]
        lines.append(f"{'':21}/translation=\"{prot_chunks[0]}")
        for chunk in prot_chunks[1:]:
            lines.append(f"{'':21}{chunk}")
        lines[-1] += '"'

    # ── ORIGIN ──────────────────────────────────────────────────
    lines.append("ORIGIN")
    seq_lc = sequence.lower()
    for pos in range(0, sl, 60):
        chunk = seq_lc[pos:pos+60]
        # Split into groups of 10
        groups = ' '.join(chunk[i:i+10] for i in range(0, len(chunk), 10))
        lines.append(f"{pos+1:>9} {groups}")
    lines.append("//")

    with open(filepath, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(lines) + '\n')
