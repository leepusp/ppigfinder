"""AlphaFold results TSV export helpers.

GUI-independent helpers extracted from ppigFinder v29.14.
"""

from __future__ import annotations

AF3_RESULTS_TSV_HEADERS = [
    'job_name', 'job_dir', 'n_chains', 'chains',
    'iptm', 'ptm', 'mean_pLDDT', 'ranking_score',
    'pae_inter_mean', 'pae_min_inter', 'contact_frac_pct',
    'hotspot_score', 'hotspot_mean_pae', 'hotspot_size',
    'pae_min_best_pair',
    'best_pair', 'contact_region',
    'anchor_seq_A', 'anchor_seq_B',
    'chain_iptm', 'chain_ptm',
    'fraction_disordered', 'has_clash',
    'n_diffusion_samples',
    'high_confidence',
    # ── v2.5 — validation/duplicate metadata ──────────────────
    'seq_status',
    'seq_fingerprint',
    'duplicate_group_id',
    'duplicate_role',
    'duplicate_count',
    'completeness',
    'chain_layout_source',
    'truncated_chains',
    'ambiguous_chains',
    'validation_warnings',
]

def _fmt(v, fmt='{:.4f}'):
    if v is None:
        return ''
    try:
        return fmt.format(v)
    except (TypeError, ValueError):
        return str(v)

def write_af3_results_tsv_file(path, results, anchor_sequence_provider=None):
    """Write AF3 analysis results using the v29.14 TSV schema."""
    with open(path, 'w', encoding='utf-8') as f:
        f.write('\t'.join(AF3_RESULTS_TSV_HEADERS) + '\n')
        for res in results:
            chains = ';'.join(
                f"{res['chain_to_orf'].get(c, c)}({res['chain_lens'][c]})"
                for c in res.get('chain_order', []))
            bp = res.get('best_pair', ('?', '?'))
            bp_str = (f"{res['chain_to_orf'].get(bp[0], bp[0])}"
                      f"-{res['chain_to_orf'].get(bp[1], bp[1])}"
                      if isinstance(bp, tuple) and len(bp) == 2
                      else '')
            iptm = res.get('iptm')
            paei = res.get('pae_inter')
            hc = ('yes' if (iptm is not None and paei is not None
                            and iptm > 0.75 and paei < 8.0)
                  else 'no')
            # Rich per-pair metrics for best pair
            bp_pm = res.get('pair_metrics', {}).get(bp, {})
            pae_min_s     = _fmt(bp_pm.get('pair_pae_min'), '{:.3f}')
            contact_frac_s= _fmt((res.get('contact_frac') or 0) * 100, '{:.1f}')
            _hs           = res.get('hotspot') or {}
            hotspot_sc_s  = _fmt(_hs.get('hotspot_score'),    '{:.3f}')
            hotspot_pae_s = _fmt(_hs.get('hotspot_mean_pae'), '{:.2f}')
            _hsa, _hsb    = _hs.get('hotspot_size_A'), _hs.get('hotspot_size_B')
            hotspot_sz_s  = f"{_hsa}x{_hsb}" if (_hsa and _hsb) else ''
            anchor_seq_a, anchor_seq_b = (anchor_sequence_provider(res) if anchor_sequence_provider is not None else ('', ''))
            anchor_seq_a = anchor_seq_a.replace('\t', ' ').replace('\n', ' ')
            anchor_seq_b = anchor_seq_b.replace('\t', ' ').replace('\n', ' ')
            # Per-chain metric stringification
            ci_s = (';'.join(_fmt(x, '{:.3f}')
                             for x in res['chain_iptm'])
                    if res.get('chain_iptm') else '')
            cp_s = (';'.join(_fmt(x, '{:.3f}')
                             for x in res['chain_ptm'])
                    if res.get('chain_ptm') else '')
            has_clash = res.get('has_clash')
            clash_s = ('yes' if (has_clash is True
                                 or (isinstance(has_clash,
                                                 (int, float))
                                     and has_clash >= 0.5))
                       else ('no' if has_clash is not None else ''))
            row = [
                res.get('job_name', ''),
                res.get('job_dir', ''),
                str(res.get('n_chains', '')),
                chains,
                _fmt(iptm),
                _fmt(res.get('ptm')),
                _fmt(res.get('mean_plddt'), '{:.2f}'),
                _fmt(res.get('ranking_score')),
                _fmt(paei, '{:.3f}'),
                pae_min_s,
                contact_frac_s,
                hotspot_sc_s,
                hotspot_pae_s,
                hotspot_sz_s,
                pae_min_s,
                bp_str,
                (res.get('contact_region') or '')
                    .replace('\t', ' ').replace('\n', ' '),
                anchor_seq_a,
                anchor_seq_b,
                ci_s,
                cp_s,
                _fmt(res.get('fraction_disordered'), '{:.4f}'),
                clash_s,
                str(len(res.get('ranking_samples', []))),
                hc,
                # ── v2.5 fields ────────────────────────────────
                res.get('seq_status', '') or '',
                res.get('seq_fingerprint', '') or '',
                str(res.get('duplicate_group_id', -1)),
                res.get('duplicate_role', '') or '',
                str(res.get('duplicate_count', 1)),
                res.get('completeness', '') or '',
                res.get('chain_layout_source', '') or '',
                # Truncation info → semicolon-joined "chain:modeled/expected"
                ';'.join(
                    f"{t['chain']}:{t['modeled']}/{t['expected']}"
                    for t in (res.get('truncation_info') or [])),
                # Ambiguity info → semicolon-joined "chain:cands"
                ';'.join(
                    f"{chr(65 + a['chain_index'])}:"
                    f"{'|'.join(a['candidates'])}"
                    for a in (res.get('ambiguity_info') or [])),
                # Validation warnings → joined with " || " (no tabs)
                ' || '.join(
                    (res.get('validation_warnings') or []))
                    .replace('\t', ' ').replace('\n', ' '),
            ]
            f.write('\t'.join(row) + '\n')
