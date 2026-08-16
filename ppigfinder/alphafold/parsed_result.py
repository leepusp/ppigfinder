"""AlphaFold parsed-result builders.

This module contains GUI-independent result-dictionary construction extracted
from ppigFinder v29.14.
"""

from __future__ import annotations

from typing import Any, Mapping


def build_af3_job_result(ctx: Mapping[str, Any]) -> dict:
    """Build the canonical v29.14 AF3 analysis result dictionary.

    The input mapping follows the local variable names produced by
    legacy_v29_14._af3a_parse_job. Keeping this as a separate builder makes
    the final result schema testable while leaving file parsing in the GUI
    layer for now.
    """
    job_dir = ctx.get('job_dir')
    sum_path = ctx.get('sum_path')
    conf_path = ctx.get('conf_path')
    model_cif = ctx.get('model_cif')
    ranking_csv = ctx.get('ranking_csv')
    input_json = ctx.get('input_json')
    lazy = ctx.get('lazy')
    orf_names = ctx.get('orf_names')
    chain_order = ctx.get('chain_order')
    chain_lens = ctx.get('chain_lens')
    chain_to_orf = ctx.get('chain_to_orf')
    n_chains = ctx.get('n_chains')
    iptm = ctx.get('iptm')
    ptm = ctx.get('ptm')
    mean_plddt = ctx.get('mean_plddt')
    ranking_score = ctx.get('ranking_score')
    fraction_disordered = ctx.get('fraction_disordered')
    has_clash = ctx.get('has_clash')
    chain_iptm = ctx.get('chain_iptm')
    chain_ptm = ctx.get('chain_ptm')
    pae_matrix = ctx.get('pae_matrix')
    contact_probs = ctx.get('contact_probs')
    plddt_arr = ctx.get('plddt_arr')
    token_res_ids = ctx.get('token_res_ids')
    pair_metrics = ctx.get('pair_metrics')
    best_pae_inter = ctx.get('best_pae_inter')
    best_pair = ctx.get('best_pair')
    best_cr = ctx.get('best_cr')
    motifs = ctx.get('motifs')
    seq_status = ctx.get('seq_status')
    seq_status_legacy = ctx.get('seq_status_legacy')
    seq_chains = ctx.get('seq_chains')
    chain_layout_source = ctx.get('chain_layout_source')
    seq_fingerprint = ctx.get('seq_fingerprint')
    seq_seed_fingerprint = ctx.get('seq_seed_fingerprint')
    model_seeds = ctx.get('model_seeds')
    completeness = ctx.get('completeness')
    truncation_info = ctx.get('truncation_info')
    ambiguity_info = ctx.get('ambiguity_info')
    validation_warnings = ctx.get('validation_warnings')
    ranking_samples = ctx.get('ranking_samples')

    return {
                'job_name':      job_dir.name,
                'job_dir':       str(job_dir),
                # v2.8.1 lazy-loading file paths.  The scan stores only paths
                # and lightweight metrics; heavy matrices are loaded on demand.
                'summary_path':   str(sum_path) if sum_path else None,
                'conf_path':      str(conf_path) if conf_path else None,
                'model_path':     str(model_cif) if model_cif else None,
                'ranking_csv_path': str(ranking_csv) if ranking_csv else None,
                'input_json_path': str(input_json) if input_json else None,
                'heavy_loaded':   (not lazy),
                'orf_names':     orf_names,
                'chain_order':   chain_order,
                'chain_lens':    chain_lens,
                'chain_to_orf':  chain_to_orf,
                'n_chains':      n_chains,
                # Global metrics
                'iptm':          iptm,
                'ptm':           ptm,
                'mean_plddt':    mean_plddt,
                'ranking_score': ranking_score,
                'fraction_disordered': fraction_disordered,
                'has_clash':     has_clash,
                # Per-chain metrics (v1.15)
                'chain_iptm':    chain_iptm,
                'chain_ptm':     chain_ptm,
                # Arrays for plotting
                'pae_matrix':    pae_matrix,
                'contact_probs': contact_probs,
                'plddt_arr':     plddt_arr,
                'token_res_ids': token_res_ids,
                # Pairwise interface metrics
                'pair_metrics':  pair_metrics,
                'pae_inter':     best_pae_inter,
                'best_pair':     best_pair,
                'contact_region': best_cr,
                # Focal metrics (v2.0) — correctly identify domain-limited interactions
                'pae_min_inter':  (pair_metrics[best_pair].get('pair_pae_min') or
                                   pair_metrics[best_pair].get('pae_min'))
                                  if best_pair in pair_metrics else None,
                'cp_iptm_inter':  pair_metrics[best_pair].get('pair_iptm')
                                  if best_pair in pair_metrics else None,
                'contact_frac':   pair_metrics[best_pair].get('contact_frac')
                                  if best_pair in pair_metrics else None,
                # Hotspot metrics (v2.7)
                'hotspot':        pair_metrics[best_pair].get('hotspot', {})
                                  if best_pair in pair_metrics else {},
                # Motifs (v1.16)
                'motifs':        motifs,
                # Sequence verification (v2.0)
                'seq_status':    seq_status,
                'seq_status_legacy': seq_status_legacy,
                'seq_chains':    seq_chains,
                # v2.7 — diagnostic: which fallback tier produced the chain
                # layout ('confidences', 'pae_matrix', 'input_json',
                # 'summary_chain_iptm', 'cif', or 'none').  Useful when
                # debugging why some folders show n_chains=0.
                'chain_layout_source': chain_layout_source,
                # ── v2.5 Tier 1: fingerprinting & duplicate detection ─────
                'seq_fingerprint':       seq_fingerprint,
                'seq_seed_fingerprint':  seq_seed_fingerprint,
                'model_seeds':           list(model_seeds),
                # Filled in by _af3a_assign_duplicate_groups after scan:
                'duplicate_group_id':    -1,    # -1 = singleton
                'duplicate_role':        None,  # 'canonical' | 'replicate' | None
                'duplicate_count':       1,
                # ── v2.5 Tier 2: validation taxonomy ──────────────────────
                'completeness':          completeness,
                'truncation_info':       truncation_info,
                'ambiguity_info':        ambiguity_info,
                'validation_warnings':   list(validation_warnings),
                # Extras
                'ranking_samples': ranking_samples,
                'partner_name':  (orf_names[1] if len(orf_names) > 1
                                  else chain_order[1] if len(chain_order) > 1
                                  else '-'),
            }
