"""
Sequence Search Module.
Handles homology searches (BLAST) and domain identification (HMM).
"""
import logging
from typing import List, Dict, Any

logger = logging.getLogger(__name__)

class SequenceSearcher:
    def __init__(self):
        pass

    def kmer_blast(self, query, subjects, params=None):
        """BLAST rápido: k-mer index + diagonal filter + banded SW."""
        if params is None: params = {}
        threshold = params.get('threshold', 30)
        gap_open = params.get('gap_open', -11)
        gap_extend = params.get('gap_extend', -1)
        word_size = params.get('word_size', 4)
        min_diag_hits = params.get('min_diag_hits', 2)
        evalue_max = params.get('evalue', 0.05)

        query = query.upper().replace('*', '').strip()
        if not query: return []
        q_len = len(query)
        db_size = sum(len(s.replace('*','')) for s in subjects)

        # Fase 1: indexar k-mers do query
        q_index = defaultdict(list)
        for i in range(q_len - word_size + 1):
            kmer = query[i:i+word_size]
            q_index[kmer].append(i)

        hits = []
        band_width = 15

        for idx, subject in enumerate(subjects):
            subj = subject.upper().replace('*', '').strip()
            if not subj or len(subj) < word_size: continue
            s_len = len(subj)

            # Fase 2: contar hits por diagonal
            diag_hits = defaultdict(list)
            for j in range(s_len - word_size + 1):
                kmer = subj[j:j+word_size]
                if kmer in q_index:
                    for i in q_index[kmer]:
                        diag = i - j
                        diag_hits[diag].append((i, j))

            # Keep only diagonals with enough seed hits
            hot_diags = [(d, pts) for d, pts in diag_hits.items()
                         if len(pts) >= min_diag_hits]
            if not hot_diags: continue

            # Phase 3: banded SW on the best diagonals
            best_diag = max(hot_diags, key=lambda x: len(x[1]))
            diag_val, seed_pts = best_diag

            # Determine the band around the best diagonal
            min_i = min(p[0] for p in seed_pts)
            max_i = max(p[0] for p in seed_pts)
            min_j = min(p[1] for p in seed_pts)
            max_j = max(p[1] for p in seed_pts)

            # Expand the alignment window
            q_start = max(0, min_i - band_width)
            q_end = min(q_len, max_i + word_size + band_width)
            s_start = max(0, min_j - band_width)
            s_end = min(s_len, max_j + word_size + band_width)

            q_sub = query[q_start:q_end]
            s_sub = subj[s_start:s_end]
            ql = len(q_sub); sl = len(s_sub)

            if ql < 3 or sl < 3: continue

            # Banded Smith-Waterman with affine gaps
            H = [[0]*(sl+1) for _ in range(ql+1)]
            tb = [[0]*(sl+1) for _ in range(ql+1)]
            max_score = 0; mi = mj = 0

            for i in range(1, ql+1):
                j_center = i - diag_val + s_start - q_start
                j_lo = max(1, j_center - band_width)
                j_hi = min(sl, j_center + band_width)
                for j in range(j_lo, j_hi + 1):
                    diag_s = H[i-1][j-1] + self._blosum_score(q_sub[i-1], s_sub[j-1])
                    # Affine gap: open+extend for new gap, extend for continuation
                    up_open = H[i-1][j] + gap_open + gap_extend
                    left_open = H[i][j-1] + gap_open + gap_extend
                    best = max(0, diag_s, up_open, left_open)
                    H[i][j] = best
                    if best <= 0: tb[i][j] = 0
                    elif best == diag_s: tb[i][j] = 1
                    elif best == up_open: tb[i][j] = 2
                    else: tb[i][j] = 3
                    if best > max_score:
                        max_score = best; mi = i; mj = j

            if max_score <= 0: continue

            # Traceback
            aln_q, aln_s, aln_m = [], [], []
            i, j = mi, mj
            ids = pos = gaps = al = 0
            while i > 0 and j > 0 and H[i][j] > 0:
                t = tb[i][j]
                if t == 1:
                    qa, sa = q_sub[i-1], s_sub[j-1]
                    aln_q.append(qa); aln_s.append(sa)
                    sc = self._blosum_score(qa, sa)
                    if qa == sa: aln_m.append(qa); ids += 1; pos += 1
                    elif sc > 0: aln_m.append('+'); pos += 1
                    else: aln_m.append(' ')
                    al += 1; i -= 1; j -= 1
                elif t == 2:
                    aln_q.append(q_sub[i-1]); aln_s.append('-'); aln_m.append(' ')
                    al += 1; gaps += 1; i -= 1
                elif t == 3:
                    aln_q.append('-'); aln_s.append(s_sub[j-1]); aln_m.append(' ')
                    al += 1; gaps += 1; j -= 1
                else: break

            aln_q.reverse(); aln_s.reverse(); aln_m.reverse()
            id_pct = (ids/al*100) if al > 0 else 0
            pos_pct = (pos/al*100) if al > 0 else 0
            cov = (al/q_len*100) if q_len > 0 else 0
            evalue = self.calc_evalue(max_score, q_len, db_size)

            if id_pct >= threshold and evalue <= evalue_max:
                hits.append({
                    'orf_index': idx, 'identity': round(id_pct, 1),
                    'positives': round(pos_pct, 1), 'score': max_score,
                    'aln_length': al, 'identities_count': ids,
                    'positives_count': pos, 'gaps': gaps,
                    'q_start': q_start + i + 1, 'q_end': q_start + mi,
                    's_start': s_start + j + 1, 's_end': s_start + mj,
                    'coverage': round(cov, 1), 'evalue': evalue,
                    'aln_query': ''.join(aln_q), 'aln_midline': ''.join(aln_m),
                    'aln_subject': ''.join(aln_s),
                })
        return sorted(hits, key=lambda x: x['score'], reverse=True)

    # ═══════ METHOD 2: NCBI BLAST+ (EXTERNAL) ═══════

    def run_ncbi_blast(self, query_protein, orfs, params=None):
        """Run NCBI BLAST+ blastp via subprocess."""
        if not BACKENDS.get('blast+', {}).get('available'):
            return None  # fallback para Python
        if params is None: params = {}
        evalue_thresh = params.get('evalue', 0.05)
        matrix = params.get('matrix', 'BLOSUM62')
        word_size = params.get('word_size', 5)
        max_targets = params.get('max_targets', 100)
        gap_open = abs(params.get('gap_open', -11))
        gap_extend = abs(params.get('gap_extend', -1))
        low_complexity = params.get('low_complexity', True)

        tmpdir = tempfile.mkdtemp(prefix='blast_')
        try:
            # Write query FASTA
            qfile = os.path.join(tmpdir, 'query.fasta')
            with open(qfile, 'w') as f:
                f.write(">query\n")
                for i in range(0, len(query_protein), 80):
                    f.write(query_protein[i:i+80] + "\n")

            # Write subject FASTA (ORF proteins)
            sfile = os.path.join(tmpdir, 'subjects.fasta')
            with open(sfile, 'w') as f:
                for i, orf in enumerate(orfs):
                    prot = orf['protein'].rstrip('*')
                    if prot:
                        f.write(f">ORF{i+1}\n")
                        for j in range(0, len(prot), 80):
                            f.write(prot[j:j+80] + "\n")

            # Run blastp (protein vs ORF proteins)
            outfile = os.path.join(tmpdir, 'results.txt')
            cmd = [
                'blastp',
                '-query', qfile,
                '-subject', sfile,
                '-outfmt', '6 sseqid score pident positive length qstart qend sstart send evalue gaps',
                '-evalue', str(evalue_thresh),
                '-matrix', matrix,
                '-word_size', str(word_size),
                '-gapopen', str(gap_open),
                '-gapextend', str(gap_extend),
                '-max_target_seqs', str(max_targets),
                '-seg', 'yes' if low_complexity else 'no',
                '-out', outfile,
            ]
            subprocess.run(cmd, capture_output=True, text=True, timeout=120)

            # Parse tabular output
            hits = []
            if os.path.exists(outfile):
                with open(outfile) as f:
                    for line in f:
                        parts = line.strip().split('\t')
                        if len(parts) >= 11:
                            orf_name = parts[0]
                            orf_idx = int(orf_name.replace('ORF','')) - 1
                            hits.append({
                                'orf_index': orf_idx,
                                'score': int(float(parts[1])),
                                'identity': round(float(parts[2]), 1),
                                'positives': round(float(parts[3]), 1),
                                'aln_length': int(parts[4]),
                                'q_start': int(parts[5]),
                                'q_end': int(parts[6]),
                                's_start': int(parts[7]),
                                's_end': int(parts[8]),
                                'evalue': float(parts[9]),
                                'gaps': int(parts[10]),
                                'identities_count': int(round(float(parts[2])*int(parts[4])/100)),
                                'positives_count': int(round(float(parts[3])*int(parts[4])/100)),
                                'coverage': round(int(parts[4])/len(query_protein)*100, 1),
                                'aln_query': '', 'aln_midline': '', 'aln_subject': '',
                            })

            # Retrieve detailed alignments for top hits
            if hits:
                aln_file = os.path.join(tmpdir, 'results_aln.txt')
                cmd2 = [
                    'blastp', '-query', qfile, '-subject', sfile,
                    '-outfmt', '0', '-evalue', str(evalue_thresh),
                    '-matrix', matrix, '-word_size', str(word_size),
                    '-gapopen', str(gap_open), '-gapextend', str(gap_extend),
                    '-max_target_seqs', str(min(10, max_targets)),
                    '-out', aln_file,
                ]
                subprocess.run(cmd2, capture_output=True, text=True, timeout=120)
                if os.path.exists(aln_file):
                    self._parse_blast_alignments(aln_file, hits)

            return sorted(hits, key=lambda x: x['score'], reverse=True)

        except Exception:
            return None
        finally:
            import shutil as sh
            sh.rmtree(tmpdir, ignore_errors=True)

    def _parse_blast_alignments(self, aln_file, hits):
        """Parse alignments from BLAST output format 0."""
        try:
            with open(aln_file) as f:
                content = f.read()
            # blastp -subject writes "> ORF42" (with space); accept both forms
            blocks = re.split(r'>\s*ORF(\d+)', content)
            hit_dict = {h['orf_index']: h for h in hits}
            for i in range(1, len(blocks) - 1, 2):
                try:
                    orf_idx = int(blocks[i]) - 1
                except ValueError:
                    continue
                block = blocks[i + 1]
                q_seqs, s_seqs = [], []
                for qm in re.finditer(r'Query\s+\d+\s+([A-Z\-]+)\s+\d+', block):
                    q_seqs.append(qm.group(1))
                for sm in re.finditer(r'Sbjct\s+\d+\s+([A-Z\-]+)\s+\d+', block):
                    s_seqs.append(sm.group(1))
                if orf_idx in hit_dict and q_seqs and s_seqs:
                    full_q = ''.join(q_seqs)
                    full_s = ''.join(s_seqs)
                    # Extrair midline real do BLAST (linha entre Query e Sbjct)
                    mid_lines = re.findall(
                        r'Query\s+\d+\s+[A-Z\-]+\s+\d+\n([ A-Z\+]*)\nSbjct', block)
                    if mid_lines:
                        raw_mid = ''.join(mid_lines)
                        # Alinhar comprimento com a sequência
                        mid = list(raw_mid[:len(full_q)].ljust(len(full_q)))
                    else:
                        # Fallback: reconstruir via BLOSUM
                        mid = []
                        for qa, sa in zip(full_q, full_s):
                            if qa == sa:
                                mid.append(qa)
                            elif qa != '-' and sa != '-' and self._blosum_score(qa, sa) > 0:
                                mid.append('+')
                            else:
                                mid.append(' ')
                    hit_dict[orf_idx]['aln_query']   = full_q
                    hit_dict[orf_idx]['aln_subject']  = full_s
                    hit_dict[orf_idx]['aln_midline']  = ''.join(mid)
        except Exception:
            pass

    # ═══════ METHOD 3: FULL SMITH-WATERMAN ═══════

    def sw_blast(self, query, subjects, params=None):
        """Full Smith-Waterman alignment (most sensitive, slower)."""
        if params is None: params = {}
        threshold = params.get('threshold', 30)
        gap_open = params.get('gap_open', -11)
        gap_extend = params.get('gap_extend', -1)
        evalue_max = params.get('evalue', 0.05)
        query = query.upper().replace('*', '').strip()
        if not query: return []
        q_len = len(query)
        db_size = sum(len(s.replace('*','')) for s in subjects)
        hits = []
        query_kmers = set()
        for ki in range(q_len - 2):
            query_kmers.add(query[ki:ki+3])

        for idx, subject in enumerate(subjects):
            sc = subject.upper().replace('*', '').strip()
            if not sc: continue
            sl = len(sc)
            # Pre-filtro 3-mer
            found = False
            for si in range(sl - 2):
                if sc[si:si+3] in query_kmers: found = True; break
            if not found: continue

            NEG_INF = float('-inf')
            E_prev = [NEG_INF]*(sl+1)
            H = [[0]*(sl+1) for _ in range(q_len+1)]
            tb = [[0]*(sl+1) for _ in range(q_len+1)]
            ms = mi = mj = 0

            for i in range(1, q_len+1):
                F_val = NEG_INF
                E_curr = [NEG_INF]*(sl+1)
                for j in range(1, sl+1):
                    diag = H[i-1][j-1] + self._blosum_score(query[i-1], sc[j-1])
                    E_curr[j] = max(H[i-1][j]+gap_open+gap_extend, E_prev[j]+gap_extend)
                    f_o = H[i][j-1]+gap_open+gap_extend
                    F_val = max(f_o, F_val+gap_extend)
                    best = max(0, diag, E_curr[j], F_val)
                    H[i][j] = best
                    if best <= 0: tb[i][j] = 0
                    elif best == diag: tb[i][j] = 1
                    elif best == E_curr[j]: tb[i][j] = 2
                    else: tb[i][j] = 3
                    if best > ms: ms = best; mi = i; mj = j
                E_prev = E_curr

            if ms <= 0: continue
            aln_q, aln_s, aln_m = [], [], []
            i, j = mi, mj
            ids = pos = gaps = al = 0
            while i > 0 and j > 0 and H[i][j] > 0:
                t = tb[i][j]
                if t == 1:
                    qa, sa = query[i-1], sc[j-1]
                    aln_q.append(qa); aln_s.append(sa)
                    s = self._blosum_score(qa, sa)
                    if qa == sa: aln_m.append(qa); ids += 1; pos += 1
                    elif s > 0: aln_m.append('+'); pos += 1
                    else: aln_m.append(' ')
                    al += 1; i -= 1; j -= 1
                elif t == 2:
                    aln_q.append(query[i-1]); aln_s.append('-'); aln_m.append(' ')
                    al += 1; gaps += 1; i -= 1
                elif t == 3:
                    aln_q.append('-'); aln_s.append(sc[j-1]); aln_m.append(' ')
                    al += 1; gaps += 1; j -= 1
                else: break
            aln_q.reverse(); aln_s.reverse(); aln_m.reverse()
            id_p = (ids/al*100) if al > 0 else 0
            pos_p = (pos/al*100) if al > 0 else 0
            cov = (al/q_len*100) if q_len > 0 else 0
            ev = self.calc_evalue(ms, q_len, db_size)
            if id_p >= threshold and ev <= evalue_max:
                hits.append({
                    'orf_index': idx, 'identity': round(id_p,1), 'positives': round(pos_p,1),
                    'score': ms, 'aln_length': al, 'identities_count': ids,
                    'positives_count': pos, 'gaps': gaps,
                    'q_start': i+1, 'q_end': mi, 's_start': j+1, 's_end': mj,
                    'coverage': round(cov,1), 'evalue': ev,
                    'aln_query': ''.join(aln_q), 'aln_midline': ''.join(aln_m),
                    'aln_subject': ''.join(aln_s),
                })
        return sorted(hits, key=lambda x: x['score'], reverse=True)

    # ═══════ HMM SEARCH ═══════

    def hmm_scan_orfs(self, hmm_file, orfs, params=None):
        """Search HMM profile against all ORFs. Uses HMMER3 if available."""
        if params is None: params = {}
        if BACKENDS.get('hmmer3', {}).get('available'):
            return self._hmmer3_search(hmm_file, orfs, params)
        return self._pssm_scan(hmm_file, orfs, params)

    def _hmmer3_search(self, hmm_file, orfs, params=None):
        """Chama hmmsearch do HMMER3. Usa --domtblout + -A para coordenadas e alinhamentos."""
        if params is None: params = {}
        hmm_evalue  = params.get('hmm_evalue', 10.0)
        hmm_score   = params.get('hmm_score_thresh', None)
        dom_evalue  = params.get('hmm_dom_evalue', 10.0)
        use_wsl     = BACKENDS.get('hmmer3', {}).get('wsl', False)
        tmpdir      = tempfile.mkdtemp(prefix='hmm_')
        try:
            # Escrever proteínas
            sfile = os.path.join(tmpdir, 'orfs.fasta')
            n_written = 0
            with open(sfile, 'w') as f:
                for i, orf in enumerate(orfs):
                    prot = orf['protein'].rstrip('*')
                    if prot and len(prot) >= 10:
                        f.write(f">ORF{i+1}\n")
                        for j in range(0, len(prot), 80):
                            f.write(prot[j:j+80] + "\n")
                        n_written += 1

            domtbl_file = os.path.join(tmpdir, 'results.domtbl')
            aln_file    = os.path.join(tmpdir, 'results.sto')

            if use_wsl:
                import shutil as sh_copy
                hmm_copy = os.path.join(tmpdir, 'profile.hmm')
                sh_copy.copy2(hmm_file, hmm_copy)

                def to_wsl_path(p):
                    p = os.path.abspath(p).replace('\\', '/')
                    if len(p) > 1 and p[1] == ':':
                        return f"/mnt/{p[0].lower()}{p[2:]}"
                    return p

                wsl_sfile   = to_wsl_path(sfile)
                wsl_hmm     = to_wsl_path(hmm_copy)
                wsl_domtbl  = to_wsl_path(domtbl_file)
                wsl_aln     = to_wsl_path(aln_file)

                hmm_cmd  = f'hmmsearch --domtblout "{wsl_domtbl}" -A "{wsl_aln}"'
                hmm_cmd += f' -E {hmm_evalue} --domE {dom_evalue}'
                if hmm_score is not None:
                    hmm_cmd += f' -T {hmm_score}'
                hmm_cmd += f' "{wsl_hmm}" "{wsl_sfile}"'
                cmd = ['wsl', 'bash', '-c', hmm_cmd]
            else:
                cmd = ['hmmsearch',
                       '--domtblout', domtbl_file,
                       '-A', aln_file,
                       '-E', str(hmm_evalue),
                       '--domE', str(dom_evalue)]
                if hmm_score is not None:
                    cmd.extend(['-T', str(hmm_score)])
                cmd.extend([hmm_file, sfile])

            result = subprocess.run(cmd, capture_output=True, text=True, timeout=180)

            # ── Parse --domtblout (coordenadas precisas) ──────────────
            hits = []
            if os.path.exists(domtbl_file):
                with open(domtbl_file) as f:
                    for line in f:
                        if line.startswith('#'): continue
                        parts = line.split()
                        if len(parts) < 23: continue
                        orf_name = parts[0]
                        try:
                            orf_idx = int(orf_name.replace('ORF', '')) - 1
                            # domtblout columns (0-based):
                            # 0=target, 2=tlen, 3=query, 5=qlen,
                            # 11=c-evalue, 12=i-evalue, 13=score, 14=bias
                            # 15=hmm_from, 16=hmm_to
                            # 17=ali_from, 18=ali_to
                            # 19=env_from, 20=env_to
                            hmm_from = int(parts[15])
                            hmm_to   = int(parts[16])
                            ali_from = int(parts[17])
                            ali_to   = int(parts[18])
                            score    = float(parts[13])
                            evalue   = float(parts[12])  # i-evalue (independent)
                            tlen     = int(parts[2])
                            qlen     = int(parts[5])
                            hits.append({
                                'orf_index':  orf_idx,
                                'orf_name':   orf_name,
                                'hmm_name':   parts[3],
                                'score':      score,
                                'evalue':     evalue,
                                'bias':       float(parts[14]),
                                'hmm_from':   hmm_from,
                                'hmm_to':     hmm_to,
                                'ali_from':   ali_from,
                                'ali_to':     ali_to,
                                'hmm_len':    qlen,
                                'target_len': tlen,
                                # Região formatada para exibição
                                'match_region': f"HMM:{hmm_from}-{hmm_to}/{qlen}  Prot:{ali_from}-{ali_to}/{tlen}",
                            })
                        except (ValueError, IndexError):
                            continue

            # ── Parse -A Stockholm (alinhamentos) ─────────────────────
            aln_dict = {}
            if os.path.exists(aln_file) and os.path.getsize(aln_file) > 10:
                try:
                    aln_dict = self._parse_stockholm_aln(aln_file)
                except Exception:
                    pass

            # Anexar alinhamento a cada hit
            for h in hits:
                # Stockholm usa nomes como "ORF104/28-251" — normalizar
                orf_name_bare  = h['orf_name']                       # "ORF104"
                orf_name_range = f"{orf_name_bare}/{h.get('ali_from','')}-{h.get('ali_to','')}"

                aln = (aln_dict.get(orf_name_bare)
                    or aln_dict.get(orf_name_range)
                    or next((v for k, v in aln_dict.items()
                             if k.startswith(orf_name_bare + '/')), None)
                    or {})

                if aln:
                    h['aln_hmm']    = aln.get('hmm',    '')
                    h['aln_target'] = aln.get('target', '')
                    h['aln_match']  = aln.get('match',  '')
                else:
                    # Fallback: mostrar a subsequência proteica alinhada
                    # even without a Stockholm alignment file
                    oi = h.get('orf_index', -1)
                    if 0 <= oi < len(orfs):
                        prot      = orfs[oi]['protein'].rstrip('*')
                        ali_from  = h.get('ali_from', 1)
                        ali_to    = h.get('ali_to',   len(prot))
                        try:
                            subseq = prot[int(ali_from) - 1 : int(ali_to)]
                        except (TypeError, ValueError):
                            subseq = prot[:50]
                        h['aln_hmm']    = ''        # HMM consensus não disponível sem -A
                        h['aln_target'] = subseq
                        h['aln_match']  = ''
                    else:
                        h['aln_hmm']    = ''
                        h['aln_target'] = ''
                        h['aln_match']  = ''

            if not hits:
                return [{'error': f"0 hits. (ORFs: {n_written}, return: {result.returncode})"}]
            return sorted(hits, key=lambda x: x['score'], reverse=True)

        except subprocess.TimeoutExpired:
            return [{'error': 'Timeout (180s).'}]
        except Exception as e:
            return [{'error': str(e)}]
        finally:
            import shutil as sh
            sh.rmtree(tmpdir, ignore_errors=True)

    def _parse_stockholm_aln(self, sto_file: str) -> dict:
        """
        Parse Stockholm output from hmmsearch -A.
        Returns dict: orf_name → {'hmm': str, 'target': str, 'match': str}
        """
        result = {}
        with open(sto_file) as f:
            content = f.read()

        # Each alignment block starts with '# STOCKHOLM 1.0' and ends with '//'
        blocks = content.split('//')
        for block in blocks:
            lines = block.strip().splitlines()
            seqs  = {}     # name → sequence fragments
            rf    = ''
            for line in lines:
                if not line.strip() or line.startswith('#=GF'): continue
                if line.startswith('#=GC RF'):
                    rf += line.split()[-1]
                    continue
                if line.startswith('#=GC PP_cons') or line.startswith('#=GC seq_cons'):
                    continue
                if line.startswith('#=GC'):
                    continue
                if line.startswith('#=GR'):
                    # per-residue annotation — skip
                    continue
                if line.startswith('#'):
                    continue
                parts = line.split(None, 1)
                if len(parts) == 2:
                    name = parts[0]
                    seq  = parts[1].strip()
                    seqs.setdefault(name, '')
                    seqs[name] += seq

            # Identify the HMM query (not an ORF)
            hmm_seq    = ''
            target_seqs = {}
            for name, seq in seqs.items():
                if name.startswith('ORF'):
                    target_seqs[name] = seq
                else:
                    hmm_seq  = seq  # name ignored — only sequence is used

            # Build midline by comparing HMM consensus vs target
            for orf_name, tgt_seq in target_seqs.items():
                if not hmm_seq:
                    result[orf_name] = {'hmm': tgt_seq, 'target': tgt_seq, 'match': ''}
                    continue
                mid = []
                for hc, tc in zip(hmm_seq, tgt_seq):
                    if hc == '-' or tc == '-':
                        mid.append(' ')
                    elif hc.upper() == tc.upper():
                        mid.append('|')
                    elif hc.upper() in 'ACDEFGHIKLMNPQRSTVWY' and tc.upper() in 'ACDEFGHIKLMNPQRSTVWY':
                        sc = self._blosum_score(hc.upper(), tc.upper())
                        mid.append('+' if sc > 0 else '.')
                    else:
                        mid.append('.')
                result[orf_name] = {
                    'hmm':    hmm_seq,
                    'target': tgt_seq,
                    'match':  ''.join(mid),
                }
        return result

    def _pssm_scan(self, hmm_file, orfs, params=None):
        """PSSM scan derived from the HMM file (fallback when HMMER3 is unavailable)."""
        if params is None: params = {}
        score_thresh = params.get('hmm_score_thresh', None)
        evalue_thresh = params.get('hmm_evalue', 10.0)
        pssm, aa_order = self._parse_hmm_to_pssm(hmm_file)
        if not pssm: return [{'error': 'Failed to parse the HMM file'}]
        hits = []
        model_len = len(pssm)

        # Calculate null model score (average across all aa per position)
        null_score = 0
        for pos_scores in pssm:
            valid = [v for v in pos_scores.values() if v > -900]
            null_score += sum(valid) / len(valid) if valid else -4
        # Threshold: must beat null by at least 0.5 nats per position
        auto_thresh = null_score + model_len * 0.5
        min_score = score_thresh if score_thresh is not None else auto_thresh

        for i, orf in enumerate(orfs):
            prot = orf['protein'].rstrip('*').upper()
            if len(prot) < model_len: continue
            best_score = float('-inf')
            best_pos = 0
            for start in range(len(prot) - model_len + 1):
                score = 0
                for k in range(model_len):
                    aa = prot[start + k]
                    score += pssm[k].get(aa, -4)
                if score > best_score:
                    best_score = score; best_pos = start
            ev = self.calc_evalue(max(0, best_score - null_score), model_len, len(prot))
            if best_score > min_score and ev <= evalue_thresh:
                hits.append({
                    'orf_index': i, 'orf_name': f'ORF{i+1}',
                    'score': round(best_score, 1), 'position': best_pos,
                    'evalue': ev, 'hmm_name': Path(hmm_file).stem,
                    'match_region': f'{best_pos+1}-{best_pos+model_len}',
                })
        return sorted(hits, key=lambda x: x['score'], reverse=True)

    def _parse_hmm_to_pssm(self, hmm_file):
        """Extract PSSM from HMMER3 Match states. Reads AA order from the file header."""
        try:
            with open(hmm_file) as f:
                content = f.read()
            if 'HMMER3' not in content: return None, None
            lines = content.split('\n')
            in_model = False; pssm = []
            aa_order = list('ACDEFGHIKLMNPQRSTVWY')  # default

            for line in lines:
                # Parse AA order from HMM header line
                if line.strip().startswith('HMM') and not line.strip().startswith('HMMER'):
                    parts = line.split()
                    if len(parts) >= 20:
                        aa_order = [p.strip() for p in parts[1:] if len(p.strip()) == 1 and p.strip().isalpha()]
                        if len(aa_order) < 20:
                            aa_order = list('ACDEFGHIKLMNPQRSTVWY')
                    in_model = True
                    continue
                if not in_model: continue
                if line.startswith('//'): break
                parts = line.split()
                # Match emission lines start with node number
                if len(parts) >= len(aa_order) + 1 and parts[0].isdigit():
                    scores = {}
                    for j, aa in enumerate(aa_order):
                        try:
                            val = parts[j+1]
                            if val == '*':
                                scores[aa] = -999
                            else:
                                # HMMER3: scores are -ln(prob/null), lower=better match
                                # Convert: higher score = better (negate)
                                scores[aa] = -float(val)
                        except (ValueError, IndexError):
                            scores[aa] = -4
                    pssm.append(scores)
            return (pssm, aa_order) if pssm else (None, None)
        except Exception:
            return None, None

