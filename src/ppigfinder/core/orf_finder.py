"""
ORF Finder Module.
Handles the identification of Open Reading Frames in DNA sequences
using standard codon translation or Pyrodigal.
"""
import logging
from typing import List, Dict, Any, Optional

logger = logging.getLogger(__name__)

class ORFFinder:
    def __init__(self, min_aa_length: int = 30):
        self.min_aa_length = min_aa_length

    def find_orfs(self, dna_sequence, min_aa=30, start_codons=None):
        if start_codons is None: start_codons = {'ATG','GTG','TTG'}
        stop_codons = {'TAA','TAG','TGA'}
        min_len = min_aa * 3
        orfs = []
        for frame in range(3):
            for strand_seq, strand_name in [(dna_sequence, '+'),
                                            (self.reverse_complement(dna_sequence), '-')]:
                i = frame
                while i < len(strand_seq) - 2:
                    codon = strand_seq[i:i+3]
                    if codon in start_codons:
                        j = i + 3
                        while j < len(strand_seq):
                            if strand_seq[j:j+3] in stop_codons:
                                length = j + 3 - i
                                if length >= min_len:
                                    dna = strand_seq[i:j+3]
                                    protein = self.translate(dna)
                                    orfs.append({
                                        'frame': frame + (3 if strand_name == '-' else 0),
                                        'strand': strand_name,
                                        'start': i if strand_name == '+' else len(dna_sequence) - (j + 3),
                                        'end': j + 3 if strand_name == '+' else len(dna_sequence) - i,
                                        'dna': dna, 'protein': protein, 'length': length,
                                        'gc': self.gc_content(dna),
                                        'domains': [],
                                        'neighborhood': [], 'candidate_score': 0.0,
                                        'source': '6frame',
                                    })
                                i = j; break
                            j += 3
                    i += 3
        # Sort by genomic start position 5'→3' so ORF numbers increase from
        # the molecule origin toward the end (lower number = closer to 5')
        orfs.sort(key=lambda x: x['start'])
        return orfs



    def find_orfs_pyrodigal(self, dna_sequence, meta=True, min_aa=30, closed_ends=False):
        """
        Find ORFs using Pyrodigal (Python binding for Prodigal gene caller).

        Pyrodigal uses dynamic programming on GC-content, RBS motifs, and
        coding potential to predict real protein-coding genes — much more
        accurate than simple start→stop scanning.

        Parameters
        ----------
        dna_sequence : str — DNA sequence (uppercase)
        meta : bool — True for metagenomic mode (no training needed),
                      False for single-genome mode (trains on the sequence)
        min_aa : int — minimum protein length in amino acids
        closed_ends : bool — allow genes to run off edges of the sequence

        Returns
        -------
        list of ORF dicts (same format as find_orfs)
        """
        if not PYRODIGAL_AVAILABLE:
            raise ImportError(
                "Pyrodigal not installed.\n\n"
                "Install with:\n"
                "  pip install pyrodigal\n\n"
                "Or use conda:\n"
                "  conda install -c bioconda pyrodigal"
            )

        orfs = []
        seq = dna_sequence.upper()

        if meta:
            # Metagenomic mode — uses pre-trained models, no training needed
            gene_finder = pyrodigal.GeneFinder(meta=True, closed=closed_ends,
                                               min_gene=min_aa * 3)
        else:
            # Single-genome mode — trains on the input sequence
            gene_finder = pyrodigal.GeneFinder(meta=False, closed=closed_ends,
                                               min_gene=min_aa * 3)
            gene_finder.train(seq.encode() if isinstance(seq, str) else seq)

        # Run gene prediction
        genes = gene_finder.find_genes(seq.encode() if isinstance(seq, str) else seq)

        for gene in genes:
            # Pyrodigal uses 1-based coordinates
            start_0 = gene.begin - 1   # convert to 0-based
            end_0 = gene.end           # end is already exclusive-like in our format
            strand = '+' if gene.strand == 1 else '-'

            # Calculate frame (0-based)
            if strand == '+':
                frame = start_0 % 3
            else:
                frame = (len(seq) - end_0) % 3 + 3  # frames 3,4,5 for minus strand

            # Extract DNA subsequence
            if strand == '+':
                dna_sub = seq[start_0:end_0]
            else:
                dna_sub = self.reverse_complement(seq[start_0:end_0])

            protein = self.translate(dna_sub)
            gc = self.gc_content(dna_sub)

            # Pyrodigal confidence score (0-100)
            try:
                cscore = gene.confidence()
            except (AttributeError, TypeError):
                cscore = gene.cscore if hasattr(gene, 'cscore') else 0.0

            # RBS motif if available
            try:
                rbs_motif = gene.rbs_motif
            except AttributeError:
                rbs_motif = None

            orfs.append({
                'frame': frame,
                'strand': strand,
                'start': start_0,
                'end': end_0,
                'dna': dna_sub,
                'protein': protein,
                'length': end_0 - start_0,
                'gc': gc,
                'domains': [],
                'neighborhood': [],
                'candidate_score': 0.0,
                'source': 'pyrodigal',
                'pyrodigal_score': round(cscore, 2) if cscore else 0.0,
                'rbs_motif': rbs_motif or '',
                'partial': getattr(gene, 'partial_begin', False) or getattr(gene, 'partial_end', False),
            })

        orfs.sort(key=lambda x: x['start'])
        return orfs

