# Get bam segment in reference region
import sys
import pysam
import re
from collections import namedtuple

in_aln = sys.argv[1]  # bam or paf
gene_region_bed = sys.argv[2]
locus_buffer = int(sys.argv[3])
ref_chrom_seq_fname = sys.argv[4] if (len(sys.argv) > 4) else None

AlignRec = namedtuple('AlignRec', [
    'query_name', 'query_length', 'query_start', 'query_end', 'strand', 'target_name', 'reference_start', 'cigartuples'
])

CIGAR_OPS = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8, 'B': 9}

def main():
  with open(gene_region_bed) as f_reg:
    for reg in f_reg:
        reg = reg.strip().split('\t')
        if len(reg) > 3:
            chrom, locus_pos, locus_end, locus_name = reg
        else:
            locus_name = None
            chrom, locus_pos, locus_end = reg
        locus_pos = int(locus_pos)
        locus_end = int(locus_end)
        locus_pos_buffered = locus_pos - locus_buffer
        locus_end_buffered = locus_end + locus_buffer

        ref_chrom_seq = None
        if ref_chrom_seq_fname:
            with open(ref_chrom_seq_fname) as f:
                ref_chrom_seq = ''.join([z.strip() for z in f.readlines() if not z.startswith('>')]).upper()

        is_paf = False
        if in_aln.endswith('.bam') or in_aln.endswith('.cram'):
            s = pysam.AlignmentFile(in_aln, "rb")
            rec_iter = s.fetch(reference=chrom, start=locus_pos_buffered, end=locus_end_buffered)
        else:
            assert in_aln.endswith('.paf')
            rec_iter = open(in_aln) 
            is_paf = True

        for r in rec_iter:
            if is_paf:
                line = r.strip().split('\t')
                # make sure overlaps locus
                if (line[5] != chrom) or (int(line[8]) < locus_pos) or (int(line[7]) > locus_end):
                    continue
                cigar_cols = [z for z in line if z.startswith('cg:Z:')]
                assert(len(cigar_cols)) == 1
                cigar_col = cigar_cols[0]
                r = AlignRec(
                    line[0], int(line[1]), int(line[2]), int(line[3]), line[4], line[5], int(line[7]),
                    [(CIGAR_OPS[op], int(length)) for length, op in re.findall(r'(\d+)([MIDNSHP=XB])', cigar_col.replace('cg:Z:', ''))]
                )
            is_rev_strand = r.strand == '-' if (is_paf) else r.flag & 0x10 != 0
            aligned_pairs, hardclip_len = get_aligned_pairs_eqx(r, is_rev_strand, is_paf, ref_chrom_seq)
            if is_paf:
                hardclip_len = (r.query_end - r.query_length) if (is_rev_strand) else r.query_start
            #print(hardclip_len, aligned_pairs[:10], aligned_pairs[-10:])
            all_pairs = [p for p in aligned_pairs if (
                (p[3] is not None) and (p[3] >= locus_pos_buffered) and (p[3] <= locus_end_buffered)
            )]  # ctg, ref pairs where ref. is within the locus range
            if not all_pairs:
                continue
            all_ref = [z[3] for z in all_pairs]
            all_ctg = [z[4] for z in all_pairs if z[4] is not None]
            if not all_ctg:
                continue
            #if locus_pos in all_ref and locus_end in all_ref:  # locus contained in contig
            # first and last contig coord in locus + buffer
            ctg_start = all_ctg[0]
            ctg_end = all_ctg[-1]
            #print(r.reference_name, r.qname, r.cigarstring[:15], r.cigarstring[-15:],  r.flag, r.pos, r.query_alignment_start, ctg_start, ctg_end, all_ref[0], all_ref[-1], len(r.query_alignment_sequence))

            # correct contig name
            ctg_name = r.query_name.replace('-adjusted', '')
            if '#' in ctg_name:
                ctg_name = ctg_name[:ctg_name.index('#')]
            out_name = f"\t{locus_name}" if (locus_name) else ""
            if is_rev_strand:
                print(f"{r.query_name}/rc\t{r.query_length - ctg_end - 1 + hardclip_len}\t{r.query_length - ctg_start + hardclip_len}{out_name}")
            else:
                print(f"{r.query_name}\t{ctg_start + hardclip_len}\t{ctg_end + hardclip_len}{out_name}")


def get_aligned_pairs_eqx(rec, is_rev_strand, is_paf, ref_seq=None):
    #if rec.is_unmapped:
    #    return []
    aligned_pairs = []
    ref_pos = rec.reference_start
    query_seq = None if (is_paf) else rec.query_sequence.upper()
    query_length = rec.query_length if (is_paf) else len(query_seq)
    query_pos = 0
    hardclip_len = 0

    for op_char, length in rec.cigartuples:
        if op_char in [pysam.CEQUAL, pysam.CMATCH]:
            # Match: both reference and query advance
            for i in range(length):
                try:
                    aligned_pairs.append((
                        ref_seq[ref_pos + i] if (ref_seq) else None,
                        query_seq[query_pos + i] if (query_seq) else None,
                        '=', ref_pos + i, query_pos + i
                    ))
                except IndexError:
                    break
            ref_pos += length
            query_pos += length
            
        elif op_char == pysam.CDIFF:
            # Mismatch: both reference and query advance
            for i in range(length):
                try:
                    aligned_pairs.append((
                        ref_seq[ref_pos + i] if (ref_seq) else None,
                        query_seq[query_pos + i] if (query_seq) else None,
                        'X', ref_pos + i, query_pos + i
                    ))
                except IndexError:
                    break
            ref_pos += length
            query_pos += length
            
        elif op_char == pysam.CDEL:
            # Deletion: only reference advances
            for i in range(length):
                try:
                    aligned_pairs.append((
                        ref_seq[ref_pos + i] if (ref_seq) else None, '-', 'D',
                        ref_pos + i, None
                    ))
                except IndexError:
                    break
            ref_pos += length
            
        elif op_char == pysam.CINS:
            # Insertion: only query advances
            for i in range(length):
                try:
                    aligned_pairs.append((
                        '-', query_seq[query_pos + i] if (query_seq) else None, 'I',
                        None, query_pos + i
                    ))
                except IndexError:
                    break
            query_pos += length
            
        elif op_char == pysam.CSOFT_CLIP:
            # Soft clip: only query advances, bases in query
            query_pos += length
            
        elif op_char == pysam.CHARD_CLIP:
            # Hard clip: nothing advances
            # use begin/ end hardclip depending on strand
            if hardclip_len == 0 or is_rev_strand:
                hardclip_len = length
        else:
            raise ValueError(f'Unrecognized CIGAR operation {op_char}')
    return (aligned_pairs, hardclip_len)


if __name__ == '__main__':
    main()


