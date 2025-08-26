import sys
import pysam
import re
from collections import namedtuple

in_aln = sys.argv[1]  # bam or paf
gene_region_bed = sys.argv[2]
locus_buffer = int(sys.argv[3])

AlignRec = namedtuple('AlignRec', [
    'query_name', 'query_length', 'query_start', 'query_end', 'strand',
    'target_name', 'reference_start', 'cigartuples'
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
            locus_start_buf = locus_pos - locus_buffer
            locus_end_buf = locus_end + locus_buffer

            is_paf = False
            if in_aln.endswith('.bam') or in_aln.endswith('.cram'):
                s = pysam.AlignmentFile(in_aln, 'rb')
                rec_iter = s.fetch(reference=chrom, start=locus_start_buf, end=locus_end_buf)
            else:
                assert in_aln.endswith('.paf')
                rec_iter = open(in_aln)
                is_paf = True

            for r in rec_iter:
                if is_paf:
                    line = r.strip().split('\t')
                    if (line[5] != chrom) or (int(line[8]) < locus_pos) or (int(line[7]) > locus_end):
                        continue
                    cigar_cols = [z for z in line if z.startswith('cg:Z:')]
                    assert len(cigar_cols) == 1
                    cigar_col = cigar_cols[0]
                    r = AlignRec(
                        line[0], int(line[1]), int(line[2]), int(line[3]), line[4], line[5], int(line[7]),
                        [(CIGAR_OPS[op], int(length)) for length, op in re.findall(r'(\d+)([MIDNSHP=XB])',
                                                                                   cigar_col.replace('cg:Z:', ''))]
                    )
                is_rev_strand = r.strand == '-' if (is_paf) else (r.flag & 0x10 != 0)

                first_ctg_coord, last_ctg_coord, hardclip_len = get_first_last_contig_coords(r, is_rev_strand, is_paf, locus_start_buf, locus_end_buf)

                if first_ctg_coord is None or last_ctg_coord is None:
                    continue

                out_name = f"\t{locus_name}" if (locus_name) else ""
                if is_rev_strand:
                    print(f"{r.query_name}/rc\t{r.query_length - last_ctg_coord - 1 + hardclip_len}\t{r.query_length - first_ctg_coord + hardclip_len}{out_name}")
                else:
                    print(f"{r.query_name}\t{first_ctg_coord + hardclip_len}\t{last_ctg_coord + hardclip_len}{out_name}")


def get_first_last_contig_coords(rec, is_rev_strand, is_paf, locus_start, locus_end):
    ref_pos = rec.reference_start
    query_seq = None if (is_paf) else rec.query_sequence.upper()
    query_length = rec.query_length if (is_paf) else len(query_seq)
    query_pos = 0
    hardclip_len = 0
    
    first_ctg_coord = None
    last_ctg_coord = None

    for op_char, length in rec.cigartuples:
        if op_char in [pysam.CEQUAL, pysam.CMATCH, pysam.CDIFF]:
            # Match or mismatch: both reference and query advance
            for i in range(length):
                # Check if this position is within our locus range
                if locus_start <= ref_pos <= locus_end:
                    if first_ctg_coord is None:
                        first_ctg_coord = query_pos
                    last_ctg_coord = query_pos
                ref_pos += 1
                query_pos += 1
            
        elif op_char == pysam.CDEL:
            # Deletion: only reference advances
            # If locus range overlaps with this deletion, include flanking query positions
            if ref_pos + length > locus_start and ref_pos < locus_end:
                if first_ctg_coord is None:
                    first_ctg_coord = query_pos
                last_ctg_coord = query_pos
            ref_pos += length
            
        elif op_char == pysam.CINS:
            # Insertion: only query advances
            # If we're within the reference range, include all inserted bases
            if locus_start <= ref_pos <= locus_end:
                for i in range(length):
                    if first_ctg_coord is None:
                        first_ctg_coord = query_pos
                    last_ctg_coord = query_pos
                    query_pos += 1
            else:
                query_pos += length
            
        elif op_char == pysam.CSOFT_CLIP:
            # Soft clip: only query advances
            query_pos += length
            
        elif op_char == pysam.CHARD_CLIP:
            # Hard clip: nothing advances
            if hardclip_len == 0 or is_rev_strand:
                hardclip_len = length
        else:
            raise ValueError(f'Unrecognized CIGAR operation {op_char}')
            
    return first_ctg_coord, last_ctg_coord, hardclip_len


if __name__ == '__main__':
    main()
