# Get bam segment in reference region
import sys
import pysam

in_bam = sys.argv[1]
gene_region_bed = sys.argv[2]
ref_chrom_seq_fname = sys.argv[3]
locus_buffer = int(sys.argv[4])

def main():
    with open(gene_region_bed) as f:
        chrom, locus_pos, locus_end = f.readline().strip().split('\t')
    locus_pos = int(locus_pos)
    locus_end = int(locus_end)
    locus_pos_buffered = locus_pos - locus_buffer
    locus_end_buffered = locus_end + locus_buffer
    with open(ref_chrom_seq_fname) as f:
        ref_chrom_seq = ''.join([z.strip() for z in f.readlines() if not z.startswith('>')]).upper()

    s = pysam.AlignmentFile(in_bam, "rb")  # 1000151-asm_h1.minimap2.bam
    for r in s.fetch(reference=chrom, start=locus_pos_buffered, end=locus_end_buffered):
        is_rev_strand = r.flag & 0x10 != 0
        aligned_pairs, hardclip_len = get_aligned_pairs_eqx(r, ref_chrom_seq, is_rev_strand)
        #print(aligned_pairs[:10], aligned_pairs[-10:])
        all_pairs = [p for p in aligned_pairs if (
            (p[3] is not None) and (p[3] >= locus_pos_buffered) and (p[3] <= locus_end_buffered)
        )]  # ctg, ref pairs where ref. is within the locus range
        if not all_pairs:
            continue
        all_ref = [z[3] for z in all_pairs]
        all_ctg = [z[4] for z in all_pairs if z[4] is not None]

        #if locus_pos in all_ref and locus_end in all_ref:  # locus contained in contig
        # first and last contig coord in locus + buffer
        ctg_start = all_ctg[0]
        ctg_end = all_ctg[-1]
        #print(r.reference_name, r.qname, r.cigarstring[:15], r.cigarstring[-15:],  r.flag, r.pos, r.query_alignment_start, ctg_start, ctg_end, all_ref[0], all_ref[-1], len(r.query_alignment_sequence))
        if is_rev_strand:
            print(f"{r.query_name}/rc\t{r.query_length - ctg_end - 1 + hardclip_len}\t{r.query_length - ctg_start + hardclip_len}")
        else:
            print(f"{r.query_name}\t{ctg_start + hardclip_len}\t{ctg_end + hardclip_len}")


def get_aligned_pairs_eqx(read, ref_seq, is_rev_strand):
    if read.is_unmapped:
        return []
    aligned_pairs = []
    ref_pos = read.reference_start
    query_seq = read.query_sequence.upper()
    query_pos = 0
    hardclip_len = 0

    for op_char, length in read.cigartuples:
        if op_char == pysam.CEQUAL:
            # Match: both reference and query advance
            for i in range(length):
                if ref_pos + i < len(ref_seq) and query_pos + i < len(query_seq):
                    aligned_pairs.append((
                        ref_seq[ref_pos + i], query_seq[query_pos + i], '=',
                        ref_pos + i, query_pos + i
                    ))
            ref_pos += length
            query_pos += length
            
        elif op_char == pysam.CDIFF:
            # Mismatch: both reference and query advance
            for i in range(length):
                if ref_pos + i < len(ref_seq) and query_pos + i < len(query_seq):
                    aligned_pairs.append((
                        ref_seq[ref_pos + i], query_seq[query_pos + i],
                        'X', ref_pos + i, query_pos + i
                    ))
            ref_pos += length
            query_pos += length
            
        elif op_char == pysam.CDEL:
            # Deletion: only reference advances
            for i in range(length):
                if ref_pos + i < len(ref_seq):
                    aligned_pairs.append((
                        ref_seq[ref_pos + i], '-', 'D', ref_pos + i, None
                    ))
            ref_pos += length
            
        elif op_char == pysam.CINS:
            # Insertion: only query advances
            for i in range(length):
                if query_pos + i < len(query_seq):
                    aligned_pairs.append(('-', query_seq[query_pos + i], 'I',
                        None, query_pos + i
                    ))
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
            raise ValueError('Unrecognized CIGAR operation')
    return (aligned_pairs, hardclip_len)


if __name__ == '__main__':
    main()
