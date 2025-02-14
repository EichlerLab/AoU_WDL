# Get bam segment in reference region
import sys
import pysam

in_bam = sys.argv[1]
gene_region_bed = sys.argv[2]
locus_buffer = int(sys.argv[3])

s = pysam.AlignmentFile(in_bam, "rb")  # 1000151-asm_h1.minimap2.bam

with open(gene_region_bed) as f:  # Updated to allow multi-line bed
    for line in f:
        reg = line.strip().split('\t')
        if len(reg) > 3:
            chrom, locus_pos, locus_end, locus_name = reg
        else:
            locus_name = None
            chrom, locus_pos, locus_end = reg
        locus_pos = int(locus_pos)
        locus_end = int(locus_end)
        locus_pos_buffered = locus_pos - locus_buffer
        locus_end_buffered = locus_end + locus_buffer

        for r in s.fetch(reference=chrom, start=locus_pos_buffered, end=locus_end_buffered):
            is_rev_strand = r.flag & 0x10 != 0
            all_pairs = [p for p in r.get_aligned_pairs() if (
                (p[1] is not None) and (p[1] >= locus_pos_buffered) and (p[1] <= locus_end_buffered)
            )]  # ctg, ref pairs within locus + buffer
            if not all_pairs:
                continue
            all_ref = [z[1] for z in all_pairs]
            all_ctg = [z[0] for z in all_pairs if z[0] is not None]
            if not all_ctg:
                continue

            ctg_name = r.query_name.replace('-adjusted', '')
            if '#' in ctg_name:
                ctg_name = ctg_name[:ctg_name.index('#')]
            if locus_pos in all_ref and locus_end in all_ref:  # locus contained in contig
                # first and last contig coord in locus + buffer
                ctg_start = all_ctg[0]
                ctg_end = all_ctg[-1]
                out_name = f"\t{locus_name}" if (locus_name) else ""
                if is_rev_strand:
                    print(f"{ctg_name}/rc\t{r.query_length - ctg_end - 1}\t{r.query_length - ctg_start}{out_name}")
                else:
                    print(f"{ctg_name}\t{ctg_start}\t{ctg_end}{out_name}")
