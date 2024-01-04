# Get bam segment in reference region
import sys
import pysam

in_bam = sys.argv[1]
gene_region_bed = sys.argv[2]
ref_pos_end_buffer = int(sys.argv[3])  # find closest matching contig start/ end coord to reference coord, within this buffer

with open(gene_region_bed) as f:
    chrom, ref_pos, ref_end = f.readline().strip().split('\t')
ref_pos = int(ref_pos)
ref_end = int(ref_end)
    
s = pysam.AlignmentFile(in_bam, "rb")  # 1000151-asm_h1.minimap2.bam
for r in s.fetch(reference=chrom, start=ref_pos, end=ref_end):
    is_rev_strand = r.flag & 0x10 != 0
    ctg_start = None
    ctg_end = None
    start_min_dist = ref_pos_end_buffer
    end_min_dist = ref_pos_end_buffer

    for ctg_coord, ref_coord in r.get_aligned_pairs():
        if ref_coord is None:
            continue
        # assume sorted
        if abs(ref_coord - ref_pos) <= start_min_dist:
            ctg_start = ctg_coord  # closest position to the start position in the bed, within buffer
            start_min_dist = abs(ref_coord - ref_pos)
        elif abs(ref_coord - ref_end) <= end_min_dist:
            ctg_end = ctg_coord
            end_min_dist = abs(ref_coord - ref_end)

    if ctg_start is not None and ctg_end is not None:
        if is_rev_strand:
            print(f"{r.query_name}/rc\t{r.query_length - ctg_end - 1}\t{r.query_length - ctg_start}\n")
        else:
            print(f"{r.query_name}\t{ctg_start}\t{ctg_end}\n")
