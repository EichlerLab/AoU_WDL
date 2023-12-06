# Get bam segment in reference region
import sys
import pysam

in_bam = sys.argv[1]
gene_region_bed = sys.argv[2]

with open(gene_region_bed) as f:
    chrom, ref_pos, ref_end = f.readline().strip().split('\t')
ref_pos = int(ref_pos)
ref_end = int(ref_end)
    
s = pysam.AlignmentFile(in_bam, "rb")  # 1000151-asm_h1.minimap2.bam
for r in s.fetch(reference=chrom, start=ref_pos, end=ref_end):
    is_rev_strand = r.flag & 0x10 != 0
    ctg_start = None
    ctg_end = None
    for ctg_coord, ref_coord in r.get_aligned_pairs():
        if ref_coord is None:
            continue
        # assume sorted
        if ref_coord == ref_pos:
            ctg_start = ctg_coord
        elif ref_coord == ref_end:
            ctg_end = ctg_coord
    assert ctg_start is not None and ctg_end is not None
    if is_rev_strand:
        print(f"{r.query_name}/rc\t{r.query_length - ctg_end}\t{r.query_length - ctg_start}\n")
    else:
        print(f"{r.query_name}\t{ctg_start}\t{ctg_end}\n")
