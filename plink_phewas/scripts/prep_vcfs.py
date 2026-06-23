#!/usr/bin/env python3
"""
Create per-SV VCF files from a long-format genotype TSV.

Input TSV columns: SVID, Sample, GT
Output: one VCF per passing SV written to --output-dir/
        a list of output VCF filenames written to --sv-list-out
"""
import argparse
import os
import sys

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--geno-tsv", required=True,
                   help="Long-format genotype file: SVID\\tSample\\tGT")
    p.add_argument("--pheno-tsv", required=True,
                   help="Phenotype matrix: IID\\t[pheno cols...]")
    p.add_argument("--relatedness-flagged-samples", required=True,
                   help="One sample per line (first line is header)")
    p.add_argument("--min-genotyped-samples", type=int, default=200_000,
                   help="Minimum number of unrelated genotyped samples required")
    p.add_argument("--sv-subset", default="",
                   help="Comma-separated SV IDs to process; omit to process all passing SVs")
    p.add_argument("--output-dir", default="plink_gt",
                   help="Directory for output VCFs")
    p.add_argument("--sv-list-out", default="SV_VCF_list.txt",
                   help="Path to write the list of created VCF filenames")
    return p.parse_args()


def main():
    args = parse_args()

    print("Loading genotype TSV…", flush=True)
    geno = pd.read_csv(args.geno_tsv, dtype=str, sep="\t")

    print("Loading sample phenotype index…", flush=True)
    sample_phenos = pd.read_csv(
        args.pheno_tsv, dtype={"IID": str}, sep="\t", header=0
    )
    pheno_sample_set = set(sample_phenos["IID"])

    # Skip header line
    with open(args.relatedness_flagged_samples) as f:
        related_samples = set(line.strip() for line in f if line.strip())
    related_samples.discard(next(iter(related_samples)))  # pop header if present
    # Safer: re-read and skip first line
    with open(args.relatedness_flagged_samples) as f:
        lines = f.read().splitlines()
    related_samples = set(lines[1:])  # index 0 is header

    os.makedirs(args.output_dir, exist_ok=True)

    svs = geno["SVID"].unique().tolist()
    if args.sv_subset:
        keep = set(args.sv_subset.split(","))
        svs = [s for s in svs if s in keep]
        print(f"#SVs after subset filter: {len(svs)}", flush=True)
    else:
        print(f"#Total SVs: {len(svs)}", flush=True)

    miss_pheno = failed_geno = 0

    with open(args.sv_list_out, "w") as sv_list_fh:
        for svid in svs:
            this_geno = geno.loc[
                (geno["SVID"] == svid) & (~geno["Sample"].isin(related_samples))
            ]
            geno_samples = this_geno["Sample"].unique().tolist()

            if not set(geno_samples).issubset(pheno_sample_set):
                print(f"WARN {svid}: missing phenotype for some samples", flush=True)
                miss_pheno += 1
                continue

            if len(geno_samples) < args.min_genotyped_samples:
                print(
                    f"WARN {svid}: only {len(geno_samples)} genotyped samples "
                    f"(min {args.min_genotyped_samples})",
                    flush=True,
                )
                failed_geno += 1
                continue

            if svid.count("-") == 3:
                chrom, pos, svtype, _ = svid.split("-")
            else:
                # Non-standard ID: use chr1/pos=1 so plink2 treats it as a normal autosomal variant;
                # the ID column carries the svid and is what appears in plink2 output
                chrom, pos, svtype = "chrUn", "1", "SV"
            ref = "N"
            alt = f"<{svtype}>"

            vcf_path = os.path.join(args.output_dir, f"{svid}_GT.vcf")
            with open(vcf_path, "w") as vcf_fh:
                vcf_fh.write("##fileformat=VCFv4.2\n")
                vcf_fh.write(f"##contig=<ID={chrom}>\n")
                vcf_fh.write(
                    "##FORMAT=<ID=GT,Number=1,Type=String,"
                    "Description=\"Genotype based on diploid assembly\">\n"
                )
                vcf_fh.write(
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                    + "\t".join(geno_samples)
                    + "\n"
                )
                gt_values = this_geno["GT"].tolist()
                vcf_fh.write(
                    f"{chrom}\t{pos}\t{svid}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t"
                    + "\t".join(gt_values)
                    + "\n"
                )

            print(f"{svid}_GT.vcf", file=sv_list_fh)
            print(f"Wrote {svid} ({len(geno_samples)} samples)", flush=True)

    n_complete = len(svs) - failed_geno - miss_pheno
    print(f"\n#SV failed genotype threshold:  {failed_geno}", flush=True)
    print(f"#SV missing phenotype:          {miss_pheno}", flush=True)
    print(f"Complete SVs:                   {n_complete}", flush=True)


if __name__ == "__main__":
    main()
