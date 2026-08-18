#!/usr/bin/env python3
"""
gene_to_chrom_vcf.py

Rewrite a gene/region-based VCF's CHROM/POS into true chromosome-based
coordinates, given a pre-computed chromosome name, region start offset,
and contig length.

Background
----------
Pipelines that call variants against a small extracted FASTA region, e.g.
    >chr15:28419321-28443019
produce a VCF where:
    - CHROM = "chr15:28419321-28443019"   (the region, not a real chromosome)
    - POS   = 1-based offset WITHIN that region
    - a header line: ##contig=<ID=chr15:28419321-28443019,length=23699>

Since every sample in a run is called against the exact same gene_seq
region, the chromosome name / start offset / true chromosome length only
need to be worked out ONCE per workflow run (see the parse_gene_seq task
in paralogous_variants.wdl), not once per sample. This script therefore
takes those three values directly as arguments instead of re-parsing a
FASTA header or a chrom.sizes table itself -- important since it may be
invoked across thousands of samples/haplotypes.

This script:
    1. Rewrites every record's CHROM to --chrom.
    2. Rewrites POS to the true genomic coordinate:
           new_POS = --start + old_POS - 1
    3. Rewrites any ##contig header line to
           ##contig=<ID=--chrom,length=--length>
    4. Leaves everything else (INFO, FORMAT, genotypes, etc.) untouched.

Usage
-----
    python gene_to_chrom_vcf.py input.vcf --chrom chr15 --start 28419321 --length 99753195 -o output.vcf
"""
import sys
import argparse


def convert_vcf(in_path, out_handle, chrom, start, contig_length):
    with open(in_path) as fin:
        for raw_line in fin:
            line = raw_line.rstrip('\n')
            if not line:
                continue

            if line.startswith('##contig='):
                out_handle.write(f'##contig=<ID={chrom},length={contig_length}>\n')
                continue

            if line.startswith('#'):
                out_handle.write(line + '\n')
                continue

            fields = line.split('\t')
            old_pos = int(fields[1])
            fields[0] = chrom
            fields[1] = str(start + old_pos - 1)
            out_handle.write('\t'.join(fields) + '\n')


def main():
    parser = argparse.ArgumentParser(
        description="Rewrite a gene/region-based VCF's CHROM/POS into "
                    "chromosome-based coordinates using pre-computed "
                    "chrom/start/length values.")
    parser.add_argument('vcf', help='Input VCF file with region-based CHROM/POS')
    parser.add_argument('--chrom', required=True,
                         help='Real chromosome name, e.g. chr15')
    parser.add_argument('--start', required=True, type=int,
                         help="Region's start coordinate on the chromosome "
                              "(1-based), e.g. 28419321")
    parser.add_argument('--length', required=True, type=int,
                         help='Contig length to declare in the ##contig header')
    parser.add_argument('-o', '--output', default=None,
                         help='Output VCF path (default: stdout)')
    args = parser.parse_args()

    if args.output:
        with open(args.output, 'w') as fout:
            convert_vcf(args.vcf, fout, args.chrom, args.start, args.length)
    else:
        convert_vcf(args.vcf, sys.stdout, args.chrom, args.start, args.length)


if __name__ == '__main__':
    main()
