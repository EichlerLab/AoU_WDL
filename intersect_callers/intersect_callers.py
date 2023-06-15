#!/bin/env python

import argparse
import numpy as np
import pandas as pd
import svpoplib


def sv_only(infile):
    df_in = pd.read_csv(infile, sep='\t')
    df_in['SVLEN'] = df_in['ID'].str.split('-', expand=True)[3].astype(float).astype(int)
    df_in = df_in.loc[df_in['SVLEN'] >= 50]
    return df_in


def temp_bed(infile, outfile):
    df_out = sv_only(infile)
    df_out[['#CHROM', 'POS', 'END', 'SVTYPE', 'ID', 'SVLEN']].to_csv(outfile, sep='\t', index=False)


parser = argparse.ArgumentParser()

parser.add_argument("--a_file", "-a", type=str, required=True, help="Bed table with columns ['#CHROM', 'POS', 'END', 'SVTYPE','ID'] where ID == chr-pos-svtype-svlen")
parser.add_argument("--b_file", "-b", type=str, required=True, help="Bed table with columns ['#CHROM', 'POS', 'END', 'SVTYPE','ID'] where ID == chr-pos-svtype-svlen")
parser.add_argument("--output", "-o", type=str, required=True, help="Output file name for processed intersect file")
parser.add_argument("--threads", "-t", type=int, required=True, help="Output file name for processed intersect file")

args = parser.parse_args()


# args = pd.Series(['a.bed', 'b.bed', 'c.out', 4], index=['a_file', 'b_file', 'output', 'threads'])

temp_bed(args.a_file, "temp_a.bed")
temp_bed(args.b_file, "temp_b.bed")

df = svpoplib.svmerge.merge_variants(
    bed_list=["temp_a.bed", "temp_b.bed"],
    sample_names=['A', 'B'],
    strategy="nr::szro",
    threads=args.threads,
)

#df.loc[df['MERGE_SAMPLES'] == 'A', 'MERGE_VARIANTS'] = df.loc[df['MERGE_SRC'] == 'A', 'MERGE_VARIANTS'] + ','


#df = df.loc[df['MERGE_SAMPLES'].str.contains('A')][['MERGE_SAMPLES', 'MERGE_SRC_ID']]

#df_a = sv_only(args.a_file)

#df = pd.merge(df, df_a, left_on='MERGE_SRC_ID', right_on='ID')

#assert len(df) == len(df_a)

df.to_csv(args.output, sep='\t', index=False)


