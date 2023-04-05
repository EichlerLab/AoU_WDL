#!/bin/env python

import argparse
import numpy as np
import pandas as pd
import svpoplib

parser = argparse.ArgumentParser()
parser.add_argument("--bed_file", "-b", type=str, required=True, help="Bed table with columns")
parser.add_argument("--sample", "-s", type=str, required=True, help="sample_name")
parser.add_argument("--output", "-o", type=str, required=True, help="Output file name for processed bed file")

args = parser.parse_args()

df = pd.read_csv(args.bed_file, sep='\t')

df = svpoplib.varbed.bcftools_query_to_tsv(df, args.sample)

df = df.loc[(df['SVTYPE'] == 'INV') | (df['SVTYPE'] == 'DUP')]
df['END'] = df['END'].astype(int)
df['SVLEN'] = df['END'] - df['POS']

del(
    df['CIPOS'], df['IMPRECISE'], df['MATEDIST'], df['MATEID'],
    df['REF'], df['ALT']
)

df['PBSV_ID'] = df['ID']
df['ID'] = svpoplib.variant.get_variant_id(df)
df = svpoplib.variant.order_variant_columns(df)

df['FAIL_REASON'] = np.nan
fail_set = set(df.loc[df['FILTER'] != 'PASS'].index)

df.loc[fail_set, 'FAIL_REASON'] = 'FILTER != PASS'

df = df.loc[[val not in fail_set for val in df.index]].copy()

del(df['FAIL_REASON'])

id_count = df.groupby('ID').apply(lambda subdf: subdf.shape[0])
df['PBSV_ID_COUNT'] = df['ID'].apply(lambda id: id_count[id]).copy()

df.drop_duplicates(subset='ID', keep='first', inplace=True)

df.sort_values(['#CHROM', 'POS'], inplace=True)

no_rename = ['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'SEQ', 'PBSV_ID', 'PBSV_ID_COUNT']

df.columns = [('PBSV_' + col if col not in no_rename else col) for col in df.columns]

del(df['PBSV_FILTER'])
df.to_csv(args.output, sep='\t', na_rep='NA', index=False)


