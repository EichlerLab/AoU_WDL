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

N_PROP_CUTOFF = 0.1

df = pd.read_csv(args.bed_file, sep='\t')

df = svpoplib.varbed.bcftools_query_to_tsv(df, args.sample)

df = df.loc[df['SVTYPE'].apply(lambda val: val in {'INS','DEL'})].copy()

del(df['CIPOS'],df['MATEID'],df['MATEDIST'])

df_coord = df.apply(svpoplib.variant.vcf_fields_to_seq, axis=1)

df['POS'] = df_coord['POS']
df['END'] = df_coord['END']
df['SVLEN'] = df_coord['SVLEN']
df['SVTYPE'] = df_coord['SVTYPE']
df['SEQ'] = df_coord['SEQ']
df['FAIL_REASON'] = np.nan

df['PBSV_ID'] = df['ID']
df['ID'] = svpoplib.variant.get_variant_id(df)

df = svpoplib.variant.order_variant_columns(df)

fail_set_nseq = set(df.loc[df.apply(
    lambda row: (
            np.sum([val == 'N' for val in row['SEQ'].upper()]) / row['SVLEN']
        ) >= N_PROP_CUTOFF if not pd.isnull(row['SEQ']) else False,
    axis=1
)].index)

df.loc[fail_set_nseq,'FAIL_REASON'] = 'SEQ N proportion'
fail_set_ins_n = set(df.loc[(df['SVTYPE'] == 'INS') & (df['REF'].apply(lambda val: val.upper() == 'N'))].index)
df.loc[fail_set_ins_n, 'FAIL_REASON'] = 'INS in N'

fail_set_filter = set(df.loc[df['FILTER'] != 'PASS'].index)

df.loc[fail_set_filter, 'FAIL_REASON'] = 'FILTER != PASS'

fail_set = fail_set_nseq | fail_set_ins_n | fail_set_filter


df = df.loc[[val not in fail_set for val in df.index]].copy()
del(df['FAIL_REASON'])
id_count = df.groupby('ID').apply(lambda subdf: subdf.shape[0])

df['PBSV_ID_COUNT'] = df['ID'].apply(lambda id: id_count[id])

df.drop_duplicates(subset='ID', keep='first', inplace=True)

del(df['REF'], df['ALT'], df['FILTER'])

no_rename = ['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'SEQ', 'PBSV_ID', 'PBSV_ID_COUNT']

df.columns = [('PBSV_' + col if col not in no_rename else col) for col in df.columns]

df.sort_values(['#CHROM', 'POS'], inplace=True)

del(df['SEQ'])

df_sub = df.loc[df['SVLEN'] >= 50]

df_sub.to_csv(args.output, sep='\t', index=False)