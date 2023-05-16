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

del(df['REF'], df['ALT'])
df = df.loc[df['SVLEN'] != '.'].copy()
df['SVLEN'] = df['SVLEN'].astype(int)
df['SVLEN'] = np.abs(df['SVLEN'])

df['SVTYPE_ORG'] = df['SVTYPE']

df['SVTYPE'] = df['SVTYPE'].apply(lambda val: val.replace('/', ''))

df['END'] = df.apply(lambda row: row['POS'] + (row['SVLEN'] if row['SVTYPE'] != 'INS' else 1), axis=1)

df['ID'] = svpoplib.variant.get_variant_id(df)

df = svpoplib.variant.order_variant_columns(df)
df.sort_values(
    ['#CHROM', 'POS', 'SVTYPE', 'SVLEN', 'GT', 'DV', 'DR'],
    ascending=(True, True, True, True, False, False, True),
    inplace=True
)

df = df.loc[df['FILTER'] == 'PASS']
 
df.drop_duplicates(['ID'], inplace=True)
 
for rm_field in (
    'PRECISE', 'IMPRECISE', 'ID.1', 
    'Kurtosis_quant_start', 'Kurtosis_quant_stop', 'STD_quant_start', 'STD_quant_stop',
    'MAPQ', 'RE', 'REF_strand', 'STRANDS', 'SUPTYPE', 'AF', 'SVMETHOD', 'ZMW', 'SAMPLE',
    'FILTER', 'SVTYPE_ORG', 'UNRESOLVED'
):
    try:
        del(df[rm_field])
    except:
        continue
prepend_caller_set = {'GT','DV','DR','CHR2'}

df.columns = [('SNIF_' + col if col in prepend_caller_set else col) for col in df.columns]

df_sub = df.loc[df['SVLEN'] >= 50]
df_sub.to_csv(args.output, sep='\t', index=False)
