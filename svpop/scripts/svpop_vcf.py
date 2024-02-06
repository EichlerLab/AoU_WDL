#!/bin/env python

import pandas as pd

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter) 

parser.add_argument('--sample', '-s', type=str, required=True, help='Sample name')
parser.add_argument('--input', '-i', type=str, required=True, help='Sniffles input vcf')
parser.add_argument('--pbsv', '-p', type=str, required=True, help='PBSV VCF for sequence')
parser.add_argument('--output', '-o', type=str, required=True, help='Output file')

args = parser.parse_args()

sample= args.sample

df = pd.read_csv(args.input, sep='\t', dtype=str)

columns = ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",f"{sample}"]

pbsv_df = pd.read_csv(args.pbsv, sep='\t', comment="#", header=None, names=columns, index_col=["#CHROM", "ID"])

df['FILTER'] = 'PASS'

df['SUPPORT'] = df['CALLERSET_LIST'].apply(lambda val : str(len(val.split(","))))

df['INFO'] = "SVTYPE="+df['SVTYPE']+";CALLERSET_SRC="+df['CALLERSET_LIST']+";SVLEN="+df['SVLEN']+";CALLERSET_VARIANTS="+df['CALLERSET_VARIANTS']+";SUPPORT_COUNT="+df['SUPPORT']

df['VCF_ALT'] = df.apply(lambda row: pbsv_df.loc[(row["#CHROM"], row['PBSV_ID'])]["ALT"].values[0] if row["CALLERSET_SRC"] == "PBSV" else row['VCF_ALT'], axis=1)
df['VCF_REF'] = df.apply(lambda row: pbsv_df.loc[(row["#CHROM"], row['PBSV_ID'])]["REF"].values[0] if row["CALLERSET_SRC"] == "PBSV" else row['VCF_REF'], axis=1)
df['GT'] = df.apply(lambda row: row["PBSV_GT"] if row["CALLERSET_SRC"] == "PBSV" else row['GT'], axis=1)

df['REF'] = df['VCF_REF']
df['ALT'] = df['VCF_ALT']
df['FORMAT'] = "GT"
df[sample] = df['GT']


df[columns].to_csv(args.output, sep='\t', index=False)
