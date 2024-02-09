#!/bin/env python 

import pandas as pd
import argparse


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter) 

parser.add_argument('--input', '-i', type=str, required=True, help='Input svpop')
parser.add_argument('--output', '-o', type=str, required=True, help='Output body')
parser.add_argument('--samples', '-s', type=str, required=True, help='Sample list')

args = parser.parse_args()

with open(args.samples, "w") as infile:
    sample_list = [line.rstrip() for line in infile]

df = pd.read_csv(args.input, sep='\t', dtype=str)

df["ALT"] = df["VCF_ALT"]
df["REF"] = df["VCF_REF"]


for sample in sample_list:
    df[sample] = df.apply(lambda row: row["GT"] if row["MERGE_SRC"] == sample else ".|.", axis=1)
    df[sample] = df.apply(lambda row: "0/1" if row[sample] == ".|." and sample in row["MERGE_SAMPLES"].split(",") else row[sample], axis=1)


df['INFO'] = "SVTYPE="+df['SVTYPE']+";CALLERSET_SRC="+df['CALLERSET_SRC']+";SVLEN="+df['SVLEN']+";CALLERSET_VARIANTS="+df['CALLERSET_VARIANTS']+";SUPPORT_COUNT="+df['SUPPORT_COUNT']
df['FORMAT'] = "GT"
df["FILTER"] = "PASS"
columns = ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]+sample_list

df[columns].to_csv(args.output, sep='\t', index=False)