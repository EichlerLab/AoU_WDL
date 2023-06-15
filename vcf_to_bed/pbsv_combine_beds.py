#!/bin/env python

import argparse
import numpy as np
import pandas as pd
import svpoplib

parser = argparse.ArgumentParser()
parser.add_argument("--a_file", "-a", type=str, required=True, help="Bed table with columns")
parser.add_argument("--b_file", "-b", type=str, required=True, help="Bed table with columns")
parser.add_argument("--output", "-o", type=str, required=True, help="Output file name for processed bed file")

args = parser.parse_args()

bed_files = [args.a_file, args.b_file]

bed_df = pd.concat([pd.read_csv(bed_file_name, sep="\t") for bed_file_name in bed_files], axis=0)

bed_df.sort_values(['#CHROM', 'POS'], inplace=True)

bed_df.to_csv(args.output, sep='\t', na_rep='NA', index=False)