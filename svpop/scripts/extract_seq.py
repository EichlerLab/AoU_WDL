#!/bin/env python

import pandas as pd
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter) 

parser.add_argument('--sample', '-s', type=str, required=True, help='Sample name')
parser.add_argument('--input', '-i', type=str, required=True, help='Sniffles input vcf')
parser.add_argument('--ref', '-r', type=str, required=True, help='Input reference fasta')
parser.add_argument('--output', '-o', type=str, required=True, help='Output file')

args = parser.parse_args()

sample=args.sample

df = pd.read_csv(f"{args.input}", sep='\t', header=None, comment="#", dtype=str).dropna()

columns=range(10)

df["ALT_HOLD"] = df.apply(lambda row: row[3] if row[4] == '<DEL>' else row[4], axis=1)

df["END"] = df.apply(lambda row: int([x for x in row[7].split(";") if x.startswith("END=")][0].split("=")[1]) if "END" in row[7] and row[4] == '<DEL>' else "NA", axis=1)


# Load your FASTA file into a dictionary
fasta_file = args.ref
sequences_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

# Read your DataFrame (assuming it has a 'Sequence_ID' column)

# Function to extract sequence based on range using SeqIO dictionary
def extract_sequence_from_dict(sequences_dict, sequence_id, start_pos, end_pos):
    sequence = sequences_dict.get(sequence_id)
    if sequence:
        extracted_sequence = str(sequence.seq)[start_pos - 1:end_pos]  # Adjust indices (Python is 0-indexed)
        return extracted_sequence
    else:
        return None

# Apply the function to each row in the DataFrame
df['Extracted_Sequence'] = df.apply(lambda row: extract_sequence_from_dict(sequences_dict, row[0], int(row[1]), row['END']) if row["END"] != "NA" else row[3], axis=1)

df[3] = df['Extracted_Sequence']
df[4] = df["ALT_HOLD"]

df[columns].to_csv(args.output, sep='\t', index=False, header=False)
