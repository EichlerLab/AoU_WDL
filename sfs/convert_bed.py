import pandas as pd
import sys

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python script.py <coverage.bed> <hap> <outfile>")
        sys.exit(1)

    cover_bed_path = sys.argv[1]
    hap = sys.argv[2]
    outfile_path = sys.argv[3]

    try:
        df = pd.read_csv(cover_bed_path, sep='\t', dtype=str, header=None, names=['#CHROM', 'POS', 'END', 'ID', 'SAMPLE', 'ITEMS', 'BASES_COV', 'BASES_WND', 'PERC_COV'])
        sample = df.iloc[0]['SAMPLE']

        df[f'{sample}_{hap}'] = df.apply(lambda row: "0" if row['PERC_COV'] == 1 else ".", axis=1)
        df[['ID', f'{sample}_{hap}']].to_csv(outfile_path, sep='\t', index=False)


    except FileNotFoundError:
        print("File not found. Please provide valid file paths.")
    except pd.errors.EmptyDataError:
        print("One or both of the provided DataFrames are empty.")


