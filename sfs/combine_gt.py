import pandas as pd
import sys

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python script.py <file1.tsv> <file2.tsv> <outpath>")
        sys.exit(1)

    h1_path = sys.argv[1]
    h2_path = sys.argv[2]
    outfile_path = sys.argv[3]


    try:
        h1_df = pd.read_csv(h1_path, sep='\t', dtype=str)
        h2_df = pd.read_csv(h2_path, sep='\t', dtype=str)

        df = h1_df.merge(h2_df, on='ID')

        assert len(df) == len(h1_df)

        df.to_csv(outfile_path, sep='\t', index=False)

    except FileNotFoundError:
        print("File not found. Please provide valid file paths.")
    except pd.errors.EmptyDataError:
        print("One or both of the provided DataFrames are empty.")