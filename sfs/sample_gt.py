import pandas as pd
import sys

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python script.py <file1.tsv> <file2.tsv> <outpath>")
        sys.exit(1)

    callable_path = sys.argv[1]
    callerset_path = sys.argv[2]
    outfile_path = sys.argv[3]

    try:
        df = pd.concat(pd.read_csv(x, sep='\t', dtype=str) for x in [callable_path, callerset_path])
        df.to_csv(outfile_path, sep='\t', index=False)

    except FileNotFoundError:
        print("File not found. Please provide valid file paths.")
    except pd.errors.EmptyDataError:
        print("One or both of the provided DataFrames are empty.")