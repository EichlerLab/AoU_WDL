import pandas as pd
import sys

if __name__ == '__main__':
    output_matix_path = sys.argv[1]
    hap = sys.argv[2:]

    try:
        for i, df_in in enumerate(hap):
            if i == 0:
                df = pd.read_csv(df_in, sep='\t')
            else:
                df = df.merge(pd.read_csv(df_in, sep='\t'), on='ID', how='left')

        df.to_csv(output_matix_path, sep='\t', index=False)


    except FileNotFoundError:
        print("File not found. Please provide valid file paths.")
    except pd.errors.EmptyDataError:
        print("One or both of the provided DataFrames are empty.")
