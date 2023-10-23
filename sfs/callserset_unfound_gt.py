import pandas as pd
import sys

def merge_and_separate_dataframes(callerset, all_calls):
    # Merge DataFrames based on 'ID' column
    merged_df = pd.merge(callerset, all_calls, on='ID', how='inner')  # Shared IDs
    unique_all_calls = all_calls[~all_calls['ID'].isin(merged_df['ID'])]  # Unique IDs in all_calls

    return merged_df, unique_all_calls

if __name__ == '__main__':
    if len(sys.argv) != 6:
        print("Usage: python script.py <file1.tsv> <file2.tsv> <all_calls_unique.bed> <callerset_found.gt> <sample_name>")
        sys.exit(1)

    callerset_path = sys.argv[1]
    all_calls_path = sys.argv[2]
    outpath_unique = sys.argv[3]
    outpath_gt = sys.argv[4]
    sample = sys.argv[5]

    try:
        callerset = pd.read_csv(callerset_path, sep='\t', dtype=str)
        all_calls = pd.read_csv(all_calls_path, sep='\t', dtype=str)

        shared_df, unique_all_calls = merge_and_separate_dataframes(callerset, all_calls)

        shared_df[f'{sample}_h1'] = shared_df['GT'].str.split("|", expand=True)[0]
        shared_df[f'{sample}_h2'] = shared_df['GT'].str.split("|", expand=True)[1]

        shared_df[['ID', f'{sample}_h1', f'{sample}_h2']].to_csv(f'{outpath_gt}', index=False, sep='\t')
        unique_all_calls['SAMPLE'] = sample

        unique_all_calls[['#CHROM', 'POS', 'END', 'ID', 'SAMPLE']].to_csv(f'{outpath_unique}', index=False, sep='\t')

        print("Merged DataFrame (shared IDs):")
        print(shared_df)

        print("\nUnique IDs in DataFrame 2:")
        print(unique_all_calls)

    except FileNotFoundError:
        print("File not found. Please provide valid file paths.")
    except pd.errors.EmptyDataError:
        print("One or both of the provided DataFrames are empty.")