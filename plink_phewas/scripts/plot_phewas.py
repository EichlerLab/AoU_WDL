#!/usr/bin/env python3
"""
Generate a PheWAS Manhattan plot from a PheTK-ready results TSV.
Labels all points passing Bonferroni correction (0.05 / n_tests).
"""
import argparse

import pandas as pd
from phetk.plot import Plot


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--phewas-results", required=True,
                   help="PheTK-ready TSV produced by format_results.py")
    p.add_argument("--output", required=True,
                   help="Output image path (e.g. my_sv_manhattan.png)")
    p.add_argument("--title", default="Phenome-Wide Association Study",
                   help="Plot title")
    p.add_argument("--label-only-bonferroni", action="store_true",
                   help="Label only Bonferroni-significant points; takes precedence over --manhattan-label-count")
    p.add_argument("--manhattan-label-count", type=int, default=15,
                   help="Number of top associations to label when --label-only-bonferroni is not set")
    p.add_argument("--min-case-carrier-ct", type=int, default=5,
                   help="Minimum CASE_HET_A1_CT + CASE_HOM_A1_CT to include in plot")
    return p.parse_args()


def main():
    args = parse_args()

    df = pd.read_csv(args.phewas_results, sep="\t")
    df_plot = df.dropna(subset=["phecode_category"])
    df_plot = df_plot[
        (df_plot["CASE_HET_A1_CT"] + df_plot["CASE_HOM_A1_CT"]) >= args.min_case_carrier_ct
    ]
    n_tests = len(df_plot)
    df_plot.to_csv("phewas_results_for_plot.tsv", sep="\t", index=False)

    if args.label_only_bonferroni:
        bonferroni_threshold = 0.05 / n_tests if n_tests > 0 else 0.05
        label_count = int((df["p_value"] < bonferroni_threshold).sum())
        print(f"n_tests={n_tests}, Bonferroni threshold={bonferroni_threshold:.2e}, "
              f"significant={label_count}", flush=True)
    else:
        label_count = args.manhattan_label_count
        print(f"n_tests={n_tests}, labelling top {label_count}", flush=True)

    phewas_plot = Plot("phewas_results_for_plot.tsv")
    phewas_plot.manhattan(
        title=args.title,
        output_file_path=args.output,
        label_count=label_count,
    )
    print(f"Manhattan plot written to {args.output}", flush=True)


if __name__ == "__main__":
    main()
