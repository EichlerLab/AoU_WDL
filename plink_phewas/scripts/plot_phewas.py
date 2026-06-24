#!/usr/bin/env python3
"""
Generate a PheWAS Manhattan plot from a PheTK-ready results TSV.
Labels all points passing Bonferroni correction (0.05 / n_tests).
"""
import argparse

import numpy as np
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

    # Bonferroni uses all phenotypes with a valid category (before carrier filter)
    df_all = df.dropna(subset=["phecode_category"])
    n_tests = len(df_all)
    bonferroni_threshold = 0.05 / n_tests if n_tests > 0 else 0.05

    # Points shown in plot also require sufficient alt-allele cases
    df_plot = df_all[
        (df_all["CASE_HET_A1_CT"] + df_all["CASE_HOM_A1_CT"]) >= args.min_case_carrier_ct
    ].copy()

    # Annotate labels with OR and case/control carrier counts
    def fmt_label(row):
        try:
            or_val = f"OR={np.exp(row['beta']):.2f}" if pd.notna(row.get("beta")) else "OR=NA"
        except Exception:
            or_val = "OR=NA"
        def ct(col):
            v = row.get(col)
            return int(v) if pd.notna(v) else 0
        cases = ct("CASE_HET_A1_CT") + ct("CASE_HOM_A1_CT") + ct("CASE_NON_A1_CT")
        controls = ct("CTRL_HET_A1_CT") + ct("CTRL_HOM_A1_CT") + ct("CTRL_NON_A1_CT")
        return f"{row['phecode_string']} ({or_val}, cases={cases}, controls={controls})"

    df_plot["phecode_string"] = df_plot.apply(fmt_label, axis=1)

    df_plot.to_csv("phewas_results_for_plot.tsv", sep="\t", index=False)

    neg_log_bonferroni = -np.log10(bonferroni_threshold)
    # Ensure y-axis extends at least to the Bonferroni line (+ buffer)
    max_neg_log_p = df_plot["neg_log_p_value"].dropna().max()
    y_limit = int(max(max_neg_log_p, neg_log_bonferroni) * 1.1) + 1

    if args.label_only_bonferroni:
        sig = df_plot[df_plot["p_value"] < bonferroni_threshold]
        label_values = sig["phecode"].dropna().astype(str).tolist()
        label_count = len(label_values)
        print(f"n_tests={n_tests}, Bonferroni threshold={bonferroni_threshold:.2e}, "
              f"significant in plot={label_count}", flush=True)
    else:
        label_values = "p_value"
        label_count = args.manhattan_label_count
        print(f"n_tests={n_tests}, Bonferroni threshold={bonferroni_threshold:.2e}, "
              f"labelling top {label_count}", flush=True)

    phewas_plot = Plot("phewas_results_for_plot.tsv", bonferroni=neg_log_bonferroni)
    phewas_plot.manhattan(
        title=args.title,
        label_values=label_values,
        label_count=label_count,
        y_limit=y_limit,
        output_file_path=args.output,
        save_plot=True,
    )
    print(f"Manhattan plot written to {args.output}", flush=True)


if __name__ == "__main__":
    main()