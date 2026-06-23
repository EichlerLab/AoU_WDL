#!/usr/bin/env python3
"""
Merge raw plink2 GLM output with phecode definitions to produce a
PheTK-ready results TSV.

Input:  aggregated plink results TSV (PHENOTYPE prepended by awk)
        phecode_definitions1.2.csv
Output: all original plink columns plus phecode, phecode_string, phecode_category,
        neg_log_p_value, beta, and converged
"""
import argparse

import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--plink-results", required=True,
                   help="Aggregated plink2 Firth GLM results (with PHENOTYPE column)")
    p.add_argument("--phecode-definitions", required=True,
                   help="phecode_definitions1.2.csv")
    p.add_argument("--output", required=True,
                   help="Path for formatted PheTK-ready TSV")
    return p.parse_args()


def main():
    args = parse_args()

    df = pd.read_csv(args.plink_results, sep="\t", na_values=["NA", "nan"],
                     keep_default_na=False)

    # Derived columns
    df["phecode"]         = pd.to_numeric(df["PHENOTYPE"].str.replace("phecode_", "", regex=False), errors="coerce")
    df["converged"]       = df["ERRCODE"].isin([".", "NM_ERR_NONE"])
    p                     = pd.to_numeric(df["P"], errors="coerce")
    df["p_value"]         = p
    df["neg_log_p_value"] = (-np.log10(p)).replace([np.inf, -np.inf], np.nan)
    df["beta"]            = np.log(pd.to_numeric(df["OR"], errors="coerce")).replace([np.inf, -np.inf], np.nan)

    phecode_map = pd.read_csv(args.phecode_definitions)
    phecode_map = phecode_map.rename(columns={
        "category":  "phecode_category",
        "phenotype": "phecode_string",
    })[["phecode", "phecode_string", "phecode_category"]]

    df_final = df.merge(phecode_map, on="phecode", how="inner")

    df_final.to_csv(args.output, sep="\t", index=False)
    print(f"Formatted {len(df_final)} phenotypes → {args.output}", flush=True)


if __name__ == "__main__":
    main()
