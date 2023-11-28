import pandas as pd
import argparse


parser = argparse.ArgumentParser()

parser.add_argument("--support", "-s", type=int, required=True, help="Support count")
parser.add_argument(
    "--output", "-o", type=str, required=False, help="Output file to write to"
)
parser.add_argument("--matrix", "-m", type=str, required=True, help="GT Matrix file")
parser.add_argument("--callerset", "-c", type=str, required=True, help="All callerset")

args = parser.parse_args()


df = pd.read_csv(args.matrix, sep="\t", dtype=str).fillna(".")
print("Done reading Matrix")

callerset_df = pd.read_csv(
    args.callerset,
    sep="\t",
    usecols=["ID", "SUPPORT"],
    dtype={"ID": str, "SUPPORT": int},
)
print("Done reading callerset")

callerset_df = callerset_df.loc[callerset_df["SUPPORT"] >= args.support].copy()

df = pd.merge(callerset_df, df).set_index("ID")

print("Done merging df")

sample_list = [x for x in df.columns if x not in ["VAR", "OBS", "FREQ", "SUPPORT"]][2:]

out_df = pd.DataFrame()

df["VAR"] = 0
df["OBS"] = 0

for sample in sample_list:
    df["VAR"] = df.apply(
        lambda row: row["VAR"] + 1 if row[sample] == "1" else row["VAR"], axis=1
    )
    df["OBS"] = df.apply(
        lambda row: row["OBS"] + 1 if row[sample] != "." else row["OBS"], axis=1
    )
    df["FREQ"] = df["VAR"] / df["OBS"]
    df_var = df.loc[df["VAR"] > 0].copy()
    out_df = pd.concat(
        [
            out_df,
            pd.DataFrame.from_dict(
                {
                    "SINGLETON": [
                        len(df_var.loc[(df_var["VAR"] == 1) & (df_var["FREQ"] != 1)])
                    ],
                    "POLY": [
                        len(df_var.loc[(df_var["VAR"] > 1) & (df_var["FREQ"] < 0.5)])
                    ],
                    "MAJOR": [
                        len(
                            df_var.loc[
                                (df_var["VAR"] > 1)
                                & (df_var["FREQ"] >= 0.5)
                                & (df_var["FREQ"] < 1)
                            ]
                        )
                    ],
                    "FIXED": [len(df_var.loc[df_var["FREQ"] == 1])],
                }
            ),
        ]
    ).reset_index(drop=True)
    print(f"{sample}")

out_df.to_csv(args.output, sep="\t", index=False)
