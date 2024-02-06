import polars as pl
import argparse


parser = argparse.ArgumentParser()

parser.add_argument("--support", "-s", type=int, required=True, help="Support count")
parser.add_argument(
    "--output", "-o", type=str, required=False, help="Output file to write to"
)
parser.add_argument("--matrix", "-m", type=str, required=True, help="GT Matrix file")
parser.add_argument("--callerset", "-c", type=str, required=True, help="All callerset")

args = parser.parse_args()


df = pl.read_csv(args.matrix, sep="\t").fill_nulls(".")

print("Done reading Matrix")

callerset_df = pl.read_csv(args.callerset, sep="\t", columns=["ID", "SUPPORT"], dtype={"ID": str, "SUPPORT": int})
callerset_df = callerset_df.filter(callerset_df["SUPPORT"] >= args.support).select(["ID"]).set_index("ID")

df = df.join(callerset_df, on="ID", how="inner")

callerset_df = pl.DataFrame()

print("Done merging df")

sample_list = [x for x in df.columns if x not in ["VAR", "OBS", "FREQ", "SUPPORT"]][2:]

out_df = pl.DataFrame()

df = df.with_column(pl.col("VAR", pl.lit(0))).with_column(pl.col("OBS", pl.lit(0)))

for sample in sample_list:
    df = df.with_column(
        "VAR",
        pl.when(df[sample] == "1").then(df["VAR"] + 1).otherwise(df["VAR"]),
    )
    df = df.with_column(
        "OBS",
        pl.when(df[sample] != ".").then(df["OBS"] + 1).otherwise(df["OBS"]),
    )
    df = df.with_column("FREQ", df["VAR"] / df["OBS"])
    df_var = df.filter(df["VAR"] > 0).clone()
    out_df = pl.concat(
        [
            out_df,
            pl.DataFrame(
                {
                    "SINGLETON": [len(df_var.filter((df_var["VAR"] == 1) & (df_var["FREQ"] != 1)))],
                    "POLY": [len(df_var.filter((df_var["VAR"] > 1) & (df_var["FREQ"] < 0.5)))],
                    "MAJOR": [len(df_var.filter((df_var["VAR"] > 1) & (df_var["FREQ"] >= 0.5) & (df_var["FREQ"] < 1)))],
                    "FIXED": [len(df_var.filter(df_var["FREQ"] == 1))],
                }
            ),
        ],
        ignore_index=True,
    )
    print(f"{sample}")

out_df.write_csv(args.output, delimiter="\t", with_header=True)