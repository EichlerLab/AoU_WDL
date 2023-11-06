import sys
import os
import pandas as pd
import glob

ctg_region_file = sys.argv[1]
out_hap1 = sys.argv[2]
out_hap2 = sys.argv[3]

with open(out_hap1, 'w'), open(out_hap2, 'w'):
    pass
        
df = pd.read_csv(
    ctg_region_file, sep=' ', index_col=False,
    names=[
        'CONTIG', 'CONTIG_POS', 'CONTIG_END', 'ORIGCONTIGS',
        'ORIGPOS', 'ORIGEND', 'REF', 'ALT', 'SAMPLE'
    ]
)
if df.empty:
    sys.exit()
sample = df.loc[0]['SAMPLE']
assert df['CONTIG'].str.startswith(('h1', 'h2')).all()
df['OUT_STR'] = df.apply(
    lambda x: f"{x['CONTIG']}:{x['CONTIG_POS']}-{x['CONTIG_END']}", axis=1
)
df_h1 = df[df['CONTIG'].str.startswith('h1')]
df_h2 = df[df['CONTIG'].str.startswith('h2')]
if not df_h1.empty:
    df_h1.to_csv(f"{out_hap1}", index=False, header=False, columns=['OUT_STR'])
if not df_h2.empty:
    df_h2.to_csv(f"{out_hap2}", index=False, header=False, columns=['OUT_STR'])
