import pandas as pd
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--pav_tab", type=str, required=True)
parser.add_argument("--pbsv_tab", type=str, required=True)
parser.add_argument("--sniffles_tab", type=str, required=True)
parser.add_argument("--pav_bg", type=str, required=True)
parser.add_argument("--pav_pbsv", type=str, required=True)
parser.add_argument("--pav_sniffles", type=str, required=True)
parser.add_argument("--pbsv_bg", type=str, required=True)
parser.add_argument("--pbsv_sniffles", type=str, required=True)
parser.add_argument("--sniffles_bg", type=str, required=True)
parser.add_argument("--output", type=str, required=True)

args = parser.parse_args()

df_pav = pd.read_csv(args.pav_tab, sep='\t')
df_pav['SVLEN'] = df_pav['ID'].str.split('-', expand=True)[3].astype(float).astype(int)
df_pav = df_pav.loc[df_pav['SVLEN'] >= 50]
df_sniffles = pd.read_csv(args.sniffles_tab, sep='\t')
df_pbsv = pd.read_csv(args.pbsv_tab, sep='\t')
int_pav_bg = pd.read_csv(args.pav_bg, sep='\t')
int_pav_bg = int_pav_bg.loc[int_pav_bg['MERGE_SAMPLES'] == 'A,B']
int_pav_pbsv = pd.read_csv(args.pav_pbsv, sep='\t')
int_pav_sniffles = pd.read_csv(args.pav_sniffles, sep='\t')
int_pbsv_bg = pd.read_csv(args.pbsv_bg, sep='\t')
int_pbsv_bg = int_pbsv_bg.loc[int_pbsv_bg['MERGE_SAMPLES'] == 'A,B']
int_pbsv_sniffles = pd.read_csv(args.pbsv_sniffles, sep='\t')
int_sniffles_bg = pd.read_csv(args.sniffles_bg, sep='\t')
int_sniffles_bg = int_sniffles_bg.loc[int_sniffles_bg['MERGE_SAMPLES'] == 'A,B']

df_pav['ID_PAV'] = df_pav['ID']

pav_pbsv_shared = int_pav_pbsv.loc[int_pav_pbsv['MERGE_SAMPLES'] == 'A,B'].copy()
pbsv_only = int_pav_pbsv.loc[int_pav_pbsv['MERGE_SAMPLES'] == 'B'].copy()
pav_pbsv_shared['ID_PBSV'] = pav_pbsv_shared['MERGE_VARIANTS'].str.split(',', expand=True)[1]
pav_pbsv_shared['PBSV'] = 'PBSV'


out_df = df_pav.merge(pav_pbsv_shared[['ID_PBSV', 'PBSV', 'ID']], how='outer')

pav_sniffles_shared = int_pav_sniffles.loc[int_pav_sniffles['MERGE_SAMPLES'] == 'A,B'].copy()
sniffles_only = int_pav_sniffles.loc[int_pav_sniffles['MERGE_SAMPLES'] == 'B'].copy()
pav_sniffles_shared['ID_SNIFFLES'] = pav_sniffles_shared['MERGE_VARIANTS'].str.split(',', expand=True)[1]
pav_sniffles_shared['SNIFFLES'] = 'SNIFFLES'

out_df = out_df.merge(pav_sniffles_shared[['ID_SNIFFLES', 'SNIFFLES', 'ID']], how='outer')

pbsv_sniffles_shared = int_pbsv_sniffles.loc[int_pbsv_sniffles['MERGE_SAMPLES'] == 'A,B'].copy()
pbsv_sniffles_shared['ID_SNIFFLES'] = pbsv_sniffles_shared['MERGE_VARIANTS'].str.split(',', expand=True)[1]
pbsv_sniffles_shared['ID_PBSV'] = pbsv_sniffles_shared['MERGE_VARIANTS'].str.split(',', expand=True)[0]

df_pbsv['ID_PBSV'] = df_pbsv['ID']
df_sniffles['ID_SNIFFLES'] = df_sniffles['ID']

pbsv_sniffles = df_pbsv.merge(pbsv_sniffles_shared[['ID', 'ID_PBSV', 'ID_SNIFFLES']])

out_df = out_df.append(pbsv_sniffles.loc[(~pbsv_sniffles['ID'].isin(out_df['ID_PBSV'])) & (~pbsv_sniffles['ID_SNIFFLES'].isin(out_df['ID_SNIFFLES']))][['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'ID_PBSV', 'ID_SNIFFLES']])

out_df = out_df.append(df_pbsv.loc[~df_pbsv['ID'].isin(out_df['ID_PBSV'])][['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'ID_PBSV']])
out_df = out_df.append(df_sniffles.loc[~df_sniffles['ID'].isin(out_df['ID_SNIFFLES'])][['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'ID_SNIFFLES']])


out_df['PAV'] = out_df.apply(lambda row: True if type(row['ID_PAV']) == str else False, axis=1)
out_df['PBSV'] = out_df.apply(lambda row: True if type(row['ID_PBSV']) == str else False, axis=1)
out_df['SNIFFLES'] = out_df.apply(lambda row: True if type(row['ID_SNIFFLES']) == str else False, axis=1)

assert len(out_df.loc[out_df['PBSV'] == True]) == len(df_pbsv)
assert len(out_df.loc[out_df['PAV'] == True]) == len(df_pav)
assert len(out_df.loc[out_df['SNIFFLES'] == True]) == len(df_sniffles)

out_df['PAV_BG'] = out_df.apply(lambda row : True if row['ID_PAV'] in int_pav_bg['MERGE_SRC_ID'].values else False, axis=1 )
out_df['PBSV_BG'] = out_df.apply(lambda row : True if row['ID_PBSV'] in int_pbsv_bg['MERGE_SRC_ID'].values else False, axis=1 )
out_df['SNIFFLES_BG'] = out_df.apply(lambda row : True if row['ID_SNIFFLES'] in int_sniffles_bg['MERGE_SRC_ID'].values else False, axis=1 )

out_df['SEEN'] = out_df.apply(lambda row : 'HGSVC_HPRC' if True in [row['PAV_BG'], row['PBSV_BG'], row['SNIFFLES_BG']] else 'NOVEL', axis=1)

out_df['CALLERSET_LIST'] = out_df.apply(lambda row: ','.join([x for x in ['PAV', 'PBSV', 'SNIFFLES'] if row[x] == True]), axis=1)

out_df.to_csv(args.output, sep='\t', index=False)


