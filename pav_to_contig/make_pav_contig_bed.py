# Write per-haplotype contig coordinates using TIG_REGION INFO field
import sys
import pandas as pd

df = pd.read_csv(
    sys.stdin, sep='\t', index_col=False,
    names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'TIG_REGION', 'GT', 'SAMPLE']
)
#assert (df['TIG_REGION'].str.count(',') == 1).all()  # assume format is h1,h2
assert (df['GT'].str.count('|') == 1).all()

# Expand hap-specific columns into separate rows
df['HAP']=[[0, 1]] * len(df)
df = df.explode('HAP', ignore_index=True)
df['TIG_REGION'] = df.apply(
    lambda x:(x['TIG_REGION'].split(',')[x['HAP']]) if (x['GT'] == '1|1') else (x['TIG_REGION']), axis=1
)
df['GT'] = df.apply(
    lambda x:x['GT'].split('|')[x['HAP']], axis=1
)
df = df[df['GT'] == '1']  # alt only
# Write coordinates 
df['CONTIG'] = df['TIG_REGION'].apply(lambda x: x.split(':')[0])
df['CONTIG_POS'] = df['TIG_REGION'].apply(lambda x: int(x.split(':')[1].split('-')[0]) - 1)
df['CONTIG_END'] = df['TIG_REGION'].apply(lambda x: x.split(':')[1].split('-')[1])
df.to_csv(sys.stdout, sep='\t', columns=[
    'CONTIG', 'CONTIG_POS', 'CONTIG_END', 'REF', 'ALT', 'SAMPLE'
], header=False, index=False)
