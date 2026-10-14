#!/usr/bin/env python

import pandas as pd
import os

netmhc_folder = 'output/netmhc/gisaid'

netmhc_files = [fn for fn in os.listdir(netmhc_folder) if fn.startswith('HLA-')]
alleles = set([fn.split('.')[0] for fn in netmhc_files])
#alleles = set(['A01:01'])

dfs = []
for al in alleles:
    allele = al.replace(':', '')
    print(allele)
    aldfs = []
    for fn_allele in [f for f in netmhc_files if al in f]:
        fn = os.path.join(netmhc_folder, fn_allele)
        tdf = pd.read_csv(fn, sep='\t', skiprows=[0], header=0, index_col=[0,1,2])
        aldfs.append(tdf)
    temp_df = pd.concat(aldfs, axis=0, sort=False, copy=False)
    dfs.append(temp_df)
    dfs[-1].columns = pd.MultiIndex.from_product([[allele], dfs[-1].columns])

df = pd.concat(dfs, axis=1, sort=False, copy=False)
#fn_out = os.path.join('output/tables', 'bs.43hla.full.tsv.gz')
#df.to_csv(fn_out, sep='\t', float_format='%g', compression='gzip')

subdf = df.loc[:, (slice(None), ['Rank'])]
subdf.columns = ['Rank_' + x for x in list(zip(*df.loc[:, (slice(None), ['Rank'])].columns))[0]]
rank = subdf.reset_index().drop(['Pos', 'ID'], axis=1).set_index('Peptide')
rank['MinRank'] = rank.min(axis=1)

subdf = df.loc[:, (slice(None), ['nM'])]
subdf.columns = ['BS_' + x for x in list(zip(*df.loc[:, (slice(None), ['nM'])].columns))[0]]
bs = subdf.reset_index().drop(['Pos', 'ID'], axis=1).set_index('Peptide')
bs['MinBS'] = bs.min(axis=1)

subdf = df.loc[:, (slice(None), ['1-log50k'])]
subdf.columns = ['LogBS_' + x for x in list(zip(*df.loc[:, (slice(None), ['1-log50k'])].columns))[0]]
logbs = subdf.reset_index().drop(['Pos', 'ID'], axis=1).set_index('Peptide')
logbs['MaxLogBS'] = logbs.max(axis=1)

df_min = pd.concat([rank, bs, logbs], axis=1, sort=False, copy=False)
fn_out = os.path.join('output/tables', 'bs.43hla.all.alleles.tsv.gz')
df_min.to_csv(fn_out, sep='\t', float_format='%g', compression='gzip')

df_slim = df_min[['MinRank', 'MinBS', 'MaxLogBS']]
fn_out = os.path.join('output/tables', 'bs.43hla.tsv')
df_slim.to_csv(fn_out, sep='\t', float_format='%g')
