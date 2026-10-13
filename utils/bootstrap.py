#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor as ppe
from functools import partial
from itertools import islice

df = pd.read_csv('virus_project/output/predictions/NC_004718.BS.tsv.gz', sep='\t', compression='gzip')
df = df[(df.Method == 'SGD') & (df.Context_size == 162)]

with open('virus_project/validation/SarsCoV1/grifoni.txt', 'r') as f:
    grifoni_peplist = set([x.strip() for x in f.readlines()])

valid = pd.read_csv('virus_project/validation/SarsCoV1/EpitopeList_SARS-COV.csv', index_col=0)
peplist = []
for pep in valid[~(valid.Tcell == 'CD4')].Sequence:
    for i in range(len(pep)):
        for peplen in [8, 9, 10, 11]:
            sub_pep = pep[i:i+peplen]
            if len(sub_pep) == peplen:
                peplist.append(sub_pep)
peplist = set(peplist)
subdf = df[df.Peptide_sequence.isin(peplist)]
subdf = subdf[subdf.Rank < 1]
peplist.update(grifoni_peplist)
subdf_wgrif = df[df.Peptide_sequence.isin(peplist)]
subdf_wgrif = subdf_wgrif[subdf_wgrif.Rank < 1]

valid = pd.read_csv('virus_project/validation/SarsCoV1/EpitopeList_SARS-COV.csv', index_col=0)
peplist_top = []
for pep in valid[~(valid.Tcell == 'CD4')].Sequence:
    this_peplist = []
    for i in range(len(pep)):
        for peplen in [8, 9, 10, 11]:
            sub_pep = pep[i:i+peplen]
            if len(sub_pep) == peplen:
                this_peplist.append(sub_pep)
    temp_df = df[df.Peptide_sequence.isin(this_peplist) & (df.Rank < 1)]
    if not temp_df.empty:
        peplist_top.append(temp_df.loc[temp_df.Score.idxmax()].Peptide_sequence)

peplist_top = set(peplist_top)
subdf_top = df[df.Peptide_sequence.isin(peplist_top)]
peplist_top.update(grifoni_peplist)
subdf_top_wgrif = df[df.Peptide_sequence.isin(peplist_top)]

# Merge pp1ab and pp1a except when it is unique to pp1a or if ANN Score is greater in pp1a
def redmax(x):
    if x.shape[0] > 1:
        if not set(x.Protein_description) - set(['orf1ab polyprotein (pp1ab)', 'orf1a polyprotein (pp1a)']):
            if x.iloc[0].Score == x.iloc[1].Score:
                return x[x.Protein_description == 'orf1ab polyprotein (pp1ab)']
            else:
                print('Taking max')
                return x.loc[x.Score.idxmax()]
        else:
            print('ERROR')
            return -1
    else:
        return x

subdf_red = subdf.groupby('Peptide_sequence').apply(redmax)
subdf_red.index = subdf_red.index.get_level_values(1)
subdf_wgrif_red = subdf_wgrif.groupby('Peptide_sequence').apply(redmax)
subdf_wgrif_red.index = subdf_wgrif_red.index.get_level_values(1)
subdf_top_red = subdf_top.groupby('Peptide_sequence').apply(redmax)
subdf_top_red.index = subdf_top_red.index.get_level_values(1)
subdf_top_wgrif_red = subdf_top_wgrif.groupby('Peptide_sequence').apply(redmax)
subdf_top_wgrif_red.index = subdf_top_wgrif_red.index.get_level_values(1)

subdf_grifoni = df[df.Peptide_sequence.isin(grifoni_peplist)]

nboot = 10000000
chunk_size = nboot//100
workers = 32

temp_dict = df[(df.Rank < 1)].groupby('Protein_accession').apply(lambda x: x.Score.to_list()).to_dict()

def get_boot(i_range):
    def calc(i):
        pd.np.random.seed(i)
        a = pd.np.concatenate([pd.np.random.choice(temp_dict[name], count) for name, count in counts.items()])
        rand = pd.np.median(a)
        return target - rand
    return [calc(x) for x in i_range]

print('N. bootsraps:', nboot)

temp_subdf = subdf_red
counts = temp_subdf.groupby('Protein_accession').count()[['Method']].rename({'Method': 'Hits'}, axis=1).Hits.to_dict()
target = temp_subdf.Score.median()

with ppe(max_workers=workers) as ex:
    iter_boot = iter(range(nboot))
    vals_ann = [y for x in ex.map(get_boot, iter(lambda: tuple(islice(iter_boot, chunk_size)), ())) for y in x]

pval = sum([x <= 0 for x in vals_ann])/nboot
print("All hits pvalue:", pval)

temp_subdf = subdf_wgrif_red
counts = temp_subdf.groupby('Protein_accession').count()[['Method']].rename({'Method': 'Hits'}, axis=1).Hits.to_dict()
target = temp_subdf.Score.median()

with ppe(max_workers=32) as ex:
    iter_boot = iter(range(nboot))
    vals_ann = [y for x in ex.map(get_boot, iter(lambda: tuple(islice(iter_boot, chunk_size)), ())) for y in x]

pval = sum([x <= 0 for x in vals_ann])/nboot
print("All hits (with grifoni) pvalue:", pval)

temp_subdf = subdf_top_red#[subdf_top_red.Rank < 1]
counts = temp_subdf.groupby('Protein_accession').count()[['Method']].rename({'Method': 'Hits'}, axis=1).Hits.to_dict()
target = temp_subdf.Score.median()

with ppe(max_workers=workers) as ex:
    iter_boot = iter(range(nboot))
    vals_ann_top = [y for x in ex.map(get_boot, iter(lambda: tuple(islice(iter_boot, chunk_size)), ())) for y in x]

pval = sum([x <= 0 for x in vals_ann_top])/nboot
print("Top hits pvalue:", pval)

temp_subdf = subdf_top_wgrif_red
counts = temp_subdf.groupby('Protein_accession').count()[['Method']].rename({'Method': 'Hits'}, axis=1).Hits.to_dict()
target = temp_subdf.Score.median()

with ppe(max_workers=workers) as ex:
    iter_boot = iter(range(nboot))
    vals_ann_top = [y for x in ex.map(get_boot, iter(lambda: tuple(islice(iter_boot, chunk_size)), ())) for y in x]

pval = sum([x <= 0 for x in vals_ann_top])/nboot
print("Top hits (with grifoni) pvalue:", pval)

temp_subdf = subdf_grifoni
counts = temp_subdf.groupby('Protein_accession').count()[['Method']].rename({'Method': 'Hits'}, axis=1).Hits.to_dict()
target = temp_subdf.Score.median()

with ppe(max_workers=workers) as ex:
    iter_boot = iter(range(nboot))
    vals_ann_top = [y for x in ex.map(get_boot, iter(lambda: tuple(islice(iter_boot, chunk_size)), ())) for y in x]

pval = sum([x <= 0 for x in vals_ann_top])/nboot
print("Grifoni hits pvalue:", pval)
