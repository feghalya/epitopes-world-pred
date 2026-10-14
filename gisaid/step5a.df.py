#!/usr/bin/env python3

import pickle as pkl
import gzip
import pandas as pd

p = pkl.load(gzip.open('output/pickles/all_peptides_scored.p.gz', 'rb'))

d = {'Peptide': [], 'Protein': [], 'Position': [], 'Length': [], 'Score': [], 'Count': [], 'Reference': []}
for prot in p:
    for pos in p[prot]:
        for pep in p[prot][pos]:
            peplen = len(pep)
            for score in p[prot][pos][pep]:
                n_strains = p[prot][pos][pep][score][0]
                gbref = p[prot][pos][pep][score][1]
                d['Peptide'].append(pep)
                d['Protein'].append(prot)
                d['Position'].append(pos)
                d['Length'].append(peplen)
                d['Score'].append(score)
                d['Count'].append(n_strains)
                d['Reference'].append(gbref)

df = pd.DataFrame.from_dict(d)
print('Number peptides:', df[['Peptide']].drop_duplicates().shape[0])
df.to_csv('output/tables/peptides.tsv.gz', sep='\t', index=False, compression='gzip')
