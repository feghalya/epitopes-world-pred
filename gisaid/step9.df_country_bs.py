#!/usr/bin/env python3

import pickle as pkl
import gzip
import pandas as pd
from collections import defaultdict
import numpy as np

p = pkl.load(gzip.open('output/pickles/all_peptides_scored_by_country.p.gz', 'rb'))

nstrains = 42109

protein_names = {
    'YP_009724389.1': 'ORF1ab',
    'YP_009725295.1': 'ORF1a',
    'YP_009724390.1': 'Surface',
    'YP_009724391.1': 'ORF3a',
    'YP_009724392.1': 'Envelope',
    'YP_009724393.1': 'Membrane',
    'YP_009724394.1': 'ORF6',
    'YP_009724395.1': 'ORF7a',
    'YP_009725318.1': 'ORF7b',
    'YP_009724396.1': 'ORF8',
    'YP_009724397.2': 'Nucleocapsid',
    'YP_009725255.1': 'ORF10',
}

netmhc = {}
with open('output/tables/bs.43hla.tsv', 'r') as f:
    for line in f:
        peptide, rank, bs, logbs = line.strip().split('\t')
        assert peptide not in netmhc
        netmhc[peptide] = [rank, bs, logbs]

stats = defaultdict(lambda: [[], []])

strain_counts = {}
with open('data/strain_country_counts.tsv', 'r') as f:
    for line in f:
        strain_count, continent, country = line.strip().split('\t')
        assert country not in strain_counts
        strain_counts[country] = int(strain_count)

d = {
    'Peptide': [],
    'Length': [],
    'Protein': [],
    'Position': [],
    'Score': [],
    'LogBindingScore': [],
    'Rank': [],
    'Country': [],
    'Count': [],
    'Fraction': [],
    'FractionCountry': [],
    'ReferenceContext': []
    }

for prot in p:
    name = protein_names[prot]
    print(prot, name)
    for pos in p[prot]:
        for pep in p[prot][pos]:
            peplen = len(pep)
            rank, bs, lbs = netmhc[pep]
            for score in p[prot][pos][pep]:
                ref = 'Reference' in p[prot][pos][pep][score]
                if ref:
                    assert len(p[prot][pos][pep][score].keys()) > 1
                for country in p[prot][pos][pep][score]:
                    if country != 'Reference':
                        count = p[prot][pos][pep][score][country]
                        fra = count / nstrains
                        frac = count / strain_counts[country]
                        d['Peptide'].append(pep)
                        d['Length'].append(peplen)
                        d['Protein'].append(name)
                        d['Position'].append(pos)
                        d['Score'].append(score)
                        d['LogBindingScore'].append(lbs)
                        d['Rank'].append(rank)
                        d['Country'].append(country)
                        d['Count'].append(count)
                        d['Fraction'].append(fra)
                        d['FractionCountry'].append(frac)
                        d['ReferenceContext'].append(ref)
                        stats[(pep, name, pos)][0].append(score)
                        stats[(pep, name, pos)][1].append(count)

ref_peptides = set([p for p, r in zip(d['Peptide'], d['ReferenceContext']) if r])
d['ReferencePeptide'] = [True if p in ref_peptides else False for p in d['Peptide']]

stats_d = {'Peptide': [], 'Protein': [], 'Position': [], 'MeanScore': [], 'StdScore': []}
for (pep, prot, pos), (vals, weights) in stats.items():
    stats_d['Peptide'].append(pep)
    stats_d['Protein'].append(prot)
    stats_d['Position'].append(pos)
    mean = np.average(vals, weights=weights)
    std = np.sqrt(np.average((vals - mean)**2, weights=weights))
    stats_d['MeanScore'].append(mean)
    stats_d['StdScore'].append(std)

stats_df = pd.DataFrame.from_dict(stats_d).set_index(['Peptide', 'Protein', 'Position'])
print('Stats dataframe constructed')
#stats_df.to_csv('stats_df.country.tsv', sep='\t', float_format='%g')

df = pd.DataFrame.from_dict(d)
df.Rank = df.Rank.astype(float)  # pandas doesn't detect automatically for some reason
print('Number peptides:', df[['Peptide']].drop_duplicates().shape[0])

df = df.merge(stats_df, left_on=['Peptide', 'Protein', 'Position'], right_index=True)
print('Stats merged')

# ==> TOO SLOW
#
#def mean_std(subdf):
#    vals = [[t.Score]*t.Count for t in subdf[['Score', 'Count']].itertuples()]
#    vals = [x for y in vals for x in y]
#    return np.mean(vals), np.std(vals)

#def mean_std(subdf):
#    mean = np.average(subdf.Score, weights=subdf.Count)
#    std = np.sqrt(np.average((subdf.Score - mean)**2, weights=subdf.Count))
#    return mean, std
#
#df['MeanScore'], df['StdScore'] = zip(*df.groupby(['Peptide', 'Protein', 'Position']).apply(mean_std))

df.to_csv('output/tables/peptides.country.tsv.gz', sep='\t', index=False, float_format='%g', compression='gzip')
