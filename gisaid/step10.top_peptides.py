#!/usr/bin/env python3

import pandas as pd

df = pd.read_csv('output/tables/peptides.country.tsv.gz', sep='\t', compression='gzip')

df['MergedScore'] = df.LogBindingScore + df.Score + df.FractionCountry

#top10 = df[(df.Rank < 1) & (df.FractionCountry > 0.25)].sort_values(
#        ['Protein', 'Country', 'Peptide', 'Score'], ascending=False).drop_duplicates(
#        ['Protein', 'Country', 'Peptide']).sort_values(
#        ['Protein', 'Country', 'Score'], ascending=False).groupby(
#        ['Protein', 'Country']).head(10).set_index('Peptide').sort_values(
#        ['Protein', 'Country', 'Score'], ascending=False)

top10 = df.sort_values(
         ['Protein', 'Country', 'Peptide', 'MergedScore'], ascending=False).drop_duplicates(
         ['Protein', 'Country', 'Peptide']).sort_values(
         ['Protein', 'Country', 'MergedScore'], ascending=False).groupby(
         ['Protein', 'Country']).head(10).set_index('Peptide').drop(
         'MergedScore', axis=1)

top10.to_csv('output/tables/top10.country.tsv', sep='\t', float_format='%g')

#top100 = df[(df.Rank < 1) & (df.FractionCountry > 0.25)].sort_values(
#        ['Protein', 'Country', 'Peptide', 'Score'], ascending=False).drop_duplicates(
#        ['Protein', 'Country', 'Peptide']).sort_values(
#        ['Protein', 'Country', 'Score'], ascending=False).groupby(
#        ['Protein', 'Country']).head(100).set_index('Peptide').sort_values(
#        ['Protein', 'Country', 'Score'], ascending=False)

top100 = df.sort_values(
         ['Protein', 'Country', 'Peptide', 'MergedScore'], ascending=False).drop_duplicates(
         ['Protein', 'Country', 'Peptide']).sort_values(
         ['Protein', 'Country', 'MergedScore'], ascending=False).groupby(
         ['Protein', 'Country']).head(100).set_index('Peptide').drop(
         'MergedScore', axis=1)

top100.to_csv('output/tables/top100.country.tsv', sep='\t', float_format='%g')

cdf = df[df.ReferenceContext].drop(['Country', 'FractionCountry'], axis=1).groupby(
         ['Peptide', 'Protein', 'Position', 'Score']).sum()[['Count', 'Fraction']]

refdf = df[df.ReferenceContext].drop(
        ['Country', 'Count', 'Fraction', 'FractionCountry', 'ReferenceContext', 'ReferencePeptide'],
        axis=1).drop_duplicates().merge(cdf, on=['Peptide', 'Protein', 'Position', 'Score'])

refdf.to_csv('output/tables/public.reference.only.tsv', sep='\t', float_format='%g', index=False)
