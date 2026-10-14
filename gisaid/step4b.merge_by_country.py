#!/usr/bin/env python3

from collections import defaultdict, Counter
from pyGeno.tools.UsefulFunctions import translateDNA
import pickle as pkl
import gzip
from concurrent.futures import ProcessPoolExecutor


country_origin = {'NC_045512': 'Reference'}
with open('data/strain_country.tsv', 'r') as f:
    for line in f:
        strain_id, continent, country = line.strip().split('\t')
        country_origin[strain_id] = country
countries = set(country_origin.values())


def undd(dct):
    for k in dct:
        if type(dct[k]) == defaultdict or type(dct[k]) == Counter:
            dct[k] = dict(dct[k])
        if type(dct[k]) == dict:
            undd(dct[k])


def fill(prot):
    pk = pkl.load(gzip.open('output/scored/contexts.%s.p.gz' % prot, 'rb'))
    index = {v: k for k, v in pk['index'].items()}
    p = pk['genes']
    d = {key: defaultdict(lambda: defaultdict(lambda: Counter({c: 0 for c in countries}))) for key in p[prot]['162']}
    for pos in p[prot]['162']:
        for c in p[prot]['162'][pos]:
            peptide = translateDNA(c[162:-162])
            score = p[prot]['162'][pos][c]['score']
            positions_in_index = set(p[prot]['162'][pos][c]['positions'])  # strains that contain context
            strains = [index[s] for s in positions_in_index]
            n_strains = len(strains)
            countries = Counter([country_origin[s] for s in strains])
            d[pos][peptide][score].update(countries)
    undd(d)
    return prot, d


def main(proteins):
    dct = {}
    with ProcessPoolExecutor() as ex:
        for prot, d in ex.map(fill, proteins):
            dct[prot] = d
    
    pkl.dump(dct, gzip.open('output/pickles/all_peptides_scored_by_country.p.gz', 'wb'))


if __name__ == '__main__':
    proteins = ['YP_009724389.1', 'YP_009725295.1', 'YP_009724390.1', 'YP_009724391.1', 'YP_009724392.1', 'YP_009724393.1', 'YP_009724394.1', 'YP_009724395.1', 'YP_009725318.1', 'YP_009724396.1', 'YP_009724397.2' ,'YP_009725255.1']
    main(proteins)
