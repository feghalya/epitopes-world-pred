#!/usr/bin/env python3

from collections import defaultdict
from pyGeno.tools.UsefulFunctions import translateDNA
import pickle as pkl
import gzip
from concurrent.futures import ProcessPoolExecutor


def undd(dct):
    for k in dct:
        if type(dct[k]) == defaultdict:
            dct[k] = dict(dct[k])
        if type(dct[k]) == dict:
            undd(dct[k])


def fill(prot):
    p = pkl.load(gzip.open('output/scored/contexts.%s.p.gz' % prot, 'rb'))['genes']
    d = {key: defaultdict(lambda: defaultdict(lambda: [0, False])) for key in p[prot]['162']}
    for pos in p[prot]['162']:
        for c in p[prot]['162'][pos]:
            peptide = translateDNA(c[162:-162])
            score = p[prot]['162'][pos][c]['score']
            positions_in_index = set(p[prot]['162'][pos][c]['positions'])  # strains containing context
            n_strains = len(positions_in_index)
            reference = 0 in positions_in_index
            d[pos][peptide][score][0] += n_strains
            d[pos][peptide][score][1] = d[pos][peptide][score][1] or reference
    undd(d)
    return prot, d


def main(proteins):
    dct = {}
    with ProcessPoolExecutor() as ex:
        for prot, d in ex.map(fill, proteins):
            dct[prot] = d
    
    pkl.dump(dct, gzip.open('output/pickles/all_peptides_scored.p.gz', 'wb'))


if __name__ == '__main__':
    proteins = ['YP_009724389.1', 'YP_009725295.1', 'YP_009724390.1', 'YP_009724391.1', 'YP_009724392.1', 'YP_009724393.1', 'YP_009724394.1', 'YP_009724395.1', 'YP_009725318.1', 'YP_009724396.1', 'YP_009724397.2' ,'YP_009725255.1']
    main(proteins)
