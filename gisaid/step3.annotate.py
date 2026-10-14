#!/usr/bin/env python3

import sys
import pickle as pkl
import gzip

from camaptools.MLP import Model


def main(protein_id, model_name):
    model = Model(model_name, 162, 'cpu')

    coding_dct = pkl.load(gzip.open('output/sequences/contexts.%s.p.gz' % protein_id, 'rb'))

    for pid in coding_dct['genes']:
        for pos in coding_dct['genes'][pid]['162']:
            coding_dct['genes'][pid]['162'][pos] = dict(coding_dct['genes'][pid]['162'][pos])
            for context in coding_dct['genes'][pid]['162'][pos]:
                score = float(model.get_seq_score_from_nt(context)[0][1])
                dct = coding_dct['genes'][pid]['162'][pos]
                dct[context] = {'positions': dct[context], 'score': score}

    pkl.dump(coding_dct, gzip.open('output/scored/contexts.%s.p.gz' % protein_id, 'wb'))


if __name__ == '__main__':
    protein = sys.argv[1]
    main(protein, './output/pytorch-sgd/p162e4000pep9t5-7/validation-bestMin-score.pytorch')
