#!/usr/bin/env python3

import sys

def discover_peptides(sequence, context, peptide_lengths):
    for peplen in peptide_lengths:
        assert len(sequence) > peplen

        ntpeplen = peplen*3
        for i in range(0, len(sequence) - ntpeplen, 3):
            nt_peptide = sequence[i:i+ntpeplen]
            before = sequence[max(0, i-context):i]
            after = sequence[min(len(sequence), i+ntpeplen):i+ntpeplen+context]
            seq = before + nt_peptide + after

            if set(seq) - set('ATCGUatcgu'):
                pass
                #sys.stderr.write('Non ATCGU encoutered in peptide or context\n')
            else:
                if len(before) < context:
                    before = ''.join(['0'] * (context - len(before))) + before
                assert len(before) == context

                if len(after) < context:
                    after = after + ''.join(['0'] * (context - len(after)))
                assert len(after) == context

                seq = before + nt_peptide + after

                #score = float(model.get_seq_score_from_nt(seq)[0][1])
                zero_base_peptide_start = i

                yield seq, zero_base_peptide_start


def main(protein_id):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.Alphabet import generic_dna
    import os
    import re
    from collections import defaultdict
    import numpy as np
    import pickle as pkl
    import gzip

    coding_dct = {}
    coding_dct['genes'] = {}
    coding_dct['index'] = {}

    fn = 'NC_045512.gb'
    gb_record = next(SeqIO.parse(fn, 'genbank'))
    # accession with version
    record_id = gb_record.id
    # accession
    name = gb_record.name
    # genome
    genome_seq = gb_record.seq.__str__()
    # metadata
    #description = gb_record.description
    #date = gb_record.annotations['date']
    #source = gb_record.annotations['source']
    #organism = gb_record.annotations['organism']
    #entity, realm, order, suborder, family, subfamily, genus, subgenus = gb_record.annotations['taxonomy']

    for feature in gb_record.features:
        if feature.type == 'CDS':
            # protein id
            pid = feature.qualifiers['protein_id'][0]
            # protein/gene name
            pro = feature.qualifiers['product'][0]
            # gene name with aliases
            try:
                gene = ','.join(feature.qualifiers['gene'])  # expecting only one entry
            except KeyError:
                gene = ''
            try:
                gene += ' (' + ','.join(feature.qualifiers['gene'] + feature.qualifiers['gene_synonym']) + ')'
            except KeyError:
                pass
            # location
            locations = []
            for part in feature.location.parts:
                locations.append((int(part.start), int(part.end)))
            loc = feature.location.__str__()
            # dna and protein sequences
            coding_seq = feature.extract(gb_record).seq.__str__()
            protein_seq = feature.translate(gb_record).seq.__str__()

            if pid == protein_id:
                coding_dct['genes'][pid] = {}
                coding_dct['genes'][pid]['info'] = (pid, gene, pro, loc)
                coding_dct['genes'][pid]['location'] = locations
                coding_dct['genes'][pid]['162'] = defaultdict(lambda: defaultdict(list))

    coding_dct['index']['NC_045512'] = 0
    for i, f in enumerate(os.listdir('data/gisaid.alignment')):
        coding_dct['index'][f.replace('.aln.fa', '')] = i+1

    for k, v in coding_dct['genes'].items():
        location = v['location']
        seq = ''.join([genome_seq[i:j] for (i, j) in location])
        protein = Seq(seq, generic_dna).translate()
        assert protein[-1] == '*'
        for tup in discover_peptides(seq, 162, [8, 9, 10, 11]):
            seq, pos = tup
            coding_dct['genes'][k]['162'][pos][seq].append(0)

    for i, f in enumerate(os.listdir('data/gisaid.alignment')):
        if not i % 100:
            print('*********\n%d done\n*********' % i)

        gisaid, ref = SeqIO.parse(os.path.join('data/gisaid.alignment', f), 'fasta')
        ref_seq = ref.seq.__str__()
        gisaid_seq = gisaid.seq.__str__()

        offset = np.zeros(len(genome_seq), dtype=int)
        if ref_seq != genome_seq:
            positions = [m.start() for m in re.finditer('-', ref_seq)]
            insertions = [gisaid_seq[p] for p in positions]
            s = re.sub('[a-zA-Z]', '0', ref_seq)
            s = re.sub('-', '1', s)
            s = np.cumsum([int(x) for x in s])
            s.dtype = int
            offset = s

        # get peptides
        for k, v in coding_dct['genes'].items():
            location = v['location']
            oriseq = ''.join([ref_seq[i+offset[i]:j+offset[j]] for (i, j) in location])
            #p = Seq(oriseq.replace('-', ''), generic_dna).translate()
            #assert p[-1] == '*' and p[0] == 'M'
            seq = ''.join([gisaid_seq[i+offset[i]:j+offset[j]] for (i, j) in location])
            seq_short = seq.replace('-', '')
            if not seq_short:
                print('NONUCLEOTIDE:', f, v['info'][:3])
            else:
                if seq[0] == '-' and seq[-1] == '-':
                    print('NOANCHORENDS:', f, v['info'][:3])
                else:
                    if seq[0] == '-':
                        seq_short = seq_short[len(seq_short) % 3:]
                    elif seq[-1] == '-':
                        seq_short = seq_short[:-(len(seq_short) % 3)] if len(seq_short) % 3 else seq_short
                    if seq_short == '' or len(seq_short) % 3:
                        print('INCORRECTLEN:', f, v['info'][:3])
                    else:
                        protein = Seq(seq_short, generic_dna).translate()
                        if (protein[0] != 'M' and protein[0] != 'X') or (protein[-1] != '*' and protein[-1] != 'X'):
                            print('TRUNCATEDEND:', f, v['info'][:3])
                        else:
                            if '*' in protein[:-1]:
                                print('NONSENSESTOP:', f, v['info'][:3])
                            else:
                                seq = seq[:-3]
                                gis_offset = [x for x, y in zip(range(len(seq)), seq) if y != '-']
                                oriseq = oriseq[:-3]
                                l = []
                                c = 0
                                c2 = 0
                                for x in oriseq:
                                    if x == '-':
                                        c2 += 1
                                        l.append(c)
                                    else:
                                        c += c2
                                        c2 = 0
                                        l.append(c)
                                ori_offset = l
                                for tup in discover_peptides(seq_short, 162, [8, 9, 10, 11]):
                                    cs, pos = tup
                                    new_pos = gis_offset[pos] - ori_offset[gis_offset[pos]]
                                    assert new_pos < len(seq)-24+1
                                    ix = coding_dct['index'][f.replace('.aln.fa', '')]
                                    coding_dct['genes'][k]['162'][new_pos][cs].append(ix)

    for pid in coding_dct['genes']:
        coding_dct['genes'][pid]['162'] = dict(coding_dct['genes'][pid]['162'])
    pkl.dump(coding_dct, gzip.open('output/sequences/contexts.%s.p.gz' % protein_id, 'wb'))


if __name__ == '__main__':
    protein = sys.argv[1]
    main(protein)
