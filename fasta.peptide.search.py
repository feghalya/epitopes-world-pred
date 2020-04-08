#!/usr/bin/env python3

import pandas as pd
from parse import parse
import sys
from model_tools import Model, available_models as am

fasta_file = 'virus_project/data/refseq/coding.fasta'
meta = pd.read_csv('virus_project/data/refseq/metadata.csv', keep_default_na=False).set_index('Accession')
meta_dct = meta.to_dict(orient='index')

# Include strain information in genotype field (for now)
meta_dct['NC_004718']['Genotype'] = 'SARS-CoV-1'
meta_dct['NC_045512']['Genotype'] = 'SARS-CoV-2'
meta_dct['NC_006213']['Genotype'] = 'Human coronavirus OC43'

peptide_lengths = [8, 9, 10, 11]
context_size_selection = [162]#None
method_selection = ['SGD', 'SGD_shuffle']
replicate_selection = [int(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else None

model_listing = am()

print('\t'.join(['Method', 'Context_size', 'Model_run', 'Accession', 'Protein_accession',
    'Position_start_one_based', 'Position_end_one_based', 'Peptide_length', 'Peptide_sequence', 'Score',
    'Protein_description', 'Sub_accession', 'Genbank_name', 'Genome_size',
    'Family', 'Genus', 'Species', 'Genotype']))

for context_size in model_listing.keys():
    for method in model_listing[context_size].keys():
        for i, model in enumerate(model_listing[context_size][method]):
            replicate = i+1
            ok_context_size = context_size_selection is None or  context_size in context_size_selection
            ok_method = method_selection is None or method in method_selection
            ok_replicate = replicate_selection is None or replicate in replicate_selection
            if ok_context_size and ok_method and ok_replicate:
                model = Model(model, context_size)
                for tup in model.discover_peptides(fasta_file, peptide_lengths, coding=True):
                    header, one_base_pos, peplen, peptide, score = tup
                    one_base_end = one_base_pos + context_size * 2 + peplen * 3
                    info = []
                    s = ''
                    stack = 0
                    for c in header:
                        if c == '(':
                            stack += 1
                        if c == ')':
                            stack -= 1
                        if not stack and c == '|':
                            info.append(s)
                            s = ''
                        else:
                            s += c
                    info.append(s)
                    if len(info) == 3:
                        sub_accession, description, protein_accession = info
                    elif len(info) == 4:
                        sub_accession, description, what, protein_accession = info
                        description = description + '|' + what
                    sub_accession = sub_accession.strip()
                    accession = sub_accession.split('.')[0].split('(')[-1]
                    protein_accession = protein_accession.strip()
                    description = parse('{p} [{n}]', description.strip())
                    protein_description = description['p']
                    genbank_name = description['n']
                    
                    species = meta_dct[accession]['Species']
                    genotype = meta_dct[accession]['Genotype']
                    genus = meta_dct[accession]['Genus']
                    family = meta_dct[accession]['Family']
                    genome_size = meta_dct[accession]['Length']
                    
                    print(method, context_size, replicate, accession, protein_accession,
                            one_base_pos, one_base_end, peplen, peptide, score,
                            protein_description, sub_accession, genbank_name, genome_size,
                            family, genus, species, genotype, sep='\t')
