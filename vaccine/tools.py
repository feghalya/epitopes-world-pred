#!/usr/bin/env python3

def genbank_to_fasta(fn):
    from Bio import SeqIO
    import pandas as pd

    base = fn.rsplit('/', 1)[0]
    genbank = SeqIO.parse(fn, "genbank")
    coding = []
    protein = []
    for record in genbank:
        record_id = record.id
        name = record.name
        # genome
        genome_seq = record.seq.__str__()
        # metadata
        description = record.description
        date = record.annotations['date']
        source = record.annotations['source']
        organism = record.annotations['organism']
        entity, realm, order, suborder, family, subfamily, genus, subgenus = record.annotations['taxonomy']
        meta_header = ['Accession', 'Version', 'Genbank_date', 'Description', 'Source', 'Organism',
                       'Taxonomic_entity', 'Realm', 'Order', 'Suborder', 'Family', 'Subfamily', 'Genus', 'Subgenus']
        meta_content = [name, record_id.split('.')[-1], date, description, source, organism,
                        entity, realm, order, suborder, family, subfamily, genus, subgenus]

        try:
            meta_df = pd.read_csv(base + '/' + 'metadata.tsv', sep='\t')
        except FileNotFoundError:
            meta_df = pd.DataFrame(columns=meta_header)
        if name not in meta_df.Accession:
            meta_df.loc[meta_df.shape[0]] = meta_content
            meta_df.to_csv(base + '/' + 'metadata.tsv', sep='\t', index=False)

        with open(base + '/' + name + '.genome.fa', 'w') as f:
            f.write('>' + record_id + '\n')
            for i in range(0, len(genome_seq), 60):
                f.write(genome_seq[i:i+60] + '\n')

        for feature in record.features:
            if feature.type == 'CDS':
                try:
                    gene = ','.join(feature.qualifiers['gene'])  # expecting only one entry
                except KeyError:
                    gene = ''
                try:
                    gene += ' (' + ','.join(feature.qualifiers['gene'] + feature.qualifiers['gene_synonym']) + ')'
                except KeyError:
                    pass
                coding_header = ' |'.join([
                    record_id + ':' + feature.location.__str__(),
                    gene,
                    feature.qualifiers['product'][0],
                    feature.qualifiers['protein_id'][0]
                    ])
                protein_header = ' |'.join([
                    feature.qualifiers['protein_id'][0],
                    gene,
                    feature.qualifiers['product'][0]
                    ])
                coding_seq = feature.extract(record).seq.__str__()
                protein_seq = feature.translate(record).seq.__str__()

                coding.append((coding_header, coding_seq))
                protein.append((protein_header, protein_seq))

        with open(base + '/' + name + '.coding.fa', 'w') as f:
            for header, sequence in coding:
                f.write('>' + header + '\n')
                for i in range(0, len(sequence), 60):
                    f.write(sequence[i:i+60] + '\n')

        with open(base + '/' + name + '.protein.fa', 'w') as f:
            for header, sequence in protein:
                f.write('>' + header + '\n')
                for i in range(0, len(sequence), 60):
                    f.write(sequence[i:i+60] + '\n')


def viral_peptides_from_fasta(accession='NC_004718', replicate_selection=7):
    import pandas as pd
    from parse import parse
    from vaccine.model_tools import Model
    from camaptools.utils import available_models as am

    if type(replicate_selection) == str or type(replicate_selection) == int:
        replicate_selection = [int(replicate_selection)]

    fasta_file = 'data/%s.coding.fa' % accession

    peptide_lengths = [8, 9, 10, 11]
    context_size_selection = [9, 27, 81, 162]
    method_selection = ['SGD', 'SGD_Shuffle']
    parameters = 'e4000'
    target = 'test-bestMin-score.pytorch'

    print('\t'.join(['Method', 'Context_size', 'Model_run', 'Accession', 'Sub_accession',
        'Peptide_start_one_based', 'Peptide_length', 'Peptide_sequence', 'Score',
        'Gene_name', 'Protein_accession', 'Protein_description', 'Location']))

    # camaptools lists one context at a time and keys on {method: {params: [files]}}
    for context_size in context_size_selection:
        model_listing = am(target, context_size, method_selection, [parameters])
        for method in model_listing.keys():
            for i, model in enumerate(model_listing[method][parameters]):
                replicate = i+1
                ok_replicate = replicate_selection is None or replicate in replicate_selection
                if ok_replicate:
                    model = Model(model, context_size)
                    for tup in model.discover_peptides(fasta_file, peptide_lengths, coding=True):
                        header, one_base_peptide_start, peptide, score = tup

                        peplen = len(peptide)
                        location, gene_name, protein_description, protein_accession = header.split(' |')

                        print(method, context_size, replicate, accession, 'CDS_of_' + protein_accession,
                                one_base_peptide_start, peplen, peptide, score,
                                gene_name, protein_accession, protein_description, location, sep='\t')


def merge_netmhc_files():
    import pandas as pd
    import os
    netmhc_folder = 'output/netmhc/sarscov'
    dfs = []
    for fn_allele in os.listdir(netmhc_folder):
        allele = fn_allele.split('.')[0].replace(':', '')
        print(allele)
        fn = os.path.join(netmhc_folder, fn_allele)
        temp_df = pd.read_csv(fn, sep='\t', skiprows=[0], header=0, index_col=[0,1,2])
        dfs.append(temp_df)
        dfs[-1].columns = pd.MultiIndex.from_product([[allele], dfs[-1].columns])
    df = pd.concat(dfs, axis=1, sort=False, copy=False)
    fn_out = 'output/tables/sarscov.peptides.NetMHCpan4.0a.43HLAs.FULL.tsv.gz'
    df.to_csv(fn_out, sep='\t', float_format='%g', compression='gzip')

    subdf = df.loc[:, (slice(None), ['Rank'])]
    subdf.columns = ['Rank_' + x for x in list(zip(*df.loc[:, (slice(None), ['Rank'])].columns))[0]]
    rank = subdf.reset_index().drop(['Pos', 'ID'], axis=1).set_index('Peptide')
    rank['Rank_min'] = rank.min(axis=1)

    subdf = df.loc[:, (slice(None), ['nM'])]
    subdf.columns = ['Binding_score_' + x for x in list(zip(*df.loc[:, (slice(None), ['nM'])].columns))[0]]
    bs = subdf.reset_index().drop(['Pos', 'ID'], axis=1).set_index('Peptide')
    bs['Binding_score_min'] = bs.min(axis=1)

    subdf = df.loc[:, (slice(None), ['1-log50k'])]
    subdf.columns = ['Log_BS_' + x for x in list(zip(*df.loc[:, (slice(None), ['1-log50k'])].columns))[0]]
    logbs = subdf.reset_index().drop(['Pos', 'ID'], axis=1).set_index('Peptide')
    logbs['Log_BS_max'] = logbs.max(axis=1)

    df_min = pd.concat([rank, bs, logbs], axis=1, sort=False, copy=False)
    fn_out = 'output/tables/sarscov.peptides.NetMHCpan4.0a.43HLAs.MAIN.tsv.gz'
    df_min.to_csv(fn_out, sep='\t', float_format='%g', compression='gzip')


def join_netmhc_to_virus(accession):
    """ NOTE: Assumes we have NetMHC file ready for merge
    """
    import pandas as pd

    fn_out = 'output/tables/sarscov.peptides.NetMHCpan4.0a.43HLAs.MAIN.tsv.gz'
    viruses = 'output/predictions/' + accession + '.tsv'

    netmhc_df = pd.read_csv(netmhc, sep='\t', header=0, index_col=0,
        usecols=['Peptide', 'Rank_min', 'Binding_score_min', 'Log_BS_max'], compression='gzip')
    netmhc_df = netmhc_df.rename({'Rank_min': 'Rank',
                                  'Binding_score_min':'Binding_score',
                                  'Log_BS_max':'Log_BS'}, axis=1)
    print(netmhc_df.head())
    netmhc_df = netmhc_df[['Rank', 'Binding_score', 'Log_BS']]

    df_cur = pd.read_csv(viruses, sep='\t', header=0, index_col=None)
    df = df_cur.merge(netmhc_df, left_on='Peptide_sequence', right_index=True, sort=False, how='left')

    df = df[['Method', 'Context_size', 'Model_run', 'Accession', 'Sub_accession',
        'Peptide_start_one_based', 'Peptide_length', 'Peptide_sequence', 'Score',
        'Rank', 'Binding_score', 'Log_BS',
        'Gene_name', 'Protein_accession', 'Protein_description', 'Location']]

    df.to_csv('output/predictions/' + accession + '.BS.tsv', sep='\t', index=False,
        float_format='%g', na_rep='nan')


def format_human_files():
    import pandas as pd
    import numpy as np

    dfs = []
    for status in ['Source', 'NonSource']:
        fn_1 = 'output/finalTables/uniqpeptides/HumanBLCL_trainDS_NetMHCpan-4.0a_affinity.' + status + '.7.tsv'
        fn_2 = 'output/finalTables/uniqpeptides/HumanBLCL_trainDS_NetMHCpan-4.0a_rank.' + status + '.7.tsv'
        df_1 = pd.read_csv(fn_1, sep='\t')
        df_2 = pd.read_csv(fn_2, sep='\t')
        for col_name in df_1.columns.intersection(df_2.columns):
            assert df_1[col_name].equals(df_2[col_name])
        df = df_1.merge(df_2[df_2.columns.difference(df_1.columns)], left_index=True, right_index=True)

        df = df.rename({
            'peptide': 'Peptide_sequence',
            'transcript': 'Sub_accession',
            'NetMHC_min_rank': 'Rank',
            'NetMHC_min_affinity': 'Binding_score',
            'logBS': 'Log_BS',
            'ANNScore': 'SGD',
            'ANNShuffleScore': 'SGD_shuffle'
            }, axis=1).drop([
            'Class', 'transcript_TPM', 'NetMHC_min_transcript_affinity', 'NetMHC_min_transcript_rank',
            'tr_is_source', 'logTPM'
            ], axis=1).melt([
            'Peptide_sequence', 'Sub_accession', 'Rank', 'Binding_score', 'Log_BS'
            ], var_name='Method', value_name='Score')

        df['Context_size'] = 162
        df['Model_run'] = 7
        df['Accession'] = 'GRCh37.75[%s]' % {"Source": "Hits", "NonSource":"Decoy"}[status]
        df['Peptide_start_one_based'] = np.nan
        df['Peptide_length'] = df['Peptide_sequence'].str.len()
        df['Gene_name'] = ''
        df['Protein_accession'] = ''
        df['Protein_description'] = ''
        df['Location'] = ''
        #df['Genome_size'] = np.nan
        #df['Genbank_name'] = np.nan
        #df['Family'] = 'Hominidae'
        #df['Genus'] = 'Homo'
        #df['Species'] = 'Homo sapiens'
        #df['Genotype'] = status

        df = df[["Method", "Context_size", "Model_run", "Accession", "Sub_accession",
            "Peptide_start_one_based", "Peptide_length",
            "Peptide_sequence", "Score", "Rank", "Binding_score", "Log_BS",
            "Gene_name", "Protein_accession", "Protein_description", "Location"]]

        dfs.append(df)

    df = pd.concat(dfs, axis=0)
    df.to_csv('output/predictions/human.7.tsv.gz', sep='\t',
        index=False, float_format='%g', compression='gzip')


def join_human_virus_predictions():
    """ NOTE: Assumes we have final human and Sars-CoV viruses ready
    """
    import pandas as pd

    human = pd.read_csv("output/predictions/human.7.tsv.gz", sep="\t")
    human = human[human.Method == "SGD"]
    human["Name"] = human.Accession
    human.Name = human.Name.str.replace("GRCh37.75", "Human")
    print("Human Hits %d Peptides" % len(human[human.Name == "Human[Hits]"].Peptide_sequence.unique()))
    print("Human Hits %d Entries" % len(human[human.Name == "Human[Hits]"].Peptide_sequence))
    print("Human Decoy %d Peptides" % len(human[human.Name == "Human[Decoy]"].Peptide_sequence.unique()))
    print("Human Decoy %d Entries" % len(human[human.Name == "Human[Decoy]"].Peptide_sequence))

    v1 = pd.read_csv("output/predictions/NC_045512.BS.tsv.gz", sep="\t")
    v1 = v1[(v1.Method == "SGD") & (v1.Context_size == 162) & (v1.Protein_description != "ORF1a polyprotein")]
    v1["Name"] = "SARS-CoV-2"
    print("SARS-CoV-2 %d Peptides" % len(v1.Peptide_sequence.unique()))
    print("SARS-CoV-2 %d Entries" % len(v1.Peptide_sequence))

    v2 = pd.read_csv("output/predictions/NC_004718.BS.tsv.gz", sep="\t")
    v2 = v2[(v2.Method == "SGD") & (v2.Context_size == 162) & (v2.Protein_description != "orf1a polyprotein (pp1a)")]
    v2["Name"] = "SARS-CoV-1"
    print("SARS-CoV-1 %d Peptides" % len(v2.Peptide_sequence.unique()))
    print("SARS-CoV-1 %d Entries" % len(v2.Peptide_sequence))

    df = pd.concat([v1, v2, human])
    df.to_csv("output/predictions/virman.filtered.tsv.gz",
        sep="\t", float_format="%g", index=False, compression="gzip")
