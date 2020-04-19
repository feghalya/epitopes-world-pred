#!/usr/bin/env python3

def viral_peptides_from_fasta(replicate_selection=None):
    import pandas as pd
    from parse import parse
    from model_tools import Model, available_models as am
    
    fasta_file = 'virus_project/refseq/coding.fasta'
    meta = pd.read_csv('virus_project/refseq/metadata.csv', keep_default_na=False).set_index('Accession')
    meta_dct = meta.to_dict(orient='index')
    
    # Include strain information in genotype field (for now)
    meta_dct['NC_004718']['Genotype'] = 'SARS-CoV-1'
    meta_dct['NC_045512']['Genotype'] = 'SARS-CoV-2'
    meta_dct['NC_006213']['Genotype'] = 'Human coronavirus OC43'
    
    peptide_lengths = [8, 9, 10, 11]#[:1]
    context_size_selection = [162]#None
    method_selection = ['SGD', 'SGD_shuffle']#[:1]
    
    model_listing = am()
    
    print('\t'.join(['Method', 'Context_size', 'Model_run', 'Accession', 'Protein_accession',
        'Position_start_one_based', 'Position_end_one_based', 'Peptide_length', 'Peptide_sequence', 'Score',
        'Protein_description', 'Sub_accession', 'Genbank_name', 'Genome_size',
        'Family', 'Genus', 'Species', 'Genotype']))
    
    for context_size in model_listing.keys():
        for method in model_listing[context_size].keys():
            for i, model in enumerate(model_listing[context_size][method]):
                replicate = i+1
                ok_context_size = context_size_selection is None or context_size in context_size_selection
                ok_method = method_selection is None or method in method_selection
                ok_replicate = replicate_selection is None or replicate in replicate_selection
                if ok_context_size and ok_method and ok_replicate:
                    model = Model(model, context_size)
                    for tup in model.discover_peptides(fasta_file, peptide_lengths,
                                                       coding=True):#False, whitelist=['NC_045512']):
                        header, one_base_pos, peplen, peptide, score = tup
                        # minus context_size ... if neg or plus context_size ... if pos
                        one_base_end = one_base_pos + (context_size * 2 + peplen * 3) * ((one_base_pos > 0) * 2 - 1)
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


def blcl9_split():
    from functools import partial
    import dataset_tools as dt
    
    class Manager(dt.RegressionManager):
        def submit(self):
            self.ex = self.executor(max_workers=self.main_workers)
            for t in self.trainers:
                t.dataset_name += '_' + t.variable
                fct = partial(t._data_preprocessing, retain_cols=True)
                p = self.ex.submit(fct)
                self.processes.append(p)
            print(self.ex)
            self.ex.run()  # takes care of threads in main
            self.trainers = [p.result() for p in self.processes]
            self.ex.shutdown()
            print('SUCCESS')

    manager = Manager(blcl=True, mpi=False, workers=3)
    manager.ann_versions = ['7']
    manager.get_datasets()
    for variable in ['affinity', 'rank']:
        manager.initialize_trainers(variable=variable, threshold=None, max_duplications=None)
    manager.main_workers += 2  # because we add 2 trainers
    manager.submit()

    for trainer in manager.trainers:
        trainer.df.Class = trainer.df.Class.astype(bool)
        df_pos = trainer.df[trainer.df.Class]
        df_neg = trainer.df[~trainer.df.Class]
        fn_base = 'virus_project/output/human/' + trainer.dataset_name
        df_pos.to_csv(fn_base + '.Source.7.tsv', sep='\t', float_format='%g')
        df_neg.to_csv(fn_base + '.NonSource.7.tsv', sep='\t', float_format='%g') 


def format_human_files():
    import pandas as pd
    import numpy as np
    
    dfs = []
    for status in ['Source', 'NonSource']:
        fn_1 = 'virus_project/output/human/HumanBLCL_trainDS_NetMHCpan-4.0a_affinity.' + status + '.7.tsv'
        fn_2 = 'virus_project/output/human/HumanBLCL_trainDS_NetMHCpan-4.0a_rank.' + status + '.7.tsv'
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
            'Class', 'transcript_TPM', 'NetMHC_min_transcript_affinity', 'NetMHC_min_transcript_affinity',
            'tr_is_source', 'logTPM'
            ], axis=1).melt([
            'Peptide_sequence', 'Sub_accession', 'Rank', 'Binding_score', 'Log_BS'
            ], var_name='Method', value_name='Score')

        df['Context_size'] = 162
        df['Model_run'] = 7
        df['Accession'] = 'GRCh37.75'
        df['Protein_accession'] = ''
        df['Position_start_one_based'] = np.nan
        df['Position_end_one_based'] = np.nan
        df['Peptide_length'] = df['Peptide_sequence'].str.len()
        df['Protein_description'] = ''
        df['Genome_size'] = np.nan
        df['Genbank_name'] = np.nan
        df['Family'] = 'Hominidae'
        df['Genus'] = 'Homo'
        df['Species'] = 'Homo sapiens'
        df['Genotype'] = status

        df = df[["Method", "Context_size", "Model_run", "Accession", "Protein_accession",
            "Position_start_one_based", "Position_end_one_based", "Peptide_length",
            "Peptide_sequence", "Score", "Rank", "Binding_score", "Log_BS", "Protein_description",
            "Sub_accession", "Genbank_name", "Genome_size", "Family", "Genus", "Species", "Genotype"]]

        dfs.append(df)
 
    df = pd.concat(dfs, axis=0)
    df.to_csv('virus_project/output/predictions/human.7.tsv', sep='\t', index=False, float_format='%g')


def merge_netmhc_files():
    import pandas as pd
    import os
    netmhc_folder = 'virus_project/output/peptides/NetMHCpan4.0a'
    dfs = []
    for fn_allele in os.listdir(netmhc_folder):
        allele = fn_allele.split('.')[0].replace(':', '')
        print(allele)
        fn = os.path.join(netmhc_folder, fn_allele)
        temp_df = pd.read_csv(fn, sep='\t', skiprows=[0], header=0, index_col=[0,1,2])
        dfs.append(temp_df)
        dfs[-1].columns = pd.MultiIndex.from_product([[allele], dfs[-1].columns])
    df = pd.concat(dfs, axis=1, sort=False, copy=False)
    fn_out = 'virus_project/output/peptides/viral.peptides.NetMHCpan4.0a.43HLAs.FULL.tsv.gz'
    df.to_csv(fn_out, sep='\t', float_format='%g', compression='gzip')
    
    subdf = df.loc[:, (slice(None), ['Rank'])]
    subdf.columns = ['Rank_' + x for x in list(zip((*df.loc[:, (slice(None), ['Rank'])].columns)))[0]]
    rank = subdf.reset_index().drop(['Pos', 'ID'], axis=1).set_index('Peptide')
    rank['Rank_min'] = rank.min(axis=1)
    
    subdf = df.loc[:, (slice(None), ['nM'])]
    subdf.columns = ['Binding_score_' + x for x in list(zip((*df.loc[:, (slice(None), ['nM'])].columns)))[0]]
    bs = subdf.reset_index().drop(['Pos', 'ID'], axis=1).set_index('Peptide')
    bs['Binding_score_min'] = bs.min(axis=1)
    
    subdf = df.loc[:, (slice(None), ['1-log50k'])]
    subdf.columns = ['Log_BS_' + x for x in list(zip((*df.loc[:, (slice(None), ['1-log50k'])].columns)))[0]]
    logbs = subdf.reset_index().drop(['Pos', 'ID'], axis=1).set_index('Peptide')
    logbs['Log_BS_max'] = logbs.max(axis=1)
    
    df_min = pd.concat([rank, bs, logbs], axis=1, sort=False, copy=False)
    fn_out = 'virus_project/output/peptides/viral.peptides.NetMHCpan4.0a.43HLAs.MAIN.tsv.gz'
    df_min.to_csv(fn_out, sep='\t', float_format='%g', compression='gzip')


def merge_human_virus(replicate):
    """ NOTE: Assumes we have NetMHC file, human.[replicate].tsv and viruses.[replicate].tsv ready for merge
    """
    import pandas as pd

    f = 'virus_project/output/peptides/viral.peptides.NetMHCpan4.0a.43HLAs_minBS.tsv'
    netmhc_df = pd.read_csv(f, sep='\t', header=0, index_col=0)
    netmhc_df = netmhc_df.rename({'Rank_min': 'Rank',
                                  'Binding_score_min':'Binding_score',
                                  'Log_BS_max':'Log_BS'}, axis=1)

    filename = 'virus_project/output/predictions/viruses.' + replicate + '.tsv'
    human = 'virus_project/output/predictions/human.' + replicate + '.tsv'

    df_cur = pd.read_csv(filename, sep='\t', header=0, index_col=None, dtype={'Genotype': str})
    df = df_cur.merge(netmhc_df, left_on='Peptide_sequence', right_index=True, sort=False)

    df = df[["Method", "Context_size", "Model_run", "Accession", "Protein_accession",
        "Position_start_one_based", "Position_end_one_based", "Peptide_length",
        "Peptide_sequence", "Score", "Rank", "Binding_score", "Log_BS", "Protein_description",
        "Sub_accession", "Genbank_name", "Genome_size", "Family", "Genus", "Species", "Genotype"]]

    df_hum = pd.read_csv(human, sep='\t', header=0, index_col=None, dtype={'Genotype': str})
    df = pd.concat([df, df_hum], axis=0)
    df.Genotype.fillna('', inplace=True)
    df.Accession.fillna('', inplace=True)
    df.Protein_accession.fillna('', inplace=True)
    df.Protein_description.fillna('', inplace=True)
    df.to_csv('virus_project/output/predictions/virman.%s.tsv' % replicate, sep='\t', index=False,
        float_format='%g', na_rep='nan')
