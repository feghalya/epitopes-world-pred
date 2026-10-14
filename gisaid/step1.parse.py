#!/usr/bin/env python3

from mpi4py.futures import MPIPoolExecutor as Executor
#from concurrent.futures import ProcessPoolExecutor as Executor

from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
import gzip
import os
import sys

def align_with_muscle(ix):
    in_file = 'data/gisaid.fasta/' + ix + '.fa'
    out_file = 'data/gisaid.alignment/' + ix + '.aln.fa'
    cmd = MuscleCommandline(input=in_file, out=out_file)
    stdout, stderr = cmd()

def main(workers):
    os.makedirs('data/gisaid.fasta', exist_ok=True)
    os.makedirs('data/gisaid.alignment', exist_ok=True)
    
    f = 'data/NC_045512.gb'
    
    gb_record = next(SeqIO.parse(f, 'genbank'))
    
    fasta = SeqIO.parse(gzip.open('data/gisaid_hcov-19_2020_06_09_03.fasta.gz', 'rt'), 'fasta')
    
    indexes = []
    for record in fasta:
        ix = record.description.split('|')[1]
        assert 'EPI' in ix
        indexes.append(ix)
        SeqIO.write([record, gb_record], 'data/gisaid.fasta/' + ix + '.fa', 'fasta')

    with Executor(max_workers=workers) as ex:
        ex.map(align_with_muscle, indexes)

if __name__ == '__main__':
    workers = int(sys.argv[1]) - 1
    main(workers)
