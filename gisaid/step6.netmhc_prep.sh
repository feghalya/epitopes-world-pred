#!/usr/bin/env bash

PEPTIDES=output/tables/peptides.tsv.gz
zcat $PEPTIDES | awk -F '\t' '$4 == 8' | cut -f 1 | uniq | sort -u > data/peptides.8.txt
zcat $PEPTIDES | awk -F '\t' '$4 == 9' | cut -f 1 | uniq | sort -u > data/peptides.9.txt
zcat $PEPTIDES | awk -F '\t' '$4 == 10' | cut -f 1 | uniq | sort -u > data/peptides.10.txt
zcat $PEPTIDES | awk -F '\t' '$4 == 11' | cut -f 1 | uniq | sort -u > data/peptides.11.txt

cat validation/alleleTargetList.tsv | cut -f 2 | tail -n +2 | tr -d '*' > data/alleles.txt
