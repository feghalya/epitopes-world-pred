#!/usr/bin/env bash

cat data/gisaid_hcov-19_2020_06_09_03.tsv |\
    tail -n +2 | cut -f 2,4 | sed 's/[ ]*\+\/[ ]*/\//g; s/ \+$//g' |\
    tr '/' '\t' | cut -f 1,2,3 | sed 's/\xe2\x80\x8e//' |
    sed 's/america/America/' |\
    sed 's/North America\tCosta Rica/Central America\tCosta Rica/' |\
    sed 's/Europe\tGeorgia/Asia\tGeorgia/' |\
    sed 's/Europe\tMorocco/Africa\tMorocco/' |\
    sed 's/Europe\tTurkey/Asia\tTurkey/' > data/strain_country.tsv

cat data/strain_country.tsv | cut -f 2,3 | sort |\
    tr ' ' '@' | uniq -c | sed 's/[ ]\+/\t/g; s/^\t//' | tr '@' ' ' > data/strain_country_counts.tsv
