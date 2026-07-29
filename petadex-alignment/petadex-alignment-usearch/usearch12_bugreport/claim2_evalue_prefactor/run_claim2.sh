#!/usr/bin/env bash
# Claim 2: identical alignments, e-values differ ONLY by a constant prefactor K.
# Controlled comparison: version is the ONLY variable.
#   same command (-usearch_local), same input FASTA, same parameters, same post-filter.
set -eu
U11=../../bin/usearch11; U12=../../bin/usearch12
FA=../inputs/v1_nodes_213.fasta; MRE=../inputs/mre_n84_n35.fasta
PARAMS="-id 0.05 -evalue 1000 -maxaccepts 0 -maxrejects 0 -fulldp -wordlength 2"
UF="query+target+id+alnlen+evalue+bits"

# --- controlled all-vs-all (query = db = 213 seqs), both versions, identical everything ---
$U11 -usearch_local $FA -db $FA $PARAMS -userout pairs_u11.tsv -userfields $UF -threads 16
$U12 -usearch_local $FA -db $FA $PARAMS -userout pairs_u12.tsv -userfields $UF -threads 16

# --- minimal reproducible example: query = the 2 MRE seqs vs the same 213-seq db ---
$U11 -usearch_local $MRE -db $FA $PARAMS -userout mre_u11.tsv -userfields $UF -threads 16
$U12 -usearch_local $MRE -db $FA $PARAMS -userout mre_u12.tsv -userfields $UF -threads 16

python3 analyze_claim2.py
