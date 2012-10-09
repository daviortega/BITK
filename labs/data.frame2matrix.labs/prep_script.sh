#!/bin/bash

echo "Copying files"
cp ~/MIST/fasta/Aug2012/chea_mist_eng_Aug2012.class.fa ./
cp ~/MIST/fasta/Aug2012/chew_mist_eng_Aug2012.class.fa ./
cp ~/MIST/fasta/Aug2012/chev_mist_eng_Aug2012.fa ./

echo "Trimming"
trimbyhmm2 ~/PhD/pfammodels/PF01584.hmm chea_mist_eng_Aug2012.class.fa -m
trimbyhmm2 ~/PhD/pfammodels/PF01584.hmm chew_mist_eng_Aug2012.class.fa -m
trimbyhmm2 ~/PhD/pfammodels/PF01584.hmm chev_mist_eng_Aug2012.fa -m

echo "Inserting tag"
IIintag chea_mist_eng_Aug2012.class.hmmtrim.fa 2 CheA
IIintag chew_mist_eng_Aug2012.class.hmmtrim.fa 2 CheW
IIintag chev_mist_eng_Aug2012.hmmtrim.fa 2 CheV

echo "Selecting random $1 sequences from full files"
sel.rand.fasta chea_mist_eng_Aug2012.class.hmmtrim.IIintag.fa $1
sel.rand.fasta chew_mist_eng_Aug2012.class.hmmtrim.IIintag.fa $1
sel.rand.fasta chev_mist_eng_Aug2012.hmmtrim.IIintag.fa $1

echo "Concatenating"
cat chea_mist_eng_Aug2012.class.hmmtrim.IIintag.sel.rand.$1.fa chew_mist_eng_Aug2012.class.hmmtrim.IIintag.sel.rand.$1.fa chev_mist_eng_Aug2012.hmmtrim.IIintag.sel.rand.$1.fa > chea.chew.chev.rand$1.fa

echo "Running BAAR"
time ../BAAR.labs/BAAR1.1 chea.chew.chev.rand$1.fa

echo "Running data.frame2matrix"
time ./data.frame2matrix1.2b chea.chew.chev.rand$1.blastpall.dat -D -sep ';' -fd0

