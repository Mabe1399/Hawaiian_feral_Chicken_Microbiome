#!/bin/sh
# Phylogenetic analysis shell

# set the currrent directory

cd /Users/matiasbeckerburgos/Desktop/Master_in_Bioinformatics/Master_Thesis/Hawaiian_feral_Chicken_Microbiome/Analysis/Phylogeny/

## Multiple sequence alignement
mafft --auto --adjustdirection $1 > ASV.mafft_auto.fasta

## Trimming highly empty columns
trimal -in ASV.mafft_auto.fasta -out ASV.mafft.trimal.fasta -gt 0.03 -st 0.001

## Build ML tree 
iqtree -s ASV.mafft.trimal.fasta -st DNA -pre ASV.mafft.trimal -bb 1000