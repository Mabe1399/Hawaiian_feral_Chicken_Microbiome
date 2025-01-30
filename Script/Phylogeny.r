# Phylogenetic Analysis pipeline

# Library
library(ape)

# Set the path
path <- "/Users/matiasbeckerburgos/Desktop/Master_in_Bioinformatics/Master_Thesis/Hawaiian_feral_Chicken_Microbiome/Analysis"
script <- "/Users/matiasbeckerburgos/Desktop/Master_in_Bioinformatics/Master_Thesis/Hawaiian_feral_Chicken_Microbiome/Script"

# Create folder for Analysis
Phylogeny <- file.path(path, "Phylogeny")
if (!dir.exists(Phylogeny)) dir.create(Phylogeny)

# Create File
file.create(file.path(Phylogeny, "ASV.fasta"))

# Write sequence inside file
Biostrings::writeXStringSet(refseq(ps_rare), file.path(Phylogeny, "ASV.fasta"))

# read bash file
system2("sh", args = c(file.path(script, "Bash_phylo.sh"), file.path(Phylogeny, "ASV.fasta")))

# Transfer to the phyloseq object
Treefile <- file.path(Phylogeny, "ASV.mafft.trimal.treefile")

Tree <- ape::read.tree(Treefile)

ps_rare2 <- merge_phyloseq(ps_rare, Tree)
