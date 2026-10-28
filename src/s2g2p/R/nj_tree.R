#!/usr/bin/env Rscript
# Neighbour-joining tree of a core-gene alignment, for phylogeny-aware folds: Hamming distances,
# then phangorn's NJ, written as Newick. This is AMR_benchmarking src/cv_folds/phylo_tree.r
# (Hu et al. 2024, MIT) without its option parser; it loaded phytools too, which it never used.
#
#   nj_tree <core_gene_alignment.aln> <out.newick>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: nj_tree <alignment.fasta> <out.newick>")
suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
})
aln <- read.phyDat(args[1], format = "fasta", type = "DNA")
dm <- dist.hamming(aln)
tree <- NJ(dm)
write.tree(tree, file = args[2])
cat(sprintf("nj_tree: %d leaves, phangorn %s, ape %s -> %s\n", length(tree$tip.label),
            packageVersion("phangorn"), packageVersion("ape"), args[2]))
