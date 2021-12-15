rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(ape)
require(ggtree)
require(homoplasyFinder)
files <- list.files("data/alignments/", pattern = "aln.fasta")
prefixes <- separate(tibble(x = files), x, into = c("prefixes"), sep = "\\.")$prefixes

experiment_type <- "mini_animal_trees"
prefixes <- list.files(str_glue("data/alignments/{experiment_type}"))
prefixes <- prefixes[grepl(".aln.fasta", prefixes)]
prefixes <- gsub(".audacity_only.v8_masked.aln.fasta", "", prefixes)

for (prefix in prefixes) {
  tree_path <- str_glue("data/trees/{experiment_type}/{prefix}.audacity_only.v8_masked.tree")
  aln_path <- str_glue("data/alignments/{experiment_type}/{prefix}.audacity_only.v8_masked.aln.fasta")
  
  tree <- read.tree(tree_path)
  
  resdir <- paste0(getwd(), str_glue("/results/mini_animal_trees/homoplasy_out/{prefix}/"))
  dir.create(resdir)
  
  # Run the HomoplasyFinder jar tool
  inconsistentPositions <- runHomoplasyFinderInJava(treeFile = tree_path, 
                                                    fastaFile = aln_path, 
                                                    path = resdir)
}
# Read in the output table
# results <- read.table("results/homoplasy_out/consistencyIndexReport_21-11-21.txt", 
#                       header=TRUE, sep="\t", stringsAsFactors=FALSE)
# 
# # Read in the annotated tree
# tree <- readAnnotatedTree(resdir)
# 
# # Plot the annotated tree
# plotAnnotatedTree(tree, inconsistentPositions, aln_path)
