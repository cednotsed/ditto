rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(ape)
require(ggtree)
require(homoplasyFinder)
files <- list.files("data/alignments/all_animals/", pattern = "v8_masked.fasta")
prefixes <- separate(tibble(x = files), x, into = c("prefix1", "prefix2"), sep = "\\.") %>%
  mutate(prefix = paste0(prefix1, ".", prefix2))
prefixes <- prefixes$prefix

# experiment_type <- "all_animals"
# prefixes <- list.files(str_glue("data/alignments/{experiment_type}"))
# prefixes <- prefixes[grepl(".aln.fasta", prefixes)]
# prefixes <- gsub(".audacity_only.v8_masked.fasta", "", prefixes)

for (i in seq(length(files))) {
  prefix <- prefixes[i]
  tree_path <- str_glue("data/trees/all_animals/{prefix}.audacity_only.v8_masked-delim.fasta.contree")
  aln_path <- str_glue("data/alignments/all_animals/{files[i]}")
  
  # Remove ref
  tree <- read.tree(tree_path)
  tree <- drop.tip(tree, "EPI_ISL_402124")
  no_ref_path <- gsub(".contree", ".no_ref.contree", tree_path)
  write.tree(tree, no_ref_path)
  
  aln <- read.FASTA(aln_path)
  aln <- aln[names(aln) != "EPI_ISL_402124"]
  no_ref_aln_path <- gsub(".fasta", ".no_ref.fasta", aln_path)
  write.FASTA(aln, no_ref_aln_path)
  
  resdir <- paste0(getwd(), str_glue("/results/homoplasy_out/{prefix}/"))
  dir.create(resdir)
  
  # Run the HomoplasyFinder jar tool
  inconsistentPositions <- runHomoplasyFinderInJava(treeFile = no_ref_path,
                                                    fastaFile = no_ref_aln_path,
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
