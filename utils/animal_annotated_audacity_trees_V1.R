rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(ape)

audacity_tree <- read.tree("data/GISAID-hCoV-19-phylogeny-2022-02-21/global.tree")

cluster_meta <- fread("results/cluster_annotation/deer_mink_raw_clusters.csv")

mink <- fread("data/metadata/human_animal_subsets/V1/all_mink.n22381.csv") %>%
  left_join(cluster_meta)

deer <- fread("data/metadata/human_animal_subsets/V1/all_deer.n7325.csv") %>%
  left_join(cluster_meta)

host_list <- list(mink = mink, deer = deer)

for (i in seq(length(host_list))) {
  host_name <- names(host_list)[i]
  host_meta <- host_list[[i]]
  
  print(host_name)
  
  # Subset tree
  to_drop <- audacity_tree$tip.label[!(audacity_tree$tip.label %in% host_meta$accession_id)]
  audacity_filt <- drop.tip(audacity_tree, to_drop)
  
  # Match metadata to tips
  meta.match <- host_meta[match(audacity_filt$tip.label, host_meta$accession_id), ]
  
  # Annotate global tree tips
  annotated_tree <- audacity_filt
  tip_annot <- annotated_tree$tip.label
  all(tip_annot == meta.match$accession_id)
  tip_annot <- paste0(tip_annot,"|", meta.match$cluster) 
  annotated_tree$tip.label <- tip_annot
  write.tree(annotated_tree, str_glue("results/cluster_annotation/{host_name}.V1.tree"))
}
