rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(ape)

audacity_tree <- read.tree("data/GISAID-hCoV-19-phylogeny-2022-02-21/global.tree")

mink <- fread("data/metadata/all_animals/mink_only.n929.csv")
deer <- fread("data/metadata/all_animals/deer_only.n96.csv")

all_meta <- fread("data/metadata/all_sequence_metadata_260322.audacity.parsed.tsv")
cluster_meta <- fread("results/cluster_annotation/deer_mink_raw_clusters.csv")

meta_filt <- all_meta %>%
  mutate(collection_date = as.Date(collection_date)) %>%
  mutate(collection_date = if_else(location == "Europe / Lithuania", 
                                   as.Date("2020-11-26"), 
                                   collection_date)) %>%
  left_join(cluster_meta) %>%
  filter(host %in% c("Neovison vison", "Odocoileus virginianus", "Human"))

date_stats <- meta_filt %>%
  group_by(host, location) %>%
  summarise(max = max(collection_date, na.rm = T))

ref <- meta_filt %>%
  filter(accession_id == "EPI_ISL_402124")

for (host_species in c("Neovison vison", "Odocoileus virginianus")) {
  date_stats_filt <-  date_stats %>%
    filter(host == host_species)
  
  animal_meta <- meta_filt %>%
    filter(host == host_species)
  
  morsels <- foreach (loc = date_stats_filt$location) %do% {
    max_date <- (date_stats_filt %>% filter(location == loc))$max
    temp_human <- meta_filt %>%
      filter(host == "Human") %>%
      filter(collection_date <= max_date) %>%
      filter(location == loc)
    
    return(temp_human)
  } 
  
  human_meta <- bind_rows(morsels)

  human_animal <- ref %>%
    bind_rows(animal_meta, human_meta)
  
  # Subset tree
  to_drop <- audacity_tree$tip.label[!(audacity_tree$tip.label %in% human_animal$accession_id)]
  audacity_filt <- drop.tip(audacity_tree, to_drop)
  
  # Match metadata to tips
  meta.match <- human_animal[match(audacity_filt$tip.label, human_animal$accession_id), ]
  
  host_name <- gsub(" ", "_", host_species)
  fwrite(meta.match, str_glue("results/cluster_annotation/{host_name}.csv"))
  
  # Annotate global tree tips
  annotated_tree <- audacity_filt
  tip_annot <- annotated_tree$tip.label
  all(tip_annot == meta.match$accession_id)
  tip_annot <- paste0(tip_annot,"|", meta.match$cluster) 
  annotated_tree$tip.label <- tip_annot
  write.tree(annotated_tree, str_glue("results/cluster_annotation/{host_name}.tree"))
}
