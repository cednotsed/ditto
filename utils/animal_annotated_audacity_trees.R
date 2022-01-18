rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(ape)

audacity_tree <- read.tree("data/GISAID-hCoV-19-phylogeny-2021-11-16/global.tree")

all_meta <- fread("data/metadata/all_sequence_metadata_231121.tsv") %>% 
  rename_all(~ tolower(gsub(" ", "_", .))) %>%
  filter(accession_id %in% audacity_tree$tip.label) %>%
  separate(location, into = c("loc1", "loc2", "loc3"), sep = " / ") %>%
  mutate(location = paste0(loc1, " / ", loc2))

cluster_meta <- fread("results/cluster_annotation/deer_mink_raw_clusters.csv")

meta_filt <- all_meta %>%
  left_join(cluster_meta) %>%
  filter(host %in% c("Neovison vison", "Odocoileus virginianus", "Human"))

date_stats <- meta_filt %>%
  group_by(host) %>%
  summarise(max = max(collection_date, na.rm = T))

for (host_species in c("Neovison vison", "Odocoileus virginianus")) {
  max_date <- date_stats[date_stats$host == host_species, ]$max
  
  animal_meta <- meta_filt %>%
    filter(host == host_species)
  
  locations <- unique(animal_meta$loc2)
  
  human_animal <- meta_filt %>%
    filter(host %in% c("Human", host_species)) %>%
    filter(collection_date <= max_date) %>%
    filter(loc2 %in% locations)
  
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
