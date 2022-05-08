rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(data.table)
require(tidyverse)
require(ape)
require(foreach)

# Load metadata
mink <- fread("data/metadata/all_animals/mink_only.n929.csv")
deer <- fread("data/metadata/all_animals/deer_only.n96.csv")
mink_deer <- bind_rows(mink, deer) %>%
  filter(host != "Human") %>%
  mutate(host_name = ifelse(host == "Neovison vison", "mink", "deer"))

cluster_stats <- mink_deer %>%
  group_by(host, host_name, loc2, pango_lineage) %>%
  summarise(n = n())

locs <- unique(cluster_stats$loc2)

morsels <- foreach (loc = locs) %do% {
    cluster_temp <- cluster_stats %>% 
      filter(loc2 == loc)
    
    cluster_temp %>%
      add_column(cluster_num = seq(nrow(cluster_temp))) %>%
      mutate(cluster = str_glue("{host_name}_{loc2}_{cluster_num}"))
}  

merged <- bind_rows(morsels)
fwrite(merged, "results/cluster_annotation/deer_mink_raw_clusters.csv")



