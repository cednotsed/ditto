rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(data.table)
require(tidyverse)
require(ape)
require(foreach)

# Load metadata
audacity <- fread("data/GISAID-hCoV-19-phylogeny-2021-11-16/metadata.csv")

all_meta <- fread("data/metadata/all_sequence_metadata_231121.tsv") %>% 
  rename_all(~ tolower(gsub(" ", "_", .))) %>%
  filter(accession_id %in% audacity$accession_id) %>%
  separate(location, into = c("loc1", "loc2", "loc3"), sep = " / ") %>%
  mutate(location = paste0(loc1, " / ", loc2))

mink_deer <- all_meta %>%
  filter(host %in% c("Neovison vison", "Odocoileus virginianus")) %>%
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

## After manual inspection ##
deer_cluster <- fread("results/cluster_annotation/deer_clusters.csv") %>%
  mutate(cluster2 = paste0("deer_USA_", cluster)) %>%
  separate(accession_id, into = c("accession_id", NA), sep = "\\|") %>%
  select(-cluster)

parsed_meta <- mink_deer %>%
  left_join(merged, by = c("loc2", "host", "pango_lineage")) %>%
  left_join(deer_cluster) %>%
  # For mink
  mutate(cluster = case_when(cluster %in% c("mink_Denmark_8", "mink_Denmark_6") ~ "mink_Denmark_1",
                             cluster == "mink_Denmark_7" ~ "mink_Denmark_6",
                             accession_id %in% c("EPI_ISL_447629", "EPI_ISL_447630", "EPI_ISL_523075",
                                                 "EPI_ISL_447632", "EPI_ISL_447628", "EPI_ISL_523085",
                                                 "EPI_ISL_523070") ~ "mink_Netherlands_7",
                             accession_id %in% c("EPI_ISL_984307", "EPI_ISL_984305") ~ "mink_Poland_1",
                             cluster %in% c("mink_USA_3", "mink_USA_4") ~ "mink_USA_1",
                             accession_id %in% c("EPI_ISL_2896210", "EPI_ISL_2896209", "EPI_ISL_2834937",
                                                 "EPI_ISL_2834938", "EPI_ISL_2896205", "EPI_ISL_2896204",
                                                 "EPI_ISL_2896206", "EPI_ISL_2896207", "EPI_ISL_2896208",
                                                 "EPI_ISL_2896203") ~ "mink_Netherlands_7",
                             TRUE ~ as.character(cluster))) %>%
  mutate(cluster = ifelse(accession_id %in% c("EPI_ISL_683023", "EPI_ISL_683024", "EPI_ISL_683025",
                                              "EPI_ISL_641422", "EPI_ISL_641423", "EPI_ISL_683005",
                                              "EPI_ISL_641421", "EPI_ISL_683010", "EPI_ISL_683009"), 
                          "mink_Denmark_7", cluster)) %>%
  # For deer
  mutate(cluster = ifelse(accession_id %in% deer_cluster$accession_id, cluster2, cluster)) %>%
  select(accession_id, pango_lineage, host, cluster)

fwrite(parsed_meta, "results/cluster_annotation/deer_mink_parsed_clusters.csv")



