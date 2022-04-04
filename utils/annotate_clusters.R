rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(data.table)
require(tidyverse)

# Load metadata
mink <- fread("data/metadata/all_animals/mink_only.n929.csv")
deer <- fread("data/metadata/all_animals/deer_only.n96.csv")
mink_deer <- bind_rows(mink, deer) %>%
  filter(host != "Human") %>%
  mutate(host_name = ifelse(host == "Neovison vison", "mink", "deer"))

## After manual inspection ##
raw_clusters <- fread("results/cluster_annotation/deer_mink_raw_clusters.csv") %>%
  select(-n, -cluster_num, -host_name)

# Utah clusters
utah_a <- as.vector(t(read.table("results/cluster_annotation/mink_utah_A.txt")))
utah_a <- str_split(utah_a, "\\|", simplify = T)[, 1]
utah_b <- as.vector(t(read.table("results/cluster_annotation/mink_utah_B.txt")))
utah_b <- str_split(utah_b, "\\|", simplify = T)[, 1]

parsed_meta <- mink_deer %>%
  left_join(raw_clusters, by = c("loc2", "host", "pango_lineage")) %>%
  rename(raw_clusters = cluster) %>%
  ## For mink ##
  mutate(cluster = case_when(
                             # Denmark 
                             raw_clusters %in% c("mink_Denmark_5") ~ "mink_Denmark_1",
                             raw_clusters %in% c("mink_Denmark_7") ~ "mink_Denmark_5",
                             # Poland
                             raw_clusters %in% c("mink_Poland_2", "mink_Poland_5") ~ "mink_Poland_1",
                             raw_clusters %in% c("mink_Poland_3", "mink_Poland_4") ~ "mink_Poland_2",
                             # USA
                             accession_id %in% utah_a ~ "mink_USA_1",
                             accession_id %in% utah_b ~ "mink_USA_2",
                             host_name == "mink" & loc3 == "Wisconsin" ~ "mink_USA_3",
                             host_name == "mink" & loc3 == "Michigan" ~ "mink_USA_4",
                             host_name == "mink" & loc3 == "Wisconsin" ~ "mink_USA_5",
                             host_name == "mink" & loc3 == "Oregon" ~ "mink_USA_6",
                             # France
                             raw_clusters %in% c("mink_France_2") ~ "mink_France_1",
                             # Netherlands
                             raw_clusters %in% c("mink_Netherlands_2", "mink_Netherlands_5") ~ "mink_Netherlands_2",
                             raw_clusters == "mink_Netherlands_4" ~ "mink_Netherlands_3",
                             raw_clusters == "mink_Netherlands_6" ~ "mink_Netherlands_4",
                             TRUE ~ as.character(raw_clusters))) %>%
  # Others
  mutate(cluster = ifelse(accession_id %in% c("EPI_ISL_641423", "EPI_ISL_683009", "EPI_ISL_683010"), "mink_Denmark_5", cluster),
         cluster = ifelse(raw_clusters %in% c("mink_Netherlands_3", "mink_Belarus_1", "mink_Denmark_6"), "Unassigned", cluster),
         cluster = ifelse(accession_id %in% c("EPI_ISL_3218557"), "Unassigned", cluster),
         cluster = ifelse(accession_id %in% c("EPI_ISL_626344", "EPI_ISL_626340", "EPI_ISL_626341", 
                                              "EPI_ISL_626348", "EPI_ISL_626346", "EPI_ISL_626345", 
                                              "EPI_ISL_626347", "EPI_ISL_626350", "EPI_ISL_626342", 
                                              "EPI_ISL_626351", "EPI_ISL_626349", "EPI_ISL_626343"), "Unassigned", cluster)) %>%
  ## For Deer ##
  mutate(cluster = case_when(
         host_name == "mink" ~ cluster,
         raw_clusters %in% c("deer_USA_21") ~ "deer_USA_1",
         raw_clusters %in% c("deer_USA_17", "deer_USA_24") ~ "deer_USA_2",
         accession_id == "EPI_ISL_5804748" ~ "deer_USA_2",
         raw_clusters %in% c("deer_USA_20") ~ "deer_USA_3",
         raw_clusters %in% c("deer_USA_18") ~ "deer_USA_4",
         raw_clusters %in% c("deer_USA_22", "deer_USA_23") ~ "deer_USA_7",
         grepl("_14|_10|_12|_11|_7|_9|_13|_15|_8", raw_clusters) ~ "deer_USA_8",
         TRUE ~ as.character(NA))) %>%
  mutate(cluster = ifelse(accession_id %in% c("EPI_ISL_5804740", "EPI_ISL_5804741", "EPI_ISL_5804737",
                                              "EPI_ISL_5804738", "EPI_ISL_5804739"), "deer_USA_5", cluster),
         cluster = ifelse(accession_id %in% c("EPI_ISL_5804735", "EPI_ISL_5804736", "EPI_ISL_5804730",
                                              "EPI_ISL_5804729"), "Unassigned", cluster)) %>%
  select(accession_id, pango_lineage, host, raw_clusters, cluster)

fwrite(parsed_meta, "results/cluster_annotation/deer_mink_parsed_clusters.csv")

# mink_deer %>%
#   left_join(raw_clusters, by = c("loc2", "host", "pango_lineage")) %>%
#   rename(raw_clusters = cluster)
# 
  parsed_meta %>%
  filter(accession_id == "EPI_ISL_3218557")
