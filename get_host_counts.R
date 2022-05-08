rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(ape)

# Anthroponosis
ant_df <- read.csv("data/metadata/netherlands_humans_anthroponoses.txt", header = F)$V1

# Animals only
fread("data/metadata/all_animals/mink_only.n929.csv") %>%
  group_by(host) %>%
  summarise(n = n())

# V1
fread("data/metadata/human_animal_subsets/V1/all_deer.n7325.csv") %>%
  group_by(host) %>%
  summarise(n = n())
fread("data/metadata/human_animal_subsets/V1/all_mink.n22381.csv") %>%
  group_by(host) %>%
  summarise(n = n())

# V5
fread("data/metadata/human_animal_subsets/V5/all_deer.n189.csv") %>%
  group_by(host) %>%
  summarise(n = n())
fread("data/metadata/human_animal_subsets/V5/all_mink.n1769.csv") %>%
  filter(!(accession_id %in% ant_df)) %>%
  group_by(host) %>%
  summarise(n = n())

# V6
fread("data/metadata/human_animal_subsets/V6/all_animals.n16911.csv") %>%
  group_by(host) %>%
  summarise(n = n()) %>%
  summarise(sum = sum(n))

fread("data/metadata/human_animal_subsets/V6/denmark_minks.n15895.csv") %>%
  group_by(host) %>%
  summarise(n = n())

# Cluster
fread("results/cluster_annotation/deer_mink_parsed_clusters.csv") %>%
  filter(cluster != "Unassigned") %>%
  group_by(host) %>%
  summarise(n = n_distinct(cluster))
fread("results/cluster_annotation/deer_mink_parsed_clusters.csv") %>%
  group_by(cluster) %>%
  summarise(n = n()) %>%
  View()

# Get alignment counts
aln <- read.FASTA("data/alignments/human_animal_subsets/V5/all_mink.n1769.audacity_only.v8_masked.fasta")
length(aln)

# Get counts of animal isolates
all_meta <- fread("data/metadata/all_sequence_metadata_260322.audacity.parsed.tsv")

all_meta %>%
  filter(!(host %in% c("Human", "Environment"))) %>%
  group_by(host) %>%
  summarise(n = n()) %>%
  View()

all_meta %>%
  filter(!(host %in% c("Human", "Environment")))

# Get counts of animal isolates > 10 in Audacity release
host_table <- all_meta %>%
  filter(!(host %in% c("Human", "Environment"))) %>%
  group_by(host) %>%
  summarise(n = n()) %>%
  filter(n > 10)

host_table

# Get no. of locations for animal isolates
all_meta %>%
  filter(!(host %in% c("Human", "Environment"))) %>%
  mutate(host = ifelse(grepl("Panthera", host), "Panthera spp.", host)) %>%
  filter(host %in% host_table$host) %>%
  group_by(host) %>%
  summarise(n = n_distinct(location))

# Count number of clusters
fread("results/cluster_annotation/deer_mink_parsed_clusters.csv") %>%
  filter(cluster != "Unassigned") %>%
  group_by(host) %>%
  summarise(n = n_distinct(cluster))

# Count no. of lineages
fread("data/metadata/all_sequence_metadata_260322.audacity.parsed.tsv") %>%
  filter(host == "Human") %>%
  summarise(n = n_distinct(pango_lineage))
fread("data/metadata/human_animal_subsets/V6/all_animals.n16911.csv") %>%
  filter(host != "Human") %>%
  summarise(n = n_distinct(pango_lineage), n_host = n_distinct(host))

fread("data/metadata/human_animal_subsets/V6/all_animals.n16911.csv") %>%
  filter(host != "Human") %>%
  distinct(variant)

# Count v6 hosts
fread("data/metadata/human_animal_subsets/V6/all_animals.n16911.csv") %>% nrow()
fread("data/metadata/human_animal_subsets/V6/all_animals.n16911.csv") %>%
  group_by(host) %>%
  summarise(n = n())

# Count v5 hosts
fread("data/metadata/human_animal_subsets/V5/all_mink.n1769.csv") %>% nrow()
fread("data/metadata/human_animal_subsets/V5/all_mink.n1769.csv") %>%
  group_by(host) %>%
  summarise(n = n())

# Count all_animal hosts
fread("data/metadata/all_animals/mink_only.n929.csv") %>% nrow()
fread("data/metadata/all_animals/mink_only.n929.csv") %>%
  group_by(host) %>%
  summarise(n = n())

# Count cluster annotation hosts
mink <- fread("results/cluster_annotation/Neovison_vison.csv")
mink %>% 
  filter(host == "Neovison vison") %>%
  group_by(location) %>%
  summarise(max = max(collection_date, na.rm = T))

mink %>% 
  filter(host == "Human") %>%
  group_by(location) %>%
  summarise(max = max(collection_date, na.rm = T))

table(mink$host)
nrow(mink)

deer <- fread("results/cluster_annotation/Odocoileus_virginianus.csv")
deer %>% 
  filter(host == "Odocoileus virginianus") %>%
  group_by(location) %>%
  summarise(max = max(collection_date, na.rm = T))

deer %>% 
  filter(host == "Human") %>%
  group_by(location) %>%
  summarise(max = max(collection_date, na.rm = T))

deer %>% 
  filter(host == "Odocoileus virginianus") %>%
  summarise(max = max(collection_date, na.rm = T))
table(deer$host)
nrow(deer)
