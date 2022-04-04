rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(ape)
require(data.table)
require(tidyverse)
require(foreach)
require(doParallel)
require(lubridate)
registerDoParallel(cores=8)

host_species <- "Neovison vison"
host_name <- "mink"

# Load metadata
all_meta <- fread("data/metadata/all_sequence_metadata_260322.audacity.parsed.tsv")
cluster_meta <- fread("results/cluster_annotation/deer_mink_parsed_clusters.csv") %>%
  select(accession_id, cluster)

all_meta_cluster <- all_meta %>%
  left_join(cluster_meta)

animal_meta <- all_meta_cluster %>%
  filter(host == host_species, 
         loc2 == "Denmark")

human_meta <- all_meta_cluster %>%
  filter(host == "Human",
         loc2 == "Denmark",
         collection_date <= ymd("2020-12-01"))

human_animal <- human_meta %>%
  bind_rows(animal_meta)

accessions <- human_animal %>% select(accession_id)

fwrite(human_animal,
       str_glue("data/metadata/human_animal_subsets/V6/denmark_minks.n{nrow(human_animal)}.csv"))

# fwrite(dates,
#        str_glue("data/metadata/human_animal_subsets/V5/all_{host_name}.n{nrow(accessions)}.dates_only.csv"))

fwrite(accessions,
       str_glue("data/metadata/human_animal_subsets/V6/denmark_minks.n{nrow(accessions)}.accessions_only.csv"))


# all_meta_cluster %>%
#   filter(loc2 == "Denmark") %>%
#   filter(pango_lineage == "B.1.1.298") %>%
#   # summarise(max = max(collection_date)) %>%
#   select(aa_substitutions) %>%
#    View()
# 
# all_meta_cluster %>%
#   # filter(pango_lineage %in% c("B.1.1.298", "B.1.536", "B.1")) %>%
#   filter(cluster == "mink_Denmark_1") %>%
#   summarise(max = max(collection_date, na.rm = T))
#   distinct(pango_lineage)
# animal_meta$cluster
