setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(lubridate)
registerDoParallel(cores = 10)

# Load metadata
audacity <- fread("data/GISAID-hCoV-19-phylogeny-2021-11-16/metadata.csv")
meta <- fread("data/metadata/all_sequence_metadata_231121.tsv", nThread = 10) %>%
  rename_all(~ tolower(gsub(" ", "_", .))) %>%
  filter(accession_id %in% audacity$accession_id) %>%
  separate(location, into = c("loc1", "loc2", "loc3"), sep = " / ") %>%
  mutate(location = paste0(loc1, " / ", loc2))

ant <- read.csv("data/metadata/netherlands_anthroponoses.txt", header = F)
ant_df <- tibble(accession_id = ant$V1) %>%
  left_join(meta) %>%
  filter(host == "Human")

fwrite(ant_df %>% select(accession_id), 
       "data/metadata/netherlands_humans_anthroponoses.txt", 
       col.names = F)

aln_meta <- fread("data/metadata/human_animal_subsets/V1/mink.Netherlands.n2213.accessions_only.csv")
aln_meta %>% filter(accession_id %in% ant_df$accession_id)
