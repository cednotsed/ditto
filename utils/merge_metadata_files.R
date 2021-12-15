setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(foreach)

# Merge metadata
metas <- foreach (file = list.files("data/metadata/all_animals/", full.names = T)) %do% {
  fread(file) %>% 
    rename_all(~ tolower(gsub(" ", "_", .))) %>%
    mutate(collection_date = as.character(collection_date))
}

meta <- bind_rows(metas)
fwrite(meta, "data/metadata/combined_metadata_201121.tsv")

# Remove entries not in audacity tree
audacity <- fread("data/GISAID-hCoV-19-phylogeny-2021-11-16/metadata.csv")
meta_filt <- meta %>% filter(accession_id %in% audacity$accession_id)
fwrite(meta_filt, "data/metadata/combined_metadata_201121.audacity_only.tsv")
