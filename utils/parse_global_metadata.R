rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(data.table)
require(tidyverse)
require(lubridate)

# Load metadata
audacity <- fread("data/GISAID-hCoV-19-phylogeny-2022-02-21/metadata.csv")

all_meta <- fread("data/metadata/all_sequence_metadata_260322.tsv") %>%
  rename_all(~ tolower(gsub(" ", "_", .))) %>%
  filter(accession_id %in% audacity$accession_id) %>%
  separate(location, into = c("loc1", "loc2", "loc3"), sep = " / ") %>%
  mutate(location = paste0(loc1, " / ", loc2),
         collection_date = as.Date(collection_date)) %>%
  select(-type, -gender, -`is_complete?`, -`gc-content`, 
         -`n-content`, -`is_high_coverage?`, -`is_reference?`, -pangolin_version,
         -patient_age, -sequence_length, -additional_location_information,
         -virus_name)

fwrite(all_meta, "data/metadata/all_sequence_metadata_260322.audacity.parsed.tsv")

all_meta
