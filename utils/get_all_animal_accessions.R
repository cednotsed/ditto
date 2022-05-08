rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(ape)
require(data.table)
require(tidyverse)
require(foreach)
require(doParallel)
require(lubridate)
registerDoParallel(cores=8)


# Load metadata
all_meta <- fread("data/metadata/all_sequence_metadata_260322.audacity.parsed.tsv") %>%
  mutate(collection_date = as.Date(collection_date))

# Run for mink and deer
host_list <- list(mink = "Neovison vison", deer = "Odocoileus virginianus")

for (i in seq(length(host_list))) {
  host_name <- names(host_list)[i]
  host_species <- host_list[[i]]
  
  ref_seq <- all_meta %>% filter(accession_id == "EPI_ISL_402124")
  animal_meta <- all_meta %>%
    filter(host == host_species)
  
  animal_meta <- ref_seq %>%
    bind_rows(animal_meta)
  
  accessions <- animal_meta %>% select(accession_id)
  
  dates <- animal_meta %>%
    rename(accession = accession_id) %>%
    select(accession, collection_date) %>%
    mutate(collection_date = decimal_date(as.Date(collection_date)))
  
  fwrite(animal_meta,
         str_glue("data/metadata/all_animals/{host_name}_only.n{nrow(accessions)}.csv"))
  
  fwrite(dates,
         str_glue("data/metadata/all_animals/{host_name}_only.n{nrow(accessions)}.dates_only.csv"))
  
  fwrite(accessions,
         str_glue("data/metadata/all_animals/{host_name}_only.n{nrow(accessions)}.accessions_only.csv"))
}
