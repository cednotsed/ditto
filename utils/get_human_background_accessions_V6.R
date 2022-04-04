rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(ape)
require(data.table)
require(tidyverse)
require(foreach)
require(doParallel)
require(lubridate)
registerDoParallel(cores=8)

hosts <- c("Neovison vison", "Odocoileus virginianus")

# Load metadata
all_meta <- fread("data/metadata/all_sequence_metadata_260322.audacity.parsed.tsv") %>%
  mutate(collection_date = as.Date(collection_date))

animal_filt <- all_meta %>% 
  filter(host %in% hosts)

animal_all <- all_meta %>% 
  filter(!(host %in% c("Human", "Environment"))) %>%
  mutate(host = ifelse(grepl("Panthera", host), "Panthera spp.", host))

animal_stats_filt <- animal_all %>%
  group_by(host) %>%
  summarise(n = n()) %>%
  filter(n > 10)

animal_stats

animal_all <- animal_all %>%
  filter(host %in% animal_stats_filt$host)

# Print no. of isolates per country
country_stats <- animal_filt %>%
  mutate(collection_date = if_else(location == "Europe / Lithuania", 
                                   ymd("2020-11-26"), 
                                   collection_date)) %>%
  group_by(location) %>% 
  summarise(n_total = n(), 
            min_date = min(collection_date, na.rm = T),
            max_date = max(collection_date, na.rm = T))

set.seed(66)
  
# Get 10 human isolates per lineage
n_samples <- 10

human_meta <- all_meta %>%
  filter(host == "Human")

lineages <- unique(human_meta$pango_lineage)

morsels <- foreach(lineage = lineages, .packages = "tidyverse") %do% {
  lineage_df <- human_meta %>% 
    filter(pango_lineage == lineage)
  
  n_accessions <- nrow(lineage_df)
  
  if (n_accessions >= n_samples) {
    return(lineage_df %>% sample_n(n_samples))
  } else {
    return(lineage_df %>% sample_n(n_accessions))
  }
}

human_filt <- bind_rows(morsels)

human_animal <- animal_filt %>%
  bind_rows(human_filt)

human_animal_all <- animal_all %>%
  bind_rows(human_filt)

accessions <- human_animal %>% select(accession_id)

dates <- human_animal %>%
  rename(accession = accession_id) %>%
  select(accession, collection_date) %>%
  mutate(collection_date = decimal_date(as.Date(collection_date)))

fwrite(human_animal,
       str_glue("data/metadata/human_animal_subsets/V6/mink_deer.n{nrow(accessions)}.csv"))
fwrite(human_animal_all,
       str_glue("data/metadata/human_animal_subsets/V6/all_animals.n{nrow(human_animal_all)}.csv"))

fwrite(dates,
       str_glue("data/metadata/human_animal_subsets/V6/mink_deer.n{nrow(accessions)}.dates_only.csv"))

fwrite(accessions,
       str_glue("data/metadata/human_animal_subsets/V6/mink_deer.n{nrow(accessions)}.accessions_only.csv"))


