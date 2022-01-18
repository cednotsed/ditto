rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(ape)
require(data.table)
require(tidyverse)
require(foreach)
require(doParallel)
registerDoParallel(cores=8)

host_species <- "Neovison vison"
host_name <- "mink"

# host_species <- "Odocoileus virginianus"
# host_name <- "deer"

# Load metadata
audacity <- fread("data/GISAID-hCoV-19-phylogeny-2021-11-16/metadata.csv")

all_meta <- fread("data/metadata/all_sequence_metadata_231121.tsv") %>% 
  rename_all(~ tolower(gsub(" ", "_", .))) %>%
  filter(accession_id %in% audacity$accession_id) %>%
  separate(location, into = c("loc1", "loc2", "loc3"), sep = " / ") %>%
  mutate(location = paste0(loc1, " / ", loc2)) %>%
  as_tibble()

# Print no. of isolates per country
animal_stats <- all_meta %>%
  filter(host == host_species) %>%
  group_by(location) %>%
  summarise(n = n())

## Allele frequency background ##
countries <- animal_stats %>% filter(n > 10)
countries <- countries$location

set.seed(66)

for (country in countries) {
  # Get animal accessions in country of preference
  animal_filt <- all_meta %>%
    filter(host == host_species) %>%
    filter(location == country)
  
  # Get human isolates that are in same country
  in_same_country <- all_meta %>% 
    filter(host == "Human") %>%
    filter(grepl(country, location))
  
  # Subsample accessions while ensuring all lineages are represented
  lineages <- unique(in_same_country$pango_lineage) 
  n_samples <- 10
  length(lineages)
  
  morsels <- foreach(lineage = lineages, .packages = "tidyverse") %do% {
    lineage_df <- in_same_country %>% 
      filter(pango_lineage == lineage)
    
    n_accessions <- nrow(lineage_df)
    
    if (n_accessions >= n_samples) {
      return(lineage_df %>% sample_n(n_samples))
    } else {
      return(lineage_df %>% sample_n(n_accessions))
    }
  }
  
  human_filt <- bind_rows(morsels)
  human_animal <- human_filt %>% 
    bind_rows(animal_filt)
  
  accessions <- human_animal %>% select(accession_id)
  
  dates <- human_animal %>%
    rename(accession = accession_id) %>%
    select(accession, collection_date) %>%
    mutate(collection_date = decimal_date(as.Date(collection_date)))
  
  country_name <- strsplit(country, " / ")[[1]][2]
  
  fwrite(human_animal,
         str_glue("data/metadata/human_animal_subsets/V1/{host_name}.{country_name}.n{nrow(accessions)}.csv"))
  
  fwrite(dates,
         str_glue("data/metadata/human_animal_subsets/V1/{host_name}.{country_name}.n{nrow(accessions)}.dates_only.csv"))
  
  fwrite(accessions,
         str_glue("data/metadata/human_animal_subsets/V1/{host_name}.{country_name}.n{nrow(accessions)}.accessions_only.csv"))
  
}



