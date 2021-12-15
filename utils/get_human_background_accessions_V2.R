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
  mutate(location = paste0(loc1, " / ", loc2))

animal_meta <- fread("data/metadata/combined_metadata_201121.audacity_only.csv") %>%
  select(accession_id, cluster)

meta <- all_meta %>%
  left_join(animal_meta) %>%
  select(accession_id, virus_name, location, 
         collection_date, host, pango_lineage,
         clade, cluster, variant, aa_substitutions)

# Print no. of isolates per country
animal_stats <- meta %>%
  filter(host == host_species) %>%
  group_by(location) %>%
  summarise(n = n())

## Allele frequency background ##
countries <- animal_stats %>% filter(n > 50)
countries <- countries$location

set.seed(66)
# country <- countries[1]
for (country in countries) {
  # Get animal accessions in country of preference
  animal <- meta %>%
    filter(host == host_species,
           location == country) 
  
  # Get human isolates that are in same country
  human <- meta %>%
    filter(host == "Human", location == country) %>%
    filter(collection_date >= min(animal$collection_date) & 
             collection_date <= max(animal$collection_date))
  
  human_animal <- human %>% 
    bind_rows(animal)
  
  accessions <- human_animal %>% select(accession_id)
  dates <- human_animal %>%
    rename(accession = accession_id) %>%
    select(accession, collection_date) %>%
    mutate(collection_date = decimal_date(as.Date(collection_date)))
  
  country_name <- strsplit(country, " / ")[[1]][2]
  
  fwrite(human_animal,
         str_glue("data/metadata/human_animal_subsets/{host_name}.{country_name}.n{nrow(accessions)}.csv"))

  fwrite(dates,
         str_glue("data/metadata/human_animal_subsets/{host_name}.{country_name}.n{nrow(accessions)}.dates_only.csv"))

  fwrite(accessions,
         str_glue("data/metadata/human_animal_subsets/{host_name}.{country_name}.n{nrow(accessions)}.accessions_only.csv"))
}



