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
animal_stats <- all_animal_meta %>%
  group_by(location) %>%
  summarise(n = n())

## Allele frequency background ##
countries <- animal_stats %>% filter(n > 50)
countries <- countries$location

set.seed(66)

for (country in countries) {
  # Get animal accessions in country of preference
  animal_filt <- all_animal_meta %>%
    filter(location == country)
  
  # Get human isolates that are in same country
  in_same_country <- all_meta %>% 
    filter(accession_id %in% audacity$accession_id,
           host %in% c("Human")) %>%
    filter(grepl(country, location))
  
  # Subsample accessions while ensuring all lineages are represented
  lineages <- unique(in_same_country$pango_lineage) 
  n_samples <- 10
  length(lineages)
  
  morsels <- foreach(lineage = lineages, .packages = "tidyverse") %dopar% {
    lineage_df <- in_same_country %>% 
      filter(pango_lineage == lineage)
    
    n_accessions <- nrow(lineage_df)
    
    if (n_accessions > n_samples) {
      return(lineage_df %>% sample_n(n_samples))
    } else {
      return(lineage_df %>% sample_n(n_accessions))
    }
  }
  
  human_filt <- bind_rows(morsels)
  human_animal <- human_filt %>% 
    bind_rows(animal_filt)
  
  accessions <- human_animal %>% select(accession_id)
  country_name <- strsplit(country, " / ")[[1]][2]
  
  fwrite(accessions, 
         str_glue("data/metadata/human_animal_subsets/{host_name}.{country_name}.n{nrow(accessions)}.accessions_only.csv"))
}



