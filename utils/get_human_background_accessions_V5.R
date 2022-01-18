rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(ape)
require(data.table)
require(tidyverse)
require(foreach)
require(doParallel)
require(lubridate)
registerDoParallel(cores=8)

# host_species <- "Neovison vison"
# host_name <- "mink"

host_species <- "Odocoileus virginianus"
host_name <- "deer"

# Load metadata
audacity <- fread("data/GISAID-hCoV-19-phylogeny-2021-11-16/metadata.csv")

cluster_meta <- fread("data/metadata/combined_metadata_201121.audacity_only.csv") %>%
  select(accession_id, cluster)

all_meta <- fread("data/metadata/all_sequence_metadata_231121.tsv") %>% 
  rename_all(~ tolower(gsub(" ", "_", .))) %>%
  filter(accession_id %in% audacity$accession_id) %>%
  separate(location, into = c("loc1", "loc2", "loc3"), sep = " / ") %>%
  mutate(location = paste0(loc1, " / ", loc2),
         collection_date = as.Date(collection_date)) %>%
  left_join(cluster_meta) %>%
  as_tibble()

animal_meta <- all_meta %>%
  filter(host == host_species)

human_meta <- all_meta %>%
  filter(host == "Human")

# Print no. of isolates per country
country_stats <- animal_meta %>%
  mutate(collection_date = if_else(location == "Europe / Lithuania", 
                            ymd("2020-11-26"), 
                            collection_date)) %>%
  group_by(location) %>% 
  
  summarise(n_total = n(), 
            min_date = min(collection_date, na.rm = T),
            max_date = max(collection_date, na.rm = T))

animal_stats <- animal_meta %>%
  filter(host == host_species) %>%
  group_by(location, pango_lineage) %>%
  summarise(n = n()) %>%
  left_join(country_stats) %>%
  filter(n_total > 10)

animal_meta_filt <- animal_meta %>%
  filter(location %in% unique(animal_stats$location))

set.seed(66)

morsels <- foreach(i = seq(nrow(animal_stats)), .packages = "tidyverse") %do% {
  print(i)
  row <- animal_stats[i, ]
  
  animal_filt <- all_meta %>% 
    filter(host == host_species) %>%
    filter(location == row$location,
           pango_lineage == row$pango_lineage)
 
   # Get human isolates that are in same country, lineage and time span
  human_filt <- human_meta %>% 
    filter(location == row$location,
           collection_date >= row$min_date - 30 &
             collection_date <= row$max_date + 30,
           pango_lineage == row$pango_lineage) 
  
  if (nrow(human_filt) >= nrow(animal_filt) * 10) {
    human_filt <- human_filt %>% 
      sample_n(row$n * 10)
  }
  
  human_filt
}

human_animal <- animal_meta_filt %>%
  bind_rows(morsels)

# Check
test_human <- human_animal %>% filter(host == "Human")
test_animal <- human_animal %>% filter(host == "Neovison vison")

accessions <- human_animal %>% select(accession_id)

dates <- human_animal %>%
  rename(accession = accession_id) %>%
  select(accession, collection_date) %>%
  mutate(collection_date = decimal_date(collection_date))

fwrite(human_animal,
       str_glue("data/metadata/human_animal_subsets/V5/all_{host_name}.n{nrow(accessions)}.csv"))

fwrite(dates,
       str_glue("data/metadata/human_animal_subsets/V5/all_{host_name}.n{nrow(accessions)}.dates_only.csv"))

fwrite(accessions,
       str_glue("data/metadata/human_animal_subsets/V5/all_{host_name}.n{nrow(accessions)}.accessions_only.csv"))


