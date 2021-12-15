rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(data.table)
require(tidyverse)
require(foreach)
require(doParallel)
registerDoParallel(cores=8)
# 
# file <- "mink_n789_201121.tsv"
# host_name <- "mink"
file <- "deer_n76_201121.tsv"
host_name <- "deer"

set.seed(66)

# Load metadata
audacity <- fread("data/GISAID-hCoV-19-phylogeny-2021-11-16/metadata.csv")

all_meta <- fread("data/metadata/all_sequence_metadata_231121.tsv") %>% 
  rename_all(~ tolower(gsub(" ", "_", .))) %>%
  filter(accession_id %in% audacity$accession_id) %>%
  separate(location, into = c("loc1", "loc2", "loc3"), sep = " / ") %>%
  mutate(location = paste0(loc1, " / ", loc2))

all_animal_meta <- fread("data/metadata/combined_metadata_201121.audacity_only.tsv") %>%
  filter(host == host_species) %>%
  separate(location, into = c("loc1", "loc2", "loc3"), sep = " / ") %>%
  mutate(location = paste0(loc1, " / ", loc2))

animal_meta <- fread(str_glue("data/metadata/all_animals/{file}")) %>%
  rename_all(~ tolower(gsub(" ", "_", .))) %>%
  mutate(collection_date = as.Date(collection_date, "%d/%m/%Y")) %>%
  filter(!is.na(collection_date), 
         accession_id %in% audacity$accession_id) %>%
  separate(location, into = c("loc1", "loc2", "loc3"), sep = " / ") %>%
  mutate(location = paste0(loc1, " / ", loc2))

loc_meta <- animal_meta %>%
  group_by(location) %>%
  summarise(min_date = min(collection_date, na.rm = T),
            max_date = max(collection_date, na.rm = T)) %>%
  mutate(span = max_date - min_date) %>%
  filter(span > 0)

locations <- loc_meta$location

lin_meta <- animal_meta %>%
  filter(location %in% locations) %>%
  group_by(location, lineage) %>%
  summarise(n_isolates = n())

date_meta <- lin_meta %>%
  left_join(loc_meta)

locations <- unique(date_meta$location)

humans_in_countries <- all_meta %>%
  filter(host == "Human",
         location %in% locations) %>%
  mutate(collection_date = as.Date(collection_date)) %>%
  filter(!is.na(collection_date))

human_samples <- foreach(i = seq(nrow(date_meta)), .packages = "tidyverse") %dopar% {
  row <- date_meta[i, ]
  human_filt <- humans_in_countries %>%
    filter(host == "Human") %>% 
    filter(location == row$location) %>% 
    filter(pango_lineage == row$lineage) %>% 
    filter(collection_date >= row$min_date - 30) %>%
    filter(collection_date <= row$max_date + 30)
  
  if (nrow(human_filt) < row$n_isolates) {
    human_filt
  } else {
    human_filt %>% sample_n(row$n_isolates)
  }
  
}

human_df <- bind_rows(human_samples)

# Check if lineage and location matching is successful
generated_df <- human_df %>%
  select(pango_lineage, location) %>%
  rename(lineage = pango_lineage) %>%
  group_by(location, lineage) %>%
  summarise(n_isolates = n())

generated_df <- tibble(location = names(generated_df), sampled = as.numeric(generated_df))

lin_meta %>%
  rename(animal_isolates = n_isolates) %>%
  left_join(generated_df)

# Save background
n_animals <- nrow(animal_meta)
n_humans <- nrow(human_df)
fwrite(animal_meta %>% select(accession_id), 
       str_glue("data/genomes/matched/{host_name}_{n_animals}_201121.accessions.txt"))
fwrite(human_df %>% select(accession_id), 
       str_glue("data/genomes/matched/{host_name}_matched_humans_{n_humans}_201121.accessions.txt"))
