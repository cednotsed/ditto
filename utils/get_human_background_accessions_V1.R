rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(ape)
require(data.table)
require(tidyverse)
require(foreach)
require(doParallel)
require(lubridate)
# cl <- makeCluster(4)
# registerDoParallel(cl)

# Load metadata
all_meta <- fread("data/metadata/all_sequence_metadata_260322.audacity.parsed.tsv") %>%
  mutate(collection_date = as.Date(collection_date))

host_species <- "Neovison vison"
host_name <- "mink"

# host_species <- "Odocoileus virginianus"
# host_name <- "deer"

n_samples <- 10

animal_meta <- all_meta %>%
  filter(host == host_species)

human_meta <- all_meta %>%
  filter(host == "Human")

# Print no. of animal isolates per country
country_stats <- animal_meta %>%
  mutate(collection_date = if_else(location == "Europe / Lithuania", 
                                   ymd("2020-11-26"), 
                                   collection_date)) %>%
  group_by(location) %>% 
  summarise(n_total = n())

# Get no. of isolates per lineage per country 
lineage_stats <- human_meta %>%
  filter(location %in% country_stats$location) %>%
  group_by(location, pango_lineage) %>%
  summarise(n_total = n())

animal_meta_filt <- animal_meta %>%
  filter(location %in% unique(country_stats$location))

human_meta_filt <- human_meta %>%
  filter(location %in% unique(country_stats$location))

# rm(list = c("human_meta", "animal_meta", "all_meta"))

set.seed(66)

morsels <- foreach(i = seq(nrow(lineage_stats)), .packages = "tidyverse") %do% {
  row <- lineage_stats[i, ]
  
  # Get human isolates with matching pango lineage and country
  human_filt <- human_meta_filt %>% 
    filter(location == row$location, pango_lineage == row$pango_lineage)
  
  # Get sample of 10 isolates
  if (nrow(human_filt) >= n_samples) {
    return(human_filt %>% sample_n(n_samples))
  } else {
    return(human_filt)
  }

}

human_animal <- animal_meta_filt %>%
  bind_rows(morsels)

# Check
test_human <- human_animal %>% filter(host == "Human")
test_human %>% 
    group_by(location, pango_lineage) %>%
    summarise(n_total = n()) %>%
    arrange(desc(n_total))

accessions <- human_animal %>% select(accession_id)

dates <- human_animal %>%
  rename(accession = accession_id) %>%
  select(accession, collection_date) %>%
  mutate(collection_date = decimal_date(collection_date))

fwrite(human_animal,
       str_glue("data/metadata/human_animal_subsets/V1/all_{host_name}.n{nrow(accessions)}.csv"))

fwrite(dates,
       str_glue("data/metadata/human_animal_subsets/V1/all_{host_name}.n{nrow(accessions)}.dates_only.csv"))

fwrite(accessions,
       str_glue("data/metadata/human_animal_subsets/V1/all_{host_name}.n{nrow(accessions)}.accessions_only.csv"))

# stopCluster(cl)




