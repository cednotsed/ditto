setwd("../Desktop/git_repos/ditto")
require(tidyverse)
require(data.table)
require(ape)

prefixes <- c("deer_n73_201121", "mink_n789_201121")
meta <- fread("data/metadata/combined_metadata_201121.audacity_only.csv")

for (prefix in prefixes) {
  animal_name <- str_split(prefix, "_")[[1]][1]
  fasta <- read.FASTA(str_glue("data/genomes/all_animals/{prefix}.audacity_only.fasta"))
  
  # Rename fasta
  # fasta_names <- separate(tibble(accessions = names(fasta)), accessions,
  #                         into = c(NA, "accessions", NA), 
  #                         sep = "\\|")$accessions
  # names(fasta) <- fasta_names
  
  meta_filt <- meta %>% filter(accession_id %in% names(fasta))
  to_keep <- table(meta_filt$cluster)
  to_keep <- to_keep[to_keep > 10]
  to_keep <- to_keep[names(to_keep) != ""]
  
  for (i in seq(nrow(to_keep))) {
    cluster <- to_keep[i]
    to_keep_names <- names(cluster)
    to_keep_freq <- as.numeric(cluster)
    
    # Get accessions in cluster
    meta_filt2 <- meta_filt %>% 
      filter(cluster == to_keep_names)
    accessions <- meta_filt2$accession_id
    
    # Get location
    loc <- separate(meta_filt2, location, into = c(NA, "location"), sep = " / ")
    location <- unique(loc$location)
    
    fasta_filt <- fasta[names(fasta) %in% accessions]
    write.FASTA(fasta_filt, 
                str_glue("data/genomes/mini_animal_trees/{animal_name}.{location}.{to_keep_names}.n{to_keep_freq}.audacity_only.fasta"))
  }
}
