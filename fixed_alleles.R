setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(ape)
require(ggplot2)
require(foreach)

host_name <- "Neovison vison"
save_name <- "mink"
min_freq_change <- 0.1
aln <- read.dna(str_glue("data/alignments/human_animal_subsets/V1/all_mink.audacity_only.v8_masked.aln.fasta"),
                format = "fasta",
                as.matrix = T)

# Remove wuhan-hu-1
aln <- aln[rownames(aln) != "NC_045512.2", ]
  
# Load metadata
meta <- fread(str_glue("data/metadata/human_animal_subsets/V1/all_mink.tsv")) %>%
  rename_all(~ tolower(gsub(" ", "_", .))) %>%
  select(accession_id, host, clade, lineage, location, collection_date) %>%
  separate(location, into = c("loc1", "loc2", "loc3"), sep = " / ") %>%
  mutate(location = paste0(loc1, " / ", loc2, " / ", loc3)) %>%
  rename(Accession = accession_id) %>%
  as_tibble()
  
  animal <- aln[meta$host == host_name, ]
  human <- aln[meta$host != host_name, ]
  
  animal_morsels <- foreach(pos = seq(ncol(aln)), .packages = c("ape", "tidyverse")) %do% {
    base_freq <- base.freq(animal[, pos])
    
    tibble(pos = pos,
           base = names(base_freq),
           base_freq = as.numeric(base_freq))
  }
  
  human_morsels <- foreach(pos = seq(ncol(aln)), .packages = c("ape", "tidyverse")) %do% {
    base_freq <- base.freq(human[, pos])
    
    tibble(pos = pos,
           base = names(base_freq),
           base_freq = as.numeric(base_freq))
  }
  
  animal_df <- bind_rows(animal_morsels)
  human_df <- bind_rows(human_morsels)
  
  merged_df <- animal_df %>% 
    left_join(human_df, by = c("pos", "base"), suffix = c("_animal", "_human"))
  
  # Get sites variant between humans and animals
  unique_df <- merged_df %>%
    filter(base_freq_animal != base_freq_human)
  
  # Calculate change in frequency
  unique_df <- unique_df %>%
    mutate(abs_base_diff = abs(base_freq_animal - base_freq_human),
           base_diff = base_freq_animal - base_freq_human) %>%
    arrange(desc(abs_base_diff))
  
  unique_df
}

df <- bind_rows(morsels)
df
