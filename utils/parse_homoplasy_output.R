setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(seqinr)
require(foreach)

result_dir <- "results/all_animals/homoplasy_out"
# result_dir <- "results/mini_animal_trees/homoplasy_out"
hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.csv") %>% as_tibble()
to_remove <- c("coronavirus frameshifting stimulation element stem-loop 1 (NSP12a)",
               "coronavirus frameshifting stimulation element stem-loop 2 (NSP12a)", 
               "3'UTR (pseudoknot stem-loop 1)", "3'UTR (pseudoknot stem-loop 2)",
               "3'UTR  (stem-loop II-like motif)")

hookup <- hookup %>% 
  filter(!(protein_name %in% to_remove))

hookup_unique <- hookup %>% 
  distinct(nucleotide_pos, .keep_all = T) %>%
  select(nucleotide_pos, ref_nuc) %>%
  mutate(ref_nuc = paste0("n", ref_nuc))

nuc_count_cols <- c("nA", "nC", "nG", "nT")
prefixes <- list.files(result_dir)
prefixes <- prefixes[!grepl("csv", prefixes)]

for (prefix in prefixes) {
  df <- fread(str_glue("{result_dir}/{prefix}/consistencyIndexReport_21-11-21.txt")) %>%
    as_tibble()
  
  if(nrow(df) > 0) {
    merged_df <- df %>% 
      separate(CountsACGT, into = nuc_count_cols) %>%
      mutate_at(nuc_count_cols, as.numeric) %>%
      rename(nucleotide_pos = Position) %>%
      left_join(hookup_unique)
       
    merged_df
    
    morsels <- foreach(i = seq(nrow(merged_df))) %do% {
      ref <- merged_df[i, ]$ref_nuc
      alts <- nuc_count_cols
      alts <- alts[!(alts %in% ref)]
      counts <- merged_df[i, nuc_count_cols]
      n_alleles <- sum(counts > 0)
      alt_counts <- merged_df[i, alts]
      n_alt <- max(alt_counts)
      alt <- colnames(alt_counts)[alt_counts == n_alt]
      n_ref <- as.numeric(merged_df[i, ref])
      
      merged_df[i, ] %>% 
        add_column(var_nuc = gsub("n", "", alt), 
                   n_alleles = n_alleles,
                   n_ref = n_ref, 
                   n_alt = n_alt) %>%
        mutate(ref_nuc = gsub("n", "", ref_nuc))
    }
    
    parsed_df <- bind_rows(morsels) %>%
      left_join(hookup, by = c("nucleotide_pos", "ref_nuc", "var_nuc"))
    
    fwrite(parsed_df, str_glue("{result_dir}/{prefix}_homoplasies.csv"))
  } else {
    fwrite(tibble(), str_glue("{result_dir}/{prefix}/{prefix}_homoplasies.csv"))
  }
}
 




  
  

