rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(ape)
require(ggplot2)
require(foreach)

human_background <- "V5"
host_name <- "Mink"
species_name <- "Neovison vison"
file <- "all_mink.n1769.audacity_only.v8_masked.fasta"
meta_file <- "all_mink.n1769.csv"

mut_df <- fread("results/allele_frequency/mink_homoplasy_alele_frequency_V5.csv") %>%
  filter(variant_name %in% c("Y453F", "F486L", "N501T"))
mut_df

aln <- read.dna(str_glue("data/alignments/all_animals/mink_only.n929.audacity_only.v8_masked.fasta"),
                format = "fasta",
                as.matrix = T)

mut_df$nucleotide_pos

morsels <- foreach(i = seq(nrow(aln))) %do% {
  aln_filt <- aln[i, mut_df$nucleotide_pos]
  acc <- rownames(aln_filt)
  as.data.frame(as.character(aln_filt))
}

df <- bind_rows(morsels)
colnames(df) <- mut_df$variant_name
parsed <- df %>%
  mutate(Y453F = ifelse(Y453F == tolower(deframe(mut_df[1, "var_nuc"])), T, F),
         F486L = ifelse(F486L == tolower(deframe(mut_df[2, "var_nuc"])), T, F),
         N501T = ifelse(N501T == tolower(deframe(mut_df[3, "var_nuc"])), T, F)) %>%
  mutate(n_mut = Y453F + F486L + N501T) %>%
  rownames_to_column("accession_id") %>%
  arrange(desc(n_mut))

parsed 
parsed %>%
  pivot_longer(!c(accession_id, n_mut), names_to = "mutation", values_to = "presence") %>%
  ggplot(aes(y = accession_id, x = mutation, fill = presence)) +
  geom_tile()
parsed$n_mut
table(parsed$n_mut)

parsed %>%
  mutate(type = case_when(Y453F & !F486L & !N501T ~ "Y453F",
                          Y453F & !F486L & N501T ~ "Y453F + N501T", 
                          Y453F & F486L & !N501T ~ "Y453F + F486L",
                          !Y453F & F486L & !N501T ~ "F486L",
                          !Y453F & F486L & N501T ~ "F486L + N501T",
                          !Y453F & !F486L & N501T ~ "N501T")) %>%
  View()
sum(parsed$n_mut == 2)
