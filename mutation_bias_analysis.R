rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(ape)
require(ggplot2)
require(ggrepel)
require(foreach)
require(doParallel)
registerDoParallel(cores = 6)

human_background <- "V5"
host_name <- "Mink"
species_name <- "Neovison vison"
file <- "all_mink.n1769.audacity_only.v8_masked.fasta"
meta_file <- "all_mink.n1769.csv"

human_background <- "V5"
host_name <- "Deer"
species_name <- "Odocoileus virginianus"
file <- "all_deer.n189.audacity_only.v8_masked.fasta"
meta_file <- "all_deer.n189.csv"

aln <- read.dna(str_glue("data/alignments/human_animal_subsets/{human_background}/{file}"),
                format = "fasta",
                as.matrix = T)

ref <- read.dna("data/genomes/WIV04_EPI_ISL_402124.fasta", 
                format = "fasta",
                as.matrix = T, 
                as.character = T)
ref_df <- tibble(pos = seq(ncol(ref)), ref_nuc = as.vector(ref))

ant_df <- read.csv("data/metadata/netherlands_humans_anthroponoses.txt", header = F)

# Remove anthroponoses
aln <- aln[!(rownames(aln) %in% ant_df$V1), ]

# Load metadata
meta <- fread(str_glue("data/metadata/human_animal_subsets/{human_background}/{meta_file}"))

meta.match <- meta[match(rownames(aln), meta$accession_id), ]

get_mutation_freqs <- function(aln, species_name, meta.match) {
  aln_filt <- aln[meta.match$host == species_name, ]
  n_isolates <- nrow(aln_filt)
  
  morsels <- foreach(pos = seq(ncol(aln)), .packages = c("ape", "tidyverse")) %dopar% {
    base_freq <- base.freq(aln_filt[, pos], freq = T)
    tibble(pos = pos,
           base = names(base_freq),
           base_freq = as.numeric(base_freq),
           norm_base_freq = as.numeric(base_freq) / n_isolates)
  }
  
  morsel <- bind_rows(morsels)
  
  ref_df %>%
    right_join(morsel) %>%
    filter(base_freq != 0, ref_nuc != base) %>%
    mutate(ref_nuc = toupper(ref_nuc),
           ref_nuc = ifelse(ref_nuc == "T", "U", ref_nuc),
           base = toupper(base),
           base = ifelse(base == "T", "U", base),
           nuc_trans = paste0(ref_nuc, "->", base)) %>%
    group_by(nuc_trans) %>%
    summarise(norm_count = sum(norm_base_freq),
              count = sum(base_freq)) %>%
    add_column(host = species_name)
}

animal_freqs <- get_mutation_freqs(aln, species_name, meta.match)
human_freqs <- get_mutation_freqs(aln, "Human", meta.match)

# Permutation test
# # Observed
# animal_prop <- animal_freqs %>%
#   mutate(prop = count / sum(count))
# human_prop <- human_freqs %>%
#   mutate(prop = count / sum(count))
# 
# obs <- animal_prop %>%
#   inner_join(human_prop, by = "nuc_trans", suffix = c("_animal", "_human")) %>%
#   select(nuc_trans, prop_animal, prop_human) %>%
#   mutate(diff = prop_animal - prop_human,
#          log_diff = log(prop_animal / prop_human, base = 10))
# 
# perm_morsels <- foreach(i = seq(1000), .packages = "tidyverse") %do% {
#   meta.perm <- tibble(accession_id = meta.match$accession_id,
#                       host = sample(meta.match$host))
#   animal_perm <- get_mutation_freqs(aln, species_name, meta.perm) %>%
#     mutate(prop = count / sum(count))
#   human_perm <- get_mutation_freqs(aln, "Human", meta.perm) %>%
#     mutate(prop = count / sum(count))
# 
#   merged_df_perm <- animal_perm %>%
#     inner_join(human_perm, by = "nuc_trans", suffix = c("_animal", "_human")) %>%
#     select(nuc_trans, prop_animal, prop_human) %>%
#     mutate(log_diff = log(prop_animal / prop_human, base = 10)) %>%
#     add_column(iter_n = i)
# 
#   merged_df_perm
# }
# 
# perm_df <- bind_rows(perm_morsels)
# 
# mut <- "G->A"
# mut_obs <- obs %>% filter(nuc_trans == mut)
# 
# perm_df %>%
#   filter(nuc_trans == mut) %>%
#   mutate(diff = prop_animal - prop_human) %>%
#   summarise(sum(diff > mut_obs$diff) / 1000)
# mut_obs$log_diff
# 
# test_morsels <- foreach (mut = obs$nuc_trans) %do% {
# # mut <- "G->A"
#   mut_obs <- obs %>% filter(nuc_trans == mut)
# 
#   if (mut_obs$log_diff > 0) {
#     n_ext <- perm_df %>%
#         filter(nuc_trans == mut) %>%
#         summarise(n = sum(log_diff > mut_obs$log_diff))
#     tibble(n_extreme = n_ext$n / 1000, mut = mut) 
#   } else {
#     n_ext <- perm_df %>%
#       filter(nuc_trans == mut) %>%
#       summarise(n = sum(log_diff < mut_obs$log_diff))
#     tibble(n_extreme = n_ext$n / 1000, mut = mut)
#   }
# }
# 
# bind_rows(test_morsels)
# # fwrite(perm_df, "results/permutation/permutation_df.n1000.csv")
# 
# perm_df
# 

## Chiseq.test ##
merged_df <- animal_freqs %>%
  inner_join(human_freqs, by = "nuc_trans", suffix = c("_animal", "_human")) %>%
  column_to_rownames("nuc_trans") %>%
  select(count_animal, count_human)
test <- chisq.test(merged_df, simulate.p.value = T)
pval <- test$p.value
chisq <- as.numeric(test$statistic)
print(str_glue("For {host_name}, Chisq = {chisq}, pval = {pval}"))

plot_df <- bind_rows(animal_freqs, human_freqs)

# Order bars
order_df <- plot_df %>% 
  filter(host == "Human") %>%
  arrange(desc(count))

plot_df %>%
  mutate(nuc_trans = factor(nuc_trans, order_df$nuc_trans)) %>%
  ggplot(aes(x = nuc_trans, y = norm_count, fill = nuc_trans)) +
  facet_grid(vars(host)) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none") +
  labs(x = "Mutation", y = "Avg. mutation per isolate")

ggsave(str_glue("results/mutation_bias/nucleotide_transitions_{host_name}.png"), 
       width = 6, height = 4)
############## S/NS mutations ##################################
# hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.csv") %>%
#   filter(!grepl("coronavirus|3'UTR ", protein_name))
# 
# # Consider only coding regions
# coding_regions <- hookup %>%
#   filter(region_type == "coding")
# coding_regions <- unique(coding_regions$nucleotide_pos)
# 
# get_mutations_by_pos <- function(aln, species_name) {
#   aln_filt <- aln[meta.match$host == species_name, ]
#   pos_of_interest <- seq(ncol(aln_filt))
#   pos_of_interest <- pos_of_interest[pos_of_interest %in% coding_regions]
#   
#   morsels <- foreach(pos = pos_of_interest, .packages = c("ape", "tidyverse")) %dopar% {
#     base_freq <- base.freq(aln_filt[, pos], freq = T)
#     
#     tibble(pos = pos,
#            base = names(base_freq),
#            base_freq = as.numeric(base_freq))
#   }
#   
#   morsel <- bind_rows(morsels)
#   
#   ref_df %>%
#     right_join(morsel) %>%
#     filter(base_freq != 0) %>%
#     filter(ref_nuc != base) %>%
#     mutate(ref_nuc = toupper(ref_nuc),
#            var_nuc = toupper(base)) %>%
#     rename(nucleotide_pos = pos) %>%
#     left_join(hookup)
# }
# 
# animal_pos <- get_mutations_by_pos(aln, species_name)
# human_pos <- get_mutations_by_pos(aln, "Human")
# 
## DNDS ##
# Get no. of isolates per host
# n_animal <- nrow(aln[meta.match$host == species_name, ])
# n_human <- nrow(aln[meta.match$host == "Human", ])
# n_mut_animal <- as.numeric(animal_pos %>% summarise(sum(base_freq) / n_animal))
# n_mut_human <- as.numeric(human_pos %>% summarise(sum(base_freq) / n_human))
# 
# animal_df <- animal_pos %>%
#   group_by(mutation_type) %>%
#   summarise(value = sum(base_freq) / n_animal) %>%
#   rename(stat = mutation_type)
# 
# animal_dnds <- animal_df[1, ]$value / animal_df[2, ]$value
# animal_df <- animal_df %>%
#   bind_rows(tibble(stat = c("NS/S", "NS + S", "Chi2", "pval"),
#                    value = c(animal_dnds, n_mut_animal, chisq, pval)))
# 
# human_df <- human_pos %>%
#   group_by(mutation_type) %>%
#   summarise(value = sum(base_freq) / n_human) %>%
#   rename(stat = mutation_type)
# 
# human_dnds <- human_df[1, ]$value / human_df[2, ]$value
# human_df <- human_df %>%
#   bind_rows(tibble(stat = c("NS/S", "NS + S", "Chi2", "pval"),
#                    value = c(human_dnds, n_mut_human, chisq, pval)))
# final_df <- human_df %>% 
#   inner_join(animal_df, by = "stat", suffix = c("_human", "_animal"))
# 
# final_df 
# fwrite(final_df, str_glue("results/mutation_bias/{host_name}_{human_background}_dnds.csv"))
