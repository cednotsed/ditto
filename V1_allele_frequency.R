rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(jcolors)
require(ggpubr)
require(ggrepel)
require(foreach)
require(ggnewscale)

deer_cluster <- fread("results/founder_alleles/deer_allele_frequency_by_cluster.csv")
mink_cluster <- fread("results/founder_alleles/mink_allele_frequency_by_cluster.csv")

# Get final list of mutations
deer_mutations <- fread("results/allele_frequency/deer_homoplasy_alele_frequency_V5.csv") %>% 
  mutate(mut = paste(protein_name, mutation)) %>%
  left_join(deer_cluster) %>%
  mutate(high_conf = ifelse(fold > 2 &
                              base_freq_animal > 0.1 &
                              MinimumNumberChangesOnTree > 2 &
                              n > 2, T, F)) %>%
  select(mutation, high_conf)

mink_mutations <- fread("results/allele_frequency/mink_homoplasy_alele_frequency_V5.csv") %>% 
  mutate(mut = paste(protein_name, mutation)) %>%
  left_join(mink_cluster) %>%
  mutate(high_conf = ifelse(fold > 2 &
                              base_freq_animal > 0.1 &
                              MinimumNumberChangesOnTree > 2 &
                              n > 2, T, F)) %>%
  select(mutation, high_conf)

high_conf_deer <- (deer_mutations %>% filter(high_conf))$mutation
high_conf_mink <- (mink_mutations %>% filter(high_conf))$mutation
high_conf_mink <- c("A22920T (Y453F)", high_conf_mink)

# V1 frequencies
deer <- fread("results/allele_frequency/deer_homoplasy_alele_frequency_V1.csv") %>%
  mutate(host = "deer") %>%
  filter(mutation %in% deer_mutations$mutation)

mink <- fread("results/allele_frequency/mink_homoplasy_alele_frequency_V1.csv") %>%
  mutate(host = "mink") %>%
  filter(mutation %in% mink_mutations$mutation)

plot_df <- bind_rows(deer, mink) %>%
  mutate(mutation_annot = ifelse(mutation_type == "S", mutation, mutation_annot)) %>%
  mutate(high_conf_annot = ifelse(mutation %in% c(high_conf_deer, high_conf_mink), mutation_annot, NA),
         high_conf_shape = ifelse(mutation %in% c(high_conf_deer, high_conf_mink), mutation_type, NA))
  
nrow(deer) + nrow(mink) == nrow(plot_df)

plt <- plot_df %>% 
  ggplot(aes(x = base_freq_animal, y = base_freq_human)) +
  geom_point(aes(shape = host, fill = host),
             alpha = 0.8) +
  scale_shape_manual(values = c(21, 24)) +
  geom_text_repel(aes(label = high_conf_annot),
                  color = "black",
                  size = 3,
                  max.overlaps = 25) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  geom_abline(slope = 0.5, intercept = 0, lty = "dotted", color = "red") +
  # theme(legend.position = "none"
  #       # legend.title = element_text(size = 9),
  #       # legend.text = element_text(size = 9)
  # ) +
  labs(x = "Animal allele freq.", y = "Human allele freq.", 
       fill = "Host", shape = "Host") +
  ylim(-0.05, 0.3) +
  xlim(0, 1)
plt

ggsave("results/allele_frequency/V1_allele_freq.png", 
       dpi = 300,
       width = 8, 
       height = 2)



