rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(jcolors)
require(ggpubr)
require(ggrepel)
require(foreach)
require(ggnewscale)

parse_results <- function(host_name, hp_path, freq_path) {
  hp <- fread(hp_path) %>%
    select(mutation, MinimumNumberChangesOnTree)
  
  freq <- fread(freq_path)
  
  freq_filt <- freq %>%
    mutate(fold_change = base_freq_animal / base_freq_human) %>%
    filter(fold_change >= 2,
           base_freq_animal > 0.1)
  
  merged_df <- freq %>%
    left_join(hp) %>% 
    rename(min_changes = MinimumNumberChangesOnTree) %>%
    filter(min_changes > 3 | mutation %in% freq_filt$mutation) %>%
    mutate(min_changes = case_when(is.na(min_changes) ~ "<2",
                                   min_changes == 2 ~ "2",
                                   min_changes >= 3 & min_changes <= 5 ~ "3-5",
                                   min_changes >= 6 & min_changes <= 10 ~ "6-10",
                                   min_changes > 10 ~ ">10")) %>% 
    mutate(mutation_annot = paste0(ref_AA, codon_number, var_AA),
           min_changes = factor(min_changes, levels = c("<2", "2", "3-5", "6-10", ">10")),
           region = factor(region, levels = c("ORF1ab", "S", "ORF3a", "M", "ORF8", "N")),
           mutation_type = factor(mutation_type, levels = c("NS", "S"))) %>%
    mutate(mutation_annot = ifelse(mutation_annot == "-1-1-1", "Non-coding", mutation_annot)) %>%
    add_column(host = host_name)
  
  return(merged_df)
}

mink_df <- parse_results("mink",
                         "results/all_animals/mink_n789_201121_parsed.csv",
                         "results/human_animal_subsets/allele_frequency/V1/Mink.global_allele_frequency.csv")

deer_df <- parse_results("deer", 
                         "results/all_animals/deer_n73_201121_parsed.csv",
                         "results/human_animal_subsets/allele_frequency/V1/Deer.global_allele_frequency.csv")

# All mutations
v5_freq <- fread("results/mink_homoplasy_alele_frequency_V5.csv") %>%
  bind_rows(fread("results/deer_homoplasy_alele_frequency_V5.csv"))

plot_df <- mink_df %>% 
  bind_rows(deer_df) %>%
  filter(nucleotide_pos %in% v5_freq$nucleotide_pos) %>%
  mutate(high_conf = ifelse(mutation_annot %in% c("Y453F", "N501T", "F486L", "T229I", "L219V",
                                                  "L1035F", "I82T", "D377Y"),
                            mutation_annot, NA)) %>%
  mutate(is_high_conf = ifelse(is.na(high_conf) | high_conf %in% c("I82T", "D377Y"), NA, "red"))

# High conf mutations
plot_df %>%
  ggplot(aes(x = base_freq_animal, y = base_freq_human, color = host, shape = host)) +
  geom_point() +
  labs(x = "Animal allele freq.", y = "Human allele freq.",
       shape = "Host", color = "Host") +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  geom_abline(slope = 0.5, intercept = 0, lty = "dotted", color = "red") +
  geom_text_repel(aes(label = high_conf), min.segment.length = 15, 
                  color = "black") +
  new_scale_color() +
  geom_point(aes(x = base_freq_animal, y = base_freq_human, color = is_high_conf),
             shape = 0) +
  scale_color_manual(values = "red", na.translate = F) +
  ylim(-0.05, 0.25)

ggsave(str_glue("results/mink_deer_homoplasy_allele_frequency_V1_results.png"), 
       dpi = 300,
       width = 8, 
       height = 2)


# combined <- ggpubr::ggarrange(plotlist = plots, ncol = 1, nrow = 2, 
#                               common.legend = T,
#                               legend = "bottom")
# combined
# ggsave(str_glue("results/homoplasy_allele_frequency_{human_background}_results.png"), combined, dpi = 300, width = 6, height = 8)

############ PROVEAN #################
# # Check duplicates
# positions <- c(merged_df_deer$nucleotide_pos, merged_df_mink$nucleotide_pos)
# any(duplicated(positions))
# 
# merged_df_deer$host <- "Deer"
# merged_df_mink$host <- "Mink"
# 
# protein_meta <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.csv")
# protein_order <- unique(protein_meta$protein_name)
# 
# plot_df <- merged_df_deer %>% 
#   bind_rows(merged_df_mink) %>%
#   filter(base_freq_human < 0.1)
#   
# prot_df <- plot_df %>%
#   group_by(host, protein_name) %>%
#   summarise(n = n()) %>% 
#   mutate(protein_name = factor(protein_name,
#                                levels = protein_order))
# 
# type_df <- plot_df %>%
#   group_by(host, mutation_type) %>%
#   summarise(n = n())
# 
# type_plt <- type_df %>%
#   ggplot(aes(x = mutation_type, y = n, fill = host)) +
#   geom_bar(stat = "identity",
#            position = position_dodge2(preserve = "single")) +
#   labs(x = "Type", y = "No. of mutations", fill = "Host") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# bar_plt <- prot_df %>%
#   ggplot(aes(x = protein_name, y = n, fill = host)) +
#     geom_bar(stat = "identity",
#              position = position_dodge2(preserve = "single")) +
#   labs(x = "Protein", y = "No. of mutations", fill = "Host") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# combined <- ggarrange(type_plt, bar_plt, nrow = 2, common.legend = T)
# ggsave("results/mutation_breakdown.png", combined, dpi = 300, width = 5, height = 4)


# merged_df_deer %>% 
#   bind_rows(merged_df_mink) %>%
#   # ggplot(aes(x = protein_name)) +
#   # geom_histogram(stat = "count")
#   filter(base_freq_human < 0.1) %>% 
#   filter(region == "S", ref_AA != var_AA) %>% 
#   mutate(variant_name = paste0(ref_AA, codon_from_gene_start, var_AA)) %>%
#   select(variant_name) %>%
#   fwrite("results/mutations.txt")
