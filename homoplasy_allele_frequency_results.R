rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(jcolors)
require(ggpubr)
require(ggrepel)
require(foreach)

human_background <- "V1"
iter_list <- list(mink = c("mink_n789_201121_parsed.csv", 
                           "Mink.global_allele_frequency.csv",
                          str_glue( "mink_homoplasy_alele_frequency_{human_background}.csv")),
                  deer = c("deer_n73_201121_parsed.csv", 
                           "Deer.global_allele_frequency.csv",
                           str_glue("deer_homoplasy_alele_frequency_{human_background}.csv")))

plots <- foreach(i = seq(length(iter_list))) %do% {
  hp <- fread(str_glue("results/all_animals/{iter_list[[i]][1]}")) %>%
    select(mutation, MinimumNumberChangesOnTree)
  
  freq <- fread(str_glue("results/human_animal_subsets/allele_frequency/{human_background}/{iter_list[[i]][2]}"))
  
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
    mutate(mutation_annot = ifelse(mutation_annot == "-1-1-1", "Non-coding", mutation_annot))
  
  nrow(merged_df)
  
  fwrite(merged_df, str_glue("results/{iter_list[[i]][3]}"))
  
  pal <- jcolors("pal5")
  names(pal)[1:length(levels(merged_df$region))] <- levels(merged_df$region)
  
  plt <- merged_df %>% 
    filter(mutation_type == "NS") %>%
    ggplot(aes(x = base_freq_animal, y = base_freq_human)) +
    geom_point(aes(size = min_changes,
                   shape = mutation_type,
                   color = region),
               alpha = 0.8) +
    geom_text_repel(aes(label = mutation_annot),
                    size = 3,
                    color = "black",
                    max.overlaps = 25) +
    geom_abline(slope = 1, intercept = 0, color = "black") +
    geom_abline(slope = 0.5, intercept = 0, lty = "dotted", color = "red") +
    theme(legend.position = "none"
          # legend.title = element_text(size = 9),
          # legend.text = element_text(size = 9)
          ) +
    labs(x = "Animal allele freq.", y = "Human allele freq.",
         shape = "Mutation type", color = "Region",
         size = "Min. changes") +
    scale_size_manual(values = c(1, 2, 3, 4, 6)) +
    scale_color_manual(values = pal) +
    ylim(-0.1, 0.6)
  plt
  
  ggsave(str_glue("results/{names(iter_list)[i]}_homoplasy_allele_frequency_{human_background}_results.png"), 
         dpi = 300,
         width = 5, 
         height = 4)
}

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
