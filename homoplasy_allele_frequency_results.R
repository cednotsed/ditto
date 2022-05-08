rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(jcolors)
require(ggpubr)
require(ggrepel)
require(foreach)

human_background <- "V5"
iter_list <- list(mink = c("mink_only.n929_homoplasies.csv", 
                           str_glue("mink.{human_background}.global_allele_frequency.csv"),
                           str_glue("mink_homoplasy_alele_frequency_{human_background}.csv")),
                  deer = c("deer_only.n96_homoplasies.csv", 
                           str_glue("deer.{human_background}.global_allele_frequency.csv"),
                           str_glue("deer_homoplasy_alele_frequency_{human_background}.csv")))

plots <- foreach(i = seq(length(iter_list))) %do% {
  host_name <- Hmisc::capitalize(names(iter_list[i]))
  hp <- fread(str_glue("results/homoplasy_out/{iter_list[[i]][1]}")) %>%
    select(mutation, MinimumNumberChangesOnTree)
  
  freq <- fread(str_glue("results/allele_frequency/{iter_list[[i]][2]}"))
  
  freq_filt <- freq %>%
    mutate(fold_change = base_freq_animal / base_freq_human) %>%
    filter(fold_change >= 2,
           base_freq_animal > 0.1)
  
  merged_df <- freq %>%
    left_join(hp) %>% 
    mutate(min_changes = MinimumNumberChangesOnTree) %>%
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
  
  fwrite(merged_df, str_glue("results/allele_frequency/{iter_list[[i]][3]}"))
  
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
    labs(x = str_glue("{host_name} allele freq."), y = "Human allele freq.",
         shape = "Mutation type", color = "Region",
         size = "Min. changes") +
    scale_size_manual(values = c(1, 2, 3, 4, 6)) +
    scale_color_manual(values = pal) +
    ylim(-0.1, 0.6)
  
  ggsave(str_glue("results/allele_frequency/{names(iter_list)[i]}_homoplasy_allele_frequency_{human_background}_results.png"), 
         plot = plt,
         dpi = 300,
         width = 5, 
         height = 4)
  
  plt_legend <- plt + theme(legend.position = "top",
                            legend.title = element_text(size = 9),
                            legend.text = element_text(size = 9))
  
  ggsave(str_glue("results/allele_frequency/legend.homoplasy_allele_frequency_{human_background}_results.png"), 
         plot = plt_legend,
         dpi = 300,
         width = 10, 
         height = 4)
}

