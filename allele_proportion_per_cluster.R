rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(ape)
require(jcolors)
require(ggpubr)
require(ggrepel)
require(foreach)

human_background <- "V5"
host_name <- "deer"

plots <- foreach(host_name = c("mink", "deer")) %do% {
  mut_df <- fread(str_glue("results/{host_name}_homoplasy_alele_frequency_{human_background}.csv")) %>%
    mutate(mutation_annot = paste(protein_name, mutation))
  
  cluster_meta <- fread("results/cluster_annotation/deer_mink_parsed_clusters.csv")

  aln_file <- list.files("data/alignments/all_animals/", pattern = host_name, full.names = T)
  aln_file <- aln_file[grepl("audacity_only.v8_masked.aln.fasta", aln_file)]
  
  aln <- read.dna(aln_file, 
                  format = "fasta", 
                  as.matrix = T,
                  as.character = T)
  aln <- aln[rownames(aln) != "NC_045512.2", ]
  aln_filt <- as_tibble(aln)
  aln_filt <- aln_filt[, mut_df$nucleotide_pos]
  aln_filt <- aln_filt %>% 
    mutate_all(toupper) %>%
    add_column(accession_id = rownames(aln)) %>%
    left_join(cluster_meta)
  
  pos_of_interest <- 245
  
  morsels <- foreach (pos_of_interest = mut_df$nucleotide_pos) %do% {
    pos_meta <- mut_df %>% filter(nucleotide_pos == pos_of_interest)
    
    aln_filt %>%
      group_by(cluster, get(paste0("V", pos_of_interest))) %>% 
      summarise(n = n() ) %>%
      rename(allele = `get(paste0("V", pos_of_interest))`) %>%
      mutate(prop = n / sum(n)) %>%
      filter(allele == pos_meta$var_nuc) %>%
      add_column(mut = pos_meta$mutation_annot)
  }
  
  plot_df <- bind_rows(morsels) %>% 
    mutate(mut = gsub(" \\(stem-loop II-like motif\\) ", "", mut))
  
  order_df <- plot_df %>%
    group_by(mut) %>%
    summarise(n = n()) %>%
    arrange(n)
  
  plot_df %>%
    mutate(mut = factor(mut, levels = order_df$mut),
           cluster = gsub(paste0(host_name, "_"), "", cluster)) %>%
    ggplot(aes(x = cluster, y = mut, fill = prop)) +
    geom_tile() +
    scale_fill_gradient(low = "blue", high = "salmon") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Cluster", 
          # title = str_to_title(host_name),
         y = "Mutation", 
         fill = "Prop. of isolates")
  
  ggsave(str_glue("results/{host_name}_founder_effect_plot.png"), 
         dpi = 600,
         width = 6, 
         height = 4.5)
}

# ggpubr::ggarrange(plotlist = plots,
#                   ncol = 1,
#                   common.legend = T)
# ggsave("results/founder_effect_plot.png", 
#        dpi = 300,
#        width = 5, 
#        height = 9)
# meta.match <- cluster_meta[match(rownames(aln), cluster_meta$accession_id), ]
# 
# all(meta.match$accession_id == rownames(aln))


