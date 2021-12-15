setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(foreach)
require(jcolors)

prefixes <- list.files("results/homoplasy_out/")
prefixes <- prefixes[!grepl(".csv", prefixes)]

# Get protein list
morsels <- foreach(prefix = prefixes) %do% {
  fread(str_glue("results/homoplasy_out/{prefix}_homoplasies.csv"))
}

protein_list <- bind_rows(morsels) %>% distinct(region)
protein_list <- protein_list$region

plots <- foreach (prefix = prefixes) %do% {
  df <- fread(str_glue("results/homoplasy_out/{prefix}_homoplasies.csv"))
  
  df %>%
    mutate(n_alleles = factor(n_alleles),
           region = factor(region, levels = protein_list),
           mutation_name = paste0("(", ref_AA, codon_number, var_AA, ")")) %>%
    mutate(mutation_name = ifelse(mutation_type == "NS", mutation_name, "")) %>%
    mutate(to_annotate = ifelse(MinimumNumberChangesOnTree > 2, 
                                str_glue("{nucleotide_pos} {mutation_name}"), "")) %>%
    # filter(MinimumNumberChangesOnTree > 8) %>%
    ggplot(aes(x = nucleotide_pos, 
               y = MinimumNumberChangesOnTree, 
               shape = mutation_type,
               color = region,
               size = n_alleles)) +
    geom_point(alpha = 0.5) +
    scale_color_jcolors(palette = "pal8") +
    geom_text(aes(label = to_annotate), 
              color = "black", size = 2, 
              vjust = 0, hjust = 0,
              position = position_jitter(width = 0.5, height = 0.5, seed = 66)) +
    labs(x = "Position", y = "Min. no. of changes", 
         shape = "Mutation type", 
         color = "Genomic region",
         size = "No. of alleles",
         title = prefix)
}

combined <- ggpubr::ggarrange(plotlist = plots, common.legend = T, ncol = 2, nrow = 3) 
ggsave("results/homoplasy_plots_all.png", plot = combined, dpi = 300, width = 15, height = 10)
