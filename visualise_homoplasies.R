setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(foreach)
require(jcolors)
require(ggpubr)
require(ggrepel)

prefixes <- list.files("results/all_animals/homoplasy_out/")
prefixes <- prefixes[!grepl(".csv", prefixes)]
# prefixes <- prefixes[grepl("mink|deer", prefixes)]

# # Get protein list
# morsels <- foreach(prefix = prefixes) %do% {
#   fread(str_glue("results/all_animals/homoplasy_out/{prefix}_homoplasies.csv"))
# }
# 
# protein_list <- bind_rows(morsels) %>% distinct(region)
# protein_list <- protein_list$region

plots <- foreach (prefix = prefixes) %do% {
  df <- fread(str_glue("results/all_animals/homoplasy_out/{prefix}_homoplasies.csv"))
  host_name <- tools::toTitleCase(str_split(prefix, "_")[[1]][1])
  n <- gsub("n", "", str_split(prefix, "_")[[1]][2])

  parsed_df <- df %>%
    mutate(mutation_name = paste0(ref_nuc, nucleotide_pos, var_nuc), 
           variant_name = paste0(ref_AA, codon_number, var_AA)) %>%
    mutate(variant_name = case_when(mutation_type == "S" & region_type != "non-coding" ~ "Syn",
                                    region_type == "non-coding" ~ "non-coding",
                                    TRUE ~ variant_name)) %>%
    mutate(mutation = str_glue("{mutation_name} ({variant_name})")) %>%
    arrange(nucleotide_pos) %>%
    filter(protein_name != "3'UTR") %>% # Remove 3'UTR entry since separate entry is present
    select(nucleotide_pos, region, region_type, mutation_type, protein_name, mutation, MinimumNumberChangesOnTree) %>%
    mutate(region = ifelse(region_type == "non-coding", "Non-coding", region))
  
  fwrite(parsed_df, str_glue("results/all_animals/{prefix}_parsed.csv"))
  
  if (prefix == "deer_n73_201121") {
    threshold <- 3
  } else(
    threshold <- 2
  )
  
  parsed_df %>%
    mutate(to_annotate = ifelse(MinimumNumberChangesOnTree > threshold, mutation, NA)) %>%
    ggplot(aes(x = factor(nucleotide_pos), 
               y = MinimumNumberChangesOnTree, 
               shape = mutation_type,
               color = region)) +
    geom_point(alpha = 0.8, size = 4) +
    scale_color_jcolors(palette = "pal8") +
    geom_text_repel(aes(label = to_annotate), 
              color = "black", size = 4) +
    labs(x = "Position", y = "Min. no. of changes", 
         shape = "Mutation type", 
         color = "Genomic region",
         size = "No. of alleles",
         title = str_glue("{host_name} (n = {n})")) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

  }

combined <- ggpubr::ggarrange(plotlist = plots, common.legend = T, ncol = 2, nrow = 1)
combined <- ggpubr::ggarrange(plotlist = plots, common.legend = F) 
ggsave("results/all_animals/homoplasy_plots_all_animals.png", plot = combined, dpi = 300, width = 10, height = 5)

