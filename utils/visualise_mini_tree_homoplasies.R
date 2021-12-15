setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(foreach)
require(jcolors)

result_dir <- "results/mini_animal_trees/homoplasy_out"
prefixes <- list.files(result_dir)

# Add variant names
all_meta <- fread("data/metadata/all_sequence_metadata_231121.tsv") %>%
  rename_all(~ tolower(gsub(" ", "_", .))) %>%
  rename(Accession = accession_id)

View(all_meta)

df_morsels <- foreach(prefix = prefixes) %do% {
  df <- fread(str_glue("{result_dir}/{prefix}/{prefix}_homoplasies.csv")) %>%
      add_column(prefix = prefix) %>%
      separate(prefix, into = c("animal", "country", "cluster", "n_samples"))
  
  if(nrow(df) > 0) {
    return(df)
  }
}

df <- bind_rows(df_morsels)

cluster_counts <- df %>% 
  group_by(nucleotide_pos, animal) %>%
  summarise(n_clusters = n()) 

# Get protein list
morsels <- foreach(prefix = prefixes) %do% {
  fread(str_glue("{result_dir}/{prefix}/{prefix}_homoplasies.csv"))
}

protein_list <- bind_rows(morsels) %>% distinct(region)
protein_list <- protein_list$region

plot_df <- df %>%
  filter(MinimumNumberChangesOnTree >= 2) %>% 
  group_by(nucleotide_pos) %>%
  slice_max(order_by = MinimumNumberChangesOnTree, n = 1) %>% # Get highest for each position
  sample_n(1) %>% # Randomly select 1 row if tie in min_changes
  left_join(cluster_counts) %>%
  mutate(region = factor(region, levels = protein_list),
         mutation_name = case_when(mutation_type == "NS" ~ paste0(ref_AA, codon_number, var_AA),
                                   mutation_type == "S" ~ paste0(ref_nuc, nucleotide_pos, var_nuc))) %>%
  mutate(to_annotate = ifelse(mutation_type == "NS", mutation_name, ""))

  
plot_df %>%
    ggplot(aes(x = factor(nucleotide_pos), 
               y = factor(MinimumNumberChangesOnTree), 
               shape = mutation_type,
               color = region,
               size = factor(n_clusters))) +
    theme_bw() +
    facet_grid(rows = vars(animal)) +
    geom_point(alpha = 0.5) +
    scale_color_jcolors(palette = "pal8") +
    geom_text(aes(label = to_annotate),
              color = "black", 
              size = 2, 
              angle = 45, 
              vjust = 0.5, hjust = -0.1) +
    labs(x = "Position", y = "Min. no. of changes", 
         shape = "Mutation type", 
         color = "Genomic region",
         size = "No. of spillover events") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) 

ggsave("results/mini_animal_trees/homoplasy_plots_mini_trees.png", dpi = 300, width = 8, height = 5)


mut <- "G142D"
mut <- "Y453F"

nrow(all_meta[grepl(mut, all_meta$aa_substitutions), ])
all_meta[grepl(mut, all_meta$aa_substitutions), ]
