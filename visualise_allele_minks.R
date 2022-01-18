setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(ape)
require(ggtree)
require(ggplot2)
require(ggtreeExtra)
require(ggnewscale)
require(jcolors)

prefix <- "mink_n789_201121"
# prefix <- "deer_n73_201121"
tree <- read.tree(str_glue("data/trees/all_animals/{prefix}.audacity_only.v8_masked.tree"))

animal_meta <- fread("data/metadata/combined_metadata_201121.audacity_only.csv") %>%
  separate(location, into = c("loc1", "loc2"), sep = " / ") %>%
  mutate(location = paste0(loc2)) %>%
  select(accession_id, host, clade, lineage, location, -cluster)

cluster_meta <- fread("results/cluster_annotation/deer_mink_parsed_clusters.csv") %>%
  select(accession_id, cluster)

animal_meta <- animal_meta %>% left_join(cluster_meta)

aln <- read.dna(str_glue("data/alignments/all_animals/{prefix}.audacity_only.v8_masked.aln.fasta"),
                format = "fasta",
                as.matrix = T)

mut_df <- fread("results/mink_homoplasy_alele_frequency_V5.csv") %>%
  filter(mutation_annot %in% c("G37E","F486L", "N501T", "T229I", "L219V", "D178Y", "A192V", "I258V",
                               "T18A", "N507I", "D377Y", "I82T"))
mut_df

for (i in seq(nrow(mut_df))) {
  row <- mut_df[i, ]  
  
  aln_df <- tibble(accession_id = rownames(aln), allele = toupper(as.vector(as.character(aln[, row$nucleotide_pos]))))
  
  animal_parsed <- animal_meta %>%
    left_join(aln_df) %>%
    mutate(allele_annot = case_when(allele == row$ref_nuc ~ "WT",
                                    allele == row$var_nuc ~ paste0(row$protein_name, 
                                                                   "_", 
                                                                   row$mutation_annot),
                                    TRUE ~ as.character(NA)))
  
  # Match metadata
  animal_parsed <- animal_parsed[match(tree$tip.label, animal_parsed$accession_id), ]
  
  dd <- data.frame(Accession = tree$tip.label,
                   lineage = animal_parsed$lineage,
                   cluster = animal_parsed$cluster,
                   clade = animal_parsed$clade,
                   location = animal_parsed$location,
                   allele = animal_parsed$allele_annot)
  
  # Plot tree
  p <- ggtree(tree, color = "darkgrey")
  p <- p %<+% dd
  p <- p + 
    geom_fruit(geom = geom_tile, aes(fill = location), width=0.0018, offset=-0.6, alpha = 0.3) +
    labs(fill = "Location") +
    scale_fill_jcolors("pal8", na.translate = F) +
    new_scale_fill() +
    geom_fruit(geom = geom_tile, aes(fill = cluster), width=0.0001, offset=0.5, alpha = 1) +
    geom_tippoint(aes(color = allele), size = 3) +
    scale_color_discrete(na.translate = F) +
    labs(fill = "Cluster", color = str_glue("Position {row$nucleotide_pos}")) +
    theme_tree2()
  
  p
  ggsave(str_glue("results/all_animals/{prefix}_allele_{row$mutation_annot}_annotation.png"), 
         plot = p, 
         dpi = 300, 
         width = 10, 
         height = 10)

}

# gen_dist <- cophenetic(tree)[, "NC_045512.2"]
# df <- tibble(accession_id = names(gen_dist), root.to.tip = gen_dist) %>%
#   inner_join(meta) %>%
#   mutate(collection_date = as.Date(collection_date))
# 
# df %>%
#   ggplot(aes(x = collection_date, y = root.to.tip)) +
#   geom_point()
# 
# lr <- lm(df$root.to.tip ~ df$collection_date)
# summary(lr)
