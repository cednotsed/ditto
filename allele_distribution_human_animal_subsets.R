setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(ape)
require(ggtree)
require(ggplot2)
require(ggtreeExtra)
require(ggnewscale)
require(jcolors)

# prefix <- "all_mink.n1769"
# host_name <- "mink"
prefix <- "all_deer.n189"
host_name <- "deer"
tree <- read.tree(str_glue("data/trees/human_animal_subsets/V5/{prefix}.audacity_only.v8_masked.unambiguous.dedup.tree"))

animal_meta <- fread(str_glue("data/metadata/human_animal_subsets/V5/{prefix}.csv"))

cluster_meta <- fread("results/cluster_annotation/deer_mink_parsed_clusters.csv") %>%
  select(accession_id, cluster)

animal_meta <- animal_meta %>% left_join(cluster_meta)

aln <- read.dna(str_glue("data/alignments/human_animal_subsets/V5/{prefix}.audacity_only.v8_masked.fasta"),
                format = "fasta",
                as.matrix = T)

mut_df <- fread(str_glue("results/allele_frequency/{host_name}_homoplasy_alele_frequency_V5.csv"))
  # filter(mutation_annot %in% c("G37E","F486L", "N501T", "T229I", "L219V", "D178Y", "A192V", "I258V",
  #                              "T18A", "N507I", "D377Y", "I82T", "L1035F",
  #                              "Y453F", "H182Y", "P199L", 
  #                              # syn changes
  #                              "P181P", "I1528I", "I292I", "L45L"))
mut_df

df %>%
  filter(fold > 2, base_freq_animal > 0.1,
         !(min_changes %in% c("<2", "2"))) %>%
  mutate(full_annot = paste0(protein_name, "_", mutation))

for (i in seq(nrow(mut_df))) {
  row <- mut_df[i, ]  
  
  aln_df <- tibble(accession_id = rownames(aln), allele = toupper(as.vector(as.character(aln[, row$nucleotide_pos]))))
  
  animal_parsed <- animal_meta %>%
    left_join(aln_df) %>%
    mutate(allele_annot = case_when(allele == row$var_nuc ~ paste0(row$protein_name, 
                                                                   "_", 
                                                                   row$mutation_annot),
                                    TRUE ~ as.character(NA)))
  
  # Match metadata
  animal_parsed <- animal_parsed[match(tree$tip.label, animal_parsed$accession_id), ]
  
  dd <- data.frame(Accession = tree$tip.label,
                   lineage = animal_parsed$pango_lineage,
                   host = animal_parsed$host,
                   cluster = animal_parsed$cluster,
                   clade = animal_parsed$clade,
                   location = animal_parsed$location,
                   allele = animal_parsed$allele_annot)
  
  # Plot tree
  p <- ggtree(tree, aes(color = host)) %<+% dd +
    geom_fruit(geom = geom_tile, aes(fill = location), width=0.0018, offset=-0.6, alpha = 0.3) +
    labs(fill = "Location") +
    scale_fill_jcolors("pal8", na.translate = F) +
    new_scale_fill() +
    geom_tippoint(aes(shape = host)) +
    labs(shape = "Host", color = "host") +
    new_scale_color() +
    geom_fruit(geom = geom_tile, aes(fill = cluster), width=0.0001, offset=0.5, alpha = 1) +
    geom_tippoint(aes(color = allele), shape = 0) +
    scale_color_manual(values = "red", na.translate = F) +
    labs(fill = "Cluster", color = str_glue("Position {row$nucleotide_pos}")) +
    theme_tree2()
  
  # p_zoom <- ggtree(tree, aes(color = host)) %<+% dd +
  #   labs(fill = "Location") +
  #   scale_fill_jcolors("pal8", na.translate = F) +
  #   new_scale_fill() +
  #   geom_tippoint(aes(shape = host)) +
  #   labs(shape = "Host", color = "host") +
  #   new_scale_color() +
  #   geom_tippoint(aes(color = allele), shape = 0) +
  #   scale_color_manual(values = "red", na.translate = F) +
  #   labs(fill = "Cluster", color = str_glue("Position {row$nucleotide_pos}")) +
  #   theme_tree2()
  # ggsave(str_glue("results/all_animals/founder_alleles/D178Y_subtree.png"),
  #        plot = p_zoom,
  #        dpi = 600,
  #        width = 10,
  #        height = 10)
  ggsave(str_glue("results/founder_alleles/founder_trees/{host_name}/{prefix}_allele_{row$mutation_annot}_annotation.png"),
         plot = p,
         dpi = 100,
         width = 10,
         height = 10)
  
}


