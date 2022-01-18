setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(ape)
require(ggtree)
require(ggplot2)
require(ggtreeExtra)
require(ggnewscale)
require(jcolors)

prefix <- "denmark_minks.n15879"
host_name <- "mink"
# prefix <- "all_deer.n145"
# host_name <- "deer"
audacity <- read.tree(str_glue("data/GISAID-hCoV-19-phylogeny-2021-11-16/global.tree"))

animal_meta <- fread(str_glue("data/metadata/human_animal_subsets/V6/{prefix}.csv"))

# Drop tips
to_drop <- audacity$tip.label[!(audacity$tip.label %in% animal_meta$accession_id)]
audacity_filt <- drop.tip(audacity, to_drop)

aln <- read.dna(str_glue("data/alignments/human_animal_subsets/V6/{prefix}.audacity_only.v8_masked.fasta"),
                format = "fasta",
                as.matrix = T)

mut_df <- fread(str_glue("results/{host_name}_homoplasy_alele_frequency_V5.csv")) %>%
  filter(mutation_annot %in% c("Y453F"))
mut_df

for (i in seq(nrow(mut_df))) {
  row <- mut_df[i, ]  
  
  aln_df <- tibble(accession_id = rownames(aln), allele = toupper(as.vector(as.character(aln[, row$nucleotide_pos]))))
  
  animal_parsed <- animal_meta %>%
    left_join(aln_df) %>%
    mutate(allele_annot = case_when(allele == row$var_nuc ~ paste0(row$protein_name, 
                                                                   "_", 
                                                                   row$mutation_annot),
                                    TRUE ~ as.character(NA)),
           cluster = ifelse(cluster == "", NA, cluster))
  
  # Match metadata
  animal_parsed <- animal_parsed[match(audacity_filt$tip.label, animal_parsed$accession_id), ]
  
  dd <- data.frame(Accession = audacity_filt$tip.label,
                   lineage = animal_parsed$pango_lineage,
                   host = animal_parsed$host,
                   cluster = animal_parsed$cluster,
                   clade = animal_parsed$clade,
                   location = animal_parsed$location,
                   allele = animal_parsed$allele_annot)
  
  # Plot tree
  p <- ggtree(audacity_filt, aes(color = host)) %<+% dd +
    geom_tippoint(aes(shape = host)) +
    labs(shape = "Host", color = "host") +
    new_scale_color() +
    geom_fruit(geom = geom_tile, aes(fill = cluster), width=10, offset=0, alpha = 1) +
    scale_fill_jcolors("pal8", na.translate = F) +
    geom_tippoint(aes(color = allele), shape = 0) +
    scale_color_manual(values = "red", na.translate = F) +
    labs(fill = "Cluster", color = str_glue("Position {row$nucleotide_pos}")) +
    theme_tree2()
  
  p_zoom <- p +
    geom_tiplab(aes(size = host)) +
    scale_size_manual(values = c(0, 3))
  
  # ggutils::ggzoom(p_zoom)
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
  ggsave(str_glue("results/all_animals/founder_alleles/{prefix}_allele_{row$mutation_annot}_annotation.png"),
         plot = p,
         dpi = 600,
         width = 8,
         height = 10)
  
}


