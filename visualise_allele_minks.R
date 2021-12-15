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
tree <- read.tree(str_glue("data/trees/{prefix}.audacity_only.v8_masked.tree"))

animal_meta <- fread("data/metadata/combined_metadata_201121.audacity_only.tsv") %>%
  separate(location, into = c("loc1", "loc2"), sep = " / ") %>%
  mutate(location = paste0(loc2)) %>%
  select(accession_id, host, clade, lineage, location) %>%
  rename(Accession = accession_id)

# Get allele
allele_pos <- 12795
aln <- read.dna(str_glue("data/alignments/all_animals/{prefix}.audacity_only.v8_masked.aln.fasta"),
                format = "fasta",
                as.matrix = T)
aln_df <- tibble(Accession = rownames(aln), allele = toupper(as.vector(as.character(aln[, allele_pos]))))

animal_meta <- animal_meta %>%
  left_join(aln_df)

# Match metadata
animal_meta <- animal_meta[match(tree$tip.label, animal_meta$Accession), ]
# animal_meta[is.na(animal_meta$Accession), ] <- "root"

dd <- data.frame(Accession = tree$tip.label,
                 lineage = animal_meta$lineage,
                 clade = animal_meta$clade,
                 location = animal_meta$location,
                 allele = animal_meta$allele)

# Plot tree
p <- ggtree(tree, color = "darkgrey")
p <- p %<+% dd
p <- p + 
  geom_fruit(geom = geom_tile, aes(fill = location), width=0.0018, offset=-0.6, alpha = 0.3) +
  labs(fill = "Location") +
  scale_fill_jcolors("pal8", na.translate = F) +
  new_scale_fill() +
  geom_fruit(geom = geom_tile, aes(fill = lineage), width=0.0001, offset=0.5, alpha = 1) +
  geom_tippoint(aes(color = allele), size = 3) +
  scale_color_discrete(na.translate = F) +
  labs(fill = "Lineage", color = str_glue("Alllele at pos. {allele_pos}")) +
  theme_tree2()

p
ggsave(str_glue("results/all_animals/{prefix}_allele_{allele_pos}_annotation.png"), 
       plot = p, 
       dpi = 300, 
       width = 10, 
       height = 10)

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
