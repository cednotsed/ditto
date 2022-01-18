rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(ape)
require(ggtree)
require(ggplot2)
require(ggtreeExtra)
require(ggnewscale)
require(jcolors)

host_species <- "Neovison vison"
prefix <- "mink"

# Load Audacity tree
audacity_tree <- read.tree("data/GISAID-hCoV-19-phylogeny-2021-11-16/global.tree")

# Load metadata
cluster_meta <- fread("results/cluster_annotation/deer_mink_parsed_clusters.csv")

all_meta <- fread("data/metadata/all_sequence_metadata_231121.tsv") %>% 
  rename_all(~ tolower(gsub(" ", "_", .))) %>%
  filter(accession_id %in% audacity_tree$tip.label) %>%
  separate(location, into = c("loc1", "loc2", "loc3"), sep = " / ") %>%
  mutate(location = paste0(loc1, " / ", loc2)) %>%
  as_tibble()

animal_meta <- all_meta %>%
  filter(host == host_species) %>%
  left_join(cluster_meta)

# Get subtree
to_drop <- audacity_tree$tip.label[!(audacity_tree$tip.label %in% animal_meta$accession_id)]
animal_tree <- drop.tip(audacity_tree, to_drop)

# Match metadata
meta.match <- animal_meta[match(animal_tree$tip.label, animal_meta$accession_id), ]

dd <- data.frame(Accession = animal_tree$tip.label,
                 lineage = meta.match$pango_lineage,
                 clade = meta.match$clade,
                 country = meta.match$loc2,
                 cluster = meta.match$cluster)

# Plot tree
p <- ggtree(animal_tree, color = "darkgrey") %<+% dd + 
  geom_tippoint(aes(color = cluster), size = 3) +
  scale_color_discrete(na.translate = F) +
  geom_fruit(geom = geom_tile, aes(fill = country), width=2, offset=0.1, alpha = 1) +
  scale_fill_jcolors("pal8", na.translate = F) +
  labs(color = "Cluster", fill = "Country")
  # geom_fruit(geom = geom_tile, aes(fill = lineage), width=1, offset=0.1, alpha = 1) +
  # scale_fill_discrete(na.translate = F) +
  # new_scale_fill() +


p
ggsave(str_glue("results/all_animals/{prefix}_lineages.png"), 
       plot = p, 
       dpi = 300, 
       width = 10, 
       height = 10)

