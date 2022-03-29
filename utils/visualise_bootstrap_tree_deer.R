rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(ape)
require(ggtree)
require(ggplot2)
require(ggtreeExtra)
require(ggnewscale)
require(treeio)
require(randomcoloR)

# Root to WIV04
unrooted_path <- "data/trees/all_animals/deer_only.n96.audacity_only.v8_masked.treefile"
tree <- read.iqtree(unrooted_path)
tree <- ape::root(tree, "EPI_ISL_402124", resolve.root = TRUE)
animal_meta <- fread(str_glue("data/metadata/all_animals/deer_only.n96.csv"))

cluster_meta <- fread(str_glue("results/cluster_annotation/deer_mink_raw_clusters.csv")) %>%
  select(-n, -cluster_num, -host_name)

animal_meta <- animal_meta %>%
  left_join(cluster_meta) %>%
  mutate(cluster = ifelse(accession_id == "EPI_ISL_402124", "WIV04_Root", cluster))

# Match metadata
tips <- as.phylo(tree)$tip.label
animal_parsed <- animal_meta[match(tips, animal_meta$accession_id), ]

dd <- data.frame(Accession = tips,
                 lineage = animal_parsed$pango_lineage,
                 clade = animal_parsed$clade,
                 cluster = animal_parsed$cluster,
                 location = animal_parsed$loc2)

# Get color pal
colordict <- as.list(distinctColorPalette(length(unique(dd$cluster))))

# Plot tree
p <- ggtree(tree, aes(color = cluster), branch.length = "none") %<+% dd +
  geom_tippoint(aes(color = cluster)) +
  geom_tiplab(aes(label = cluster, color = cluster)) +
  scale_color_manual(values = colordict) +
  labs(shape = "Lineage", color = "Lineage") +
  geom_text(aes(x = branch, label=UFboot), color = "blue") +
  geom_text(aes(x = branch, label=SH_aLRT), color = "red", hjust = 2) +
  theme_tree2()


ggsave(str_glue("results/all_animals/all_deer_bootstrap_tree.sh-aLRT.pdf"),
       plot = p,
       dpi = 600,
       width = 30,
       height = 30)


