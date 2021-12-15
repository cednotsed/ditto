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

prefix <- "mink_n789_201121"
# prefix <- "deer_n73_201121"
tree <- read.tree(str_glue("data/trees/all_animals/{prefix}.audacity_only.v8_masked.tree"))

animal_meta <- fread("data/metadata/combined_metadata_201121.audacity_only.csv") %>%
  separate(location, into = c("loc1", "loc2"), sep = " / ") %>%
  mutate(location = paste0(loc1, " / ", loc2)) %>%
  select(accession_id, host, clade, lineage, location, cluster) %>%
  rename(Accession = accession_id) %>%
  as_tibble()

# Match metadata
animal_meta <- animal_meta[match(tree$tip.label, animal_meta$Accession), ]
# animal_meta[is.na(animal_meta$Accession), ] <- "root"

dd <- data.frame(Accession = tree$tip.label,
                 lineage = animal_meta$lineage,
                 clade = animal_meta$clade,
                 location = animal_meta$location,
                 cluster = animal_meta$cluster)

# Plot tree
p <- ggtree(tree, color = "darkgrey")
p <- p %<+% dd
p <- p + 
  geom_tippoint(aes(color = cluster), size = 3) +
  geom_fruit(geom = geom_tile, aes(fill = location), width=0.0001, offset=0.1, alpha = 1) +
  scale_fill_jcolors("pal8", na.translate = F) +
  scale_color_discrete(na.translate = F) +
  labs(fill = "Location", color = "Cluster") +
  theme_tree2()
 
p
ggsave(str_glue("results/all_animals/{prefix}_lineages.png"), 
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
