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

# prefix <- "deer.USA.n6670"
prefix <- "mink.Netherlands.n3750"
prefix <- "mink.Denmark.n10512"
prefix <- "all_mink.n1487"
# prefix <- "mink.USA.n6726"
animal_name <- str_split(prefix, "\\.")[[1]][1]
meta <- fread(str_glue("data/metadata/human_animal_subsets/V5/{prefix}.unambiguous.dedup.csv"))

tree <- read.nexus(str_glue("results/human_animal_subsets/V5/dating_out/{prefix}.unambiguous.dedup/divergence_tree.nexus"))
Ntip(tree) == nrow(meta)
sum(meta$accession_id %in% tree$tip.label)

# Match metadata
meta.match <- meta[match(tree$tip.label, meta$accession_id), ]

dd <- data.frame(Accession = tree$tip.label,
                 host = meta.match$host,
                 lineage = meta.match$pango_lineage,
                 clade = meta.match$clade,
                 location = meta.match$loc2,
                 variant = meta.match$variant,
                 cluster = meta.match$cluster)

# dd[dd$host == "Human", "host"] <- NA
dd[dd$cluster == "", "cluster"] <- NA

p <- ggtree(tree, color = "darkgrey") %<+% dd +
  geom_fruit(geom = geom_tile, aes(fill = location), width=0.0018, offset=-0.6, alpha = 0.3) +
  labs(fill = "Location") +
  scale_fill_jcolors("pal8", na.translate = F) +
  geom_tippoint(aes(color = host), size = 3) +
  labs(color = "Host") +
  scale_color_discrete(na.translate = F) +
  labs(color = "")
# scale_size_manual(values = c(3), na.translate = F, guide = "none") +
# theme(legend.position = "top")
p

ggsave(str_glue("results/human_animal_subsets/V5/dating_out/{prefix}.unambiguous.dedup/{prefix}.divergence.png"),
       plot = p,
       dpi = 300,
       width = 8,
       height = 8)

# ggutils::ggzoom(p)
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
