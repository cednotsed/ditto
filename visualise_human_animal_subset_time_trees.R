rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(ape)
require(ggtree)
require(ggplot2)
require(ggtreeExtra)
require(ggnewscale)

# prefix <- "deer.USA.n6670"
prefix <- "mink.Netherlands.n3750"
prefix <- "mink.Denmark.n10512"
# prefix <- "mink.USA.n6726"
animal_name <- str_split(prefix, "\\.")[[1]][1]
meta <- fread(str_glue("data/metadata/human_animal_subsets/{prefix}.unambiguous.dedup.csv"))

node_meta <- fread(str_glue("results/human_animal_subsets/dating_out/{prefix}.unambiguous.dedup/dates.tsv"),
                   skip = 1) %>%
  rename_all(~ tolower(gsub(" ", "_", .))) %>%
  rename(node = "#node") %>%
  filter(numeric_date != "--")
root_meta <- node_meta %>%
  filter(node == "NODE_0000001")
root_date <- as.numeric(root_meta$numeric_date)

tree <- read.nexus(str_glue("results/human_animal_subsets/dating_out/{prefix}.unambiguous.dedup/timetree.nexus"))
Ntip(tree) == nrow(meta)
sum(meta$accession_id %in% tree$tip.label)

# Match metadata
meta.match <- meta[match(tree$tip.label, meta$accession_id), ]

dd <- data.frame(Accession = tree$tip.label,
                 host = meta.match$host,
                 lineage = meta.match$pango_lineage,
                 clade = meta.match$clade,
                 location = meta.match$location,
                 variant = meta.match$variant,
                 cluster = meta.match$cluster)

dd[dd$host == "Human", "host"] <- NA
dd[dd$cluster == "", "cluster"] <- NA

# # Annotate big clusters
# freqs <- table(dd$cluster)
# no_annot <- names(freqs[freqs < 20])
# dd[dd$cluster %in% no_annot, "cluster"] <- NA

# # timescale
# ts <- seq(as.numeric(min(node_meta$numeric_date)),
#          as.numeric(max(node_meta$numeric_date)), by = 1/12)
# ts_text <- as.yearmon(ts, "%Y %m")

# Plot tree
options(ignore.negative.edge=TRUE)

p <- ggtree(tree, color = "darkgrey", root.position = root_date, right = T) %<+% dd +
  theme_tree2() + 
  theme(panel.grid.major.x = element_line(colour = "grey50"),
        panel.grid.minor.x = element_line()) +
  geom_tippoint(aes(color = cluster), size = 3) +
  labs(color = "Host") +
  scale_color_discrete(na.translate = F) +
  labs(color = "Mink cluster")
  # scale_size_manual(values = c(3), na.translate = F, guide = "none") +
  # theme(legend.position = "top")
p

ggsave(str_glue("results/human_animal_subsets/dating_out/{prefix}.unambiguous.dedup/{prefix}.timetree.png"),
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
