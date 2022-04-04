rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(ape)
require(ggtree)
require(ggtreeExtra)

# Load Audacity tree
audacity_tree <- read.tree("data/GISAID-hCoV-19-phylogeny-2022-02-21/global.tree")

# All human and animal accessions
meta <- fread("data/metadata/human_animal_subsets/V6/all_animals.n16911.csv")

# Cluster meta
cluster_meta <- fread("results/cluster_annotation/deer_mink_parsed_clusters.csv")

meta_filt <- meta %>%
  left_join(cluster_meta)

# Subset tree
to_drop <- audacity_tree$tip.label[!(audacity_tree$tip.label %in% meta_filt$accession_id)]
to_drop <- to_drop[to_drop != "EPI_ISL_402124"]
audacity_filt <- drop.tip(audacity_tree, to_drop)

# Match metadata to tips
meta.match <- meta_filt[match(audacity_filt$tip.label, meta_filt$accession_id), ]

all(audacity_filt$tip.label == meta.match$accession_id)

# Remove humans from host metadata
meta.match <- meta.match %>%
  mutate(host = ifelse(grepl("Human", ignore.case = T, host), NA, host),
         variant = case_when(grepl("Alpha", variant) ~ "Alpha",
                             grepl("Beta", variant) ~ "Beta",
                             grepl("Gamma", variant) ~ "Gamma",
                             grepl("Delta", variant) ~ "Delta",
                             grepl("Omicron", variant) ~ "Omicron",
                             TRUE ~ as.character(NA)))

meta_df <- data.frame(Accession = meta.match$accession_id,
                      host = meta.match$host,
                      voc = meta.match$variant,
                      cluster = meta.match$cluster)

# Plot tree
p <- ggtree(audacity_filt, size = 0.001, color = "darkgrey") %<+% meta_df +
  geom_tippoint(aes(hjust = 0.5, color = host), alpha = 1, size = 3) +
  scale_color_manual(values = c("cyan4", "darkolivegreen4", "darkgoldenrod3", 
                                "#2e4057", "#d1495b", "cornflowerblue"), 
                     na.value = NA,
                     na.translate = F) +
  geom_fruit(geom = geom_tile, 
             aes(fill = voc),
             width = 10) +
  scale_fill_discrete(na.translate = F) +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 15)) +
  # guides(colour = guide_legend(override.aes = list(size = 5)),
  #        fill = guide_legend(override.aes = list(size = 5))) +
  labs(color = "Host", fill = "VoC")

p

ggsave("results/global_tree_all_animals.png", 
       plot = p, 
       dpi = 300, 
       width = 10, 
       height = 10)

p_no_legend <- p + theme(legend.position = "none")

ggsave("results/global_tree_all_animals.nolegend.png", 
       plot = p_no_legend, 
       dpi = 300, 
       width = 10, 
       height = 10)

