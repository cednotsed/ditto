setwd("../Desktop/git_repos/ditto/")
setwd("/Users/Cedric/Desktop/git_repos/ditto/")
require(doParallel)
require(foreach)
require(tidyverse)
require(data.table)
require(ape)
require(ggtree)
require(ggtreeExtra)
registerDoParallel(cores=8)

audacity_tree <- read.tree("data/GISAID-hCoV-19-phylogeny-2021-11-16/global.tree")
# animal_meta <- data.frame(fread("data/metadata/combined_metadata_201121.audacity_only.tsv")) %>%
#   select(accession_id, host, clade, lineage) %>%
#   rename(Accession = accession_id)

all_meta <- fread("data/metadata/all_sequence_metadata_231121.tsv") %>%
  rename_all(~ tolower(gsub(" ", "_", .))) %>%
  select(accession_id, host, clade, pango_lineage) %>%
  rename(Accession = accession_id)
  filter(Accession %in% audacity_tree$tip.label) # Keep only accessions in audacity tree
# audacity_tree$tip.label[!(audacity_tree$tip.label %in% all_meta$Accession)]

# Get animal accessions
animals_to_keep <- c("Odocoileus virginianus", "Felis catus", "Neovison vison")
animal_meta <- all_meta %>%
  filter(host %in% animals_to_keep)

all_human_meta <- all_meta %>%
  filter(grepl("human", ignore.case = T, host))

# Subsample tree while ensuring all lineages are represented
set.seed(66)
lineages <- unique(all_meta$pango_lineage) 
n_samples <- 100
length(lineages)  # There are 1508 lineages

morsels <- foreach(lineage = lineages, .packages = "tidyverse") %dopar% {
  lineage_df <- all_human_meta %>% 
    filter(pango_lineage == lineage)
  
  n_accessions <- nrow(lineage_df)
  
  if (n_accessions > n_samples) {
    return(lineage_df %>% sample_n(n_samples))
  } else {
    return(lineage_df %>% sample_n(n_accessions))
  }
}

meta_to_keep <- bind_rows(morsels) %>% bind_rows(animal_meta)

# Check that all lineages are represented
length(unique(meta_to_keep$pango_lineage))

# Drop tips
# to_remove <- sample(tree$tip.label, length(tree$tip.label) * (1 - 1/10000))
to_remove <- tree$tip.label
to_remove <- to_remove[!(to_remove %in% c("EPI_ISL_406801", meta_to_keep$Accession))]
filt_tree <- drop.tip(audacity_tree, to_remove)
rooted <- root(filt_tree, outgroup = "EPI_ISL_406801", resolve.root = T)

# Match metadata to tips
tree_meta <- meta_to_keep[match(rooted$tip.label, meta_to_keep$Accession), ]

# Remove humans from host metadata
tree_meta <- tree_meta %>%
  mutate(host = ifelse(grepl("human", ignore.case = T, host), NA, host))

# Set metadata for root
# tree_meta[tree_meta$Accession == "EPI_ISL_406801", colnames(tree_meta) != "Accession"] <- "root"

# Plot tree
p <- ggtree(rooted, layout="circular", size = 0.001, color = "darkgrey") %<+% tree_meta +
  geom_tippoint(aes(hjust = 0.5, color = host), alpha = 1, size = 1) +
  scale_color_manual(values = c("#2e4057", "#66a182", "#d1495b"), 
                     na.value = NA,
                     na.translate = F) +
  geom_fruit(geom = geom_tile, 
             aes(fill = clade),
             width = 10) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10)) +
  guides(colour = guide_legend(override.aes = list(size = 10)),
         fill = guide_legend(override.aes = list(size = 5))) +
  labs(color = "Host", fill = "Clade")

p
ggsave("results/test_tree_ggtree_all_lineages.n100.png", 
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
