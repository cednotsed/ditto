setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(ape)
require(lubridate)
require(BactDating)
require(ggtree)

# Parse alignment names
aln <- read.FASTA("data/alignments/mink_paper/ClusterA_mink.audacity_only.v8_masked.aln.fasta")
names(aln) <- separate(tibble(x = names(aln)), x, into = c(NA, "x"), sep = "\\|")$x
aln <- aln[2:length(aln)]
write.FASTA(aln, "data/alignments/mink_paper/parsed_ClusterA_mink.audacity_only.v8_masked.aln.fasta")

tree <- read.tree("data/trees/mink_paper/ClusterA_mink.audacity_only.v8_masked.tree")
meta <- fread("data/metadata/mink_paper/mink_paper_meta.csv") %>%
  filter(host == "mink", Cluster =="ClusterA") %>%
  rename(accession_id = GISAID_AccessionID)

tree <- drop.tip(tree, "NC_045512.2")

# Parse tree tips
tip_labels <- tree$tip.label
tip_labels <- separate(tibble(x = tree$tip.label), x, into = c(NA, "x", NA), sep = "\\|")$x
tree$tip.label <- tip_labels
write.tree(tree, "data/trees/mink_paper/parsed_ClusterA_mink.audacity_only.v8_masked.tree")

# Parse tree edges
tree_scaled <- tree
tree_scaled$edge.length <- tree_scaled$edge.length * 29903

# Check duplicates
meta$accession_id[!(meta$accession_id %in% tree$tip.label)]
meta %>%
  filter(accession_id %in% c("EPI_ISL_632385", "EPI_ISL_632479"))

meta_filt <- meta %>% 
  distinct(accession_id, .keep_all = T) %>%
  mutate(collection_date = decimal_date(as.Date(collectDate, "%d/%m/%y"))) %>%
  select(accession_id, collection_date)

dates <- tibble(accession_id = tree$tip.label) %>%
  left_join(meta_filt)

to_save <- dates %>% rename(accession = accession_id)
# to_save <- tibble(accession = "NC_045512.2", 
#                   collection_date = decimal_date(as.Date("2019-12-26"))) %>%
#   bind_rows(to_save)
# fwrite(to_save, "data/metadata/mink_paper/parsed_mink_paper_dates.tsv", sep = "\t", 
#        col.names = T)

dates <- dates$collection_date
# Match metadata
roottotip(tree_scaled, dates)

bactdate(tree_scaled, dates, 
         model = "strictgamma",
         updateRoot = T,
         useCoalPrior = T)


dd <- data.frame(Accession = meta_filt$accession_id,
                 dates = dates)

ggtree(tree) %<+% dd +
  geom_tippoint(aes(color = dates))


tt <- read.nexus("results/mink_paper/cluster_A_clock/timetree.nexus")
ggtree(tt, mrsd = date_decimal(max(dates$collection_date))) +
  theme_tree2() 

ggtree(tt, mrsd = date_decimal(2020.99)) +
  theme_tree2()

plot(tt)

