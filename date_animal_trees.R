rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(ape)
require(tidyverse)
require(BactDating)
require(foreach)
require(treedater)
require(lubridate)
require(data.table)

# Load metadata
audacity <- fread("data/GISAID-hCoV-19-phylogeny-2021-11-16/metadata.csv")

all_meta <- fread("data/metadata/all_sequence_metadata_231121.tsv", nThread = 10) %>% 
  rename_all(~ tolower(gsub(" ", "_", .))) %>%
  filter(accession_id %in% audacity$accession_id) %>%
  separate(location, into = c("loc1", "loc2", "loc3"), sep = " / ") %>%
  mutate(location = paste0(loc1, " / ", loc2))

# Parse tree names
parse_tip_labels <- function(tree) {
  tip_names <- separate(tibble(accession = tree$tip.label),
                        accession, 
                        into = c(NA, "accession", NA), sep = "\\|")$accession
  tip_names[1] <- tree$tip.label[1]
  tree$tip.label <- tip_names
  return(tree)
}

parse_tree <- function(tree) {
  new_tree <- tree
  date_df <- all_meta[match(new_tree$tip.label, all_meta$accession_id), ]
  # date_df[1, "collection_date"] <- "2019-12-26"
  # date_df[1, "accession_id"] <- "NC_045512.2"
  dates <- decimal_date(as.Date(date_df$collection_date))
  names(dates) <- date_df$accession_id
  
  new_tree$edge.length <- new_tree$edge.length * 29903
  return(list(new_tree = new_tree, dates = dates))
}

# animal_aln <- read.FASTA("data/alignments/matched/deer_73_201121.audacity_only.v8_masked.aln.fasta")
# animal_tree <- read.tree("data/trees/matched/deer_73_201121.audacity_only.v8_masked.tree")
# human_aln <- read.FASTA("data/alignments/matched/deer_matched_humans_72_201121.audacity_only.v8_masked.aln.fasta")
# human_tree <- read.tree("data/trees/matched/deer_matched_humans_72_201121.audacity_only.v8_masked.tree")
# animal_aln <- read.FASTA("data/alignments/matched/mink_723_201121.audacity_only.v8_masked.aln.fasta")
# animal_tree <- read.tree("data/trees/matched/mink_723_201121.audacity_only.v8_masked.tree")
# human_aln <- read.FASTA("data/alignments/matched/mink_matched_humans_707_201121.audacity_only.v8_masked.aln.fasta")
# human_tree <- read.tree("data/trees/matched/mink_matched_humans_707_201121.audacity_only.v8_masked.tree")
prefixes <- list.files("data/genomes/mini_animal_trees/")
prefixes <- gsub(".audacity_only.fasta", "", prefixes)
# prefixes <- prefixes[grepl("N2|US2|US3|US4|", prefixes)]
prefix <- prefixes[8]
prefixes[4]

foreach(prefix = prefixes) %do% {
  animal_aln <- read.FASTA(str_glue("data/alignments/mini_animal_trees/{prefix}.audacity_only.v8_masked.aln.fasta"))
  animal_tree <- read.tree(str_glue("data/trees/mini_animal_trees/{prefix}.audacity_only.v8_masked.tree"))
  
  # Drop wuhan-hu-1
  animal_tree <- drop.tip(animal_tree, "NC_045512.2")
  n_tips <- Ntip(animal_tree)
  
  # animal_tree <- parse_tip_labels(animal_tree)
  # human_tree <- parse_tip_labels(human_tree)
  parsed_dat <- parse_tree(animal_tree)
  animal_scaled_tree <- parsed_dat[["new_tree"]]
  dates <- parsed_dat[["dates"]]
  
  # Drop tips with incomplete dates
  to_drop <- names(dates)[is.na(dates)]
  dates_filt <- dates[!is.na(dates)]
  animal_tree <- drop.tip(animal_tree, to_drop)
  animal_scaled_tree <- drop.tip(animal_scaled_tree, to_drop)
  
  print(str_glue("Original n_tips: {n_tips}
  Tips dropped due to incomplete dates: {length(to_drop)}
  Final n_tips: {Ntip(animal_tree)}"))
  
  r <- roottotip(animal_scaled_tree, dates_filt)
  
  png(str_glue("results/mini_animal_trees/dating_out/{prefix}.roottotip.png"))
    r <- roottotip(animal_scaled_tree, dates_filt)
  dev.off()
  
}

## Date particular clusters
clusters <- c("D1", "N1", "N3")
trees <- list.files("data/trees/mini_animal_trees/")
trees <- trees[grepl(paste0(clusters, collapse = "|"), trees)]

animal_tree <- read.tree(str_glue("data/trees/mini_animal_trees/{trees[2]}"))
# animal_tree <- read.tree(str_glue("data/trees/human_animal_subsets/mink.Denmark.n3219.audacity_only.v8_masked.tree"))
# animal_tree <- parse_tip_labels(animal_tree)
ggtree(animal_tree)
animal_tree <- drop.tip(animal_tree, "NC_045512.2")
animal_tree <- drop.tip(animal_tree, "EPI_ISL_683041")
ggtree(animal_tree) + geom_tiplab()
ggsave("results/test_tree.png", height = 15, width = 15)
parsed_dat <- parse_tree(animal_tree)
animal_scaled_tree <- parsed_dat[["new_tree"]]
dates <- parsed_dat[["dates"]]

# Drop tips with incomplete dates
to_drop <- names(dates)[is.na(dates)]
dates_filt <- dates[!is.na(dates)]
animal_tree <- drop.tip(animal_tree, to_drop)
animal_scaled_tree <- drop.tip(animal_scaled_tree, to_drop)

# Date tree
animal_dated <- dater(animal_tree, dates_filt, 
                      s = 29903,
                      quiet = F,
                      clock = "strict",
                      temporalConstraints = T,
                      maxit = 1000)
animal_bactdated <- bactdate(unroot(animal_scaled_tree), dates_filt)

plot(animal_bactdated,'treeCI')

ggtree(animal_bactdated$tree, mrsd=date_decimal(animal_bactdated$rootdate)) +
  theme_tree2()

  # group_by(pango_lineage) %>%
  # summarise(n = n())
