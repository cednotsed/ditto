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
  date_df[1, "collection_date"] <- "2019-12-26"
  date_df[1, "accession_id"] <- "NC_045512.2"
  # dates <- decimal_date(as.Date(date_df$collection_date))
  dates <- date_df$collection_date
  names(dates) <- date_df$accession_id
  
  new_tree$edge.length <- new_tree$edge.length * 29903
  return(list(new_tree = new_tree, dates = dates))
}

prefix <- "mink.Denmark.n3219"
prefix <- "mink.Netherlands.n2213"

animal_aln <- read.FASTA(str_glue("data/alignments/human_animal_subsets/{prefix}.audacity_only.v8_masked.aln.fasta"))
animal_tree <- read.tree(str_glue("data/trees/human_animal_subsets/{prefix}.audacity_only.v8_masked.tree"))

# Parse tip labels
# animal_tree <- parse_tip_labels(animal_tree)

# Scale tree
parsed_dat <- parse_tree(animal_tree)
animal_scaled_tree <- parsed_dat[["new_tree"]]

# Drop tips with incomplete dates
n_tips <- Ntip(animal_tree)
dates <- parsed_dat[["dates"]]
to_drop <- names(dates)[is.na(dates)]
dates_filt <- dates[!is.na(dates)]
animal_tree <- drop.tip(animal_tree, to_drop)
animal_scaled_tree <- drop.tip(animal_scaled_tree, to_drop)

print(str_glue("Original n_tips: {n_tips}
Tips dropped due to incomplete dates: {length(to_drop)}
Final n_tips: {Ntip(animal_tree)}"))

# fwrite(tibble(accession_id = names(dates_filt), dates = dates_filt),
#        str_glue("data/metadata/human_animal_subsets/{prefix}.date_file.tsv"),
#        sep = "\t",
#        col.names = F)

png(str_glue("results/human_animal_subsets/dating_out/{prefix}.roottotip.png"))
r <- roottotip(animal_scaled_tree, dates_filt)
dev.off()


# Date tree
animal_dated <- dater(animal_tree, dates_filt, 29903,
                      quiet = F,
                      clock = "uncorrelated",
                      maxit = 1000,
                      parallel_foreach = T)

parboot(animal_dated, ncpu = 10, parallel_foreach = T, quiet = F)

# 
# animal_bactdated <- bactdate(animal_scaled_tree, 
#                              dates_filt,
#                              showProgress = T)

# plot(animal_dated,'treeCI')


