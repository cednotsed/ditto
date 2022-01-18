rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(ape)
require(jcolors)
require(ggpubr)
require(foreach)
require(Biostrings)
require(seqinr)
require(see)

parse <- function(fasta, prefix) {
  dn <- dinucleotideFrequency(fasta, step=1,
                              as.prob=T, as.matrix=FALSE,
                              fast.moving.side="right", with.labels=TRUE)
  dn <- as_tibble(dn) %>%
    add_column(host = prefix)
  
  return(dn)
}

# Files
human_background <- "V5"
host_name <- "Mink"
aln_path <- "all_mink.n1487.audacity_only.v8_masked.fasta"
meta_file <- "all_mink.n1487.csv"
host_species <- "Neovison vison"

human_background <- "V5"
host_name <- "Deer"
aln_path <- "all_deer.n145.audacity_only.v8_masked.fasta"
meta_file <- "all_deer.n145.csv"
host_species <- "Odocoileus virginianus"

# Load alignment
aln <- readDNAStringSet(file = str_glue("data/alignments/human_animal_subsets/{human_background}/{aln_path}"))

# Remove anthroponoses
ant_df <- read.csv("data/metadata/netherlands_humans_anthroponoses.txt", header = F)
aln <- aln[!(names(aln) %in% ant_df$V1)]

# Load metadata
meta <- fread(str_glue("data/metadata/human_animal_subsets/{human_background}/{meta_file}")) %>%
  rename_all(~ tolower(gsub(" ", "_", .))) %>%
  select(accession_id, host, clade, pango_lineage, location, collection_date) %>%
  separate(location, into = c("loc1", "loc2", "loc3"), sep = " / ") %>%
  mutate(location = paste0(loc1, " / ", loc2, " / ", loc3)) %>%
  as_tibble()

# Get host alignments
human_meta <- meta %>% filter(host == "Human")
human <- aln[names(aln) %in% human_meta$accession_id]
animal_meta <- meta %>% filter(host == host_species)
animal <- aln[names(aln) %in% animal_meta$accession_id]


animal_dn <- parse(animal, host_species)
human_dn <- parse(human, "Human")

plot_df <- bind_rows(animal_dn, human_dn)

pca <- prcomp(plot_df %>% select(-host), center = T, scale. = T, retx = T)
eig <- pca$sdev^2
prop <- round(eig / sum(eig) * 100, 1)

plt <- as_tibble(pca$x) %>%
  bind_cols(tibble(host = plot_df$host)) %>%
  mutate(host = factor(host, c("Human", 
                               "Odocoileus virginianus", 
                               "Neovison vison"))) %>%
  ggplot(aes(x = PC1, y = PC2, color = host)) +
  theme_bw() +
  geom_point() +
  labs(x = str_glue("PC1 ({prop[1]}%)"), y = str_glue("PC2 ({prop[2]}%)")) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

ggsave(str_glue("results/mutation_bias/{host_name}_{human_background}_PCA_dinucleotide.png"), 
       plot = plt, 
       dpi = 300, 
       width = 6, 
       height = 4)

