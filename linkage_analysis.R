rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(ape)
require(foreach)
require(doParallel)
registerDoParallel(cores = 10)

# Load metadata
audacity <- fread("data/GISAID-hCoV-19-phylogeny-2021-11-16/metadata.csv")
meta <- fread("data/metadata/all_sequence_metadata_231121.tsv", nThread = 10) %>%
  rename_all(~ tolower(gsub(" ", "_", .))) %>%
  filter(accession_id %in% audacity$accession_id) %>%
  separate(location, into = c("loc1", "loc2", "loc3"), sep = " / ") %>%
  mutate(location = paste0(loc1, " / ", loc2))

host_name <- "mink"

mut_df <- fread(str_glue("results/{host_name}_homoplasy_alele_frequency_{human_background}.csv")) %>%
  mutate(mutation_annot = ifelse(mutation_type == "NS", 
                                 paste0(protein_name, "_", variant_name),
                                 paste0(protein_name, "_", mutation_name, "(syn)")))

cluster_meta <- fread("results/cluster_annotation/deer_mink_parsed_clusters.csv")
cluster_meta
table(cluster_meta$host)
aln_file <- list.files("data/alignments/all_animals/", pattern = host_name, full.names = T)
aln_file <- aln_file[grepl("audacity_only.v8_masked.aln.fasta", aln_file)]

aln <- read.dna(aln_file, 
                format = "fasta", 
                as.matrix = T,
                as.character = T)
aln <- aln[rownames(aln) != "NC_045512.2", ]
aln_filt <- as_tibble(aln)
aln_filt <- aln_filt[, mut_df$nucleotide_pos]
aln_filt <- aln_filt %>% 
  mutate_all(toupper) %>%
  add_column(accession_id = rownames(aln)) %>%
  left_join(cluster_meta)

# host_species <- "Neovison vison"

aln_mat <- as.matrix(aln)

morsels <- foreach(i = seq(nrow(mut_df))) %do% {
  row <- mut_df[i, ]
  aln_col <- toupper(as.vector(aln_mat[, row$nucleotide_pos]))
  aln_col <- ifelse(aln_col == row$var_nuc, 1, 0)
  temp <- tibble(pos = aln_col)
  colnames(temp) <- row$mutation_annot
  temp
}

corr_mat <- bind_cols(morsels)

png(str_glue("results/all_animals/linkage_analysis/{host_name}_correlation_heatmap.png"), 
    res = 300,
    units = "cm",
    width = 15, 
    height = 15)
corrplot::corrplot(cor(corr_mat, method = ), 
                   type="lower",
                   diag = F, 
                   tl.col = "black")
dev.off()

r_df <- cor(corr_mat)
mut1_name <- rownames(r_df)

corr_df <- as_tibble(r_df) %>%
  mutate(mut1 = mut1_name, .before = 1) %>%
  pivot_longer(!mut1, names_to = "mut2", values_to = "r2") %>%
  filter(mut1 != mut2)

morsels <- foreach(i = seq(nrow(corr_df))) %do% {
  row <- corr_df[i, ]
  mut1_hookup <- mut_df %>% filter(mutation_annot == row$mut1)
  mut2_hookup <- mut_df %>% filter(mutation_annot == row$mut2)
  
  test <- tibble(accession_id = rownames(aln_mat), 
                 mut1 = toupper(as.vector(aln_mat[, mut1_hookup$nucleotide_pos])), 
                 mut2 = toupper(as.vector(aln_mat[, mut2_hookup$nucleotide_pos])))
  morsel <- test %>%
    left_join(cluster_meta) %>% 
    mutate(cooccur = ifelse(mut1 == mut1_hookup$var_nuc & mut2 == mut2_hookup$var_nuc, 1, 0)) %>%
    group_by(cluster) %>%
    summarise(prop = sum(cooccur) / n()) %>%
    add_column(mut_pair = paste0(row$mut1, "|", row$mut2))
  
  morsel
}

plot_df <- bind_rows(morsels)
View(clust_df)

# At least some coccurence
max_df <- plot_df %>%
  group_by(mut_pair) %>%
  summarise(max = max(prop)) %>%
  filter(max > 0)

# At least 3 clusters
n_df <- plot_df %>%
  filter(prop > 0) %>%
  group_by(mut_pair) %>%
  summarise(n_clusters = n()) %>%
  filter(n_clusters > 2)

to_keep <- intersect(max_df$mut_pair, n_df$mut_pair)

plot_df %>%
  filter(mut_pair %in% to_keep) %>%
  ggplot(aes(y = mut_pair, x = cluster, fill = prop)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

human <- meta %>% 
  filter(host == "Human")

animal <- meta %>% 
  filter(host == host_species)

human$aa_substitutions

mutations <- c("NS3_H182Y")
mutation
mutation <- "NS3_G172V"

human %>%
  summarise(prop = sum(grepl(mutation, aa_substitutions)) / nrow(human))
