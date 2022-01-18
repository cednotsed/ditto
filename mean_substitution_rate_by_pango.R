rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(ape)
require(ggplot2)
require(see)
require(foreach)

prefixes <- list.files("results/human_animal_subsets/V5/dating_out/")
ant_df <- read.csv("data/metadata/netherlands_humans_anthroponoses.txt", header = F)
prefixes <- prefixes[!grepl("deer|all_|USA|.png", prefixes)]
human_background <- "V5"
cluster_meta <- fread("results/cluster_annotation/deer_mink_parsed_clusters.csv") %>%
  select(accession_id, cluster)

morsels <- foreach (prefix = prefixes) %do% {
  country_name <- str_split(prefix, "\\.")[[1]][2]
  time_tree <- read.nexus(str_glue("results/human_animal_subsets/{human_background}/dating_out/{prefix}/timetree.nexus"))
  div_tree <- read.nexus(str_glue("results/human_animal_subsets/{human_background}/dating_out/{prefix}/divergence_tree.nexus"))
  
  meta <- fread(str_glue("data/metadata/human_animal_subsets/{human_background}/{prefix}.csv")) %>%
    select(-cluster) %>%
    filter(accession_id %in% time_tree$tip.label) %>%
    left_join(cluster_meta)
  
  # Divide time tree by divergence tree
  rate_tree <- div_tree
  rate_tree$edge.length <- rate_tree$edge.length / time_tree$edge.length
  
  # Drop anthroponotic tips
  rate_tree <- drop.tip(rate_tree, ant_df$V1)
  
  n <- length(rate_tree$tip.label)
  rate_df <- tibble(accession_id = rate_tree$tip.label,
                    terminal_rate = rate_tree$edge.length[sapply(1:n, function(x,y) which(y == x),
                                                                 y = rate_tree$edge[,2])])
  rate_df %>% 
    left_join(meta) %>%
    add_column(country = country_name) %>%
    select(host, terminal_rate, cluster, country, pango_lineage)
}

plot_df <- bind_rows(morsels)

plot_df %>% 
  group_by(host, pango_lineage) %>%
  summarise(n = n()) %>%
  right_join(plot_df, by = c("host", "pango_lineage")) %>%
  ggplot(aes(x = host, y = log(terminal_rate, base = 10), fill = host)) +
  facet_grid(rows = vars(pango_lineage), scales = "free", space = "free") +
  theme_bw() +
  geom_point(color = "darkgrey", position = position_jitter(width = 0.08), 
             size = 0.5, 
             alpha = 0.8) +
  geom_boxplot(position = position_nudge(x = -0.2, y = 0), 
               width = 0.2, 
               outlier.shape = NA,
               alpha = 0.3) +
  geom_violinhalf(position = position_nudge(x = 0.2, y = 0), alpha = 1) +
  geom_text(aes(x = host, y = -1.5, label = paste0("n = ", n))) +
  coord_flip() +
  labs(y = "lg(subst. rate)", x = "Host") +
  theme(legend.position = "none")

# ggsave(str_glue("results/human_animal_subsets/{human_background}/dating_out/mutations_rates_terminal.png"), 
#        width = 7, 
#        height = 5)

# Plot by cluster
plot_df %>% 
  group_by(cluster) %>%
  summarise(n = n()) %>%
  right_join(plot_df, by = c("cluster")) %>%
  filter(host == "Neovison vison") %>%
  ggplot(aes(x = cluster, y = log(terminal_rate, base = 10), fill = cluster)) +
  facet_grid(rows = vars(country), scales = "free", space = "free") +
  theme_bw() +
  geom_point(color = "darkgrey", position = position_jitter(width = 0.08), 
             size = 0.5, 
             alpha = 0.8) +
  geom_boxplot(position = position_nudge(x = -0.2, y = 0), 
               width = 0.2, 
               outlier.shape = NA,
               alpha = 0.3) +
  geom_violinhalf(position = position_nudge(x = 0.2, y = 0), alpha = 1) +
  geom_text(aes(x = cluster, y = -1.5, label = paste0("n = ", n))) +
  coord_flip() +
  labs(y = "lg(subst. rate)", x = "Host") +
  theme(legend.position = "none")

ggsave(str_glue("results/human_animal_subsets/{human_background}/dating_out/mutations_rates_terminal_by_cluster.png"), 
       width = 7, 
       height = 5)

#### ALL branches ###################
morsels2 <- foreach (prefix = prefixes) %do% {
  country_name <- str_split(prefix, "\\.")[[1]][2]
  time_tree <- read.nexus(str_glue("results/human_animal_subsets/{human_background}/dating_out/{prefix}/timetree.nexus"))
  div_tree <- read.nexus(str_glue("results/human_animal_subsets/{human_background}/dating_out/{prefix}/divergence_tree.nexus"))
  
  cluster_meta <- fread("results/cluster_annotation/deer_mink_parsed_clusters.csv") %>%
    select(accession_id, cluster)
  
  meta <- fread(str_glue("data/metadata/human_animal_subsets/{human_background}/{prefix}.csv")) %>%
    filter(accession_id %in% time_tree$tip.label) %>%
    left_join(cluster_meta)
  
  human <- meta %>% filter(host == "Human")
  animal <- meta %>% filter(host == "Neovison vison")
  
  # Divide time tree by divergence tree
  rate_tree <- div_tree
  rate_tree$edge.length <- rate_tree$edge.length / time_tree$edge.length
  
  # Drop anthroponotic tips
  rate_tree <- drop.tip(rate_tree, ant_df$V1)
  
  # Drop host tips
  animal_tree <- drop.tip(rate_tree, human$accession_id)
  human_tree <- drop.tip(rate_tree, animal$accession_id)
  tibble(rate = animal_tree$edge.length, host = "Neovison vison") %>%
    bind_rows(tibble(rate = human_tree$edge.length, host = "Human")) %>%
    add_column(country = country_name) %>%
    select(host, country, rate)
}

plot_df2 <- bind_rows(morsels2)

plot_df2 %>% 
  group_by(host, country) %>%
  summarise(n = n()) %>%
  right_join(plot_df2, by = c("host", "country")) %>%
  ggplot(aes(x = host, y = log(rate, base = 10), fill = host)) +
  facet_grid(rows = vars(country), scales = "free", space = "free") +
  theme_bw() +
  geom_point(color = "darkgrey", position = position_jitter(width = 0.08), 
             size = 0.5, 
             alpha = 0.8) +
  geom_boxplot(position = position_nudge(x = -0.2, y = 0), 
               width = 0.2, 
               outlier.shape = NA,
               alpha = 0.3) +
  geom_violinhalf(position = position_nudge(x = 0.2, y = 0), alpha = 1) +
  geom_text(aes(x = host, y = -1.5, label = paste0("n = ", n))) +
  coord_flip() +
  labs(y = "lg(subst. rate)", x = "Host") +
  theme(legend.position = "none")

ggsave(str_glue("results/human_animal_subsets/{human_background}/dating_out/mutations_rates_internal.png"), 
       width = 7, 
       height = 5)


#### Specific mutations #####
morsels3 <- foreach (prefix = prefixes) %do% {
  country_name <- str_split(prefix, "\\.")[[1]][2]
  time_tree <- read.nexus(str_glue("results/human_animal_subsets/{human_background}/dating_out/{prefix}/timetree.nexus"))
  div_tree <- read.nexus(str_glue("results/human_animal_subsets/{human_background}/dating_out/{prefix}/divergence_tree.nexus"))
  
  cluster_meta <- fread("results/cluster_annotation/deer_mink_parsed_clusters.csv") %>%
    select(accession_id, cluster)
  
  meta <- fread(str_glue("data/metadata/human_animal_subsets/{human_background}/{prefix}.csv")) %>%
    filter(accession_id %in% time_tree$tip.label) %>%
    select(-cluster) %>%
    left_join(cluster_meta)
  
  # Divide time tree by divergence tree
  rate_tree <- div_tree
  rate_tree$edge.length <- rate_tree$edge.length / time_tree$edge.length
  
  # Drop anthroponotic tips
  rate_tree <- drop.tip(rate_tree, ant_df$V1)
  
  n <- length(rate_tree$tip.label)
  rate_df <- tibble(accession_id = rate_tree$tip.label,
                    terminal_rate = rate_tree$edge.length[sapply(1:n, function(x,y) which(y == x),
                                                                 y = rate_tree$edge[,2])])
  rate_df %>% 
    left_join(meta) %>%
    add_column(country = country_name) %>%
    select(host, terminal_rate, cluster, country)
}

plot_df3 <- bind_rows(morsels3) %>%
  filter(host == "Neovison vison")

plot_df3 %>% 
  group_by(cluster) %>%
  summarise(n = n()) %>%
  right_join(plot_df3, by = "cluster") %>% 
  ggplot(aes(x = cluster, y = log(terminal_rate, base = 10), fill = cluster)) +
  theme_bw() +
  geom_point(color = "darkgrey", position = position_jitter(width = 0.08), 
             size = 0.5, 
             alpha = 0.8) +
  geom_boxplot(position = position_nudge(x = -0.2, y = 0), 
               width = 0.2, 
               outlier.shape = NA,
               alpha = 0.3) +
  geom_violinhalf(position = position_nudge(x = 0.2, y = 0), alpha = 1) +
  geom_text(aes(x = cluster, y = -1.5, label = paste0("n = ", n))) +
  coord_flip() +
  labs(y = "lg(subst. rate)", x = "Cluster") +
  theme(legend.position = "none")

ggsave(str_glue("results/human_animal_subsets/{human_background}/dating_out/mutations_rates_terminal.png"), 
       width = 7, 
       height = 5)

##### Compute by mutation #####
mut_df <- fread("results/mink_homoplasy_alele_frequency_V5.csv") %>%
  filter(mutation_annot %in% c("G37E", "Y486L", "N501T", "T229I", "L219V", "Y453F"))

# Load trees
prefix <- "all_mink.n1487.unambiguous.dedup"
time_tree <- read.nexus(str_glue("results/human_animal_subsets/{human_background}/dating_out/{prefix}/timetree.nexus"))
div_tree <- read.nexus(str_glue("results/human_animal_subsets/{human_background}/dating_out/{prefix}/divergence_tree.nexus"))

meta <- fread(str_glue("data/metadata/human_animal_subsets/{human_background}/{prefix}.csv")) %>%
  filter(accession_id %in% time_tree$tip.label)

aln_prefix <- gsub(".unambiguous", ".audacity_only.v8_masked.unambiguous", prefix)

aln <- read.dna(str_glue("data/alignments/human_animal_subsets/{human_background}/{aln_prefix}.fasta"),
                format = "fasta",
                as.matrix = T)

morsels <- foreach (i = seq(nrow(mut_df))) %do% {
  row <- mut_df[i, ]
  # Get allele by position
  aln_df <- tibble(accession_id = rownames(aln), 
                   allele = toupper(as.vector(as.character(aln[, row$nucleotide_pos]))))
  
  # Divide time tree by divergence tree
  rate_tree <- div_tree
  rate_tree$edge.length <- rate_tree$edge.length / time_tree$edge.length
  
  # Drop anthroponotic tips
  rate_tree <- drop.tip(rate_tree, ant_df$V1)
  
  n <- length(rate_tree$tip.label)
  
  with_df <- tibble(accession_id = rate_tree$tip.label,
                    terminal_rate = rate_tree$edge.length[sapply(1:n, function(x,y) which(y == x),
                                                                 y = rate_tree$edge[,2])]) %>%
    left_join(aln_df) %>%
    filter(allele == row$var_nuc) %>%
    add_column(mutation = row$mutation_annot,
               mutation_present = "Present")
  
  without_df <- tibble(accession_id = rate_tree$tip.label,
                       terminal_rate = rate_tree$edge.length[sapply(1:n, function(x,y) which(y == x),
                                                                    y = rate_tree$edge[,2])]) %>%
    left_join(aln_df) %>%
    filter(allele != row$var_nuc) %>%
    add_column(mutation = row$mutation_annot,
               mutation_present = "Absent")
  
  
  with_df %>%
    bind_rows(without_df)
}

bind_rows(morsels) %>% 
  mutate(log_rate = log(terminal_rate, base = 10)) %>%
  ggplot(aes(x = as.factor(mutation), y = log_rate, fill = mutation_present)) +
  geom_violinhalf(position = position_dodge(width = 1), alpha = 1) +
  geom_point(position = position_jitterdodge(dodge.width = 1, jitter.width = 0.07), 
             size = 0.5, 
             alpha = 0.3,
             color = "black") +
  geom_boxplot(position = position_dodge2(width = 1), 
               width = 0.2, 
               outlier.shape = NA, 
               alpha = 1)

ggpubr::ggarrange(plotlist = mut_plots)

pos_1 <- position_jitterdodge(
  jitter.width  = 0.25,
  jitter.height = 0,
  dodge.width   = 0.9
)

bind_rows(morsels) %>% 
  ggplot(aes(x = mutation, y = terminal_rate, color = mutation_present)) +
  # geom_jitter(alpha = 0.4, position = pos_1) +
  stat_summary(position = position_dodge(width = 0.5), 
               fun.y = "mean", 
               geom = "point", 
               size = 3) +
  stat_summary(position = position_dodge(width = 0.5), 
               fun.data = "mean_cl_normal",
               geom = "errorbar",
               width = 0.05,
               lwd = 1,
               fun.args = list(conf.int = 0.95)) +
  labs(y = "Subst. rate", x = "Mutation", color = "")

ggsave(str_glue("results/human_animal_subsets/{human_background}/dating_out/mutations_rates_terminal_by_mutation.png"), 
       width = 7, 
       height = 5)
