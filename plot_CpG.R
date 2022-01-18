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
                        as.prob=F, as.matrix=FALSE,
                        fast.moving.side="right", with.labels=TRUE)
  cg <- dn[, "CG"]
  
  return(tibble(cg_freq = cg, host = prefix))
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

# human_background <- "V1"
# host_name <- "Mink"
# aln_path <- "all_mink.audacity_only.v8_masked.aln.fasta"
# meta_file <- "all_mink.csv"
# host_species <- "Neovison vison"



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


animal_cg <- parse(animal, host_species)
human_cg <- parse(human, "Human")

plot_df <- bind_rows(animal_cg, human_cg)

# Get statistics
stat_df <- plot_df %>%
  group_by(host) %>%
  summarise(mean = mean(cg_freq), median = median(cg_freq), n = n())

plt <- stat_df %>%
  right_join(plot_df) %>%
  mutate(host = factor(host, levels = c(host_species, "Human"))) %>%
  ggplot(aes(x = host, y = cg_freq, fill = host)) +
  geom_violinhalf(position = position_nudge(x = 0.2, y = 0), alpha = 1) +
  geom_point(position = position_jitter(width = 0.08), 
             size = 0.5, 
             alpha = 0.3,
             color = "black") +
  geom_boxplot(position = position_nudge(x = -0.2, y = 0), 
               width = 0.2, 
               outlier.shape = NA, 
               alpha = 1) +
  geom_text(aes(x = host, y = 405, label = paste0("n = ", n))) +
      labs(y = "CpG Frequency", x = "Host") +
      coord_flip() +
      theme(legend.position = "none",
      plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      axis.text.y = element_blank()) +
  ylim(400, 440)
  # stat_compare_means(method = "wilcox.test", 
  #                    label.x = 0.5, label.y = 410)

plt
ggsave(str_glue("results/mutation_bias/{host_name}_{human_background}_CpG_frequencies.png"), 
       plot = plt, 
       dpi = 300, 
       width = 6, 
       height = 4)

plot_df %>%
  group_by(host) %>%
  summarise(med = median(cg_freq))
