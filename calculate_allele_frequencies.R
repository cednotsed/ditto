rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(ape)
require(ggplot2)
require(ggrepel)
require(foreach)
require(doParallel)
# registerDoParallel(cores = 10)

human_background <- "V1"
host_name <- "Mink"
species_name <- "Neovison vison"
file <- "all_mink.audacity_only.v8_masked.fasta"
meta_file <- "all_mink.csv"

human_background <- "V1"
host_name <- "Deer"
species_name <- "Odocoileus virginianus"
file <- "deer.USA.n6670.audacity_only.v8_masked.fasta"
meta_file <- "deer.USA.n6670.csv"

# human_background <- "V2"
# host_name <- "Mink"
# species_name <- "Neovison vison"
# file <- "all_mink.audacity_only.v8_masked.fasta"
# meta_file <- "all_mink.csv"

# human_background <- "V2"
# host_name <- "Deer"
# species_name <- "Odocoileus virginianus"
# file <- "deer.USA.n10073.audacity_only.v8_masked.fasta"
# meta_file <- "deer.USA.n10073.csv"

# human_background <- "V4"
# host_name <- "Mink"
# species_name <- "Neovison vison"
# file <- "all_mink.audacity_only.v8_masked.fasta"
# meta_file <- "all_mink.csv"

# human_background <- "V4"
# host_name <- "Deer"
# species_name <- "Odocoileus virginianus"
# file <- "deer.USA.n10073.audacity_only.v8_masked.fasta"
# meta_file <- "deer.USA.n10073.csv"

# human_background <- "V5"
# host_name <- "Mink"
# species_name <- "Neovison vison"
# file <- "all_mink.n1487.audacity_only.v8_masked.fasta"
# meta_file <- "all_mink.n1487.csv"
# 
# human_background <- "V5"
# host_name <- "Deer"
# species_name <- "Odocoileus virginianus"
# file <- "all_deer.n145.audacity_only.v8_masked.fasta"
# meta_file <- "all_deer.n145.csv"

aln <- read.dna(str_glue("data/alignments/human_animal_subsets/{human_background}/{file}"),
                format = "fasta",
                as.matrix = T)

ant_df <- read.csv("data/metadata/netherlands_humans_anthroponoses.txt", header = F)

# Remove anthroponoses
aln <- aln[!(rownames(aln) %in% ant_df$V1), ]
  
# Load metadata
meta <- fread(str_glue("data/metadata/human_animal_subsets/{human_background}/{meta_file}")) %>%
  rename_all(~ tolower(gsub(" ", "_", .))) %>%
  select(accession_id, host, clade, pango_lineage, location, collection_date) %>%
  separate(location, into = c("loc1", "loc2", "loc3"), sep = " / ") %>%
  mutate(location = paste0(loc1, " / ", loc2, " / ", loc3)) %>%
  rename(Accession = accession_id) %>%
  as_tibble()

meta.match <- meta[match(rownames(aln), meta$Accession), ]

animal <- aln[meta.match$host == species_name, ]
# write.FASTA(animal, str_glue("results/human_animal_subsets/allele_frequency/{human_background}/{host_name}.fasta"))

human <- aln[meta.match$host != species_name, ]
# write.FASTA(human, str_glue("results/human_animal_subsets/allele_frequency/{human_background}/human_for_{host_name}.fasta"))

animal_morsels <- foreach(pos = seq(ncol(aln)), .packages = c("ape", "tidyverse")) %dopar% {
  base_freq <- base.freq(animal[, pos])
  
  tibble(pos = pos,
         base = names(base_freq),
         base_freq = as.numeric(base_freq))
}
  
human_morsels <- foreach(pos = seq(ncol(aln)), .packages = c("ape", "tidyverse")) %dopar% {
  base_freq <- base.freq(human[, pos])
  
  tibble(pos = pos,
         base = names(base_freq),
         base_freq = as.numeric(base_freq))
}

animal_df <- bind_rows(animal_morsels)
human_df <- bind_rows(human_morsels)

merged_df <- animal_df %>% 
  left_join(human_df, by = c("pos", "base"), suffix = c("_animal", "_human"))

# Annotate mutations
ref <- read.dna("data/genomes/Wuhan-Hu-1_NC_045512.2.fasta", 
                format = "fasta",
                as.matrix = T)
ref_positions <- as.vector(t(as.character(ref[, merged_df$pos])))

hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.csv") %>%
  filter(!grepl("coronavirus|3'UTR ", protein_name))


parsed_df <- merged_df %>%
  add_column(ref_nuc = ref_positions, .after = 1) %>%
  rename(var_nuc = base, nucleotide_pos = pos) %>%
  mutate(var_nuc = toupper(var_nuc), ref_nuc = toupper(ref_nuc)) %>%
  left_join(hookup, by = c("nucleotide_pos", "ref_nuc", "var_nuc")) %>% 
  mutate(mutation_name = paste0(ref_nuc, nucleotide_pos, var_nuc), 
         variant_name = paste0(ref_AA, codon_number, var_AA)) %>%
  mutate(variant_name = case_when(mutation_type == "S" & region_type != "non-coding" ~ "Syn",
                                  region_type == "non-coding" ~ "non-coding",
                                  TRUE ~ variant_name)) %>%
  mutate(mutation = str_glue("{mutation_name} ({variant_name})"),
         fold = base_freq_animal / base_freq_human) %>%
  mutate(annot = ifelse(base_freq_animal > 0.1 & fold >= 2,
                        mutation, NA)) %>%
  arrange(nucleotide_pos) %>%
  filter(protein_name != "3'UTR") # Remove 3'UTR entry since separate entry is present

plt <- parsed_df %>%
  ggplot(aes(x = base_freq_animal, y = base_freq_human, color = nucleotide_pos)) +
  geom_point() +
  # geom_text_repel(aes(label = annot),
  #           color = "black",
  #           size = 2,
  #           min.segment.length = 1,
  #           direction = "y",
  #           max.overlaps = 25) +
  geom_abline(slope = 2, intercept = 0, color = "red", lty = "dotted") +
  geom_abline(slope = 0.5, intercept = 0, color = "red", lty = "dotted") +
  scale_color_gradient(low = "blue", high = "salmon") +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  labs(x = str_glue("{host_name} allele freq."), 
       y = "Human allele freq.",
       color = "Ref. pos.")
plt
ggsave(str_glue("results/human_animal_subsets/allele_frequency/{human_background}/{host_name}.global_allele_frequency.png"),
       plt,
       dpi = 300, 
       width = 6,
       height = 4)

fwrite(parsed_df, str_glue("results/human_animal_subsets/allele_frequency/{human_background}/{host_name}.global_allele_frequency.csv"))
