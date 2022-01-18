rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(ape)
require(ggplot2)
require(foreach)
require(doParallel)
registerDoParallel(cores = 10)

# prefixes <- c("deer.USA.n6670")
# host_name <- "Odocoileus virginianus"
# save_name <- "deer"

prefixes <- c("mink.Netherlands.n2213", "mink.Denmark.n3219", "mink.USA.n6726")
host_name <- "Neovison vison"
save_name <- "mink"
min_freq_change <- 0.1
# prefix <- prefixes[1]
morsels <- foreach(prefix = prefixes) %do% {
  aln <- read.dna(str_glue("data/alignments/human_animal_subsets/V5/{prefix}.audacity_only.v8_masked.aln.fasta"),
                  format = "fasta",
                  as.matrix = T)
  
  # # Parse tree names
  # aln_names <- separate(tibble(accession = rownames(aln)),
  #                       accession, 
  #                       into = c(NA, "accession", NA), sep = "\\|")$accession
  # 
  # aln_names[1] <- rownames(aln)[1]
  # rownames(aln) <- aln_names
  
  # Load metadata
  meta <- fread(str_glue("data/metadata/human_animal_subsets/V1/{prefix}.tsv")) %>%
    rename_all(~ tolower(gsub(" ", "_", .))) %>%
    select(accession_id, host, clade, lineage, location, collection_date) %>%
    separate(location, into = c("loc1", "loc2", "loc3"), sep = " / ") %>%
    mutate(location = paste0(loc1, " / ", loc2, " / ", loc3)) %>%
    rename(Accession = accession_id) %>%
    as_tibble()
  
  full_meta <- tibble(Accession = rownames(aln)) %>%
    left_join(meta) %>%
    mutate(collection_date = ifelse(Accession == "NC_045512.2", "2019-12-26", collection_date)) %>%
    mutate(host = ifelse(Accession == "NC_045512.2", "Human", host))
  
  animal <- aln[full_meta$host == host_name, ]
  human <- aln[full_meta$host != host_name, ]
  
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
  
  # Get sites variant between humans and animals
  unique_df <- merged_df %>%
    filter(base_freq_animal != base_freq_human)
  
  # Calculate change in frequency
  unique_df <- unique_df %>%
    mutate(abs_base_diff = abs(base_freq_animal - base_freq_human),
           base_diff = base_freq_animal - base_freq_human) %>%
    filter(abs_base_diff > min_freq_change) %>%
    arrange(desc(abs_base_diff))
  
  # Check if all changes are biallelic
  num_alleles <- unique_df %>% group_by(pos) %>%
    summarise(n_alleles = n())
  # assertthat::assert_that(sum(num_alleles$n_alleles != 2) == 0, msg = prefix)
  pos_to_keep <- num_alleles$pos[num_alleles$n_alleles == 2]
  
  # Get full SNP and SNV names
  mut_morsels <- foreach(position = pos_to_keep) %dopar% {
    pos_df <- unique_df %>% filter(pos == position)
    ref_nuc <- (pos_df %>% filter(base_diff < 0))$base
    var_nuc <- (pos_df %>% filter(base_diff > 0))$base
    to_keep <- pos_df %>% 
      filter(base_diff > 0) %>%
      select(-base) %>%
      rename(nucleotide_pos = pos) %>%
      add_column(ref_nuc = ref_nuc, var_nuc = var_nuc)
    
    to_keep
  }
  
  mut_df <- bind_rows(mut_morsels)
  
  # Hookup to protein meta
  hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V2.csv")
  
  plot_df <- mut_df %>%
    mutate(var_nuc = toupper(var_nuc), ref_nuc = toupper(ref_nuc)) %>%
    left_join(hookup, by = c("nucleotide_pos", "ref_nuc", "var_nuc")) %>%
    mutate(mutation_name = paste0(ref_nuc, nucleotide_pos, var_nuc), 
           variant_name = paste0(ref_AA, codon_number, var_AA)) %>%
    mutate(variant_name = case_when(mutation_type == "S" & region_type != "non-coding" ~ "Syn",
                                     region_type == "non-coding" ~ "non-coding",
                                     TRUE ~ variant_name)) %>%
    mutate(mutation = str_glue("{mutation_name} ({variant_name})")) %>%
    arrange(nucleotide_pos)
  
  # Add prefix column for facet_grid
  plot_df %>% add_column(prefix = prefix)
}

all_df <- bind_rows(morsels) %>% 
  filter(protein_name != "3'UTR  (stem-loop II-like motif)") %>%
  arrange(nucleotide_pos) %>%
  separate(prefix, into = c(NA, "prefix1", "prefix2"), sep = "\\.") %>%
  mutate(prefix = str_glue("{prefix1}.{prefix2}"))

# Order bars by position
all_df <- all_df %>% 
  mutate(mutation = factor(mutation, levels = unique(all_df$mutation)),
         protein_name = factor(protein_name, levels = unique(all_df$protein_name)))

# Get vertical line info
v_df <- all_df %>% 
  group_by(protein_name) %>%
  summarise(x_int = n_distinct(mutation))
  
v_df <- v_df[match(levels(v_df$protein_name), v_df$protein_name), ]

for (i in seq(2, nrow(v_df))) {
  v_df$x_int[i] <- v_df$x_int[i] + v_df$x_int[i - 1]
}

# Plot mutations
all_df %>% 
  ggplot(aes(x = mutation, y = prefix, fill = base_diff)) +
  geom_tile(color = "black") +
  theme_classic() +
  # facet_grid(rows = vars(prefix)) +
  # geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = ),
        plot.margin = unit(c(1,1,1,1), "cm")) +
  labs(x = "Mutation", 
       y = "Dataset", 
       fill = "Allele frequency change") +
  geom_vline(xintercept = v_df$x_int + 0.5, lty = "dotted") +
  scale_fill_gradient(low = "green", high = "red")

ggsave(str_glue("results/human_animal_subsets/allele_frequency/{save_name}.allele_frequency_by_country.png"), dpi = 600, height = 6, width = 14)

