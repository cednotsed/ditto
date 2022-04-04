rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(ape)
require(ggplot2)
require(ggrepel)
require(foreach)
require(doParallel)
registerDoParallel(cores = 10)

host_list <- list(mink = "Neovison vison", deer = "Odocoileus virginianus")
bg_list <- c("V1", "V5")

for (i in seq(length(host_list))) {
  for (human_background in bg_list) {
    host_name <- names(host_list)[i]
    species_name <- host_list[[i]]
    meta_path <- str_glue("data/metadata/human_animal_subsets/{human_background}")
    meta_path <- list.files(meta_path, host_name, full.names = T)
    meta_path <- meta_path[!grepl("acces|date", meta_path)]
    meta_path <- meta_path[!grepl("unambi|iqtree|delim", meta_path) & grepl("all", meta_path)]
    aln_path <- str_glue("data/alignments/human_animal_subsets/{human_background}")
    aln_path <- list.files(aln_path, host_name, full.names = T)
    aln_path <- aln_path[!grepl("unambi|iqtree|delim", aln_path) & grepl("all", aln_path)]
    
    # Load data
    meta <- fread(meta_path)
    aln <- read.dna(aln_path,
                    format = "fasta",
                    as.matrix = T)
    
    ant_df <- read.csv("data/metadata/netherlands_humans_anthroponoses.txt", header = F)
    
    # Remove anthroponoses
    aln <- aln[!(rownames(aln) %in% ant_df$V1), ]
    
    # Match aln and meta
    meta.match <- meta[match(rownames(aln), meta$accession_id), ]
    
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
    ref <- read.dna("data/genomes/WIV04_EPI_ISL_402124.fasta", 
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
           color = "Ref. pos.") +
      xlim(0, 1) +
      ylim(0, 1)
    plt
    ggsave(str_glue("results/allele_frequency/{host_name}.{human_background}.global_allele_frequency.png"),
           plt,
           dpi = 300, 
           width = 6,
           height = 4)
    
    fwrite(parsed_df, str_glue("results/allele_frequency/{host_name}.{human_background}.global_allele_frequency.csv"))
    
  }
}

stopImplicitCluster()