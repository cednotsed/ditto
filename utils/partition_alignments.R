setwd("../Desktop/git_repos/ditto/")
require(ape)
require(data.table)
require(tidyverse)
require(foreach)

meta <- fread("data/metadata/wuhan-hu-1_genome_annotations_V2.csv") %>%
  as_tibble() %>%
  mutate(protein_name = case_when(protein_name == "NSP12a (part 1)" ~ "NSP12a_1",
                                  protein_name == "NSP12a (part 2)" ~ "NSP12a_2",
                                  TRUE ~ protein_name))

prefixes <- list.files("data/alignments/human_animal_subsets/V2/")
prefixes <- prefixes[grepl("unambiguous.dedup.fasta", prefixes)]
prefixes <- gsub(".audacity_only.v8_masked.unambiguous.dedup.fasta", "", prefixes)
prefixes

for (prefix in prefixes) {
# prefix <- "mink_n789_201121"
  dir.create(str_glue("data/alignments/human_animal_subsets/codon_V2/{prefix}"))
  aln <- read.dna(str_glue("data/alignments/human_animal_subsets/V2/{prefix}.audacity_only.v8_masked.unambiguous.dedup.fasta"),
                  format="fasta",
                  as.matrix = T)
  
  assertthat::assert_that(ncol(aln) == 29891)
  
  orf1ab <- c("ORF1ab")
  spike <- c("S")
  structural <- c("M", "E", "N")
  accessory <- c("ORF3a", "ORF6", "ORF7a", "ORF7b", "ORF8", "ORF10")
  
  partition_list <- list(orf1ab = orf1ab, spike = spike, 
                     structural = structural, 
                     accessory = accessory)
  
  for (i in seq(length(partition_list))) {
    partition <- partition_list[[i]]
    partition_name <- names(partition_list)[i]
    
    # Retrieve alignments for each protein
    protein_alns <- foreach (protein = partition) %do% {
      meta_row <- meta[meta$region == protein, ]
      genome_start <- min(meta_row$genome_start)
      genome_end <- max(meta_row$genome_end)
      protein_aln <- aln[, genome_start:genome_end]
      protein_aln
    }
    
    partition_aln <- do.call(cbind, protein_alns)
    
    write.dna(partition_aln, format = "fasta", 
              file = str_glue("data/alignments/human_animal_subsets/codon_V2/{prefix}/{prefix}.{partition_name}.audacity_only.v8_masked.unambiguous.dedup.fasta"))
  }
}
