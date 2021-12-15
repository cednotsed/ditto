setwd("../Desktop/git_repos/ditto/")
require(ape)
require(data.table)
require(tidyverse)

meta <- fread("data/metadata/wuhan-hu-1_genome_annotations_V2.csv") %>%
  as_tibble() %>%
  mutate(protein_name = case_when(protein_name == "NSP12a (part 1)" ~ "NSP12a_1",
                                  protein_name == "NSP12a (part 2)" ~ "NSP12a_2",
                                  TRUE ~ protein_name))

prefixes <- c("mink_n789_201121", "cat_n65_201121", "deer_n73_201121")

for (prefix in prefixes) {
# prefix <- "mink_n789_201121"
  dir.create(str_glue("data/alignments/codon_alignments/{prefix}"))
  aln <- read.dna("data/alignments/mink_n789_201121.audacity_only.v8_masked.aln.fasta",
                  format="fasta",
                  as.matrix = T)
  
  assertthat::assert_that(ncol(aln) == 29903)
  
  proteins <- meta$protein_name
  to_remove <- proteins[grepl("corona|\\'|\\-", proteins)]
  proteins <- proteins[!(proteins %in% to_remove)]
  
  for (protein in proteins) {
    meta_row <- meta[meta$protein_name == protein, ]
    genome_start <- meta_row$genome_start
    genome_end <- meta_row$genome_end
    codon_aln <- aln[ , genome_start:genome_end]
    write.dna(codon_aln, format = "fasta", append = F, 
              file = str_glue("data/alignments/codon_alignments/{prefix}/{prefix}.{protein}.audacity_only.v8_masked.aln.fasta"))
  }  
}
