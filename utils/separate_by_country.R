setwd("../Desktop/git_repos/ditto/")
require(ape)
require(tidyverse)
require(data.table)

meta <- fread("data/metadata/human_animal_subsets/V5/all_mink.n1487.csv")
aln <- read.FASTA("data/alignments/human_animal_subsets/V5/all_mink.n1487.audacity_only.v8_masked.fasta")

for (loc in unique(meta$location)) {
  meta_filt <- meta %>% 
    filter(location == loc)
  
  aln_filt <- aln[names(aln) %in% meta_filt$accession_id]
  country_name <- strsplit(loc, " / ")[[1]][2]
  n_isolates <- nrow(meta_filt)
  fwrite(meta_filt, str_glue("data/metadata/human_animal_subsets/V5/mink.{country_name}.n{n_isolates}.csv"))
  write.FASTA(aln_filt, str_glue("data/alignments/human_animal_subsets/V5/mink.{country_name}.n{n_isolates}.audacity_only.v8_masked.fasta"))
}

