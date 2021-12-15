setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(ape)
require(lubridate)

prefix <- "mink.Netherlands.n3750"
prefixes <- c("mink.Netherlands.n3750", "mink.Denmark.n10512", "mink.USA.n35777")

# Remove ambiguous dates
for (prefix in prefixes) {
  meta <- fread(str_glue("data/metadata/human_animal_subsets/{prefix}.csv"))
  aln <- read.FASTA(str_glue("data/alignments/human_animal_subsets/{prefix}.audacity_only.v8_masked.aln.fasta"))
  
  # Remove ambiguous dates
  meta_filt <- meta %>%
    mutate(collection_date = date_decimal(as.Date(collection_date))) %>%
    filter(!is.na(collection_date))
  
  aln_filt <- aln[names(aln) %in% meta_filt$accession_id]
  
  write.FASTA(aln_filt, 
              str_glue("data/alignments/human_animal_subsets/{prefix}.audacity_only.v8_masked.aln.unambiguous.fasta"))
}

# Filter metadata
for (prefix in prefixes) {
  meta <- fread(str_glue("data/metadata/human_animal_subsets/{prefix}.csv"))
  aln <- read.FASTA(str_glue("data/alignments/human_animal_subsets/{prefix}.audacity_only.v8_masked.aln.unambiguous.dedup.fasta"))
  
  meta_filt <- meta %>%
    filter(accession_id %in% names(aln))
  
  dates <- meta_filt %>%
    rename(accession = accession_id) %>%
    select(accession, collection_date)
  
  fwrite(meta_filt, str_glue("data/metadata/human_animal_subsets/{prefix}.unambiguous.dedup.csv"))
  fwrite(dates, str_glue("data/metadata/human_animal_subsets/{prefix}.unambiguous.dedup.dates_only.csv"))
}

