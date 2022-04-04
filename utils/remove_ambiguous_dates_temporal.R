rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(ape)
require(lubridate)

meta_dir <- "data/metadata/temporal"
aln_dir <- "data/alignments/temporal"
prefixes <- list.files(paste0(meta_dir))
prefixes <- prefixes[!grepl("accessions|dates|unamb", prefixes)]
prefixes <- gsub(".csv", "", prefixes)
prefixes
# prefixes <- c("mink.Netherlands.n3750", "mink.Denmark.n10512", "mink.USA.n35777")

# Remove ambiguous dates
for (prefix in prefixes) {
  meta <- fread(str_glue("{meta_dir}/{prefix}.csv"))
  aln <- read.FASTA(str_glue("{aln_dir}/{prefix}.audacity_only.v8_masked.fasta"))
  
  # Remove ambiguous dates
  meta_filt <- meta %>%
    mutate(collection_date = date_decimal(as.Date(collection_date))) %>%
    filter(!is.na(collection_date))
  
  aln_filt <- aln[names(aln) %in% meta_filt$accession_id]
  
  write.FASTA(aln_filt, 
              str_glue("{aln_dir}/{prefix}.audacity_only.v8_masked.unambiguous.fasta"))
}

# Filter metadata
for (prefix in prefixes) {
  meta <- fread(str_glue("{meta_dir}/{prefix}.csv"))
  aln <- read.FASTA(str_glue("{aln_dir}/{prefix}.audacity_only.v8_masked.unambiguous.dedup.fasta"))
  
  meta_filt <- meta %>%
    filter(accession_id %in% names(aln))
  
  dates <- meta_filt %>%
    rename(accession = accession_id) %>%
    select(accession, collection_date)
  
  fwrite(meta_filt, str_glue("{meta_dir}/{prefix}.unambiguous.dedup.csv"))
  fwrite(dates, str_glue("{meta_dir}/{prefix}.unambiguous.dedup.dates_only.csv"))
}

