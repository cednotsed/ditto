setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(foreach)
require(seqinr)

# Identify missing accessions
meta <- fread("data/metadata/combined_metadata_201121.tsv")
audacity <- fread("data/GISAID-hCoV-19-phylogeny-2021-11-16/metadata.csv")

to_remove <- meta[!(meta$accession_id %in% audacity$accession_id), ]

table(to_remove$host)

# Make regex string
to_search <- to_remove$accession_id
to_search <- paste0(to_search, collapse = "|")

# Remove sequences from fasta
for (file in list.files("data/genomes/", pattern = "201121.fasta", full.names = T)) {
  fasta <- read.fasta(file)
  fasta_filt <- fasta[!grepl(to_search, names(fasta))]
  n_seqs <- length(fasta)
  new_n_seqs <- length(fasta_filt)
  
  # Parse file name
  new_file <- gsub("201121.fasta", "201121.audacity_only.fasta", file)
  new_file <- gsub(paste0("_n", n_seqs), paste0("_n", new_n_seqs), new_file)
  
  # Parse fasta names
  accessions <- separate(tibble(accession = names(fasta_filt)),
           accession, 
           into = c(NA, "accession", NA), sep = "\\|")
  
  # Save fasta
  write.fasta(sequences = fasta_filt, 
              names = accessions$accession, 
              file.out = new_file)
}

