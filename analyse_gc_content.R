setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(ape)
require(ggplot2)
require(foreach)
require(doParallel)
registerDoParallel(cores=8)

prefix <- "deer.USA.n6670"
prefix <- "mink.Netherlands.n2213"
prefix <- "mink.Denmark.n3219"
prefix <- "mink.USA.n6726"
host_name <- "Neovison vison"
aln <- read.dna(str_glue("data/alignments/human_animal_subsets/{prefix}.audacity_only.v8_masked.aln.fasta"),
                format = "fasta",
                as.matrix = T)

# Parse tree names
aln_names <- separate(tibble(accession = rownames(aln)),
                      accession, 
                      into = c(NA, "accession", NA), sep = "\\|")$accession

aln_names[1] <- rownames(aln)[1]
rownames(aln) <- aln_names

# Load metadata
meta <- fread(str_glue("data/metadata/human_animal_subsets/{prefix}.tsv")) %>%
  rename_all(~ tolower(gsub(" ", "_", .))) %>%
  select(accession_id, host, clade, lineage, location, collection_date) %>%
  separate(location, into = c("loc1", "loc2", "loc3"), sep = " / ") %>%
  mutate(location = paste0(loc1, " / ", loc2, " / ", loc3)) %>%
  rename(Accession = accession_id) %>%
  as_tibble()

full_meta <- tibble(Accession = aln_names) %>%
  left_join(meta) %>%
  mutate(collection_date = ifelse(Accession == "NC_045512.2", "2019-12-26", collection_date))

# animal <- aln[meta$host == "Neovison vison", ]
# human <- aln[meta$host != "Neovison vison", ]
# 
# animal_morsels <- foreach (i = seq(nrow(animal))) %do% {
#   seq <- animal[i, ]
#   base_freq <- base.freq(seq)
#   tibble(Accession = rownames(seq), 
#          base = names(base_freq), 
#          base_freq = as.numeric(base_freq))
# }
# 
# human_morsels <- foreach (i = seq(nrow(human))) %do% {
#   seq <- human[i, ]
#   base_freq <- base.freq(seq)
#   tibble(Accession = rownames(seq), 
#          base = names(base_freq), 
#          base_freq = as.numeric(base_freq))
# }
# 
# bind_rows(human_morsels) %>%
#   left_join(meta) %>%
#   ggplot(aes(x = as.Date(collection_date), y = base_freq * 100)) +
#   facet_grid(rows = vars(base)) +
#   geom_point()

## GC content ##
morsels <- foreach (i = seq(nrow(aln))) %dopar% {
  seq <- aln[i, ]
  gc <- GC.content(seq)
  tibble(Accession = rownames(seq), 
         gc = gc)
}

bind_rows(morsels) %>%
  left_join(meta) %>%
  filter(host == host_name) %>%
  ggplot(aes(x = as.Date(collection_date), y = gc * 100, color = host)) +
  geom_point()
