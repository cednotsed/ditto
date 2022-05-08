setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(lubridate)
registerDoParallel(cores = 10)

# Load metadata
audacity <- fread("data/GISAID-hCoV-19-phylogeny-2022-02-21/metadata.csv")
meta <- fread("data/metadata/all_sequence_metadata_260322.audacity.parsed.tsv", nThread = 10)

human <- meta %>% 
  filter(host == "Human")

mink <- meta %>%
  filter(host == "Neovison vison")

table((meta %>% filter(host != "Human"))$host)

mutation <- "M_I82T"
mutation <- "NSP3_L1035F"

# Get proportions
human %>%
  summarise(prop = sum(grepl(mutation, aa_substitutions)) / nrow(human))

# Earliest date for mutations
mink_df <- fread("results/allele_frequency/mink_homoplasy_alele_frequency_V5.csv") %>%
  mutate(host = "mink")
deer_df <- fread("results/allele_frequency/deer_homoplasy_alele_frequency_V5.csv") %>%
  mutate(host = "deer")
mut_df <- bind_rows(mink_df, deer_df)
mut_df
mutation_list <- c("NSP9_G37E", "NS3_L219V", 
                   "Spike_F486L", "Spike_N501T", "Spike_Y453F", 
                   "NSP3_L1035F",
                   "N_D377Y", "M_I82T"
                   # "NS3_H182Y", "NSP12_N507I"
                   )

morsels <- foreach (mutation = mutation_list) %do% {
  human %>% 
    mutate(collection_date = as.Date(collection_date)) %>%
    filter(grepl(mutation, aa_substitutions), 
           !is.na(collection_date)) %>%
    summarise(earliest = min(collection_date, na.rm = T)) %>%
    add_column(mut = mutation)
}

time_df <- bind_rows(morsels)


morsels2 <- foreach (i = seq(nrow(time_df))) %do% {
  human %>% 
    filter(collection_date == time_df$earliest[i]) %>%
    filter(grepl(time_df$mut[i], aa_substitutions)) %>%
    group_by(loc2) %>%
    summarise(n = n()) %>%
    add_column(mut = time_df$mut[i])
}

loc_df <- bind_rows(morsels2)
time_loc_df <- time_df %>%
  left_join(loc_df)
fwrite(time_loc_df, "results/temporal_distribution/mutation_times.csv")

# Get temporal distribution
# time_plots <- foreach(mutation = mutation_list) %do% {
#   date_df <- human %>%
#     filter(grepl(mutation, aa_substitutions)) %>%
#     mutate(collection_date = as.Date(collection_date, "%Y-%m-%d")) %>%
#     filter(!is.na(collection_date))
#   
#   plt <- date_df %>% 
#     ggplot(aes(y = loc2, x = collection_date, color = loc2)) +
#     geom_point() +
#     labs(x = "Sampling date", y = "Country", title = mutation) +
#     scale_x_date(breaks = seq(min(date_df$collection_date), 
#                               max(date_df$collection_date), 
#                               "months"),
#                  date_labels = "%b-%y")
#   
#   if (length(unique(date_df$loc2)) > 40) {
#     plt <- plt + theme(axis.text.y = element_blank(),
#                        axis.ticks.y = element_blank(),
#                        axis.text.x = element_text(angle = 45, hjust = 1),
#                        legend.position = "none")
#   } else {
#     plt <- plt +
#       theme(axis.text.x = element_text(angle = 45, hjust = 1),
#             legend.position = "none")
#   }
#   
#   plt
# }
# 
# combined <- ggpubr::ggarrange(plotlist = time_plots)
# ggsave("results/temporal_distribution_mutations.png", 
#        combined, 
#        dpi = 300,
#        width = 12, 
#        height = 10)
# mink %>% 
#   mutate(collection_date = as.Date(collection_date)) %>%
#   summarise(earliest = min(collection_date, na.rm = T))
# 
# human_filt <- human %>%
#     mutate(collection_date = as.Date(collection_date)) %>%
#     filter(collection_date < ymd("2020-04-01"))
# 
# morsels_filt <- foreach (mutation = mutation_list) %do% {
#   human_filt %>% 
#     mutate(collection_date = as.Date(collection_date)) %>%
#     filter(grepl(mutation, aa_substitutions), 
#            !is.na(collection_date)) %>%
#     summarise(n = n()) %>%
#     add_column(mut = mutation)
# }
# 
# bind_rows(morsels_filt)
# 
# mutation <- "Spike_F486L"
# 
# human %>% 
#   mutate(collection_date = as.Date(collection_date)) %>%
#   filter(grepl(mutation, aa_substitutions), 
#          !is.na(collection_date))
# 
# 
# human_background <- "V5"
# host_name <- "Deer"
# species_name <- "Odocoileus virginianus"
# file <- "all_deer.n145.audacity_only.v8_masked.fasta"
# meta_file <- "all_deer.n145.csv"
# us_meta <- fread(str_glue("data/metadata/human_animal_subsets/{human_background}/{meta_file}")) %>%
#   rename_all(~ tolower(gsub(" ", "_", .))) %>%
#   select(accession_id, host, clade, pango_lineage, location, collection_date) %>%
#   separate(location, into = c("loc1", "loc2", "loc3"), sep = " / ") %>%
#   mutate(location = paste0(loc1, " / ", loc2, " / ", loc3)) %>%
#   as_tibble()
# 
# us_human <- us_meta %>% filter(host == "Human")
# meta %>%
#   filter(accession_id %in% us_human$accession_id) %>% 
#   select(collection_date, host, aa_substitutions) %>%
#   View()
