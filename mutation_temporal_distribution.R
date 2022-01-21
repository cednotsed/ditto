setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(lubridate)
require(randomcoloR)
registerDoParallel(cores = 10)

# Load metadata
audacity <- fread("data/GISAID-hCoV-19-phylogeny-2021-11-16/metadata.csv")
meta <- fread("data/metadata/all_sequence_metadata_231121.tsv", nThread = 10) %>%
  rename_all(~ tolower(gsub(" ", "_", .))) %>%
  filter(accession_id %in% audacity$accession_id) %>%
  separate(location, into = c("loc1", "loc2", "loc3"), sep = " / ") %>%
  mutate(location = paste0(loc1, " / ", loc2)) %>%
  as_tibble()

human <- meta %>% 
  filter(host == "Human")

mink <- meta %>%
  filter(host == "Neovison vison")

# Earliest date for mutations
mut_df <- fread("results/mink_deer_homoplasy_allele_frequency_V5.csv")
mut_df
mutation_list <- c("NS3_T229I", "Spike_N501T", "M_I82T",
                   "N_D377Y", "NSP3_L1035F",
                   "NSP9_G37E", "NS3_L219V", 
                   "Spike_F486L", "Spike_Y453F"
                   # "NS3_H182Y", "NSP12_N507I"
)

# Get discrete colors
country_mut <- human %>%
  filter(grepl("Spike_F486L|NSP9_G37E|NS3_L219V|Spike_Y453F", aa_substitutions)) %>%
  mutate(collection_date = as.Date(collection_date, "%Y-%m-%d")) %>%
  filter(!is.na(collection_date))
countries <- unique(country_mut$loc2)
colordict <- as.list(distinctColorPalette(length(countries)))
names(colordict) <- countries

# Get temporal distribution
time_plots <- foreach(mutation = mutation_list) %do% {
  date_df <- human %>%
    filter(grepl(mutation, aa_substitutions)) %>%
    mutate(collection_date = as.Date(collection_date, "%Y-%m-%d")) %>%
    filter(!is.na(collection_date))
  
  # # Get first 30 countries
  # to_keep <- unique(date_df$loc2)[1:30]
  # date_df <- date_df %>%
  #   filter(loc2 %in% to_keep)
  
  plt <- date_df %>% 
    ggplot(aes(y = loc2, x = collection_date, color = loc2)) +
    geom_point() +
    labs(x = "Sampling date", y = "Country", title = mutation) +
    # scale_x_date(breaks = seq(min(date_df$collection_date), 
    #                           max(date_df$collection_date), 
    #                           "months"),
    #              date_labels = "%b-%y")
    scale_x_date(breaks = seq(ymd("2020-02-01"), 
                              ymd("2021-12-01"), 
                              "2 months"),
                 limits = c(ymd("2020-02-01"), 
                            ymd("2021-12-01")),
                 date_labels = "%b-%y") +
    theme_bw()
  
  if (length(unique(date_df$loc2)) > 40) {
    plt <- plt + theme(axis.text.y = element_blank(),
                       axis.ticks.y = element_blank(),
                       axis.text.x = element_text(angle = 45, hjust = 1),
                       panel.grid.minor.y = element_blank(),
                       panel.grid.major.y = element_blank(),
                       legend.position = "none") +
      scale_color_manual(values = distinctColorPalette(length(unique(date_df$loc2))))
    n_countries <- length(unique(date_df$loc2))
    print(str_glue("{mutation}: No. of countries = {n_countries}"))
  } else {
    plt <- plt +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") +
      scale_color_manual(values = colordict[unique(date_df$loc2)])
  }
  
  plt
}

combined <- ggpubr::ggarrange(plotlist = time_plots)
ggsave("results/temporal_distribution_mutations_distinct.png", 
       combined, 
       dpi = 100,
       width = 12, 
       height = 10)

combined1 <- ggpubr::ggarrange(plotlist = time_plots[1:5], nrow = 2, ncol = 3)
combined2 <- ggpubr::ggarrange(plotlist = time_plots[6:9], nrow = 1, ncol = 4)
ggsave("results/temporal_distribution_mutations_distinct_1.png", 
       combined1, 
       dpi = 300,
       width = 10, 
       height = 6)
ggsave("results/temporal_distribution_mutations_distinct_2.png", 
       combined2, 
       dpi = 600,
       width = 13, 
       height = 4)

