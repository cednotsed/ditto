setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(foreach)
require(doParallel)
require(lubridate)
require(randomcoloR)
registerDoParallel(cores = 10)

# Load metadata
meta <- fread("data/metadata/all_sequence_metadata_260322.audacity.parsed.tsv", nThread = 10)

human <- meta %>% 
  filter(host == "Human")

mink <- meta %>%
  filter(host == "Neovison vison")

deer <- meta %>%
  filter(host == "Odocoileus virginianus")

# Earliest date for mutations
mutation_list <- c("Spike_N501T", "M_I82T",
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
set.seed(66)

time_plots <- foreach(mutation = mutation_list) %do% {
  # mutation <- "NSP3_L1035F"
  date_df <- human %>%
    filter(grepl(mutation, aa_substitutions)) %>%
    mutate(collection_date = as.Date(collection_date, "%Y-%m-%d")) %>%
    filter(!is.na(collection_date))
  
  # Get countries of isolates before first outbreaks
  print(date_df %>% 
    filter(collection_date < min(deer$collection_date, na.rm = T)) %>%
    group_by(location) %>%
    summarise(n = n()))
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
    scale_x_date(breaks = seq(ymd("2020-03-01"), 
                              ymd("2022-03-01"), 
                              "4 months"),
                 limits = c(ymd("2020-03-01"), 
                            ymd("2022-03-01")),
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
  
  if (mutation != "NSP3_L1035F") {
    plt + geom_vline(xintercept = min(mink$collection_date, na.rm = T),
                     color = "red", lty = "dashed")
  } else {
    plt + geom_vline(xintercept = min(deer$collection_date, na.rm = T),
                     color = "steelblue4", lty = "dashed")
  }
}

combined <- ggpubr::ggarrange(plotlist = time_plots)
# ggsave("results/temporal_distribution_mutations_distinct.png", 
#        combined, 
#        dpi = 100,
#        width = 12, 
#        height = 10)

combined1 <- ggpubr::ggarrange(plotlist = time_plots[c(1, 4)], nrow = 1, ncol = 2)
combined2 <- ggpubr::ggarrange(plotlist = time_plots[5:8], nrow = 1, ncol = 4)
ggsave("results/temporal_distribution/temporal_distribution_mutations_distinct_1.png", 
       combined1, 
       dpi = 300,
       width = 13, 
       height = 4)
# ggsave("results/temporal_distribution/temporal_distribution_mutations_distinct_2.png", 
#        combined2, 
#        dpi = 300,
#        width = 13, 
#        height = 4)

