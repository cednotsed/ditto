rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(lubridate)
require(vistime)

df <- fread("data/metadata/human_animal_subsets/V6/mink_deer.n16688.csv") %>%
  filter(host != "Human") %>%
  as_tibble()

country_stats <- df %>%
  mutate(collection_date = if_else(location == "Europe / Lithuania", 
                                   ymd("2020-11-26"), 
                                   as.Date(collection_date, "%Y-%m-%d"))) %>%
  group_by(location, host) %>% 
  summarise(n_total = n(), 
            start = min(collection_date, na.rm = T),
            end = max(collection_date, na.rm = T)) %>%
  separate(location, into = c(NA, "country"), sep = " / ")

plot_df <- country_stats %>%
  bind_rows(tibble(country = c("Pandemic declared", "Emergence in humans", "Audacity release in present study",
                               "Earliest Alpha", "Earliest Beta", "Earliest Gamma", "Earliest Delta", "Earliest Omicron"), 
                   host = "Pandemic timeline",
                   start = as.Date(c("2020-03-11", "2019-10-06", "2022-03-17",
                                     "2020-09-15", "2020-05-15", "2020-11-15", 
                                     "2020-10-15", "2021-11-11")),
                   end = as.Date(c("2020-03-11", "2019-12-11", "2022-03-17",
                                   "2020-09-15", "2020-05-15", "2020-11-15", 
                                   "2020-10-15", "2021-11-11"))))

# Annotate mutation times
time_meta <- fread("results/temporal_distribution/mutation_times.csv") %>%
  as_tibble() %>%
  filter(!(mut %in% c("N_D377Y", "M_I82T")))

time_plt <- gg_vistime(plot_df, 
           col.event = "country",
           col.group = "host") +
  scale_x_datetime(breaks = seq(min(as.POSIXct(plot_df$start)), 
                                max(as.POSIXct(plot_df$end)), 
                                "2 months"), 
                   date_labels = "%b-%y") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

for (i in seq(nrow(time_meta))) {
  time_plt <- time_plt +
    geom_point(x = as.POSIXct(time_meta[i, ]$earliest), y = 4, color = "black")
    # geom_text(x = as.POSIXct(time_meta[i, ]$earliest),
    #           y = 4,
    #           label = time_meta[i, ]$mut,
    #           angle = 45,
    #           size = 1,
    #           color = "black")
}

time_plt
ggsave("results/temporal_distribution/timeline.nolab.png", dpi = 300, width = 11, height = 5)
