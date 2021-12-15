rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(tidyverse)
require(data.table)
require(ape)
require(ggplot2)
require(see)
require(foreach)

# prefix <- "deer.USA.n6670"
prefix <- "mink.Netherlands.n3750"
prefix <- "mink.Denmark.n10512"
# prefix <- "mink.USA.n6726"

prefixes <- c("mink.Netherlands.n3750", "mink.Denmark.n10512")

morsels <- foreach (prefix = prefixes) %do% {
  country_name <- str_split(prefix, "\\.")[[1]][2]
  meta <- fread(str_glue("data/metadata/human_animal_subsets/{prefix}.unambiguous.dedup.csv"))
  
  time_tree <- read.nexus(str_glue("results/human_animal_subsets/dating_out/{prefix}.unambiguous.dedup/timetree.nexus"))
  div_tree <- read.nexus(str_glue("results/human_animal_subsets/dating_out/{prefix}.unambiguous.dedup/divergence_tree.nexus"))
  
  # Divide time tree by divergence tree
  rate_tree <- div_tree
  rate_tree$edge.length <- rate_tree$edge.length / time_tree$edge.length
  
  n <- length(rate_tree$tip.label)
  rate_df <- tibble(accession_id = rate_tree$tip.label,
                    terminal_rate = rate_tree$edge.length[sapply(1:n, function(x,y) which(y == x),
                                            y = rate_tree$edge[,2])])
  rate_df %>% 
    left_join(meta) %>%
    add_column(country = country_name) %>%
    select(host, terminal_rate, cluster, country)
}

plot_df <- bind_rows(morsels)

plot_df %>% 
  group_by(host, country) %>%
  summarise(n = n()) %>%
  right_join(plot_df, by = c("host", "country")) %>%
  ggplot(aes(x = host, y = log(terminal_rate, base = 10), fill = host)) +
    facet_grid(rows = vars(country)) +
    theme_bw() +
    geom_point(color = "darkgrey", position = position_jitter(width = 0.08), 
               size = 0.5, 
               alpha = 0.8) +
    geom_boxplot(position = position_nudge(x = -0.2, y = 0), 
                 width = 0.2, 
                 outlier.shape = NA,
                 alpha = 0.3) +
    geom_violinhalf(position = position_nudge(x = 0.2, y = 0), alpha = 1) +
    geom_text(aes(x = host, y = -1, label = paste0("n = ", n))) +
    coord_flip() +
    labs(x = "lg(subst. rate)", y = "Host") +
    theme(legend.position = "none")

