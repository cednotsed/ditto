rm(list = ls())
setwd("../Desktop/git_repos/ditto/")
require(ape)
require(tidyverse)
require(BactDating)
require(foreach)
require(treedater)
require(lubridate)
require(data.table)
require(ggtree)
require(doParallel)
registerDoParallel(cores=10)

prefix <- "mink.Denmark.n10512"
prefix <- "mink.Netherlands.n3750"

# Load tree and metadata
meta <- fread(str_glue("data/metadata/human_animal_subsets/{prefix}.csv")) %>%
  mutate(cluster = ifelse(cluster == "N3", "N2", cluster)) %>% 
  mutate(cluster = case_when(cluster == "N4" ~ "N3", 
                             cluster == "N5" ~ "N4",
                             cluster == "N6" ~ "N5",
                             cluster == "N7" ~ "N6",
                             TRUE ~ cluster))

tree <- read.tree(str_glue("data/trees/human_animal_subsets/{prefix}.audacity_only.v8_masked.tree"))

# Parse dates
meta_filt <- meta %>% 
  mutate(collection_date = decimal_date(as.Date(collection_date))) %>%
  filter(!is.na(collection_date))

# Drop tips with missing dates
n_original <- Ntip(tree)
missing_dates <- tree$tip.label[!(tree$tip.label %in% meta_filt$accession_id)]
tree_filt <- drop.tip(tree, missing_dates)

print(str_glue("Original n_tips: {n_original}
  Tips dropped due to incomplete dates: {length(missing_dates)}
  Final n_tips: {Ntip(tree_filt)}"))

# Functions
drop_tips <- function(tre, metadat) {
  to_drop <- tre$tip.label[!(tre$tip.label %in% metadat$accession_id)]
  tre_filt <- drop.tip(tre, to_drop)
  return(tre_filt)
}

get_dates <- function(tre, metadat) {
  date_df <- metadat[match(tre$tip.label, metadat$accession_id), ]
  dates <- date_df$collection_date
  names(dates) <- date_df$accession_id
  
  return(dates)
}

# Get animal mini trees
cluster_names <- c("N1", "N2", "N4")
cluster_name <- c("N1")
cluster_name <- c("D1")

morsels <- foreach(cluster_name = cluster_names) %do% {
  animal <- meta_filt %>%
    filter(host == "Neovison vison", cluster == cluster_name)
  
  lineage_name <- unique(animal$pango_lineage)
  
  human <- meta_filt %>%
    filter(host == "Human", 
           pango_lineage == lineage_name)
  
  human_animal <- human %>% bind_rows(animal)
  
  # Get human, animal, human-animal trees
  human_tree <- drop_tips(tree_filt, human)
  animal_tree <- drop_tips(tree_filt, animal)
  
  # get dates
  human_dates <- get_dates(human_tree, human)
  animal_dates <- get_dates(animal_tree, animal)
  
  # Scale trees
  human_scaled_tree <- human_tree
  human_scaled_tree$edge.length <- human_scaled_tree$edge.length * 29872
  
  animal_scaled_tree <- animal_tree
  animal_scaled_tree$edge.length <- animal_scaled_tree$edge.length * 29872
  
  # Plot roottotip
  pdf(str_glue("results/human_animal_subsets/dating_out/{prefix}.{cluster_name}.roottotip.pdf"))
    r_animal <- roottotip(animal_scaled_tree, animal_dates)
    title(sub = str_glue("Animal {cluster_name}"))
    r_human <- roottotip(human_scaled_tree, human_dates)
    title(sub = str_glue("Human {lineage_name}"))
  dev.off()
  
  
  # test <- BactDating::bactdate(animal_scaled_tree, animal_dates, 
  #                              model = "strictgamma", 
  #                              showProgress = T,
  #                              nbIts = 10000, 
  #                              updateRoot=F)
  # test
  # plot(test,'trace')
  # TreeDater
  animal_dated <- dater(animal_tree, animal_dates,
                        s = 29872,
                        quiet = F,
                        clock = "strict",
                        maxit = 100)
  animal_pb <- parboot(animal_dated, ncpu = 10, parallel_foreach = T, quiet = F)
  
  human_dated <- dater(human_tree, human_dates,
                        s = 29872,
                        quiet = F,
                        clock = "strict",
                        maxit = 1000)
  
  human_pb <- parboot(human_dated, ncpu = 10, parallel_foreach = T, quiet = F)
  
  res <- tibble(cluster = cluster_name, host = c("Mink", "Human"),
                    tmrca = c(animal_dated$timeOfMRCA, human_dated$timeOfMRCA), 
                    tmrca_low = c(animal_pb$timeOfMRCA_CI[1], human_pb$timeOfMRCA_CI[1]), 
                    tmrca_high = c(animal_pb$timeOfMRCA_CI[2], human_pb$timeOfMRCA_CI[2]), 
                    rate = c(animal_dated$mean.rate, human_dated$mean.rate),
                    rate_low = c(animal_pb$meanRate_CI[1], human_pb$meanRate_CI[1]), 
                    rate_high = c(animal_pb$meanRate_CI[2], human_pb$meanRate_CI[2])) %>%
    mutate(across(starts_with("tmrca"), ~ date_decimal (.x)))
  
  # Visualise on tree
  meta.match <- meta[match(tree$tip.label, meta$accession_id)]
  meta.match <- meta.match %>%
    mutate(annot = case_when(accession_id %in% human$accession_id ~ "Human background", 
                             accession_id %in% animal$accession_id ~ "Mink cluster",
                             T ~ as.character(NA)))
  
  dd <- data.frame(Accession = meta.match$accession_id, 
                   annot = meta.match$annot)
  rownames(dd) <- meta.match$accession_id
  
  plt <- ggtree(tree, size = 1, color = "grey") %<+% dd +
    geom_tippoint(aes(color = annot), size = 2) +
    scale_color_discrete(na.translate = F) +
    theme(legend.title = element_blank())
  
  ggsave(str_glue("results/human_animal_subsets/dating_out/tree.{prefix}.{cluster_name}.png"),
         plot = plt)
  
  fwrite(res, str_glue("results/human_animal_subsets/dating_out/results.{prefix}.{cluster_name}.csv"))
  res
}

results <- bind_rows(morsels)

results %>% 
  ggplot(aes(y = decimal_date(tmrca), x = cluster, color = host)) +
  geom_errorbar(aes(ymin = rate_low, ymax = rate_high), position = position_dodge(width = 0.90)) +
  geom_point(size = 3, position = position_dodge(width = 0.90))

results
# animal_bactdated <- bactdate(animal_scaled_tree, animal_dates, model = "strictgamma")
  

results %>% 
  mutate(across(starts_with("tmrca"), ~ as.Date(.x))) %>%
  ggplot(aes(y = decimal_date(tmrca), x = cluster, color = host)) +
  geom_errorbar(aes(ymin = tmrca_low, ymax = tmrca_high), position = position_dodge(width = 0.90)) +
  geom_point(size = 3, position = position_dodge(width = 0.90))

as.Date(results$tmrca_low)
results
