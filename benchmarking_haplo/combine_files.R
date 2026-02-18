library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

data_list <- lapply(args, readRDS)

combined_df <- bind_rows(data_list)

saveRDS(combined_df, "results/combined_results.rds")