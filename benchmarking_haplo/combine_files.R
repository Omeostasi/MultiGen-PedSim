library(tidyverse)

args <- c("benchmarking_haplo/bench_nInd_Ne_2_t22.rds","benchmarking_haplo/bench_nInd_Ne_2_t10.rds", "benchmarking_haplo/bench_nInd_Ne_20000.rds", "benchmarking_haplo/bench_nInd_Ne_1_t10.rds") 

data_list <- lapply(args, readRDS)

combined_df <- bind_rows(data_list)

saveRDS(combined_df, "benchmarking_haplo/combined_results_nInd_Ne_final.rds")
