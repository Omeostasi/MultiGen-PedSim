# Libraries

library(bench)
library(AlphaSimR)
library(tidyverse)

# Parameters

# nInd (Number of Individuals/Founders)

nInd_range <- c(100)

# nChr (Number of Chromosomes)

nChr_range <- c(8, 10, 16, 23)

# segSites (Number of Segregating Sites)

segSites_range <- c(100)

# threads ranges

threads_range <- c(4, 8, 10, 16, 23)

# effective pop size
Ne_range <- c(1000)

# Create Parameter Grid


param_grid <- expand.grid(
  nInd = nInd_range,
  nChr = nChr_range,
  segSites = segSites_range,
  nThreads = threads_range,
  Ne = Ne_range,
  stringsAsFactors = FALSE
)

# benchmark function
run_benchmark <- function(nInd, nChr, segSites, Ne, nThreads) {
  
  # Benchmark both methods
  bm <- bench::mark(
    quickHaplo = {
      quickHaplo(
        nInd = nInd,
        nChr = nChr,
        segSites = segSites
      )
    },
    runMacs2 = {
      runMacs2(
        nInd = nInd,
        nChr = nChr,
        segSites = segSites,
        Ne = Ne,
        nThreads = nThreads # adding multiple cores
      )
    },
    iterations = 3,  
    check = FALSE,     # Don't check if results are identical
    memory = TRUE      # Track memory usage
  )
  
  # Add parameter info
  bm$expression <- as.character(bm$expression)
  bm$nInd <- nInd
  bm$nChr <- nChr
  bm$segSites <- segSites
  bm$Ne <- Ne
  bm$nThreads <- nThreads
  
  return(bm)
}


# Loop

# Initialize results list
benchmark_results <- list()

# Progress tracking
total_runs <- nrow(param_grid)
cat("Starting benchmark with", total_runs, "parameter combinations\n")

# Run benchmarks for each parameter combination
for (i in 1:nrow(param_grid)) {
  cat(sprintf("Progress: %d/%d (%.1f%%) - nInd=%d, nChr=%d, segSites=%d, Ne=%d, nThreads=%d\n",
              i, total_runs, (i/total_runs)*100,
              param_grid$nInd[i], param_grid$nChr[i], param_grid$segSites[i], param_grid$Ne[i],param_grid$nThreads[i]))
  
  tryCatch({
    benchmark_results[[i]] <- run_benchmark(
      nInd = param_grid$nInd[i],
      nChr = param_grid$nChr[i],
      segSites = param_grid$segSites[i],
      Ne = param_grid$Ne[i],
      nThreads = param_grid$nThreads[i]
    )
  },
  # If an error occurs, print a message and store NULL for this iteration
  error = function(e) {
    cat("Error in iteration", i, ":", conditionMessage(e), "\n")
    benchmark_results[[i]] <- NULL
  })
}


# Combine all benchmark results into single dataframe
results_df <- bind_rows(benchmark_results)

# Convert time to numeric (seconds)
results_df <- results_df %>%
  mutate(
    time_sec = as.numeric(median),
    mem_alloc_mb = as.numeric(mem_alloc) / 1024^2
  )




comparison_by_param <- results_df %>%
  select(expression, nInd, nChr, segSites, time_sec, mem_alloc_mb, nThreads, Ne) %>%
  tidyr::pivot_wider(
    names_from = expression,
    values_from = c(time_sec, mem_alloc_mb)
  ) %>%
  mutate(
    speedup_Macs2 = time_sec_runMacs2 / time_sec_quickHaplo,
    mem_ratio_Macs2 = mem_alloc_mb_runMacs2 / mem_alloc_mb_quickHaplo
  )

saveRDS(comparison_by_param, file="results/bench_threads_chr.rds")

