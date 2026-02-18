# Libraries

library(bench)
library(AlphaSimR)
library(tidyverse)


# Parameters

# nInd (Number of Individuals/Founders)

nInd_range <- c(1, 2, 3)

# nChr (Number of Chromosomes)

nChr_range <- c(1, 5, 10)

# segSites (Number of Segregating Sites)

segSites_range <- c(100, 500, 1000)


# Create Parameter Grid


param_grid <- expand.grid(
  nInd = nInd_range,
  nChr = nChr_range,
  segSites = segSites_range,
  stringsAsFactors = FALSE
)


# benchmark function
run_benchmark <- function(nInd, nChr, segSites, Ne = 10000) {
  
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
        Ne = Ne
      )
    },
    iterations = 10,  # Adjust based on computation time
    check = FALSE,     # Don't check if results are identical
    memory = TRUE      # Track memory usage
  )
  
  # Add parameter info
  bm$expression <- as.character(bm$expression)
  bm$nInd <- nInd
  bm$nChr <- nChr
  bm$segSites <- segSites
  bm$Ne <- Ne
  
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
  cat(sprintf("Progress: %d/%d (%.1f%%) - nInd=%d, nChr=%d, segSites=%d\n",
              i, total_runs, (i/total_runs)*100,
              param_grid$nInd[i], param_grid$nChr[i], param_grid$segSites[i]))
  
  tryCatch({
    benchmark_results[[i]] <- run_benchmark(
      nInd = param_grid$nInd[i],
      nChr = param_grid$nChr[i],
      segSites = param_grid$segSites[i],
      Ne = 10000
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
  select(expression, nInd, nChr, segSites, time_sec, mem_alloc_mb) %>%
  tidyr::pivot_wider(
    names_from = expression,
    values_from = c(time_sec, mem_alloc_mb)
  ) %>%
  mutate(
    speedup_Macs2 = time_sec_runMacs2 / time_sec_quickHaplo,
    mem_ratio_Macs2 = mem_alloc_mb_runMacs2 / mem_alloc_mb_quickHaplo
  )

