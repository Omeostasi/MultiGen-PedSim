## founderMemory_benchmark.R
## Benchmarks MapPop-class memory usage across method × nInd × nChr × segSites
## Outputs:
##   results/memory_results.rds   – tidy data frame with all scenarios
##   results/mappops/             – individual MapPop objects (one .rds per row)

library(AlphaSimR)
library(lobstr)      # obj_size() — modern replacement for pryr::object_size()
library(tidyverse)

# ── Parameters ────────────────────────────────────────────────────────────────
# nInd_range  <- c(1000, 10000, 1e5, 1e6)
nInd_range  <- c(1000, 10000, 1e5, 1e6, 5e6)

# nChr_range  <- c(10, 22)
nChr_range  <- c(10)
segSites_range <- c(1000)
methods     <- c("quick")

# Output directories
dir.create("results/mappops", recursive = TRUE, showWarnings = FALSE)

# ── Helper: measure object size in MB ─────────────────────────────────────────
mem_mb <- function(obj) {
  as.numeric(lobstr::obj_size(obj)) / 1024^2
}

# ── MapPop generator for "quick" method ───────────────────────────────────────
generate_quick <- function(nInd, nChr, segSites) {
  quickHaplo(nInd = nInd, nChr = nChr, segSites = segSites)
}

# ── MapPop generator for "macs" method ────────────────────────────────────────
# Starts from 1000 founders with runMacs2, then propagates via random crosses
# until nInd_target is reached, then reconstructs a MapPop from haplotypes.
generate_macs <- function(nInd_target, nChr, segSites,
                          Ne = 10000, nGenerations = 10) {

  message(sprintf("  [macs] runMacs2: nInd=1000, nChr=%d, segSites=%d", nChr, segSites))
  founderHap <- runMacs2(
    nInd     = 1000,
    nChr     = nChr,
    segSites = segSites,
    nThreads = nChr,
    Ne       = Ne
  )

  SP_founder <- SimParam$new(founderHap)
  SP_founder$setTrackPed(TRUE)
  SP_founder$setSexes(sexes = "yes_rand")

  # Assign to global SP so AlphaSimR internals that rely on SP work
  SP <<- SP_founder

  founderPop <- newPop(founderHap, simParam = SP_founder)

  message(sprintf("  [macs] Crossing to nInd=%d ...", nInd_target))
  for (generation in seq_len(nGenerations)) {
    current_n <- nInd(founderPop)

    if (current_n >= nInd_target) {
      founderPop <- randCross(
        pop      = founderPop,
        nCrosses = nInd_target,
        simParam = SP_founder
      )
      break   # we have exactly nInd_target individuals → done
    } else {
      founderPop <- randCross(
        pop      = founderPop,
        nCrosses = current_n,
        nProgeny = 2,
        simParam = SP_founder
      )
    }
  }

  # If still under target after all generations, do a final cross to hit target
  if (nInd(founderPop) != nInd_target) {
    founderPop <- randCross(
      pop      = founderPop,
      nCrosses = nInd_target,
      simParam = SP_founder
    )
  }

  message(sprintf("  [macs] Extracting haplotypes for %d individuals ...", nInd(founderPop)))

  # Rebuild MapPop from pulled haplotypes (mirrors the original RMD logic)
  pulledHaplo <- pullSegSiteHaplo(founderPop, simParam = SP_founder)
  haplotypes  <- vector("list", length = nChr)

  start_col <- 1
  for (chr in seq_len(nChr)) {
    end_col          <- start_col + segSites - 1
    haplotypes[[chr]] <- pulledHaplo[, start_col:end_col, drop = FALSE]
    start_col        <- end_col + 1
  }

  mapPop <- newMapPop(genMap = SP_founder$genMap, haplotypes = haplotypes)
  return(mapPop)
}

# ── Build scenario grid ────────────────────────────────────────────────────────
scenarios <- expand.grid(
  method    = methods,
  nInd      = nInd_range,
  nChr      = nChr_range,
  segSites  = segSites_range,
  stringsAsFactors = FALSE
)

message(sprintf("Total scenarios: %d", nrow(scenarios)))

# ── Main loop ─────────────────────────────────────────────────────────────────
results_list <- vector("list", nrow(scenarios))

for (i in seq_len(nrow(scenarios))) {

  m  <- scenarios$method[i]
  n  <- scenarios$nInd[i]
  nc <- scenarios$nChr[i]
  ss <- scenarios$segSites[i]

  message(sprintf("\n[%d/%d] method=%s | nInd=%d | nChr=%d | segSites=%d",
                  i, nrow(scenarios), m, n, nc, ss))

  mapPop <- tryCatch({
    if (m == "quick") {
      generate_quick(nInd = n, nChr = nc, segSites = ss)
    } else {
      generate_macs(nInd_target = n, nChr = nc, segSites = ss)
    }
  }, error = function(e) {
    message("  ERROR: ", conditionMessage(e))
    NULL
  })

  mb <- if (!is.null(mapPop)) mem_mb(mapPop) else NA_real_

  results_list[[i]] <- tibble(
    method   = m,
    nInd     = n,
    nChr     = nc,
    segSites = ss,
    mem_mb   = mb
  )

  # Optionally save the MapPop object
  if (!is.null(mapPop)) {
    rds_name <- sprintf("memory_usage/foundersHaplo/results/mappops/mappop_%s_nInd%d_nChr%d_seg%d.rds",
                        m, n, nc, ss)
    saveRDS(mapPop, rds_name)
    message(sprintf("  Saved MapPop → %s  (%.1f MB)", rds_name, mb))
  }

  # Free memory between iterations
  rm(mapPop)
  gc(verbose = FALSE)
}

# ── Combine and save ───────────────────────────────────────────────────────────
results_df <- bind_rows(results_list)

print(results_df, n = Inf)

saveRDS(results_df, "results/memory_results.rds")
message("\nDone! Results saved to results/memory_results.rds")

