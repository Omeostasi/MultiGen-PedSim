# =============================================================================
# time_to_generate_founders.R
# Benchmarks founderHaplotypesGenerator() across Ne values.
# Records wall-clock time (start / end) and peak RSS memory for every run.
# Output: results/time_to_generate_founders.rds  (data.frame, one row per run)
# =============================================================================

# ── Libraries ─────────────────────────────────────────────────────────────────
library(AlphaSimR)
library(tidyverse)

# ── Fixed simulation constants ─────────────────────────────────────────────────
CAP_IND      <- 1000000L   # hard ceiling on number of founders
nSnpPerChr   <- 500L       # SNPs per chromosome on the chip (≤ segSites)
nGenerations <- 10L        # max crossing generations

# ── Parameters ────────────────────────────────────────────────────────────────

nInd     <- 20000L
cat("Running with", nInd, "as starting nInd\n")

nInd_range     <- nInd
nChr_range     <- 10L
segSites_range <- 1000L
Ne_range       <- c(5000L, 10000L, 15000L, 20000L, 50000L)

# ── SLURM / environment metadata ──────────────────────────────────────────────
n_threads      <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
slurm_job_id   <- Sys.getenv("SLURM_JOB_ID",    unset = NA_character_)
slurm_nodename <- Sys.getenv("SLURMD_NODENAME",  unset = NA_character_)
run_timestamp  <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
r_version      <- paste(R.version$major, R.version$minor, sep = ".")
total_ram_gb   <- as.numeric(
  system("awk '/MemTotal/ {print $2}' /proc/meminfo", intern = TRUE)
) / 1024^2
os_info        <- utils::osVersion

cat(sprintf(
  "Threads: %d | Node: %s | Job: %s | R: %s | RAM: %.1f GB\n",
  n_threads, slurm_nodename, slurm_job_id, r_version, total_ram_gb
))

threads_range <- n_threads

# ── Parameter grid ─────────────────────────────────────────────────────────────
param_grid <- expand.grid(
  nInd     = nInd_range,
  nChr     = nChr_range,
  segSites = segSites_range,
  nThreads = threads_range,
  Ne       = Ne_range,
  stringsAsFactors = FALSE
)

# ── Memory-monitoring helper ───────────────────────────────────────────────────
# Spawns a background shell loop that appends VmRSS from /proc/<pid>/status
# every 0.5 s into a temp file.  Returns the path to that file.
start_mem_monitor <- function(pid = Sys.getpid()) {
  tmp <- tempfile(pattern = "mem_monitor_", fileext = ".txt")
  cmd <- sprintf(
    "while kill -0 %d 2>/dev/null; do grep VmRSS /proc/%d/status >> %s 2>/dev/null; sleep 0.5; done &",
    pid, pid, tmp
  )
  system(cmd)
  tmp
}

# Reads the temp file produced by start_mem_monitor() and returns peak RSS in GB.
# Waits a short moment first to capture the last poll after the workload finishes.
read_peak_mem_gb <- function(tmp_path, wait_sec = 1) {
  Sys.sleep(wait_sec)
  if (!file.exists(tmp_path) || file.size(tmp_path) == 0L) {
    warning("Memory monitor file empty or missing: ", tmp_path)
    return(NA_real_)
  }
  lines   <- readLines(tmp_path, warn = FALSE)
  kb_vals <- suppressWarnings(
    as.numeric(gsub("[^0-9]", "", lines))
  )
  kb_vals <- kb_vals[!is.na(kb_vals) & kb_vals > 0]
  if (length(kb_vals) == 0L) return(NA_real_)
  max(kb_vals) / 1024^2          # KB → GB
}

# ── Founder-haplotype generator ────────────────────────────────────────────────
# Fixed vs. original:
#   • nSnpPerChr added as formal argument (was undefined global)
#   • Cap logic: break out of generation loop once ≥ CAP_IND founders exist
#     (original kept crossing for all nGenerations even after capping)
#   • SP assigned globally only after full initialisation
founderHaplotypesGenerator <- function(
    nInd,
    nChr        = 22L,
    segSites    = 1000L,
    nThreads    = nChr,
    Ne          = 10000L,
    nGenerations = 10L,
    nSnpPerChr  = 500L,
    cap         = CAP_IND,
    method      = "runMacs2"
) {

  # ── 1. Generate base haplotypes ─────────────────────────────────────────────
  if (method == "runMacs2" || method == "run") {
    founderHap <- runMacs2(
      nInd     = nInd,
      nChr     = nChr,
      segSites = segSites,
      nThreads = nThreads,
      Ne       = Ne
    )
  } else {
    founderHap <- quickHaplo(nInd = nInd, nChr = nChr, segSites = segSites)
  }

  # ── 2. Simulation parameters ────────────────────────────────────────────────
  SP_founder <- SimParam$new(founderHap)
  SP_founder$setTrackPed(TRUE)
  SP_founder$setSexes(sexes = "yes_rand")
  SP_founder$addSnpChip(nSnpPerChr = nSnpPerChr, minSnpFreq = 0.05)

  SP <<- SP_founder   # expose globally for downstream scripts

  # ── 3. Founding population ──────────────────────────────────────────────────
  founderPop <- newPop(founderHap, simParam = SP_founder)

  # ── 4. Crossing loop ────────────────────────────────────────────────────────
  # Runs for all nGenerations without early exit.
  # The cap is a size ceiling only: once the population reaches 1M we keep
  # crossing every generation but hold it at exactly cap (no further growth).
  # Below cap we double the population each generation.
  for (generation in seq_len(nGenerations)) {

    current_n <- nInd(founderPop)
    cat(sprintf("  Generation %d: %d individuals\n", generation, current_n))

    if (current_n >= cap) {
      # At or above cap → maintain at exactly cap, keep looping
      founderPop <- randCross(
        pop      = founderPop,
        nCrosses = cap,
        simParam = SP_founder
      )
    } else {
      # Below cap → double the population
      founderPop <- randCross(
        pop      = founderPop,
        nCrosses = current_n,
        nProgeny = 2L,
        simParam = SP_founder
      )
    }
  }

  cat(sprintf("  Final population size: %d\n", nInd(founderPop)))

  # ── 5. Extract haplotype matrices per chromosome ────────────────────────────
  haplotypes  <- vector("list", length = nChr)
  pulledHaplo <- pullSegSiteHaplo(founderPop)

  start_col <- 1L
  end_col   <- segSites

  for (chr in seq_len(nChr)) {
    haplotypes[[chr]] <- pulledHaplo[, start_col:end_col, drop = FALSE]
    start_col <- end_col + 1L
    end_col   <- start_col + segSites - 1L
  }

  # ── 6. Rebuild MapPop from final founders ───────────────────────────────────
  founderHap_final <- newMapPop(
    genMap     = SP_founder$genMap,
    haplotypes = haplotypes
  )

  return(founderHap_final)
}

# ── Main benchmarking loop ─────────────────────────────────────────────────────
results_list <- vector("list", nrow(param_grid))

for (i in seq_len(nrow(param_grid))) {

  p <- param_grid[i, ]

  cat(sprintf(
    "\n[Run %d / %d]  Ne = %d | nInd = %d | nChr = %d | segSites = %d | nThreads = %d\n",
    i, nrow(param_grid), p$Ne, p$nInd, p$nChr, p$segSites, p$nThreads
  ))

  # Start background memory monitor ──────────────────────────────────────────
  mem_tmp <- start_mem_monitor()

  # Timing ───────────────────────────────────────────────────────────────────
  time_start <- Sys.time()

  founderHap_result <- founderHaplotypesGenerator(
    nInd        = p$nInd,
    nChr        = p$nChr,
    segSites    = p$segSites,
    nThreads    = p$nThreads,
    Ne          = p$Ne,
    nGenerations = nGenerations,
    nSnpPerChr  = nSnpPerChr,
    cap         = CAP_IND
  )

  time_end <- Sys.time()

  # Read peak memory (waits 1 s for last poll to land) ──────────────────────
  peak_mem_gb <- read_peak_mem_gb(mem_tmp)
  unlink(mem_tmp)   # clean up temp file

  duration_sec <- as.numeric(difftime(time_end, time_start, units = "secs"))

  cat(sprintf(
    "  → Done in %.1f s  |  Peak RSS: %.2f GB\n",
    duration_sec, peak_mem_gb
  ))

  # Collect row ──────────────────────────────────────────────────────────────
  results_list[[i]] <- tibble(
    # Parameter grid
    nInd          = p$nInd,
    nChr          = p$nChr,
    segSites      = p$segSites,
    nThreads      = p$nThreads,
    Ne            = p$Ne,
    nSnpPerChr    = nSnpPerChr,
    nGenerations  = nGenerations,
    cap_ind       = CAP_IND,
    # Timing
    time_start    = time_start,
    time_end      = time_end,
    duration_sec  = duration_sec,
    # Memory
    peak_mem_gb   = peak_mem_gb,
    # System metadata
    slurm_job_id   = slurm_job_id,
    slurm_nodename = slurm_nodename,
    run_timestamp  = run_timestamp,
    r_version      = r_version,
    total_ram_gb   = total_ram_gb,
    os_info        = as.character(os_info)
  )

  # Optional: remove large object to free memory before next run
  rm(founderHap_result)
  gc()
}

# ── Combine and save ───────────────────────────────────────────────────────────
comparison_by_param <- bind_rows(results_list)

cat("\n── Summary ──────────────────────────────────────────────────────────\n")
print(comparison_by_param |> select(Ne, duration_sec, peak_mem_gb))

mkdir_cmd <- "mkdir -p results"
system(mkdir_cmd)

saveRDS(comparison_by_param, file = "results/time_to_generate_founders.rds")
cat("\nSaved: results/time_to_generate_founders.rds\n")
