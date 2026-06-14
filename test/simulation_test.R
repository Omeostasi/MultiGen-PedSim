##--------------------------------------------------------------------------------------------------------------------------------------------------------------
## QUICKTEST CAN BE RUN LOCALLY
## -------------------------------------------------------------------------------------------------------------------------------------------------------------


## ============================================================
## SANITY CHECK — packages + source files
## ============================================================
cat("==========================================\n")
cat("SANITY CHECK\n")
cat("==========================================\n")

required_packages <- c("AlphaSimR", "MASS", "Matrix", "data.table")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("[FAIL] Required package not found: ", pkg,
                "\n       Install with: install.packages('", pkg, "')"))
  }
  cat(sprintf("  [OK] Package %-12s v%s\n", paste0(pkg, ":"), as.character(packageVersion(pkg))))
}


required_src <- c(
  "src/founderHaplotypesGenerator.R",
  "src/helpers.R",
  "src/pop_generation.R",
  "src/subsetting.R",
  "src/my_writePlink.R"
)
for (f in required_src) {
  if (!file.exists(f)) {
    stop(paste0("[FAIL] Source file not found: ", f,
                "\n       Make sure the 'src/' folder is in your working directory: ", getwd()))
  }
  cat(sprintf("  [OK] Source file found: %s\n", f))
}

cat("==========================================\n")
cat("SANITY CHECK PASSED\n")
cat("==========================================\n\n")

## ============================================================
## Load libraries + source files
## ============================================================
library(AlphaSimR)
library(MASS)
library(Matrix)

source("src/founderHaplotypesGenerator.R")   # wrapper for founder haplotype generation
source("src/helpers.R")                      # make_pd_corr, calib_alpha, alloc_kids_blocks
source("src/pop_generation.R")               # pop pedigree generator
source("src/subsetting.R")                   # creates subsets of final population and genotyped population
source("src/my_writePlink.R")                # correct writePlink function for pop objects that do not have founders included


## ----------------------------
## EXTRACTING THE ARGUMENTS
## ---------------------------
args <- commandArgs(trailingOnly = TRUE)
cat(sprintf("[DEBUG] length(args)=%d | raw args: %s\n",
            length(args), paste(args, collapse=", "))) 
# Each run_id creates its own subfolder
run_id <- as.integer(if (length(args) >= 1) args[1] else 1L)  # to count the runs
set.seed(run_id)                                              # the seed changes with the run_id

## ============================================================
## Output directory  —  one sub-folder per run
## ============================================================
OUT_DIR_POP   <- file.path("results")
OUT_DIR_PHENO <- file.path("results", sprintf("run%d", run_id))
dir.create(OUT_DIR_POP, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR_PHENO, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("[INFO] Output directory for POP: %s/\n\n", normalizePath(OUT_DIR_POP)))
cat(sprintf("[INFO] Output directory for PHENO: %s/\n\n", normalizePath(OUT_DIR_PHENO)))

## ============================================================
## Start total timer
## ============================================================
t_total_start <- proc.time()

## -------------------------
## 0) Parameters
## -------------------------

## Genome simulation
nFounder              <- 50   # initial number of founders for runMacs2
nChr                  <- 2          # number of chromosomes
segSites              <- 1500        # segregating sites per chromosome
nSnpPerChr            <- 20         # SNPs per chromosome on chip (observed)
nGenerations_founders <- 10          # generations to build the final founders
Ne                    <- 100   # effective pop size
nFounders_final       <- 20 * 10^1    # final number of founders generated through random crossing
minSnpFreq            <- 0.005       # minor allele frequency to be considered when adding snps on chip (for both founders and pop)

## Population / study design
# NOTE: simulation uses per-generation vectors for lambdaKids, p_new_partner AND mothers_fraction
mothers_fraction     <- c(0.99, 0.99, 0.98, 0.82, 0.71)          # fraction of female pop that become mothers. Each value is a different generation
lambdaKids           <- c(3.00, 2.53, 2.11, 1.64, 1.34)          # mean kids per mother per generation (zero-truncated Poisson). Estimated from the real population
nGenerations_pop     <- 5                                        # number of generations with data
overlapping_fraction <- 0.00                                     # If more than 0 allows for a fraction of previous generation to cross during current generation
rm_older_generations <- TRUE                                     # Default it's TRUE,  it removes 80% first generation and founders after generating them and computing their genetic components
# Union / partnership structure
maxPartners_mom  <- 3
# NOTE: simulation uses per-generation vector for p_new_partner and mean_unions_dad
p_new_partner    <- c(0.02, 0.12, 0.13, 0.12, 0.10)   # Probability of new partner per generation for a mother. Estimated from the real population.
mean_unions_dad  <- c(1.25, 1.17, 1.15, 1.13, 1.09)   # higher => more paternal half-sibs + higher dad reuse. Estimated from real population

## -----------------------------------------------------------------
## Phenotype and genetic components
## These parameters are used after the generation of the population
## -----------------------------------------------------------------

## Prevalences
prevCMC_mother <- as.numeric(if (length(args) >= 2) args[2] else 0.10)  # lifetime prevalence of cardiometabolic condition in mothers
prevPregCMC    <- as.numeric(if (length(args) >= 3) args[3] else 0.16)  # prevalence of pregnancy-related CMC episodes, calibrated per mother
prevASC_child  <- as.numeric(if (length(args) >= 4) args[4] else 0.10)  # population prevalence of ASC in children (both sexes combined)

# Prevalences ASC by gender
ASC_male_ratio <- as.numeric(if (length(args) >= 5) args[5] else 5/1)  # male-to-female ratio for ASC diagnosis. Default: 4 males every 1 female
prevASC_multiplier <- 1/ASC_male_ratio*2  
prevASC_male   <- prevASC_child*prevASC_multiplier*(ASC_male_ratio-1)
prevASC_female <- prevASC_child*prevASC_multiplier


## Pregnancy effect on ASC liability
betaPreg <- as.numeric(if (length(args) >= 6) args[6] else 0.35)  # effect of pregnancy CMC exposure on child ASC liability (liability scale)

## ASC liability variance targets (excluding pregnancy term)
var_d <- as.numeric(if (length(args) >= 7) args[7] else 0.20)  # child's direct additive genetic effects
var_m <- as.numeric(if (length(args) >= 8) args[8] else 0.10)  # maternal additive genetic effects on child (genetic nurture)
var_c <- as.numeric(if (length(args) >= 9) args[9] else 0.10)  # shared maternal environment (full-siblings share this)
var_e <- 1 - (var_d + var_m + var_c)                           # residual variance; must be positive
stopifnot(var_e > 0)

## Correlations among causal SNP effects across components
rho_CMC_d <- as.numeric(if (length(args) >= 10) args[10] else 0.30)  # correlation between SNP effects on CMC and child direct ASC liability
rho_CMC_m <- as.numeric(if (length(args) >= 11) args[11] else 0.20)  # correlation between SNP effects on CMC and maternal indirect ASC liability
rho_d_m   <- as.numeric(if (length(args) >= 12) args[12] else 0.10)  # correlation between child direct and maternal indirect ASC genetic components

## Causal architecture
nCausalPerChr <- 10   # causal seg sites per chromosome


## Pregnancy CMC model
kappa     <- 1.0
delta_CMC <- 1.0

## Genotyped subset targets (iPSYCH-like design)
n_geno_ASC    <- 20    # target number of ASC cases to genotype
n_geno_random <- 250   # target number of random non-ASC individuals to genotype
## ---------------------
## Parameter validation
## ---------------------
if (
  length(lambdaKids) != nGenerations_pop ||
  length(p_new_partner) != nGenerations_pop ||
  length(mothers_fraction) != nGenerations_pop ||
  length(mean_unions_dad) != nGenerations_pop
) {
  stop("Mismatch between number of generations and parameter vector lengths")
}




cat("Reporting assigned parameters to the log file: \n")
cat(sprintf("[INFO] Run ID: %d | seed: %d\n", run_id, run_id))
cat(sprintf("[INFO] prevCMC_mother=%.3f | prevPregCMC=%.3f | prevASC_child=%.3f\n",
            prevCMC_mother, prevPregCMC, prevASC_child))
cat(sprintf("[INFO] var_d=%.2f | var_m=%.2f | var_c=%.2f | var_e=%.2f\n",
            var_d, var_m, var_c, var_e))
cat(sprintf("[INFO] rho_CMC_d=%.2f | rho_CMC_m=%.2f | rho_d_m=%.2f\n",
            rho_CMC_d, rho_CMC_m, rho_d_m))
cat(sprintf("[INFO] nCausalPerChr=%d | kappa=%.2f | delta_CMC=%.2f\n",
            nCausalPerChr, kappa, delta_CMC))
cat("Reporting default parameters: \n")
cat(sprintf("segSites: %g\n", segSites))
cat(sprintf("nSnpPerChr: %g\n", nSnpPerChr))
cat(sprintf("Ne: %g\n", Ne))
cat(sprintf("minSnpFreq: %g\n", minSnpFreq))
cat(sprintf("nGenerations_pop: %g\n", nGenerations_pop))
cat(sprintf("nFounders_final: %g\n", nFounders_final))
cat("run \n")

## Filenames
## pop_description  : identifies the simulated population object (genome/demography params)
pop_description   <- sprintf("Ne%g_seg%g_mSnpF%g_Founders%g_nGen%g",
                             Ne, segSites, minSnpFreq, nFounders_final, nGenerations_pop)
out_path_pop <- file.path(OUT_DIR_POP, sprintf("pop_%s.rds", pop_description))


if (!file.exists(out_path_pop)){
  ## -------------------------
  ## 1) Generate founder haplotypes (AlphaSimR)
  ## -------------------------
  cat("------------------------------------------\n")
  cat("[STEP 1] Generating founder haplotypes ...\n")
  cat("------------------------------------------\n")
  
  t_founder_start <- proc.time()
  
  founderHap <- founderHaplotypesGenerator(
    nInd            = nFounder,
    nChr            = nChr,
    segSites        = segSites,
    nThreads        = nChr,
    Ne              = Ne,
    nGenerations    = nGenerations_founders,
    method          = "run",
    nFounders_final = nFounders_final,
    nSnpPerChr      = nSnpPerChr,
    minSnpFreq      = minSnpFreq
  )
  
  t_founder_elapsed <- proc.time() - t_founder_start
  cat(sprintf("\n[TIME] founderHaplotypesGenerator : %.1f sec (%.2f min)\n\n",
              t_founder_elapsed["elapsed"],
              t_founder_elapsed["elapsed"] / 60))
  
  ## Memory: founderHap (MapPop)
  cat("[MEMORY] Objects after Step 1:\n")
  mem_report(founderHap, "founderHap (MapPop)")
  
  SP <- SimParam$new(founderHap)
  SP$setTrackPed(TRUE)
  SP$setSexes(sexes = "yes_rand")
  SP$addSnpChip(nSnpPerChr = nSnpPerChr, minSnpFreq = minSnpFreq)
  
  
  ## -------------------------
  ## 2) Build family / pedigree structure
  ## -------------------------
  cat("------------------------------------------\n")
  cat("[STEP 2] Building population pedigree ...\n")
  cat("------------------------------------------\n")
  
  t_pop_start <- proc.time()
  
  pop_all <- pop_generation(
    founderHap,
    overlapping_fraction = overlapping_fraction,
    nGenerations_pop     = nGenerations_pop,
    mothers_fraction     = mothers_fraction,
    lambdaKids           = lambdaKids,
    maxPartners_mom      = maxPartners_mom,
    p_new_partner        = p_new_partner,
    mean_unions_dad      = mean_unions_dad
  )
  
  t_pop_elapsed <- proc.time() - t_pop_start
  cat(sprintf("\n[TIME] pop_generation             : %.1f sec (%.2f min)\n\n",
              t_pop_elapsed["elapsed"],
              t_pop_elapsed["elapsed"] / 60))
  
  
  
  cat("------------------------------------------\n")
  cat("[STEP 2b] Saving population object ...\n")
  cat("------------------------------------------\n")
  cat(sprintf("  [INFO] File : %s\n", pop_description))
  cat(sprintf("  [INFO] Path : %s\n", normalizePath(out_path_pop, mustWork = FALSE)))
  
  # Save SP and founderHap using the same naming convention as pop_all
  out_path_SP         <- file.path(OUT_DIR_POP, sprintf("SP_%s.rds",         pop_description))
  out_path_founderHap <- file.path(OUT_DIR_POP, sprintf("founderHap_%s.rds", pop_description))
  
  saveRDS(SP,         out_path_SP)
  saveRDS(founderHap, out_path_founderHap)
  
  cat(sprintf("  [OK] SP saved to         : %s\n", normalizePath(out_path_SP,         mustWork = FALSE)))
  cat(sprintf("  [OK] founderHap saved to : %s\n", normalizePath(out_path_founderHap, mustWork = FALSE)))
  
  saveRDS(pop_all, out_path_pop)
} else {
  cat(sprintf("pop_%s.rds already exists!\n", pop_description))
  cat("Avoiding pop_generation steps\n")
  cat("Unpacking .rds and proceeding with phenotype and genotyping simulation\n")
  
  pop_all <- readRDS(out_path_pop)
  
  # Restore SP and founderHap
  out_path_SP         <- file.path(OUT_DIR_POP, sprintf("SP_%s.rds",         pop_description))
  
  if (!file.exists(out_path_SP)) {
    stop(sprintf(
      "[FAIL] Population .rds found but SP or founderHap is missing.\n  Expected: %s \n",
      out_path_SP
    ))
  }
  
  SP         <- readRDS(out_path_SP)
  
  cat(sprintf("  [OK] SP restored from         : %s\n", normalizePath(out_path_SP)))
}

## Unpack pop_all
pop              <- pop_all$pop
unions_all       <- pop_all$unions_all
children_ids_all <- unique(na.omit(unlist(pop_all$children_ids_list)))

## Subsets
mothers_ids_all  <- unique(na.omit(unions_all$mother_id))
mothers          <- pop[mothers_ids_all]
children         <- pop[children_ids_all]

## Memory: pop objects
cat("[MEMORY] Objects after Step 2:\n")
mem_report(pop,      "pop (Pop — full population - before removal)")
mem_report(mothers,  "mothers (Pop subset - before removal)")
mem_report(children, "children (Pop subset - before removal)")

cat(sprintf("nInd (pop - before removal): %d\n", nInd(pop)))
cat(sprintf("nInd (mothers - before removal): %d\n", nInd(mothers)))
cat(sprintf("nInd (children - before removal): %d\n", nInd(children)))

## -------------------------
## 3) Pull genotypes
## -------------------------
t_genetic_start <- proc.time()

cat("------------------------------------------\n")
cat("[STEP 3] Pulling genotypes ...\n")
cat("------------------------------------------\n")
cat("Pulling segSites, Genotypes and choosing causal variants from segSites\n")



Gchip_mom <- pullSnpGeno(mothers,  snpChip = 1, simParam = SP)
Gchip_kid <- pullSnpGeno(children, snpChip = 1, simParam = SP)

Gseg_mom  <- pullSegSiteGeno(mothers,  simParam = SP)
Gseg_kid  <- pullSegSiteGeno(children, simParam = SP)

## Causal variants
nCausal <- nChr * nCausalPerChr
if (nCausal > ncol(Gseg_mom))
  stop("nCausal exceeds available seg sites. Reduce nCausalPerChr or increase segSites.")
causal_idx   <- sample(seq_len(ncol(Gseg_mom)), nCausal)
Gcausal_mom  <- Gseg_mom[, causal_idx, drop = FALSE]
Gcausal_kid  <- Gseg_kid[, causal_idx, drop = FALSE]
m_causal     <- ncol(Gcausal_mom)

cat("[MEMORY] Genotype matrices:\n")
mem_report(Gchip_mom,   "Gchip_mom (chip SNPs — mothers)")
mem_report(Gchip_kid,   "Gchip_kid (chip SNPs — children)")
mem_report(Gseg_mom,    "Gseg_mom  (seg sites — mothers)")
mem_report(Gseg_kid,    "Gseg_kid  (seg sites — children)")
mem_report(Gcausal_mom, "Gcausal_mom (causal — mothers)")
mem_report(Gcausal_kid, "Gcausal_kid (causal — children)")

rm(Gseg_mom, Gseg_kid) # to save memory
## -------------------------
## 4) Correlated causal SNP effects
## -------------------------
Sigma <- matrix(
  c(1,         rho_CMC_d, rho_CMC_m,
    rho_CMC_d, 1,         rho_d_m,
    rho_CMC_m, rho_d_m,   1),
  nrow = 3, byrow = TRUE
)
Sigma <- make_pd_corr(Sigma)

eff      <- MASS::mvrnorm(n = m_causal, mu = c(0, 0, 0), Sigma = Sigma)
b_cmc    <- eff[, 1]
b_d      <- eff[, 2]
b_m      <- eff[, 3]

prs_cmc_mom <- as.numeric(scale(Gcausal_mom %*% b_cmc))
prs_m_mom   <- as.numeric(scale(Gcausal_mom %*% b_m))
prs_d_kid   <- as.numeric(scale(Gcausal_kid %*% b_d))

prs_m_mom_byID <- setNames(prs_m_mom, as.character(mothers_ids_all))

# Collects mother_ids per child
child_moms      <- rep(unions_all$mother_id, times = unions_all$n_kids)
child_dads      <- rep(unions_all$father_id, times = unions_all$n_kids)

# Associate m prs to child
prs_m_for_child <- prs_m_mom_byID[as.character(child_moms)]


## -------------------------
## 5) Maternal CMC
## -------------------------
L_CMC_mom <- sqrt(0.40) * prs_cmc_mom + sqrt(0.60) * rnorm(length(mothers_ids_all))
thr_CMC   <- as.numeric(quantile(L_CMC_mom, probs = 1 - prevCMC_mother))
CMC_mom   <- as.integer(L_CMC_mom > thr_CMC)


## -------------------------
## 6) Pregnancy CMC (per child)
## -------------------------
linpred_mom <- kappa * L_CMC_mom + delta_CMC * CMC_mom
alpha       <- calib_alpha(prevPregCMC, linpred_mom)
p_mom       <- plogis(alpha + linpred_mom)
p_for_child <- p_mom[match(as.character(child_moms), as.character(mothers_ids_all))]
preg_CMC    <- rbinom(n = length(p_for_child), size = 1, prob = p_for_child)



## -------------------------
## 7) Child ASC liability
## -------------------------
c_mom      <- rnorm(length(mothers_ids_all), mean = 0, sd = 1)
c_mom_byID <- setNames(c_mom, as.character(mothers_ids_all))
c_for_child <- c_mom_byID[as.character(child_moms)]

gd_z <- as.numeric(scale(prs_d_kid))
gm_z <- as.numeric(scale(prs_m_for_child))
c_z  <- as.numeric(scale(c_for_child))

nonpreg_part <- sqrt(var_d) * gd_z +
  sqrt(var_m) * gm_z +
  sqrt(var_c) * c_z

sd_e      <- sqrt(max(1 - var(nonpreg_part, na.rm = TRUE), 1e-8))
L_ASC_kid <- betaPreg * preg_CMC +
  nonpreg_part +
  sd_e * rnorm(length(nonpreg_part))

# Sex-specific thresholds
child_sex <- children@sex

is_male   <- child_sex == "M"
is_female <- child_sex == "F"

# Each sex gets its own threshold calibrated to its target prevalence
thr_ASC_male   <- quantile(L_ASC_kid[is_male],   probs = 1 - prevASC_male,   na.rm = TRUE)
thr_ASC_female <- quantile(L_ASC_kid[is_female],  probs = 1 - prevASC_female, na.rm = TRUE)

cat(sprintf("  [INFO] ASC threshold (male)   : %.4f  [target prev = %.3f]\n", thr_ASC_male,   prevASC_male))
cat(sprintf("  [INFO] ASC threshold (female) : %.4f  [target prev = %.3f]\n", thr_ASC_female, prevASC_female))


ASC_kid <- integer(length(L_ASC_kid))
ASC_kid[is_male]   <- as.integer(L_ASC_kid[is_male]   > thr_ASC_male)
ASC_kid[is_female] <- as.integer(L_ASC_kid[is_female] > thr_ASC_female)




## -------------------------
## 8) Export outputs to results/
## -------------------------
cat("------------------------------------------\n")
cat("[STEP 8] Writing output files ...\n")
cat("------------------------------------------\n")

ped <- getPed(pop)     


pheno_mom <- data.frame(
  id       = mothers@id,
  role     = "mother",
  CMC_liab = L_CMC_mom,
  CMC      = CMC_mom
)

pheno_kid <- data.frame(
  id        = children@id,
  role      = "child",
  mother_id = child_moms,
  father_id = child_dads,
  preg_CMC  = preg_CMC,
  ASC_liab  = L_ASC_kid,
  ASC       = ASC_kid,
  sex       = as.integer(1)
)

# female get sex assigned as 2
pheno_kid[is_female,]$sex <- 2


## ---------------
## Removing older generations from the set
## -----------------

cat("------------------------------------------\n")
cat("[STEP 8b] Removing older generations ...\n")
cat("------------------------------------------\n")

subset <- remove_older_generations(pheno_mom = pheno_mom,
                                   pheno_kid = pheno_kid,
                                   ped = ped,
                                   older_generations_ids = pop_all$children_ids_list, 
                                   rm_older_generations = rm_older_generations)

ped <- subset$ped
pheno_mom <- subset$pheno_mom
pheno_kid <- subset$pheno_kid


# Filenames
pheno_description <- sprintf("vd%g_vm%g_vc%g_rcd%g_rcm%g_rdm%g_run%02d",
                             var_d, var_m, var_c, rho_CMC_d, rho_CMC_m, rho_d_m, run_id)

ped_file   <- sprintf("sim_pedigree_%s.csv",           pop_description)
mom_file   <- sprintf("sim_mom_phenotypes_%s.csv",     pheno_description)
kid_file   <- sprintf("sim_child_phenotypes_%s.csv",   pheno_description)
plink_base <- sprintf("sim_genotyped_subset_%s",       pheno_description)

# Write CSVs
write.csv(ped,       file = file.path(OUT_DIR_PHENO, ped_file), row.names = FALSE)
write.csv(pheno_mom, file = file.path(OUT_DIR_PHENO, mom_file), row.names = FALSE)
write.csv(pheno_kid, file = file.path(OUT_DIR_PHENO, kid_file), row.names = FALSE)

cat(sprintf("  [OK] %s\n", ped_file))
cat(sprintf("  [OK] %s\n", mom_file))
cat(sprintf("  [OK] %s\n", kid_file))

#  ----
# PLINK export for a genotyped subset (iPSYCH-like design)
#  ----
cat("--------------------------------------------------\n")
cat("[STEP 8c] Creating iPSYCH-like cohort to genotype \n")
cat("--------------------------------------------------\n")

# TO WORK IT REQUIRES AT LEAST 2 GENERATIONS. IT TAKES A SUBSET FROM THE LAST TWO GENERATIONS
pop_geno <- genotype_subset(pop = pop, pheno_kid = pheno_kid, generations_ids = pop_all$children_ids_list,
                            ped = ped,
                            nGenerations_pop = nGenerations_pop,
                            n_geno_random = n_geno_random, n_geno_ASC = n_geno_ASC)

## write PLINK in chunks to stay within cat() size limits
chunk_size  <- 18000   # individuals per chunk
n_geno_tot   <- nInd(pop_geno)
chunk_starts <- seq(1, n_geno_tot, by = chunk_size)
chunk_files  <- character(length(chunk_starts))

cat(sprintf("  [INFO] Writing %d PLINK chunks (chunk_size = %d) ...\n",
            length(chunk_starts), chunk_size))

for (i in seq_along(chunk_starts)) {
  idx_from <- chunk_starts[i]
  idx_to   <- min(chunk_starts[i] + chunk_size - 1L, n_geno_tot)
  
  chunk_base    <- file.path(OUT_DIR_PHENO, sprintf("%s_chunk%02d", plink_base, i))
  chunk_files[i] <- chunk_base
  
  my_writePlink(pop_geno[idx_from:idx_to], baseName = chunk_base, simParam = SP, use = "rand")
  
  cat(sprintf("    [OK] chunk %02d  (rows %d-%d)\n", i, idx_from, idx_to))
}

ped_files <- paste0(chunk_files, ".ped")
map_files <- paste0(chunk_files, ".map")
final_ped <- file.path(OUT_DIR_PHENO, paste0(plink_base, ".ped"))
final_map <- file.path(OUT_DIR_PHENO, paste0(plink_base, ".map"))

# Stream-combine chunks; patch parent IDs on the fly
con_out <- file(final_ped, open = "wt") # it also creates the file
for (pf in ped_files) {
  con_in <- file(pf, open = "rt")
  while (TRUE) {
    lines <- readLines(con_in, n = -1, warn = FALSE)
    if (length(lines) == 0L) break
    lines <- sapply(lines, function(line) {
      fields <- strsplit(line, "\\s+")[[1]]
      fields <- fields[nchar(fields) > 0]
      paste(fields, collapse = " ")
    }, USE.NAMES = FALSE)
    writeLines(lines, con_out)
  }
  close(con_in)
}
close(con_out)

file.copy(map_files[1], final_map, overwrite = TRUE)
invisible(file.remove(c(ped_files, map_files)))
cat(sprintf("  [OK] %s.{ped/map} written and chunk files removed\n\n", plink_base))


## -------------------------
## 9) Save run summary
## -------------------------
cat("------------------------------------------\n")
cat("[STEP 9] Saving run summary ...\n")
cat("------------------------------------------\n")

save_run_summary(
  OUT_DIR        = OUT_DIR_PHENO,
  run_id         = run_id,
  file_description = pheno_description,
  nChr           = nChr,
  segSites       = segSites,
  nSnpPerChr     = nSnpPerChr,
  Ne             = Ne,
  minSnpFreq     = minSnpFreq,
  prevCMC_mother = prevCMC_mother,
  prevPregCMC    = prevPregCMC,
  prevASC_child  = prevASC_child,
  ASC_male_ratio = ASC_male_ratio,
  betaPreg       = betaPreg,
  var_d          = var_d,
  var_m          = var_m,
  var_c          = var_c,
  var_e          = var_e,
  rho_CMC_d      = rho_CMC_d,
  rho_CMC_m      = rho_CMC_m,
  rho_d_m        = rho_d_m,
  nCausalPerChr  = nCausalPerChr,
  kappa          = kappa,
  delta_CMC      = delta_CMC,
  pheno_mom      = pheno_mom,
  pheno_kid      = pheno_kid,
  ASC_kid        = ASC_kid,
  is_male        = is_male,
  is_female      = is_female
)


t_genetic_components <- proc.time() - t_genetic_start
cat(sprintf("  tot_time_genetic_components : %.1f sec (%.2f min)\n",
            t_genetic_components["elapsed"], t_genetic_components["elapsed"] / 60))

## ============================================================
## Final summary
## ============================================================
cat("==========================================\n")
cat("SIMULATION SUMMARY\n")
cat("==========================================\n")
cat(sprintf("  Mothers - after removal    : %d\n",   length(pheno_mom$id)))
cat(sprintf("  Children - after removal   : %d\n",   length(pheno_kid$id)))
cat(sprintf("  CMC prevalence (mothers) - after removal  : %.4f\n", mean(pheno_mom$CMC)))
cat(sprintf("  Preg CMC prevalence (kids) - after removal : %.4f\n", mean(pheno_kid$preg_CMC)))
cat(sprintf("  ASC prevalence (overall) - after removal: %.4f\n", mean(pheno_kid$ASC)))
cat(sprintf("  ASC prevalence (overall)   : %.4f\n", mean(ASC_kid)))
cat(sprintf("  ASC prevalence — males     : %.4f\n", mean(ASC_kid[is_male])))
cat(sprintf("  ASC prevalence — females   : %.4f\n", mean(ASC_kid[is_female])))
cat("------------------------------------------\n")
cat("[TIMING SUMMARY]\n")
if (exists("t_founder_elapsed")) {
  cat(sprintf("  founderHaplotypesGenerator : %.1f sec (%.2f min)\n",
              t_founder_elapsed["elapsed"], t_founder_elapsed["elapsed"] / 60))
}
if (exists("t_pop_elapsed")) {
  cat(sprintf("  pop_generation             : %.1f sec (%.2f min)\n",
              t_pop_elapsed["elapsed"], t_pop_elapsed["elapsed"] / 60))
}
t_total_elapsed <- proc.time() - t_total_start
cat(sprintf("  TOTAL                      : %.1f sec (%.2f min)\n",
            t_total_elapsed["elapsed"], t_total_elapsed["elapsed"] / 60))
cat(sprintf("  tot_time_genetic_components: %.1f sec (%.2f min)\n",
            t_genetic_components["elapsed"], t_genetic_components["elapsed"] / 60))
cat("==========================================\n")










