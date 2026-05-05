## ------------------------------------------------------------
## Pedigree + genotypes + maternal CMC + pregnancy CMC + child ASD
## with:
##  - correct mean kids/mom (zero-truncated Poisson)
##  - union-based partner structure (divorce/remarriage; full- & half-sibs)
##  - robust PD guard for correlated SNP-effect matrix (Sigma)
##  - pregnancy CMC calibrated per mother + depends on binary CMC_mom
##  - causal variants (seg sites) separate from observed SNP chip (exported)
##  - ASD liability built from standardized predictors on child scale + residual
## ------------------------------------------------------------


## --------
##  THIS THIRD VERSION, AS THE SECOND ONE, ALLOWS FOR DIFFERENT VALUES AT DIFFERENT GENERATIONS FOR lambdaKids AND p_new_partner.
##  A NUMERIC VECTOR IS GOING TO BE ENCODED UNDER THESE VARIABLES.
##  IT CHANGES THE VALUES OF CERTAIN PARAMETERS WITH VALUES ESTIMATED FROM STATISTICS DENMARK.
##  MORE REPORT/CHECKS STEPS HAVE BEEN ADDED.
##  IT ALLOWS FOR DIFFERENT VALUES AT DIFFERENT GENERATIONS FOR mothers_fraction AND mean_unions_dad.
## --------


## ============================================================
## SANITY CHECK — packages + source files
## ============================================================
cat("==========================================\n")
cat("SANITY CHECK\n")
cat("==========================================\n")

required_packages <- c("AlphaSimR", "MASS", "Matrix")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("[FAIL] Required package not found: ", pkg,
                "\n       Install with: install.packages('", pkg, "')"))
  }
  cat(sprintf("  [OK] Package %-12s v%s\n", paste0(pkg, ":"), as.character(packageVersion(pkg))))
}

## NOTE: simulation_3 uses src_3, instead of src_2/1
required_src <- c(
  "src_3/founderHaplotypesGenerator.R",
  "src_3/helpers_juan.R",
  "src_3/pop_generation.R"
)
for (f in required_src) {
  if (!file.exists(f)) {
    stop(paste0("[FAIL] Source file not found: ", f,
                "\n       Make sure the 'src_3/' folder is in your working directory: ", getwd()))
  }
  cat(sprintf("  [OK] Source file found: %s\n", f))
}

cat("==========================================\n")
cat("SANITY CHECK PASSED\n")
cat("==========================================\n\n")


## ============================================================
## Helper: memory reporter
## ============================================================
mem_report <- function(obj, label) {
  sz <- object.size(obj)
  cat(sprintf("  [MEM] %-35s %s\n", paste0(label, ":"), format(sz, units = "auto", standard = "SI")))
}


## ============================================================
## Load libraries + source files
## ============================================================
library(AlphaSimR)
library(MASS)
library(Matrix)

source("src_3/founderHaplotypesGenerator.R")   # wrapper for founder haplotype generation
source("src_3/helpers_juan.R")                  # make_pd_corr, calib_alpha, alloc_kids_blocks
source("src_3/pop_generation.R")                # pop pedigree generator


## ============================================================
## Output directory
## ============================================================
OUT_DIR <- "results"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("[INFO] Output directory: %s/\n\n", normalizePath(OUT_DIR)))


## ============================================================
## Start total timer
## ============================================================
t_total_start <- proc.time()

set.seed(1)
options(stringsAsFactors = FALSE)


## -------------------------
## 0) Parameters
## -------------------------

## Genome simulation
nFounder              <- 50   # initial number of founders for runMacs2
nChr                  <- 2
segSites              <- 20        # segregating sites per chromosome
nSnpPerChr            <- 5        # SNPs per chromosome on chip (observed)
nGenerations_founders <- 5          # generations to build the final founders
Ne                    <- 20   # effective pop size
nFounders_final       <- 100   # final number of founders used as SimParam
minSnpFreq            <- 0.001 # minor allele frequency to be considered when adding snps on chip (for both founders and pop)

## Population / study design
## NOTE: simulation_3 uses per-generation vectors for lambdaKids, p_new_partner AND mothers_fraction
mothers_fraction     <- c(0.75, 0.75, 0.75, 0.75)          # fraction of female pop that become mothers. Each value is a different generation
lambdaKids           <- c(3.01, 2.54, 2.17, 1.7)   # mean kids per mother per generation (zero-truncated Poisson). Estimated from the real population
nGenerations_pop     <- 4            # number of generations with data
overlapping_fraction <- 0.00          # NOTE: sim_2 uses 0% (sim_1 used 5%)

## Prevalences
prevCMC_mother <- 0.10
prevPregCMC    <- 0.08
prevASD_child  <- 0.02

## Pregnancy effect on ASD liability
betaPreg <- 0.35

## ASD liability variance targets (excluding pregnancy term)
var_d <- 0.35   # child direct genetic
var_m <- 0.10   # maternal genetic on child's ASD
var_c <- 0.10   # shared maternal env (per mother)
var_e <- 1 - (var_d + var_m + var_c)
stopifnot(var_e > 0)

## Correlations among causal SNP effects across components
rho_CMC_d <- 0.30
rho_CMC_m <- 0.20
rho_d_m   <- 0.10

## Causal architecture
nCausalPerChr <- 5   # causal seg sites per chromosome

## Union / partnership structure
maxPartners_mom  <- 3
## NOTE: simulation_3 uses per-generation vector for p_new_partner and mean_unions_dad
# Last value for 4th generation needs to be fix
p_new_partner    <- c(0.02, 0.11, 0.10, 0.15)   # Probability of new partner per generation for a mother. Estimated from the real population.
mean_unions_dad  <- c(1.24, 1.17, 1.12, 1.15)  #  # higher => more paternal half-sibs + higher dad reuse. Estimated from real population

## Pregnancy CMC model
kappa     <- 1.0
delta_CMC <- 1.0

## Genotyped fraction for export
propGenotyped <- 0.30

## Parameter validation (sim_2 specific)
if (length(lambdaKids) != nGenerations_pop || length(p_new_partner) != nGenerations_pop) {
  stop("there are more generations than values for lambdaKids or p_new_partner")
}



cat("Reporting the parameters: \n")
cat(sprintf("nSnpPerChr: %g\n", nSnpPerChr))
cat(sprintf("Ne: %g\n", Ne))
cat(sprintf("minSnpFreq: %g\n", minSnpFreq))
cat("run \n")




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
## testing different values
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

## Unpack pop_all
pop              <- pop_all$pop
unions_all       <- pop_all$unions_all
children_ids_all <- unique(unlist(pop_all$chidlren_ids_list))

## Subsets
mothers_ids_all  <- unique(unions_all$mother_id)
mothers          <- pop[mothers_ids_all]
children_ids_all <- unique(unlist(pop_all$chidlren_ids_list))
children         <- pop[children_ids_all]

## Memory: pop objects
cat("[MEMORY] Objects after Step 2:\n")
mem_report(pop,      "pop (Pop — full population)")
mem_report(mothers,  "mothers (Pop subset)")
mem_report(children, "children (Pop subset)")

cat(sprintf("nInd (pop): %d\n", nInd(pop)))
cat(sprintf("nInd (mothers): %d\n", nInd(mothers)))
cat(sprintf("nInd (children): %d\n", nInd(children)))

## -------------------------
## 3) Pull genotypes
## -------------------------
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

child_moms      <- rep(unions_all$mother_id, times = unions_all$n_kids)
child_dads      <- rep(unions_all$father_id, times = unions_all$n_kids)
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
## 7) Child ASD liability
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

sd_e      <- sqrt(max(1 - var(nonpreg_part), 1e-8))
L_ASD_kid <- betaPreg * preg_CMC +
             nonpreg_part +
             sd_e * rnorm(length(nonpreg_part))

thr_ASD <- as.numeric(quantile(L_ASD_kid, probs = 1 - prevASD_child))
ASD_kid <- as.integer(L_ASD_kid > thr_ASD)


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
  ASD_liab  = L_ASD_kid,
  ASD       = ASD_kid
)

write.csv(ped,       file = file.path(OUT_DIR, "sim_pedigree_quick.csv"),         row.names = FALSE)
write.csv(pheno_mom, file = file.path(OUT_DIR, "sim_mom_phenotypes_quick.csv"),   row.names = FALSE)
write.csv(pheno_kid, file = file.path(OUT_DIR, "sim_child_phenotypes_quick.csv"), row.names = FALSE)

cat("  [OK] sim_pedigree.csv\n")
cat("  [OK] sim_mom_phenotypes.csv\n")
cat("  [OK] sim_child_phenotypes.csv\n")

## PLINK export for a genotyped subset
geno_idx <- sample(seq_len(nInd(pop)), size = floor(propGenotyped * nInd(pop)))
pop_geno  <- pop[geno_idx]

writePlink(
  pop_geno,
  baseName = file.path(OUT_DIR, "sim_genotyped_subset"),
  simParam  = SP,
  use       = "rand"
)
cat("  [OK] sim_genotyped_subset.{bed/bim/fam}\n\n")


## ============================================================
## Final summary
## ============================================================
cat("==========================================\n")
cat("SIMULATION SUMMARY\n")
cat("==========================================\n")
cat(sprintf("  Mothers                    : %d\n",   nInd(mothers)))
cat(sprintf("  Children                   : %d\n",   nInd(children)))
cat(sprintf("  CMC prevalence (mothers)   : %.4f\n", mean(CMC_mom)))
cat(sprintf("  Preg CMC prevalence (kids) : %.4f\n", mean(preg_CMC)))
cat(sprintf("  ASD prevalence (kids)      : %.4f\n", mean(ASD_kid)))
cat("------------------------------------------\n")
cat("[TIMING SUMMARY]\n")
cat(sprintf("  founderHaplotypesGenerator : %.1f sec (%.2f min)\n",
            t_founder_elapsed["elapsed"], t_founder_elapsed["elapsed"] / 60))
cat(sprintf("  pop_generation             : %.1f sec (%.2f min)\n",
            t_pop_elapsed["elapsed"],   t_pop_elapsed["elapsed"] / 60))

t_total_elapsed <- proc.time() - t_total_start
cat(sprintf("  TOTAL                      : %.1f sec (%.2f min)\n",
            t_total_elapsed["elapsed"], t_total_elapsed["elapsed"] / 60))
cat("==========================================\n")
