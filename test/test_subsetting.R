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
##  THIS FOURTH VERSION, AS THE SECOND ONE, ALLOWS FOR DIFFERENT VALUES AT DIFFERENT GENERATIONS FOR lambdaKids AND p_new_partner.
##  A NUMERIC VECTOR IS GOING TO BE ENCODED UNDER THESE VARIABLES.
##  AS THE THIRD VERSION: IT CHANGES THE VALUES OF CERTAIN PARAMETERS WITH VALUES ESTIMATED FROM STATISTICS DENMARK,
##  MORE REPORT/CHECKS STEPS HAVE BEEN ADDED, AND
##  IT ALLOWS FOR DIFFERENT VALUES AT DIFFERENT GENERATIONS FOR mothers_fraction AND mean_unions_dad.
##  UNIQUE TO THE FOURTH VERSION: IT HAS BEEN BUILT THINKING ABOUT 5 GENERATIONS (SILENT, BOOMER, GENX, MILLENNIAL, GEN Z), IT GENOTYPES 
##  A SAMPLE OF THE POPULATION TO RESEMBLE THE iPSYCH COHORT, IT DIVIDES IN CHUNKS THE PLINK FILE (BINDING IT BACK TO ONE AT THE END) TO OVERCOME
##  THE SIZE LIMIT OF THE MATRIX ALLOWED IN THE writePLINK() FUNCTION, AND, FINALLY, IT REMOVES 80% OF THE SILENT GENERATION AND THE FOUNDERS
##  FROM THE FINAL DATA.
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
  "src_4/founderHaplotypesGenerator.R",
  "src_4/helpers_juan.R",
  "src_4/pop_generation.R",
  "src_4/subsetting.R"
)
for (f in required_src) {
  if (!file.exists(f)) {
    stop(paste0("[FAIL] Source file not found: ", f,
                "\n       Make sure the 'src_4/' folder is in your working directory: ", getwd()))
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

source("src_4/founderHaplotypesGenerator.R")   # wrapper for founder haplotype generation
source("src_4/helpers_juan.R")                 # make_pd_corr, calib_alpha, alloc_kids_blocks
source("src_4/pop_generation.R")               # pop pedigree generator
source("src_4/subsetting.R")     # creates subsets of final population and genotyped population


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
nFounder              <- 20   # initial number of founders for runMacs2
nChr                  <- 2
segSites              <- 10        # segregating sites per chromosome
nSnpPerChr            <- 5        # SNPs per chromosome on chip (observed)
nGenerations_founders <- 3          # generations to build the final founders
Ne                    <- 100   # effective pop size
nFounders_final       <- 100   # final number of founders used as SimParam
minSnpFreq            <- 0.005 # minor allele frequency to be considered when adding snps on chip (for both founders and pop)

## Population / study design
## NOTE: simulation_3 uses per-generation vectors for lambdaKids, p_new_partner AND mothers_fraction
mothers_fraction     <- c(0.99, 0.99, 0.98, 0.82, 0.71)          # fraction of female pop that become mothers. Each value is a different generation
lambdaKids           <- c(2.01, 2.04, 2.00, 1.64, 1.34)   # mean kids per mother per generation (zero-truncated Poisson). Estimated from the real population
nGenerations_pop     <- 5            # number of generations with data
overlapping_fraction <- 0.00          # If more than 0 allows for a fraction of previous generation to cross during current generation
rm_older_generations <- TRUE     # Default it's TRUE,  it removes 80% first generation and founders after generating them and computing their genetic correlation


## Prevalences
prevCMC_mother <- 0.30
prevPregCMC    <- 0.40
prevASD_child  <- 0.20

# prevalences ASD by gender
prevASD_ratio <- 3/2 # Defaul: 3 males every 2 females 
prevASD_male   <- prevASD_child*prevASD_ratio
prevASD_female <- round(prevASD_child*1/prevASD_ratio, 3)

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
nCausalPerChr <- 3   # causal seg sites per chromosome

## Union / partnership structure
maxPartners_mom  <- 3
## NOTE: simulation_3 uses per-generation vector for p_new_partner and mean_unions_dad
p_new_partner    <- c(0.02, 0.12, 0.13, 0.12, 0.10)   # Probability of new partner per generation for a mother. Estimated from the real population.
mean_unions_dad  <- c(1.25, 1.17, 1.15, 1.13, 1.09)   # higher => more paternal half-sibs + higher dad reuse. Estimated from real population

## Pregnancy CMC model
kappa     <- 1.0
delta_CMC <- 1.0

## Genotyped subset targets (iPSYCH-like design)
n_geno_ASD    <- 2    # target number of ASD cases to genotype
n_geno_random <- 15   # target number of random non-ASD individuals to genotype

## Parameter validation (sim_2 specific)
if (
  length(lambdaKids) != nGenerations_pop ||
  length(p_new_partner) != nGenerations_pop ||
  length(mothers_fraction) != nGenerations_pop ||
  length(mean_unions_dad) != nGenerations_pop
) {
  stop("Mismatch between number of generations and parameter vector lengths")
}




cat("Reporting the parameters for comparison: \n")
cat(sprintf("segSites: %g\n", segSites))
cat(sprintf("nSnpPerChr: %g\n", nSnpPerChr))
cat(sprintf("Ne: %g\n", Ne))
cat(sprintf("minSnpFreq: %g\n", minSnpFreq))
cat(sprintf("nGenerations_pop: %g\n", nGenerations_pop))
cat(sprintf("nFounders_final: %g\n", nFounders_final))
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
children_ids_all <- unique(na.omit(unlist(pop_all$children_ids_list)))

## Subsets
mothers_ids_all  <- unique(na.omit(unions_all$mother_id))
mothers          <- pop[mothers_ids_all]
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

sd_e      <- sqrt(max(1 - var(nonpreg_part, na.rm = TRUE), 1e-8))
L_ASD_kid <- betaPreg * preg_CMC +
  nonpreg_part +
  sd_e * rnorm(length(nonpreg_part))

# Sex-specific thresholds
child_sex <- children@sex

is_male   <- child_sex == "M"
is_female <- child_sex == "F"

# Each sex gets its own threshold calibrated to its target prevalence
thr_ASD_male   <- quantile(L_ASD_kid[is_male],   probs = 1 - prevASD_male,   na.rm = TRUE)
thr_ASD_female <- quantile(L_ASD_kid[is_female],  probs = 1 - prevASD_female, na.rm = TRUE)

ASD_kid <- integer(length(L_ASD_kid))
ASD_kid[is_male]   <- as.integer(L_ASD_kid[is_male]   > thr_ASD_male)
ASD_kid[is_female] <- as.integer(L_ASD_kid[is_female] > thr_ASD_female)






## -------------------------
## 8) Export outputs to results/
## -------------------------
cat("------------------------------------------\n")
cat("[STEP 8] Writing output files ...\n")
cat("------------------------------------------\n")







ped_old <- getPed(pop)     


pheno_mom <- data.frame(
  id       = mothers@id,
  role     = "mother",
  CMC_liab = L_CMC_mom,
  CMC      = CMC_mom
)

pheno_kid <- data.frame(
  id        = children@id,
  role      = "child",
  sex       = child_sex,
  mother_id = child_moms,
  father_id = child_dads,
  preg_CMC  = preg_CMC,
  ASD_liab  = L_ASD_kid,
  ASD       = ASD_kid
)








## ---------------
## Removing older generations from the set
## -----------------

subset <- remove_older_generations(pheno_mom = pheno_mom,
                                   pheno_kid = pheno_kid,
                                   ped = ped_old,
                                   older_generations_ids = pop_all$children_ids_list, 
                                   rm_older_generations = rm_older_generations)

ped <- subset$ped
pheno_mom <- subset$pheno_mom
pheno_kid <- subset$pheno_kid




##
## CREATE GENOTYPED SUBSET
##



pop_geno <- genotype_subset(pop = pop, pheno_kid = pheno_kid, generations_ids = pop_all$children_ids_list,
                            ped = ped,
                            nGenerations_pop = nGenerations_pop,
                            n_geno_random = n_geno_random, n_geno_ASD = n_geno_ASD)

