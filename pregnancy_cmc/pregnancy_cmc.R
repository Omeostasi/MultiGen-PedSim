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

library(AlphaSimR)
library(MASS)
library(Matrix)

# Adding wrapper function to generate founders
source("src/founderHaplotypesGenerator.R")

# helpers: make_pd_corr, calib_alpha, alloc_kids_blocks
source("src/helpers_juan.R")

# Function to generate pop pedigree
source("pull_genotypes/pop_generation.R")

set.seed(1)
options(stringsAsFactors = FALSE)

## -------------------------
## 0) Parameters
## -------------------------

## Genome simulation
nFounder   <- 50  # initial number of founders for runMacs2 
nChr       <- 2
segSites   <- 20    # segregating sites per chromosome
nSnpPerChr <- 5     # SNPs per chromosome on chip (observed)
nGenerations_founders <- 5 # the number of generations to generate the final founders
Ne <- 20
nFounders_final <- 100 # final number of founders used as SimParam


## Population / study design
# nParents   <- 2000      # size of adult generation # this is not used
# nMothers   <- 700 # this is not used
mothers_fraction <- 0.75 # fraction of female pop that become mothers
lambdaKids <- 1.6       # mean kids per mother among mothers with >=1 child
nGenerations_pop <- 3 # number generation that we have data for
overlapping_fraction <- 0.05 # 5 % of the older pop could still have children during the current pop

## Prevalences
prevCMC_mother <- 0.10
prevPregCMC    <- 0.08  # calibrated per mother (each mother counts equally)
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
nCausalPerChr <- 10      # causal segregating sites per chromosome

## Union / partnership structure (divorce/remarriage)
maxPartners_mom <- 3     # max distinct fathers per mother
p_new_partner   <- 0.25  # higher => more mothers with multiple partners
mean_unions_dad <- 1.3   # higher => more paternal half-sibs + higher dad reuse

## Pregnancy CMC model
kappa     <- 1.0         # dependence on continuous maternal CMC liability
delta_CMC <- 1.0         # additional shift if mother is CMC case (binary)

## Genotyped fraction for export
propGenotyped <- 0.30





## -------------------------
## 1) Simulate pedigree-consistent genotypes (AlphaSimR)
## -------------------------

founderHap <- founderHaplotypesGenerator(nInd = nFounder, nChr = nChr, segSites = segSites,
                                         nThreads = nChr, Ne = Ne,
                                         nGenerations = nGenerations_founders, method = "run", nFounders_final = nFounders_final)

SP <- SimParam$new(founderHap)
SP$setTrackPed(TRUE)
SP$setSexes(sexes = "yes_rand") ### CHANGED FROM YES_SYS TO YES_RAND
SP$addSnpChip(nSnpPerChr = nSnpPerChr, minSnpFreq = 0.05)



## -------------------------
## 2) Build family structure with unions (divorce/remarriage)
##    Produces full-sib blocks + maternal and paternal half-sibs.
##    Create full pedigree for nGenerations_pop
## -------------------------

# This is object contains our pop, ids and unions
pop_all <- pop_generation(founderHap, 
                          overlapping_fraction = overlapping_fraction, 
                          nGenerations_pop = nGenerations_pop, mothers_fraction = mothers_fraction,
                          lambdaKids = lambdaKids, maxPartners_mom = maxPartners_mom,
                          p_new_partner = p_new_partner, mean_unions_dad = mean_unions_dad)

# Extracting from pop_all all the different elements
pop <- pop_all$pop # our pop
unions_all <- pop_all$unions_all # list of unions
children_ids_all <- unique(unlist(pop_all$chidlren_ids_list)) # children ids


# Extracting subset mothers and children pop
mothers_ids_all <-unique(unions_all$mother_id)
mothers <- pop[mothers_ids_all] # mothers pop subset

children_ids_all <- unique(unlist(pop_all$chidlren_ids_list))
children <- pop[children_ids_all] # children pop subset (essentially missing just the founders)

## -------------------------
## 3) Pull genotypes:
##    - SNP chip: observed/exported
##    - segregating sites: used for causal effects
## -------------------------

# Report step
cat("Pulling segSites, Genotypes and choosing causal variants from segSites \n")

## Observed chip (for export / downstream "array data")
Gchip_mom <- pullSnpGeno(mothers,  snpChip = 1, simParam = SP)
Gchip_kid <- pullSnpGeno(children, snpChip = 1, simParam = SP)

## Segregating sites (for causal architecture)
Gseg_mom <- pullSegSiteGeno(mothers,  simParam = SP)
Gseg_kid <- pullSegSiteGeno(children, simParam = SP)

## Choose causal variants from segregating sites
nCausal <- nChr * nCausalPerChr
if (nCausal > ncol(Gseg_mom)) stop("nCausal exceeds available seg sites. Reduce nCausalPerChr or increase segSites.")
causal_idx <- sample(seq_len(ncol(Gseg_mom)), nCausal)

Gcausal_mom <- Gseg_mom[, causal_idx, drop = FALSE]
Gcausal_kid <- Gseg_kid[, causal_idx, drop = FALSE]
m_causal <- ncol(Gcausal_mom)


## -------------------------
## 4) Correlated causal SNP effects (shared genetics / pleiotropy)
## -------------------------

###
### MODIFYING STEP 4 TO MAKE IT WORK WITH MODIFICATIONS
###
Sigma <- matrix(
  c(1,         rho_CMC_d, rho_CMC_m,
    rho_CMC_d, 1,         rho_d_m,
    rho_CMC_m, rho_d_m,   1),
  nrow = 3, byrow = TRUE
)
Sigma <- make_pd_corr(Sigma)

eff <- MASS::mvrnorm(n = m_causal, mu = c(0, 0, 0), Sigma = Sigma)
b_cmc <- eff[, 1]
b_d   <- eff[, 2]
b_m   <- eff[, 3]

## PRS from causal variants
prs_cmc_mom <- as.numeric(scale(Gcausal_mom %*% b_cmc))
prs_m_mom   <- as.numeric(scale(Gcausal_mom %*% b_m))
prs_d_kid   <- as.numeric(scale(Gcausal_kid %*% b_d))


## Map maternal PRS to each child by mother id

###
### THIS DOESN'T WORK, REQUIRE SOME MODIFICATION
###
prs_m_mom_byID  <- setNames(prs_m_mom, as.character(mothers_ids_all)) # this works
# prs_m_for_child <- prs_m_mom_byID[as.character(child_moms)] # problem here, we need something similar to "child_moms"

## Expand union rows into per-child parent vectors
child_moms <- rep(unions_all$mother_id, times = unions_all$n_kids)
child_dads <- rep(unions_all$father_id, times = unions_all$n_kids)

prs_m_for_child <- prs_m_mom_byID[as.character(child_moms)]

### SOLVED


## -------------------------
## 5) Maternal CMC (mother-level trait; liability threshold)
## -------------------------


### Checking this step

L_CMC_mom <- sqrt(0.40) * prs_cmc_mom + sqrt(0.60) * rnorm(length(mothers_ids_all))
thr_CMC <- as.numeric(quantile(L_CMC_mom, probs = 1 - prevCMC_mother))
CMC_mom <- as.integer(L_CMC_mom > thr_CMC)

### Perfect


## -------------------------
## 6) Pregnancy CMC (per child)
##    Fixes:
##      - calibrated per mother (each mother counts equally)
##      - depends on continuous liability AND binary CMC status
## -------------------------



### Checking this step

linpred_mom <- kappa * L_CMC_mom + delta_CMC * CMC_mom
alpha <- calib_alpha(prevPregCMC, linpred_mom)

p_mom <- plogis(alpha + linpred_mom)
p_for_child <- p_mom[match(as.character(child_moms), as.character(mothers_ids_all))]
preg_CMC <- rbinom(n = length(p_for_child), size = 1, prob = p_for_child)