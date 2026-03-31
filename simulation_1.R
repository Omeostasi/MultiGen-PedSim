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
source("src/pop_generation.R")

set.seed(1)
options(stringsAsFactors = FALSE)

## -------------------------
## 0) Parameters
## -------------------------

## Genome simulation
nFounder   <- 20*10^3  # initial number of founders for runMacs2 
nChr       <- 10
segSites   <- 1000    # segregating sites per chromosome
nSnpPerChr <- 100     # SNPs per chromosome on chip (observed)
nGenerations_founders <- 10 # the number of generations to generate the final founders
Ne <- 50*10^3
nFounders_final <- 1*10^6 # final number of founders used as SimParam


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
nCausalPerChr <- 50      # causal segregating sites per chromosome

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


prs_m_mom_byID  <- setNames(prs_m_mom, as.character(mothers_ids_all))

## Expand union rows into per-child parent vectors
child_moms <- rep(unions_all$mother_id, times = unions_all$n_kids)
child_dads <- rep(unions_all$father_id, times = unions_all$n_kids)

prs_m_for_child <- prs_m_mom_byID[as.character(child_moms)]

## -------------------------
## 5) Maternal CMC (mother-level trait; liability threshold)
## -------------------------

L_CMC_mom <- sqrt(0.40) * prs_cmc_mom + sqrt(0.60) * rnorm(length(mothers_ids_all))
thr_CMC <- as.numeric(quantile(L_CMC_mom, probs = 1 - prevCMC_mother))
CMC_mom <- as.integer(L_CMC_mom > thr_CMC)

## -------------------------
## 6) Pregnancy CMC (per child)
##    Fixes:
##      - calibrated per mother (each mother counts equally)
##      - depends on continuous liability AND binary CMC status
## -------------------------

linpred_mom <- kappa * L_CMC_mom + delta_CMC * CMC_mom
alpha <- calib_alpha(prevPregCMC, linpred_mom)

p_mom <- plogis(alpha + linpred_mom)
p_for_child <- p_mom[match(as.character(child_moms), as.character(mothers_ids_all))]
preg_CMC <- rbinom(n = length(p_for_child), size = 1, prob = p_for_child)


## -------------------------
## 7) Child ASD (liability threshold)
##    Fix: stabilize liability scaling on child sample under correlated components.
## -------------------------




## Shared maternal environment (same for siblings of a mother)
c_mom <- rnorm(length(mothers_ids_all), mean = 0, sd = 1)
c_mom_byID <- setNames(c_mom, as.character(mothers_ids_all))
c_for_child <- c_mom_byID[as.character(child_moms)]

## Standardize predictors on the CHILD scale
gd_z <- as.numeric(scale(prs_d_kid))
gm_z <- as.numeric(scale(prs_m_for_child))
c_z  <- as.numeric(scale(c_for_child))

nonpreg_part <- sqrt(var_d) * gd_z +
  sqrt(var_m) * gm_z +
  sqrt(var_c) * c_z

## Choose residual SD so Var(nonpreg_part + e) ~= 1
sd_e <- sqrt(max(1 - var(nonpreg_part), 1e-8))

L_ASD_kid <- betaPreg * preg_CMC +
  nonpreg_part +
  sd_e * rnorm(length(nonpreg_part))

thr_ASD <- as.numeric(quantile(L_ASD_kid, probs = 1 - prevASD_child))
ASD_kid <- as.integer(L_ASD_kid > thr_ASD)




## -------------------------
## 8) Export: pedigree + phenotypes + PLINK genotypes for a subset
## -------------------------



ped <- getPed(pop)  # columns: id, mother, father

pheno_mom <- data.frame(
  id = mothers@id,
  role = "mother",
  CMC_liab = L_CMC_mom,
  CMC = CMC_mom
)




pheno_kid <- data.frame(
  id = children@id,
  role = "child",
  mother_id = child_moms,
  father_id = child_dads,
  preg_CMC = preg_CMC,
  ASD_liab = L_ASD_kid,
  ASD = ASD_kid
)

write.csv(ped,       file = "sim_pedigree.csv",         row.names = FALSE)
write.csv(pheno_mom, file = "sim_mom_phenotypes.csv",   row.names = FALSE)
write.csv(pheno_kid, file = "sim_child_phenotypes.csv", row.names = FALSE)

## Export PLINK for a genotyped subset (chip SNPs)
geno_idx <- sample(seq_len(nInd(pop)), size = floor(propGenotyped * nInd(pop)))
pop_geno <- pop[geno_idx]


writePlink(
  pop_geno,
  baseName = "sim_genotyped_subset",
  simParam = SP,
  use = "rand"
  
)

## Optional: quick sanity checks
cat("Mothers:", nInd(mothers), "\n")
cat("Children:", nInd(children), "\n")
# cat("Mean kids/mom:", mean(kids_per_mom), "\n") # We don't store this in this version
cat("CMC moms prevalence:", mean(CMC_mom), "\n")
cat("Preg CMC prevalence (per child):", mean(preg_CMC), "\n")
cat("ASD prevalence:", mean(ASD_kid), "\n")


