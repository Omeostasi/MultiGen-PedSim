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

set.seed(1)
options(stringsAsFactors = FALSE)

## -------------------------
## 0) Parameters
## -------------------------

## Genome simulation
nFounder   <- 100
nChr       <- 10
segSites   <- 20       # segregating sites per chromosome
nSnpPerChr <- 5     # SNPs per chromosome on chip (observed)

## Population / study design
nParents   <- 50      # size of adult generation
nMothers   <- 20
lambdaKids <- 1.6       # mean kids per mother among mothers with >=1 child

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
## Helper functions
## -------------------------

## Zero-truncated Poisson with target mean: E[X | X>=1] = mean_target
rztpois_mean <- function(n, mean_target) {
  stopifnot(mean_target > 1)  # zero-truncated Poisson mean is always > 1
  f <- function(mu) mu / (1 - exp(-mu)) - mean_target
  mu0 <- uniroot(f, interval = c(1e-8, 1e4))$root
  
  x <- rpois(n, mu0)
  while (any(x == 0)) {
    idx <- which(x == 0)
    x[idx] <- rpois(length(idx), mu0)
  }
  x
}

## Ensure correlation matrix is positive definite (for mvrnorm)
make_pd_corr <- function(Sigma, eps = 1e-8) {
  ev <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) > eps) return(Sigma)
  as.matrix(nearPD(Sigma, corr = TRUE)$mat)
}

## Calibrate logistic intercept alpha so mean(plogis(alpha + lp)) == target
calib_alpha <- function(target, lp, interval = c(-12, 12)) {
  uniroot(function(a) mean(plogis(a + lp)) - target, interval = interval)$root
}

## Split total kids into k blocks (each >=1) for union-based full-sib clusters
alloc_kids_blocks <- function(total_kids, k_blocks) {
  if (k_blocks == 1) return(total_kids)
  x <- rep(1L, k_blocks)
  rem <- total_kids - k_blocks
  if (rem > 0) {
    x <- x + as.integer(rmultinom(1, size = rem, prob = rep(1, k_blocks))[, 1])
  }
  x
}


## -------------------------
## 1) Simulate pedigree-consistent genotypes (AlphaSimR)
## -------------------------

founderHap <- quickHaplo(nInd = nFounder, nChr = nChr, segSites = segSites)

SP <- SimParam$new(founderHap)
SP$setTrackPed(TRUE)
SP$setSexes(sexes = "yes_sys") 
SP$addSnpChip(nSnpPerChr = nSnpPerChr, minSnpFreq = 0.05)

pop0 <- newPop(founderHap)

pop1 <- randCross(
  pop0,
  nCrosses = ceiling(nParents / 2),
  nProgeny = 2,
  ignoreSexes = FALSE,
  simParam = SP
)

fem1 <- pop1[isFemale(pop1)]
mal1 <- pop1[isMale(pop1)]

if (nInd(fem1) < nMothers) stop("Not enough females to pick mothers. Increase nParents.")
mothers <- fem1[sample(seq_len(nInd(fem1)), nMothers)]
mother_ids <- mothers@id
dad_ids <- mal1@id


## -------------------------
## 2) Build family structure with unions (divorce/remarriage)
##    Produces full-sib blocks + maternal and paternal half-sibs.
## -------------------------

## Kids per mother: mean(lambdaKids) among mothers with >=1 child
kids_per_mom <- rztpois_mean(length(mother_ids), mean_target = lambdaKids)

## Number of partners per mother (1..maxPartners_mom), capped by number of kids
nPartners_mom <- 1 + rbinom(length(mother_ids), size = (maxPartners_mom - 1), prob = p_new_partner)
nPartners_mom <- pmin(nPartners_mom, kids_per_mom)
nPartners_mom <- pmax(nPartners_mom, 1)

## Union table: one row per mother-partner union, with number of kids in that union
union_list <- vector("list", length(mother_ids))
for (i in seq_along(mother_ids)) {
  blocks <- alloc_kids_blocks(kids_per_mom[i], nPartners_mom[i])
  union_list[[i]] <- data.frame(
    mother_id  = mother_ids[i],
    union_order = seq_along(blocks),
    n_kids     = blocks
  )
}
unions <- do.call(rbind, union_list)

## Allocate "union capacities" to fathers to avoid extreme super-fathers
dad_capacity <- pmax(1L, rpois(length(dad_ids), lambda = mean_unions_dad))
names(dad_capacity) <- as.character(dad_ids)

## Assign a father to each union (sampling from remaining capacity)
assigned_dads <- character(nrow(unions))
for (u in seq_len(nrow(unions))) {
  eligible <- names(dad_capacity)[dad_capacity > 0]
  if (length(eligible) == 0) {
    ## If depleted, refresh capacities (or increase mean_unions_dad)
    dad_capacity[] <- pmax(1L, rpois(length(dad_ids), lambda = mean_unions_dad))
    eligible <- names(dad_capacity)[dad_capacity > 0]
  }
  probs <- dad_capacity[eligible]
  dad <- sample(eligible, size = 1, prob = probs)
  assigned_dads[u] <- dad
  dad_capacity[dad] <- dad_capacity[dad] - 1L
}
unions$father_id <- assigned_dads

## Expand union rows into per-child parent vectors
child_moms <- rep(unions$mother_id, times = unions$n_kids)
child_dads <- rep(unions$father_id, times = unions$n_kids)

## Cross plan
crossPlan <- cbind(as.character(child_moms), as.character(child_dads))
children <- makeCross(pop1, crossPlan = crossPlan, nProgeny = 1, simParam = SP)

stopifnot(nInd(children) == length(child_moms))

pop_all <- c(pop0, pop1, children)


## we need to make more gen here

## -------------------------
## 3) Pull genotypes:
##    - SNP chip: observed/exported
##    - segregating sites: used for causal effects
## -------------------------

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
prs_m_mom_byID  <- setNames(prs_m_mom, as.character(mother_ids))
prs_m_for_child <- prs_m_mom_byID[as.character(child_moms)]


## -------------------------
## 5) Maternal CMC (mother-level trait; liability threshold)
## -------------------------

L_CMC_mom <- sqrt(0.40) * prs_cmc_mom + sqrt(0.60) * rnorm(length(mother_ids))
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
p_for_child <- p_mom[match(as.character(child_moms), as.character(mother_ids))]
preg_CMC <- rbinom(n = length(p_for_child), size = 1, prob = p_for_child)


## -------------------------
## 7) Child ASD (liability threshold)
##    Fix: stabilize liability scaling on child sample under correlated components.
## -------------------------

## Shared maternal environment (same for siblings of a mother)
c_mom <- rnorm(length(mother_ids), mean = 0, sd = 1)
c_mom_byID <- setNames(c_mom, as.character(mother_ids))
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

ped <- getPed(pop_all)  # columns: id, mother, father

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
geno_idx <- sample(seq_len(nInd(pop_all)), size = floor(propGenotyped * nInd(pop_all)))
pop_geno <- pop_all[geno_idx]


writePlink(
  pop_geno,
  baseName = "sim_genotyped_subset",
  simParam = SP,
  use = "rand"
  
)

## Optional: quick sanity checks
cat("Mothers:", nMothers, "\n")
cat("Children:", nInd(children), "\n")
cat("Mean kids/mom:", mean(kids_per_mom), "\n")
cat("CMC moms prevalence:", mean(CMC_mom), "\n")
cat("Preg CMC prevalence (per child):", mean(preg_CMC), "\n")
cat("ASD prevalence:", mean(ASD_kid), "\n")

