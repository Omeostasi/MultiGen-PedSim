library(AlphaSimR)
library(MASS)

## -------------------------
## 0) Size & “truth” settings
## -------------------------
nFounder   <- 800     # founders for the genome simulator
nChr       <- 22      # number of chromosomes for simulated genotypes
segSites   <- 400     # segregating sites per chromosome (bigger = more independent loci)
nSnpPerChr <- 200     # SNPs per chromosome on the "SNP array” 

nParents   <- 2000    # adult generation size (parents of study kids) # just a placeholder size, should be larger
nMothers   <- 700     # number of mothers with children in the study generation # just a placeholder size, should be larger
lambdaKids <- 1.6     # avg children per mother # need to get the 'real' number from population data

prevCMC_mother <- 0.10   # prevalence of “mother has CMC”
prevPregCMC    <- 0.08   # prevalence of “this pregnancy had CMC complication”
prevASD_child  <- 0.02   # prevalence of ASD in children

betaPreg <- 0.35         # pregnancy CMC effect on ASD *liability* scale

## Variance component targets for ASD liability (excluding the pregnancy indicator term)
var_d <- 0.35   # child direct genetic component
var_m <- 0.10   # maternal genetic component on child's ASD
var_c <- 0.10   # shared maternal environment (same for siblings)
var_e <- 1 - (var_d + var_m + var_c)
stopifnot(var_e > 0)

## Genetic correlation settings (correlations of SNP effects across components)
rho_CMC_d <- 0.30   # shared genetics between maternal CMC and child ASD direct genetics
rho_CMC_m <- 0.20   # shared genetics between maternal CMC and maternal genetic ASD component
rho_d_m   <- 0.10   # shared genetics between ASD direct and ASD maternal component

## -------------------------
## 1) Simulate pedigree-consistent genotypes
## -------------------------
## Founders: random haplotypes across the genome
founderHap <- quickHaplo(nInd = nFounder, nChr = nChr, segSites = segSites)

## Set simulation parameters, track pedigree, assign sexes, create a SNP array
SP <- SimParam$new(founderHap)
SP$setTrackPed(TRUE)
SP$setSexes(sexes = "yes_sys")
SP$addSnpChip(nSnpPerChr = nSnpPerChr, minSnpFreq = 0.05)

## Create founders as a population
pop0 <- newPop(founderHap)

## Make an adult generation by random mating among founders
pop1 <- randCross(pop0, nCrosses = ceiling(nParents/2), nProgeny = 2,
                  ignoreSexes = FALSE, simParam = SP)

## Pick mothers and fathers from pop1
fem1 <- pop1[isFemale(pop1)]
mal1 <- pop1[isMale(pop1)]

if (nInd(fem1) < nMothers) stop("Not enough females to pick mothers. Increase nParents.")
mothers <- fem1[sample(seq_len(nInd(fem1)), nMothers)]
mother_ids <- mothers@id

## Create a study generation with multiple children per mother
kids_per_mom <- rpois(length(mother_ids), lambda = lambdaKids) + 1
child_moms <- rep(mother_ids, times = kids_per_mom)
child_dads <- sample(mal1@id, size = length(child_moms), replace = T)

## makeCross() accepts parent IDs as *character strings* (avoid potential numeric ID errors)
crossPlan <- cbind(as.character(child_moms), as.character(child_dads))
children <- makeCross(pop1, crossPlan = crossPlan, nProgeny = 1, simParam = SP)

## Combine founders, parents and children as a single population
pop_all <- c(pop0, pop1, children)



## -------------------------
## 2) Pull SNP genotypes (0/1/2) for mothers and children
## -------------------------
G_mom  <- pullSnpGeno(mothers,  snpChip = 1, simParam = SP)
G_kid  <- pullSnpGeno(children, snpChip = 1, simParam = SP)
m      <- ncol(G_mom)

## -------------------------
## 3) Create correlated SNP effects (this is where “shared genetics” is set)
##    For each SNP j, draw effects for:
##      - maternal CMC genetic component (b_cmc[j])
##      - child ASD direct genetic component (b_d[j])
##      - maternal genetic ASD component (b_m[j])
## -------------------------
Sigma <- matrix(c(1,         rho_CMC_d, rho_CMC_m,
                  rho_CMC_d, 1,         rho_d_m,
                  rho_CMC_m, rho_d_m,   1),
                nrow = 3, byrow = T)

eff <- MASS::mvrnorm(n = m, mu = c(0,0,0), Sigma = Sigma)
b_cmc <- eff[,1]
b_d   <- eff[,2]
b_m   <- eff[,3]

## "Polygenic scores" (scaled) for each component
prs_cmc_mom <- as.numeric(scale(G_mom %*% b_cmc))
prs_m_mom   <- as.numeric(scale(G_mom %*% b_m))
prs_d_kid   <- as.numeric(scale(G_kid %*% b_d))

## Map maternal genetic ASD score to each child (by mother id)
prs_m_mom_byID <- setNames(prs_m_mom, mother_ids)
prs_m_for_child <- prs_m_mom_byID[as.character(child_moms)]

## -------------------------
## 4) Simulate maternal CMC (mother-level trait)
## -------------------------
## Liability = genetic part + noise, then threshold to desired prevalence
L_CMC_mom <- sqrt(0.40) * prs_cmc_mom + sqrt(0.60) * rnorm(length(mother_ids))

thr_CMC <- as.numeric(quantile(L_CMC_mom, probs = 1 - prevCMC_mother))
CMC_mom <- as.integer(L_CMC_mom > thr_CMC)


## -------------------------
## 5) Simulate pregnancy-specific CMC (per child)
## -------------------------
## kappa is a scaling factor in the mother's underlying CMC liability (L_CMC_mom)
## drives the probability that any given pregnancy is CMC-complicated.
## - kappa = 0  : pregnancy CMC is essentially unrelated to the mother’s liability
##               (siblings have little correlation in preg_CMC aside from chance).
## - kappa ~ 1  : moderate dependence; each unit in L_CMC_mom increases log-odds by 1
##               (odds multiply by exp(1) ≈ 2.7).
## - kappa > 1  : strong dependence; pregnancies become "sticky" within mothers:
##               high-liability mothers have CMC in many pregnancies, low-liability
##               mothers in few, so siblings correlate more in preg_CMC.
## Note: pregnancies can still differ within the same mother because preg_CMC is
## drawn as a Bernoulli(p) per pregnancy (so there is pregnancy-to-pregnancy noise).
## Make pregnancy risk depend on the mother’s CMC liability (so siblings correlate),
## but still allow within-mother variation (Bernoulli draws per pregnancy).
kappa <- 1.0
linpred <- kappa * L_CMC_mom[match(as.character(child_moms), as.character(mother_ids))]

## Calibrate intercept so mean prevalence matches prevPregCMC
calib_alpha <- function(target, lp) {
  uniroot(function(a) mean(plogis(a + lp)) - target, interval = c(-12, 12))$root
}
alpha <- calib_alpha(prevPregCMC, linpred)

p_preg <- plogis(alpha + linpred)
preg_CMC <- rbinom(n = length(p_preg), size = 1, prob = p_preg)


## -------------------------
## 6) Simulate child ASD (child-level trait)
## -------------------------
## Shared maternal environment (same value for all siblings of a mother)
c_mom <- rnorm(length(mother_ids), mean = 0, sd = 1)
c_mom_byID <- setNames(c_mom, mother_ids)
c_for_child <- c_mom_byID[as.character(child_moms)]

## ASD liability:
##   pregnancy effect + direct genetics + maternal genetics + shared maternal env + residual
L_ASD_kid <- betaPreg * preg_CMC +
  sqrt(var_d) * prs_d_kid +
  sqrt(var_m) * prs_m_for_child +
  sqrt(var_c) * c_for_child +
  sqrt(var_e) * rnorm(length(prs_d_kid))

thr_ASD <- as.numeric(quantile(L_ASD_kid, probs = 1 - prevASD_child))
ASD_kid <- as.integer(L_ASD_kid > thr_ASD)


## -------------------------
## 7) Export: pedigree + phenotype table + PLINK genotypes for a subset
## -------------------------
ped <- getPed(pop_all)  # columns: id, mother, father

## Phenotype tables
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


## Save pedigree and phenotypes
write.csv(ped,  file = "sim_pedigree.csv", row.names = F)
write.csv(pheno_mom, file = "sim_mom_phenotypes.csv", row.names = F)
write.csv(pheno_kid, file = "sim_child_phenotypes.csv", row.names = F)

## PLINK genotypes for a subset (fraction for now, but should be also individuals not related to each other)
propGenotyped <- 0.30
geno_idx <- sample(seq_len(nInd(pop_all)), size = floor(propGenotyped * nInd(pop_all)))
pop_geno <- pop_all[geno_idx]

### BOTTLE NECKS OF THIS SCRIPT RUNNING TIME
## PLINK PED/MAP (AlphaSimR writes .ped/.map, can be converted later to .fam/.bed/.bim)
writePlink(pop_geno, baseName = "sim_genotyped_subset", snpChip = 1, simParam = SP, traits = 1, use = "rand",)


