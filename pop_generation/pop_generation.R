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
library(tidyverse)

# Adding wrapper function to generate founders
source("src/founderHaplotypesGenerator.R")

# Helpers: make_pd_corr, calib_alpha, alloc_kids_blocks
source("src/helpers_juan.R")

# Function to generate pop pedigree
source("src/pop_generation.R")

set.seed(1)
options(stringsAsFactors = FALSE)

## -------------------------
## 0) Parameters
## -------------------------

## Genome simulation
nFounder   <- 20    # initial number of founders for runMacs2 
nChr       <- 2
segSites   <- 20    # segregating sites per chromosome
nSnpPerChr <- 5     # SNPs per chromosome on chip (observed)
nGenerations_founders <- 5 # the number of generations to generate the founders
Ne <- 10
nFounders_final <- 100 # final number of founders used as SimParam


## Population / study design
# nParents   <- 20      # size of adult generation  # this is also not used in this version
# nMothers   <- 10    # removed the old pure nMothers
mothers_fraction <- 0.75 # fraction of female pop that become mothers
lambdaKids <- 1.6       # mean kids per mother among mothers with >=1 child
nGenerations_pop <- 3 # number generation that we have data for
overlapping_fraction <- 0.05 # 5 %

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

pop_all <- pop_generation(founderHap, 
                          overlapping_fraction = overlapping_fraction, 
                          nGenerations_pop = nGenerations_pop, mothers_fraction = mothers_fraction,
                          lambdaKids = lambdaKids, maxPartners_mom = maxPartners_mom,
                          p_new_partner = p_new_partner, mean_unions_dad = mean_unions_dad)

