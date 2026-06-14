## ---------------------------------------------------------------------------------
##  PRE-PROCESSING FOR BLUPF90 
##  THIS PART OF THE SCRIPT PRODUCES: PEDIGREE AND PHENOTYPE FILE FOR ALL INDIVIDUALS
##  FOLLOWING BLUPF90 SOFTWARE SUITE REQUIREMENTS
## ----------------------------------------------------------------------------------
cat("\n--------------------------------\n")
cat("Pre-processing for BLUPF90...\n")

# Checking package requirements
required_packages <- c("pedigree", "MCMCglmm", "AlphaSimR", "data.table")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("[FAIL] Required package not found: ", pkg,
                "\n       Install with: install.packages('", pkg, "')"))
  }
  cat(sprintf("  [OK] Package %-12s v%s\n", paste0(pkg, ":"), as.character(packageVersion(pkg))))
}



# check run and folder
args <- commandArgs(trailingOnly = TRUE)
cat(sprintf("[DEBUG] length(args)=%d | raw args: %s\n",
            length(args), paste(args, collapse=", "))) 
# Each run_id creates its own subfolder
run_id <- as.integer(if (length(args) >= 1) args[1] else 1L)  # to count the runs

library(pedigree)
library(MCMCglmm)
library(AlphaSimR)
library(data.table)

# Create file paths
RUN_DIR <- file.path("results", sprintf("run%d", run_id))
BLUP_PED_DIR <- file.path(RUN_DIR, "blupf90_ped")
BLUP_GENO_ASC_DIR <- file.path(RUN_DIR, "blupf90_geno_asc")
BLUP_GENO_CMC_DIR <- file.path(RUN_DIR, "blupf90_geno_cmc")
BLUP_COR_DIR <- file.path(RUN_DIR, "blupf90_cor_asc_cmc")
dir.create(BLUP_GENO_ASC_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(BLUP_PED_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(BLUP_GENO_CMC_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(BLUP_COR_DIR, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("[INFO] Checking directory for run: %s/\n\n", normalizePath(RUN_DIR)))
cat(sprintf("[INFO] Output directories for BLUP90: %s/ \n\n %s/\n\n %s/ \n\n", normalizePath(BLUP_GENO_ASC_DIR), normalizePath(BLUP_PED_DIR), normalizePath(BLUP_GENO_CMC_DIR)))

# Loading pedigree
# the pedigree file is .csv file starting with sim_pedigree
# inside RUN_DIR
# Check all the file that satisfy these conditions
ped_files <- list.files(RUN_DIR, pattern = "^sim_pedigree.*\\.csv$", full.names = TRUE)

# *\\ - this is to escape the backslash in the pattern, so it is treated as a literal character

ped <- read.csv(ped_files)
# data.table
setDT(ped)
# Store ids for generations to analyze
ped[, generations := pedigree::countGen(ped)]
ped_sel_ids <- ped[generations %in% c(2,3,4)]$id

# Load pheno file: sim_child_phenotypes csv
pheno_files <- list.files(RUN_DIR, pattern = "^sim_child_phenotypes.*\\.csv$", full.names = TRUE)
pheno <- read.csv(pheno_files)
setDT(pheno)

# Converting only selected phenotypes
pheno_sel <- pheno[id %in% ped_sel_ids]

# trimming ped
ped_sel <- MCMCglmm::prunePed(ped, ped_sel_ids)
# check na, turn to 0
sum(is.na(ped_sel))

# ped exported as id, father, mother
# no headers
# separeted by white spaces, no comma, no quotes
ped_file_ped <- file.path(BLUP_PED_DIR, "pedigree_blupf90.txt")
ped_file_geno_asc <- file.path(BLUP_GENO_ASC_DIR, "pedigree_blupf90.txt")
ped_file_geno_cmc <-file.path(BLUP_GENO_CMC_DIR, "pedigree_blupf90.txt")
ped_file_cor <- file.path(BLUP_COR_DIR, "pedigree_blupf90.txt")
fwrite(ped_sel[, c(1,3,2)], ped_file_ped, sep = " ", quote = FALSE, col.names = FALSE)
fwrite(ped_sel[, c(1,3,2)], ped_file_geno_asc, sep = " ", quote = FALSE, col.names = FALSE)
fwrite(ped_sel[, c(1,3,2)], ped_file_geno_cmc, sep = " ", quote = FALSE, col.names = FALSE)
fwrite(ped_sel[, c(1,3,2)], ped_file_cor, sep = " ", quote = FALSE, col.names = FALSE)

# export pheno asc
# no headers, space separeted, no quotes, no comma
# 0 here is missing, ASC == 0 --> 1, ASC == 1 --> 2
# id, ASC, sex

pheno_sel[ASC == 1]$ASC <- 2
pheno_sel[ASC == 0]$ASC <- 1
pheno_file_ped <- file.path(BLUP_PED_DIR,"pheno_blupf90.txt")
pheno_file_geno_asc <- file.path(BLUP_GENO_ASC_DIR,"pheno_blupf90.txt")

fwrite(x = pheno_sel[, c(1, 7, 8)], pheno_file_ped, sep = " ", quote = FALSE, col.names = FALSE)
fwrite(x = pheno_sel[, c(1, 7, 8)], pheno_file_geno_asc, sep = " ", quote = FALSE, col.names = FALSE)

# mother phenotype cmc - to store mybe forcorrelation idk useless for now
pheno_files_mother <- list.files(RUN_DIR, pattern = "^sim_mom_phenotypes.*\\.csv$", full.names = TRUE)
cmc <- fread(pheno_files_mother)
cmc <- cmc[id %in% ped_sel_ids]
cmc[CMC == 1]$CMC <- 2
cmc[CMC == 0]$CMC <- 1
pheno_cmc_mother <- file.path(RUN_DIR, "pheno_mother_cmc.txt")
fwrite(cmc, pheno_cmc_mother, sep=" ", quote = FALSE, col.names = FALSE)

# pheno file pregCMC
pheno_sel[preg_CMC == 1]$preg_CMC <- 2
pheno_sel[preg_CMC == 0]$preg_CMC <- 1
pheno_pregCMC <- pheno_sel[, c(3,5)]
pheno_preg_cmc <- file.path(BLUP_GENO_CMC_DIR, "pheno_preg_cmc.txt")
fwrite(pheno_pregCMC, pheno_preg_cmc, sep=" ", quote = FALSE, col.names = FALSE )

# pheno for correlation
pheno_cor_asc_cmc <- file.path(BLUP_COR_DIR, "pheno_cor_asc_cmc.txt")
fwrite(pheno_sel[, c(3, 1, 5, 7, 8)], pheno_cor_asc_cmc, sep=" ", quote = FALSE, col.names = FALSE)

phen <- fread(pheno_cor_asc_cmc)
phen_long <- phen[, .(ID = c(V1, V2), preg_cmc = c(V3, rep(0L, .N)), child_asc = c(rep(0L, .N), V4), child_sex = c(rep(0L, .N), V5))]
# col order now should be: ID, preg_cmc, child_asc and child_sex
phen_long_txt <- file.path(BLUP_COR_DIR, "pheno_asc_cmc_LONG.txt")
fwrite(phen_long, phen_long_txt, sep = " ", col.names = F, quote = F)


# -------------------------------------------------------------------------
# Write renumf90.par for BLUP_PED (pedigree-based animal model, no SNPs)
# -------------------------------------------------------------------------
renum_ped_content <- "# parameter file for renumf90 - pedigree-based animal model (BLUP_PED)
DATAFILE
pheno_blupf90.txt
TRAITS
2
FIELDS_PASSED TO OUTPUT
1
WEIGHT(S)

RESIDUAL_VARIANCE
1.00
EFFECT
3 cross alpha #sex
EFFECT
1 cross alpha # id
RANDOM
animal  # referring to the last effect. Additive genetic variable
FILE
pedigree_blupf90.txt  # pedigree file linked to random var effect id
FILE_POS
1 2 3 0 0 # ped can have 5 cols, alternative mother for embryo transfer col 4, col 5 something else
OPTION cat 2
OPTION method VCE # variance component estimation, to add otherwise just breeding values
OPTION num_threads_pcg 8
OPTION nthreads 8
OPTION use_yams # modified method to calculate matrix, complex model
OPTION fact_once memory # avoid iterative calculation coleski factor, saves in the memory, Faster time. RAM
OPTION approx_log_like # doesnt compute exact log likelihood
OPTION EM-REML 10  # restricted max likelihood
OPTION se_covar_function h2 G_2_2_1_1/(G_2_2_1_1+R_1_1)  # calculate heritability

"

renum_par_ped <- file.path(BLUP_PED_DIR, "renumf90.par")
writeLines(renum_ped_content, con = renum_par_ped)
cat(sprintf("[INFO] renumf90.par written to: %s\n", normalizePath(renum_par_ped)))


# -----------------------------------------------------------------------------
# Write renumf90.par for BLUP_GENO - pregsf90 - create parameters for h matrix
# -----------------------------------------------------------------------------
renum_geno_content <- "# parameter file for renumf90 - single-step ssGBLUP (BLUP_GENO)
DATAFILE
pheno_blupf90.txt
TRAITS
2
FIELDS_PASSED TO OUTPUT
1
WEIGHT(S)

RESIDUAL_VARIANCE
1.00
EFFECT
3 cross alpha #sex
EFFECT
1 cross alpha # id
RANDOM
animal  # referring to the last effect. Additive genetic variable
FILE
pedigree_blupf90.txt  # pedigree file linked to random var effect id
FILE_POS
1 2 3 0 0 # ped can have 5 cols, alternative mother for embryo transfer col 4, col 5 something else
PLINK_FILE
sim_genotyped_30k
OPTION nthreads 6
OPTION use_yams # modified method to calculate matrix, complex model
"

renum_par_geno <- file.path(BLUP_GENO_ASC_DIR, "renumf90_pregibbs.par")
writeLines(renum_geno_content, con = renum_par_geno)
cat(sprintf("[INFO] renumf90.par written to: %s\n", normalizePath(renum_par_geno)))

# -------------------------------------------------------------------------
# Write renumf90_gibbs_asc.par for BLUP_GENO_ASC
# This is run AFTER preGSf90 to produce the final renf90.par for gibbsf90+
# It uses the expanded model with maternal genetic and maternal PE effects
# -------------------------------------------------------------------------
renum_gibbs_asc_content <- "# parameter file for renumf90 - expanded ASC model with maternal effects (BLUP_GENO_ASC)
DATAFILE
pheno_blupf90.txt
TRAITS
2
FIELDS_PASSED TO OUTPUT
1
WEIGHT(S)

RESIDUAL_VARIANCE
1.00
EFFECT
3 cross alpha #sex
EFFECT
1 cross alpha # id - direct genetic
RANDOM
animal
OPTIONAL
mat mpe # maternal genetic and maternal permanent environment
FILE
pedigree_blupf90.txt
FILE_POS
1 2 3 0 0
(CO)VARIANCES
0.20 0.10
0.10 0.15
(CO)VARIANCES_MPE
0.05
OPTION SNP_file sim_genotyped_30k.ped sim_genotyped_30k.ped_XrefID
OPTION readGimA22i
OPTION nthreads 8
OPTION use_yams
OPTION cat 2
OPTION method VCE
OPTION num_threads_pcg 8
OPTION fact_once memory
OPTION approx_log_like
OPTION EM-REML 10
"

renum_gibbs_asc_par <- file.path(BLUP_GENO_ASC_DIR, "renumf90_gibbs_asc.par")
writeLines(renum_gibbs_asc_content, con = renum_gibbs_asc_par)
cat(sprintf("[INFO] renumf90_gibbs_asc.par written to: %s\n", normalizePath(renum_gibbs_asc_par)))

# -------------------------------------------------------------------------
# Write renumf90_gibbs_pregcmc.par for BLUP_GENO_CMC
# This is run to produce the renf90.par for gibbsf90+ for the CMC model
# It uses the CMC model with permanent environment effect
# -------------------------------------------------------------------------
renum_gibbs_pregcmc_content <- "# parameter file for renumf90 - CMC model with permanent environment (BLUP_GENO_CMC)
DATAFILE
pheno_preg_cmc.txt
TRAITS
2
FIELDS_PASSED TO OUTPUT
1
WEIGHT(S)

RESIDUAL_VARIANCE
1.00
EFFECT
1 cross alpha # id
RANDOM
animal
OPTIONAL
pe # permanent environment for repeated measures per mother
FILE
pedigree_blupf90.txt
FILE_POS
1 2 3 0 0
(CO)VARIANCES
0.08
(CO)VARIANCES_PE
0.05
OPTION SNP_file sim_genotyped_30k.ped sim_genotyped_30k.ped_XrefID
OPTION readGimA22i
OPTION nthreads 8
OPTION use_yams
OPTION cat 2
OPTION method VCE
OPTION num_threads_pcg 8
OPTION fact_once memory
OPTION approx_log_like
OPTION EM-REML 10
"

renum_gibbs_pregcmc_par <- file.path(BLUP_GENO_CMC_DIR, "renumf90_gibbs_pregcmc.par")
writeLines(renum_gibbs_pregcmc_content, con = renum_gibbs_pregcmc_par)
cat(sprintf("[INFO] renumf90_gibbs_pregcmc.par written to: %s\n", normalizePath(renum_gibbs_pregcmc_par)))


# -------------------------------------------------------------------------
# Write renumf90_gibbs_cor.par for BLUP_GENO_COR
# This is run to produce the renf90.par for gibbsf90+ for the correlation model
# 
# -------------------------------------------------------------------------

renum_gibbs_correlation_content <- "# parameter file for renumf90
DATAFILE
pheno_asc_cmc_LONG.txt
TRAITS
2 3
FIELDS_PASSED TO OUTPUT
1
WEIGHT(S)

RESIDUAL_VARIANCE
1.0 0
0 1.0
EFFECT
0 4 cross alpha      # child_sex
EFFECT
1 1 cross alpha     # direct additive: CMC mother, ASC child
RANDOM
animal
OPTIONAL
pe mat mpe
FILE
pedigree_blupf90.txt
FILE_POS
1 2 3 0 0
(CO)VARIANCES
3.50  1.00  0.00  0.70
1.00  0.40  0.00  0.20
0.00  0.00  0.00  0.00
0.70  0.20  0.00  0.20
(CO)VARIANCES_PE
0.06  0.00
0.00  0.00
(CO)VARIANCES_MPE
0.00  0.00
0.00  0.06
OPTION SNP_file sim_genotyped_30k.ped sim_genotyped_30k.ped_XrefID
OPTION readGimA22i
OPTION cat 2 2
OPTION blksize 2
OPTION save_halfway_samples 1000
OPTION method VCE # variance component estimation, to add otherwise just breeding values 
OPTION num_threads_pcg 8
OPTION nthreads 8
OPTION use_yams # modified method to calculate matrix, complex model
OPTION fact_once memory # avoid iterative calculation Choleski factor, saves in the memory, Faster time. RAM
OPTION approx_log_like # doesnt compute exact log likelihood
OPTION EM-REML 10  # restricted max likelihood
"
renum_gibbs_correlation <- file.path(BLUP_COR_DIR, "renumf90_gibbs_correlation.par")
writeLines(renum_gibbs_correlation_content, con = renum_gibbs_correlation)
cat(sprintf("[INFO] renumf90_gibbs_correlation, written to: %s\n", normalizePath(renum_gibbs_correlation)))
