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


# Create file paths
RUN_DIR <- file.path("results", sprintf("run%02d", run_id))
BLUP_DIR <- file.path(RUN_DIR, "blupf90")
dir.create(BLUP_DIR, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("[INFO] Checking directory for run: %s/\n\n", normalizePath(RUN_DIR)))
cat(sprintf("[INFO] Output directory for BLUP90: %s/\n\n", normalizePath(BLUP_DIR)))

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
ped_file <- file.path(BLUP_DIR, "pedigree_blupf90.txt")
fwrite(ped_sel[, c(1,3,2)], ped_file, sep = " ", quote = FALSE, col.names = FALSE)



# export pheno
# no headers, space separeted, no quotes, no comma
# 0 here is missing, ASD == 0 --> 1, ASD == 1 --> 2
# id, ASD, sex

pheno_sel[ASD == 1]$ASD <- 2
pheno_sel[ASD == 0]$ASD <- 1
pheno_file <- file.path(BLUP_DIR,"pheno_blupf90.txt")

fwrite(x = pheno_sel[, c(1, 7, 8)], pheno_file, sep = " ", quote = FALSE, col.names = FALSE)

