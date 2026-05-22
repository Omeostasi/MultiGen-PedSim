###
### EXAMPLE ON HOW TO MAKE BLUPF90 FILES from ped and pheno
###



# Missing values blupf90 treated as 0. Check if there are
# if present, set value for it, e.g. 999

sum(is.na(sim_pedigree_Ne10000_seg1500_mSnpF0_005_Founders850000_nGen_5_run01))



# phenotype file
library(pedigree)
library(MCMCglmm)
library(data.table)

# Get id generations 3,4,5
ped <- sim_pedigree_Ne10000_seg1500_mSnpF0_005_Founders850000_nGen_5_run01
setDT(ped)
ped[, generations := pedigree::countGen(ped)]

ped_short_ids <- ped[generations %in% c(2,3,4)]$id

pheno <- sim_child_phenotypes_Ne10000_seg1500_mSnpF0_005_Founders850000_nGen_5_run01
setDT(pheno)
pheno_short <- pheno[id %in% ped_short_ids]

# trimming ped
ped_short <- MCMCglmm::prunePed(ped, ped_short_ids)
# check na, turn to 0
sum(is.na(ped_short))

# ped exported as id, father, mother
# no headers
# separeted by white spaces, no comma, no quotes
fwrite(ped_short[, c(1,3,2)], "pedigree_blupf90.txt", sep = " ", quote = FALSE, col.names = FALSE)

# load old pop
library(AlphaSimR)
rm(sim_child_phenotypes_Ne10000_seg1500_mSnpF0_005_Founders850000_nGen_5_run01, sim_pedigree_Ne10000_seg1500_mSnpF0_005_Founders850000_nGen_5_run01 )
pop_all <- readRDS("../simulation_final/pop_Ne10000_seg1500_mSnpF0.005_founders850000_nGen5_run01_old.rds")


pop <- pop_all$pop

# add sex column
ids_male <- as.integer(pop@id[isMale(pop)])
ids_female <- as.integer(pop@id[isFemale(pop)])
rm(pop)
rm(pop_all)

pheno_short[, sex := as.integer(1)]
pheno_short[id %in% ids_female]$sex <- 2

all(pheno_short[sex == 1]$id %in% ids_male)

# export pheno
# no headers, space separeted, no quotes, no comma
# 0 here is missing, ASD == 0 --> 1, ASD == 1 --> 2
# id, ASD, sex

pheno_short[ASD == 1]$ASD <- 2
pheno_short[ASD == 0]$ASD <- 1
fwrite(x = pheno_short[, c(1, 7, 8)], "pheno_blupf90.txt", sep = " ", quote = FALSE, col.names = FALSE)

