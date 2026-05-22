### -----------------------------
### MODIFIED VERSION OF THE writePLINK() FUNCTION
### IN AlphaSimR
### ------------------------------



library(AlphaSimR)



my_writePlink <- function(pop, baseName, traits = 1, use = "pheno", snpChip = 1, 
                          useQtl = FALSE, simParam = NULL, ...){
  
  getResponse <- AlphaSimR:::getResponse # internal helper
  if (pop@ploidy != 2L) {
    stop("writePlink() only supports ploidy=2")
  }
  if (is.null(simParam)) {
    simParam = get(x = "SP", envir = .GlobalEnv)
  }
  y = getResponse(pop = pop, trait = traits, use = use, simParam = simParam, 
                  ...) # to be exported manually
  if (useQtl) {
    H1 = pullQtlHaplo(pop = pop, trait = snpChip, haplo = 1, 
                      asRaw = TRUE, simParam = simParam)
    H2 = pullQtlHaplo(pop = pop, trait = snpChip, haplo = 2, 
                      asRaw = TRUE, simParam = simParam)
    map = getQtlMap(trait = snpChip, simParam = simParam)
  }
  else {
    H1 = pullSnpHaplo(pop = pop, snpChip = snpChip, haplo = 1, 
                      asRaw = TRUE, simParam = simParam)
    H2 = pullSnpHaplo(pop = pop, snpChip = snpChip, haplo = 2, 
                      asRaw = TRUE, simParam = simParam)
    map = getSnpMap(snpChip = snpChip, simParam = simParam)
  }
  sex = pop@sex
  sex[which(sex == "H")] = "0"
  sex[which(sex == "M")] = "1"
  sex[which(sex == "F")] = "2"
  ped <- getPed(pop)
  father = ped$father                               ### 
  father[is.na(father)] = "0"                       ###    HERE LIES THE SOLUTION 
  mother = ped$mother                               ###    REGARDING THE PARENTS
  mother[is.na(mother)] = "0"                       ###
  fam = rbind(rep("1", pop@nInd), pop@id, father, mother, sex, 
              as.character(c(y)))
  H = unname(rbind(t(H1), t(H2)))
  if (ncol(H) == 1) {                                         #
    H = H[c(matrix(1:nrow(H), nrow = 2, byrow = T)), ]        # IMPORTANT TO AVOID TRUNCATION
    H = matrix(H, ncol = 1)                                    # 
  } else {
    H = H[c(matrix(1:nrow(H), nrow = 2, byrow = T)), ]
  }
  rm(H1, H2)
  H = ifelse(H, "2", "1")
  H = rbind(fam, H)
  write(H, file = paste0(baseName, ".ped"), ncolumns = nrow(H))
  map = rbind(map$chr, map$id, as.character(map$pos * 100), 
              as.character(map$site))
  write(map, file = paste0(baseName, ".map"), ncolumns = nrow(map))
  return(invisible())
  
  
}