library(AlphaSimR)

# HOLD THE FUNCTION TO INCREASE FOUNDER POP THROUGH CROSSING

founderHaplotypesGenerator <- function(nInd, nChr = 22, segSites = 1000, nThreads = nChr, Ne = 10000, nGenerations = 10, method = "runMacs2", nFounders_final = 1*10^6, nSnpPerChr = 50, minSnpFreq = 0.001){
  
  # method should be changed only for quick tests
  if (method == "runMacs2" || method == "run"){
    # Generate the first founder pop using runMacs2
    # Report step
    cat(sprintf("Generating founder haplotypes using runMacs2 with %d nChr, %d segSites, %d nThreads and %d Ne \n", nChr, segSites, nThreads, Ne))
    founderHap <- runMacs2(nInd = nInd, nChr = nChr, segSites = segSites, nThreads = nChr, Ne = Ne)
  } else {
    # Report step
    cat(sprintf("Generating founder haplotypes using quickHaplo with %d nChr and %d segSites \n", nChr, segSites))
    founderHap <- quickHaplo(nInd = nInd, nChr = nChr, segSites = segSites)
  }
  
  ## STORE RESULTS
  filename <- "results/founderHap_Ne_20000_segSites_2000.rds"

  saveRDS(founderHap, filename)
  
  cat(sprintf("Saved founder haplotypes to %s\n", filename))
  
  
  # Set parameter for the simulation
  # Used for the crossing
  
  SP_founder <- SimParam$new(founderHap)
  SP_founder$setTrackPed(TRUE)
  SP_founder$setSexes(sexes = "yes_rand") 
  SP_founder$addSnpChip(nSnpPerChr = nSnpPerChr, minSnpFreq = minSnpFreq)
  
  SP <<- SP_founder
  
  # Create the first pop
  
  founderPop <- newPop(founderHap, simParam = SP_founder)
  
  # Loop through generations to create the final founder pop
  
  for (generation in 1:nGenerations){
    
    # Report step
    cat(sprintf("Generation %d out of %d for founder pop \n", generation, nGenerations))
    
    # Upper limit for the nInd
    # Function capped at 1milion founders or nFounders_final
    if (nInd(founderPop) == nFounders_final || nInd(founderPop) > nFounders_final) {
      
      # Report step
      cat(sprintf("At generation %d the founder pop has %d individuals, which is above the nFounders_final limit. This is not an error. \n", generation, nInd(founderPop)))
      
      founderPop <- randCross(
        pop = founderPop,
        nCrosses = nFounders_final,
        simParam = SP_founder
      )
      
      # If not enough nInd, increase through random crosses
    } else  {
      
      founderPop <- randCross(
        pop = founderPop,
        nCrosses = nInd(founderPop),
        nProgeny = 2,
        simParam = SP_founder
      )
    }
  }
  
  # From founderPop, generate founderHaplotypes
  # Required for new simulation parameters based on the actual founders
  
  haplotypes <- vector("list", length = nChr)
  pulledHaplo <- pullSegSiteHaplo(founderPop, asRaw = TRUE) # asRaw reduced the memory loda
  # Check allocation memory for haplotypes
  cat(sprintf("Pulled haplotypes matrix: %d rows x %d cols | Size: %s\n", 
              nrow(pulledHaplo), ncol(pulledHaplo), format(object.size(pulledHaplo), units = "auto")))  
  start_col <- 1
  end_col <- segSites
  
  
  for (chr in 1:nChr){
    
    # create a list of matrices for each chr
    haplotypes[[chr]] <- pulledHaplo[, start_col:end_col]
    
    start_col <- end_col + 1
    end_col <- start_col + segSites - 1
  }
  
  
  # Storing the final haplo
  founderHap <- newMapPop(genMap = SP_founder$genMap,  haplotypes = haplotypes)
  
  return(founderHap)
}
