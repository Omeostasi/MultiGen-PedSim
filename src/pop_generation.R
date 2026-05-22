library(AlphaSimR)



### ---------
###  IN VERSION 3, THESE FUNCCTIONS USE VECTORS AS INPUTS FOR lambdaKids, p_new_partner, mean_unions_dad and mothers_fraction
###  IN VERSION FINAL, THESE FUNCTIONS HAVE BEEN MODIFIED TO REMOVE THE FOUNDERS AND 80% OF THE SILENT GENERATION FROM THE FINAL OUTPUT 
###  THROUGH THE PARAMETER remove_older_generations=TRUE
### ---------




### ----------------------------------------------------
### THIS IS FUNCTION IS NOT MEANT TO BE USED DIRECTLY
### USED BY pop_generation()
### GENERATES FIRST POP FROM FOUNDERS
###  ----------------------------------------------------

pop_generation_from_founders <- function(pop_founder_haplo, mothers_fraction, lambdaKids = 2.3, maxPartners_mom = 3, p_new_partner = 0.05, mean_unions_dad = 1.3){
  
  # Select male and females from founders
  pop_founder <- newPop(pop_founder_haplo)
  fem_founder <- pop_founder[isFemale(pop_founder)]
  mal_founder <- pop_founder[isMale(pop_founder)]
  
  # Select a fraction to be a mother
  mothers <- fem_founder[sample(seq_len(nInd(fem_founder)), mothers_fraction*nInd(fem_founder))]
  mother_ids <- mothers@id
  dad_ids <- mal_founder@id
  
  ## -------------------------
  ## Build family structure with unions (divorce/remarriage)
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
  children <- makeCross(pop_founder, crossPlan = crossPlan, nProgeny = 1, simParam = SP)
  
  stopifnot(nInd(children) == length(child_moms))
  
  parents_ids <- pop_founder@id
  children_ids <- children@id
  
  pop_from_founders <- c(pop_founder, children)
  
  # Storing unions
  unions_all <- unions
  
  # Returning the pop object and Ids to distinguish different generations
  return(list(pop_from_founders = pop_from_founders, parents_ids = parents_ids, children_ids = children_ids, unions_all = unions_all))
  
}

# The pop increases of 1.6, roughly


###
### THIS FUNCTION IS GOING TO BE USED TO GENERATE POP AND STORE THEIR RESULTS
### THIS IS SUPPOSE TO GENERATE THE FINAL PEDIGREE THAT IS GOING TO BE USED
###

pop_generation <- function(pop_founder_haplo, overlapping_fraction = 0, nGenerations_pop = 3, mothers_fraction = c(0.75,0.75,0.75),
                           lambdaKids = c(2.3, 2.0, 1.6), maxPartners_mom = 3, p_new_partner = c(0.05, 0.10, 0.25),
                           mean_unions_dad = c(1.3,1.3,1.3)) {
  #
  # THIS FUNCTION GENERATES FROM THE HAPLOTYPES A CERTAIN AMOUNT OF GENERATIONS 
  # AND STORES EACH GENERATION DATA
  #
  
  # This stores the data of children_ids to allow pulling genotypes later
  children_ids_list <- vector("list", length = nGenerations_pop)
  
  
  
  
  # This generates the first pop starting from founders
  # Uses the helper function created above
  from_founders <- pop_generation_from_founders(pop_founder_haplo, mothers_fraction = mothers_fraction[1],
                                                lambdaKids = lambdaKids[1], maxPartners_mom = maxPartners_mom,
                                                p_new_partner = p_new_partner[1], mean_unions_dad = mean_unions_dad[1])
  # from_founders is an object containing ids, unions and founders
  pop_from_founders <- from_founders[["pop_from_founders"]]
  
  generations_ids <- vector(mode = "list", length = nGenerations_pop + 1)
  names(generations_ids) <- seq(nGenerations_pop +1)
  
  generations_ids[1] <- from_founders["parents_ids"]
  generations_ids[2] <- from_founders["children_ids"]
  
  # Store children ids
  children_ids_list[1] <- from_founders["children_ids"]
  
  # Store unions list 
  unions_all <- from_founders$unions_all
  
  # Store all the pop data
  pop_all <- pop_from_founders
  
  # This for loop generates all the remaining generations
  for (generation in 1:nGenerations_pop){
    
    # Skipping the first generation since it has been made
    if (generation == 1) { cat(sprintf("Generation 1 out of %d \n", nGenerations_pop)); next }
    
    # Reporting steps
    cat(sprintf("Generation %d out of %d \n", generation, nGenerations_pop))
    
    # Taking a fraction from previous generation to reuse
    
    n_to_sample <- round(overlapping_fraction * length(generations_ids[[generation - 1]]))
    previous_generation_fraction_ids <- sample(generations_ids[[generation - 1]], n_to_sample)
    fem_ids_previous_fraction <- previous_generation_fraction_ids[isFemale(pop_all[previous_generation_fraction_ids])]
    mal_ids_previous_fraction <- previous_generation_fraction_ids[isMale(pop_all[previous_generation_fraction_ids])]
    
    # Taking ids from current generation
    pop_current_generation_ids <- generations_ids[[generation]]
    fem_ids <- pop_current_generation_ids[isFemale(pop_all[pop_current_generation_ids])]
    mal_ids <- pop_current_generation_ids[isMale(pop_all[pop_current_generation_ids])]
    
    # Combine ids of fraction and current
    fem_ids <- c(fem_ids_previous_fraction, fem_ids)
    mal_ids <- c(mal_ids_previous_fraction, mal_ids)
    
    # Select a fraction to be mothers and fathers
    n_to_sample <- round(mothers_fraction[generation] * length(fem_ids))
    mother_ids <- sample(fem_ids, n_to_sample)
    dad_ids <- mal_ids
    
    ###
    ### FROM HERE ON SHOULD BE EQUAL TO THE USUAL
    ###
    
    ## -------------------------
    ## Build family structure with unions (divorce/remarriage)
    ##    Produces full-sib blocks + maternal and paternal half-sibs.
    ## -------------------------
    
    ## Kids per mother: mean(lambdaKids) among mothers with >=1 child
    kids_per_mom <- rztpois_mean(length(mother_ids), mean_target = lambdaKids[generation])
    
    ## Number of partners per mother (1..maxPartners_mom), capped by number of kids
    nPartners_mom <- 1 + rbinom(length(mother_ids), size = (maxPartners_mom - 1), prob = p_new_partner[generation])
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
    dad_capacity <- pmax(1L, rpois(length(dad_ids), lambda = mean_unions_dad[generation]))
    names(dad_capacity) <- as.character(dad_ids)
    
    
    
    ## Assign a father to each union (sampling from remaining capacity)
    assigned_dads <- character(nrow(unions))
    for (u in seq_len(nrow(unions))) {
      eligible <- names(dad_capacity)[dad_capacity > 0]
      if (length(eligible) == 0) {
        ## If depleted, refresh capacities (or increase mean_unions_dad)
        dad_capacity[] <- pmax(1L, rpois(length(dad_ids), lambda = mean_unions_dad[generation]))
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
    children <- makeCross(pop_all, crossPlan = crossPlan, nProgeny = 1, simParam = SP)
    
    stopifnot(nInd(children) == length(child_moms))
    
    generations_ids[[generation + 1]] <- children@id
    
    pop_all <- c(pop_all, children)
    
    # Storing children ids and unions
    children_ids_list[[generation]] <- children@id
    unions_all <- rbind(unions_all, unions)
  }
  
  
  return(list(pop = pop_all, unions_all = unions_all, children_ids_list = children_ids_list, older_generations_ids = c(generations_ids[1], generations_ids[2])))
}
