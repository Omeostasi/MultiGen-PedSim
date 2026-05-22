### ----------------------------------------------------------------------------------------
### Removes individuals from older generations across pedigree and phenotype data frames.
### When rm_older_generations is TRUE, removes all founders and a random 80% sample
### of the first generation, keeping the remaining data structures consistent.
### ----------------------------------------------------------------------------------------


remove_older_generations <- function(ped, pheno_mom, pheno_kid,
                                     older_generations_ids,
                                     rm_older_generations = TRUE){
  
  # IF remove_older_generations == TRUE, remove first generation of 80% and founders
  
  if (rm_older_generations == TRUE){
    
    
    # Find the individuals that belong to these generations
    
    silent <- sample(older_generations_ids[[1]], 0.80 * length(older_generations_ids[[1]]))
    older_generations <- c(seq(older_generations_ids[[1]]), silent)
    
    # Remove them while maintaining all the others unaltered
    # from ped
    ped <- ped[!(ped$id %in% older_generations), ]
    
    # from pheno_mom
    pheno_mom <- pheno_mom[!(pheno_mom$id %in% older_generations), ]
    
    # pheno_kid
    pheno_kid <- pheno_kid[!(pheno_kid$id %in% older_generations), ]
  }
  
  return(list( ped = ped, pheno_mom = pheno_mom, pheno_kid = pheno_kid))
  
  
}

###
### function to sample unrelated members
###
library(data.table)

sample_unrelated <- function(ped, target_n, subset_ids = NULL, verbose = TRUE) {
  # ped must have columns: id, father, mother
  # IDs must already be integers
  # Missing parents can be NA or 0
 
  ped <- as.data.table(ped)
  ped[, id     := as.integer(id)]
  ped[, father := as.integer(father)]
  ped[, mother := as.integer(mother)]
  
  # ped mother and father that are not in id change them to zero
  valid_ids <- ped$id
  ped[!father %in% valid_ids, father := 0L]
  ped[!mother %in% valid_ids, mother := 0L]
  setorder(ped, id)
  
  required_cols <- c("id", "father", "mother")
  if (!all(required_cols %in% names(ped))) {
    stop("ped must contain columns: id, father, mother")
  }
  
  ped[is.na(father), father := 0L]
  ped[is.na(mother), mother := 0L]
  
  if (anyDuplicated(ped$id)) {
    stop("ped$id contains duplicated IDs")
  }
  
  if (any(ped$id <= 0L, na.rm = TRUE)) {
    stop("ped$id must contain positive integer IDs")
  }

  N <- max(c(ped$id, ped$mother, ped$father), na.rm = TRUE)
  
  # -----------------------------
  # Full pedigree lookup tables
  # -----------------------------
  
  setkey(ped, id)
  
  edges <- rbindlist(list(
    ped[father != 0L, .(parent = father, child = id)],
    ped[mother != 0L, .(parent = mother, child = id)]
  ))
  
  edges <- unique(edges)
  setkey(edges, parent)
  
  parents_of <- function(x) {
    if (!length(x)) return(integer())
    
    p <- ped[J(x), .(father, mother)]
    
    out <- c(p$father, p$mother)
    out <- out[!is.na(out) & out != 0L]
    
    unique(out)
  }
  
  children_of <- function(x) {
    x <- x[!is.na(x) & x != 0L]
    if (!length(x)) return(integer())
    
    unique(edges[J(x), child, nomatch = 0L])
  }
  
  relatives_ge_0125 <- function(x) {
    strong <- integer()
    border <- integer()
    
    ancestors <- data.table(ancestor = x, depth = 0L)
    
    cur <- x
    for (d in 1:3) {
      cur <- parents_of(cur)
      if (!length(cur)) break
      
      ancestors <- rbind(
        ancestors,
        data.table(ancestor = cur, depth = d),
        use.names = TRUE
      )
    }
    
    ancestors <- unique(ancestors)
    
    for (i in seq_len(nrow(ancestors))) {
      d <- ancestors$depth[i]
      cur <- ancestors$ancestor[i]
      
      max_e <- 4L - d
      
      for (e in 0:max_e) {
        if (!length(cur)) break
        
        path_len <- d + e
        
        if (path_len <= 3L) {
          strong <- c(strong, cur)
        } else if (path_len == 4L) {
          border <- c(border, cur)
        }
        
        cur <- children_of(cur)
      }
    }
    
    strong <- strong[strong != x]
    border <- border[border != x]
    
    strong_close <- unique(strong)
    
    if (length(border)) {
      border <- sort.int(border)
      rr <- rle(border)
      border_close <- rr$values[rr$lengths >= 2L]
    } else {
      border_close <- integer()
    }
    
    unique(c(strong_close, border_close))
  }
  
  # -----------------------------
  # Candidate set
  # -----------------------------
  
  if (is.null(subset_ids)) {
    subset_ids <- ped$id
  } else {
    if (is.data.frame(subset_ids) || is.data.table(subset_ids)) {
      if (!"id" %in% names(subset_ids)) {
        stop("If subset_ids is a data.frame/data.table, it must contain column 'id'")
      }
      subset_ids <- subset_ids$id
    }
    
    subset_ids <- unique(subset_ids)
    subset_ids <- subset_ids[!is.na(subset_ids) & subset_ids != 0L]
    subset_ids <- subset_ids[subset_ids %in% ped$id]
  }
  
  if (length(subset_ids) < target_n) {
    stop("subset_ids has fewer IDs than target_n after filtering to ped$id")
  }
  
  if (verbose) {
    message("Full pedigree size: ", nrow(ped))
    message("Eligible candidate set size: ", length(subset_ids))
    message("Target sample size: ", target_n)
  }
  
  # Only these IDs are eligible for sampling.
  # Relationship checks still use the full pedigree above.
  candidates <- sample(subset_ids, length(subset_ids), replace = FALSE)
 
  # -----------------------------
  # Greedy unrelated sampling
  # -----------------------------
  forbidden <- logical(N)
  
  selected <- integer(target_n)
  n_selected <- 0L
  
  for (x in candidates) {
    if (n_selected >= target_n) break
    
    if (forbidden[x]) next
    
    n_selected <- n_selected + 1L
    selected[n_selected] <- x
    
    bad <- relatives_ge_0125(x)
    
    forbidden[x] <- TRUE
    forbidden[bad] <- TRUE
    
    if (verbose && n_selected %% 1000L == 0L) {
      message("Selected: ", n_selected)
    }
  }
  
  selected <- selected[seq_len(n_selected)]
  
  if (n_selected < target_n) {
    warning(
      "Only selected ", n_selected,
      " unrelated individuals out of requested ", target_n,
      ". The candidate set may not contain enough mutually unrelated people."
    )
  }
  
  data.table(id = selected)
}



check_related_pairs <- function(ped, ids, cutoff = 0.0625) {
  setDT(ped)
  
  # convert as integer
  ped[, id     := as.integer(id)]
  ped[, father := as.integer(father)]
  ped[, mother := as.integer(mother)]
  ids <- as.integer(ids)
  
  ped[is.na(father), father := 0L]
  ped[is.na(mother), mother := 0L]
  valid_ids <- ped$id
  ped[!father %in% valid_ids, father := 0L]
  ped[!mother %in% valid_ids, mother := 0L]
  setorder(ped, id)
  
  ids <- unique(ids)
  ids <- ids[ids %in% ped$id]
  
  N <- max(c(ped$id, ped$father, ped$mother), na.rm = TRUE)
  
  id_flag <- logical(N)
  id_flag[ids] <- TRUE
  
  setkey(ped, id)
  
  edges <- rbindlist(list(
    ped[father != 0L, .(parent = father, child = id)],
    ped[mother != 0L, .(parent = mother, child = id)]
  ))
  
  edges <- unique(edges)
  setkey(edges, parent)
  
  parents_of <- function(x) {
    if (!length(x)) return(integer())
    p <- ped[J(x), .(father, mother)]
    out <- c(p$father, p$mother)
    unique(out[!is.na(out) & out != 0L])
  }
  
  children_of <- function(x) {
    x <- x[!is.na(x) & x != 0L]
    if (!length(x)) return(integer())
    unique(edges[J(x), child, nomatch = 0L])
  }
  
  relatives_ge_cutoff <- function(x) {
    contrib <- data.table(who = integer(), r = numeric())
    
    # For cutoff 1/16, we need to search slightly wider than for 1/8.
    # max_path = 5 captures first cousins once removed and similar cases.
    max_path <- 5L
    
    ancestors <- data.table(ancestor = x, depth = 0L)
    
    cur <- x
    for (d in 1:4) {
      cur <- parents_of(cur)
      if (!length(cur)) break
      ancestors <- rbind(
        ancestors,
        data.table(ancestor = cur, depth = d),
        use.names = TRUE
      )
    }
    
    ancestors <- unique(ancestors)
    
    for (i in seq_len(nrow(ancestors))) {
      d <- ancestors$depth[i]
      cur <- ancestors$ancestor[i]
      
      max_e <- max_path - d
      if (max_e < 0L) next
      
      for (e in 0:max_e) {
        if (!length(cur)) break
        
        path_len <- d + e
        
        hit <- cur[cur <= N & id_flag[cur] & cur != x]
        
        if (length(hit)) {
          contrib <- rbind(
            contrib,
            data.table(who = hit, r = 2^(-path_len)),
            use.names = TRUE
          )
        }
        
        cur <- children_of(cur)
      }
    }
    
    if (!nrow(contrib)) {
      return(data.table(id2 = integer(), r = numeric()))
    }
    
    contrib[
      ,
      .(r = sum(r)),
      by = who
    ][r >= cutoff, .(id2 = who, r)]
  }
  
  out <- list()
  k <- 0L
  
  for (i in seq_along(ids)) {
    x <- ids[i]
    
    rel <- relatives_ge_cutoff(x)
    
    if (nrow(rel)) {
      rel <- rel[id2 > x]  # avoid duplicate pairs: x-y and y-x
      
      if (nrow(rel)) {
        k <- k + 1L
        out[[k]] <- data.table(
          id1 = x,
          id2 = rel$id2,
          r = rel$r
        )
      }
    }
    
    if (i %% 1000L == 0L) {
      message("Checked: ", i, " / ", length(ids))
    }
  }
  
  if (!k) {
    return(data.table(id1 = integer(), id2 = integer(), r = numeric()))
  }
  
  rbindlist(out)
}







###
### Selects a genotyped subset from the two most recent generations of a population.
### Prioritises ASD cases (up to n_geno_ASD) and supplements with a random sample
### of non-ASD individuals (up to n_geno_random), returning a Pop object of the
### combined selection. Warns if fewer individuals are available than requested.
### Inds are chosen avoiding close relatives using the previous 2 functions.
###
genotype_subset <- function(pop, pheno_kid, generations_ids, ped, nGenerations_pop = 5, n_geno_random = 125000, n_geno_ASD = 25000){
  
  lastGeneration <- generations_ids[[nGenerations_pop]]
  secondLastGeneration <- generations_ids[[nGenerations_pop - 1]]
  subset_genotyped <- c(secondLastGeneration, lastGeneration)
  pheno_genotyped <- pheno_kid[pheno_kid$id %in% subset_genotyped,]
  
  # select subset of ASD cases from children
  asd_child_ids <- as.integer(pheno_genotyped[pheno_genotyped$ASD == 1,]$id)
  
  
  # select subset not close-related.
   selected_asd_ids <- sample_unrelated(ped = ped, target_n = n_geno_ASD, subset_ids = asd_child_ids)
   n_asd_avail   <- length(selected_asd_ids$id)
  
  
  # select subset random individuals from the full population,
  #           excluding the already-selected ASD cases
  remaining_ids    <- as.integer(setdiff(pheno_genotyped$id, selected_asd_ids))
  n_rand_avail     <- length(remaining_ids)
  selected_rand_ids <- sample_unrelated(ped = ped, target_n = n_geno_random, subset_ids = remaining_ids)
  
  # combine and test for related
  selected_ids <- c(selected_asd_ids$id, selected_rand_ids$id)
  related_ids <- check_related_pairs(ped = ped, ids = selected_ids, cutoff = 0.125)
  
  
  # remove related ids that do not have autism
  to_remove <- c()
  if(nrow(related_ids) > 0){
    for (i in 1:nrow(related_ids)){
      # check if id1 has autism
      if (related_ids[i,]$id1 %in% asd_child_ids) {
        to_remove <- c(to_remove, related_ids[i,]$id2)
      } else { # remove id1 if doesn't have autism
        to_remove <- c(to_remove, related_ids[i,]$id1)
      }
    } 
    selected_ids[-to_remove]}
  
  
  # create pop object that is going to be genotyped
  pop_geno     <- pop[selected_ids]
  
  cat(sprintf("  [INFO] Genotyped subset:  %d total\n", length(selected_ids)))
  
  return(pop_geno)
  
}