## ---------------------------------------------------------------------------------
## THIS SCRIPT CONTAINS GENERIC HELPER FUNCTIONS USED FOR MEMORY REPORTING
## AND FOR SOME OF THE MATH BEHIND THE SIMULATION.
## make_pd_corr, calb_alpha, alloc_kids_blocks, rztpois_mean ARE THE CONTRIBUTION OF
## JUAN CORDERO, Postdoc
## Department of Biomedicine, Aarhus University
## ---------------------------------------------------------------------------------


## ============================================================
## Helper: memory reporter
## ============================================================
mem_report <- function(obj, label) {
  sz <- object.size(obj)
  cat(sprintf("  [MEM] %-35s %s\n", paste0(label, ":"), format(sz, units = "auto", standard = "SI")))
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


## ----------------------
## STORES A SUMMARY AS CSV
## ---------------------
save_run_summary <- function(
    OUT_DIR, run_id, file_description,
    ## Genome
    nChr, segSites, nSnpPerChr, Ne, minSnpFreq,
    ## Phenotype / genetic components
    prevCMC_mother, prevPregCMC, prevASC_child,
    ASC_male_ratio, betaPreg,
    var_d, var_m, var_c, var_e,
    rho_CMC_d, rho_CMC_m, rho_d_m,
    nCausalPerChr, kappa, delta_CMC,
    ## Observed outcomes
    pheno_mom, pheno_kid, ASC_kid, is_male, is_female
) {
  
  summary_df <- data.frame(
    ## Run info
    run_id                = run_id,
    seed                  = run_id,
    
    ## Genome parameters
    nChr                  = nChr,
    segSites              = segSites,
    nSnpPerChr            = nSnpPerChr,
    Ne                    = Ne,
    minSnpFreq            = minSnpFreq,
    
    ## Prevalence parameters
    prevCMC_mother        = prevCMC_mother,
    prevPregCMC           = prevPregCMC,
    prevASC_child         = prevASC_child,
    ASC_male_ratio        = ASC_male_ratio,
    
    ## Pregnancy effect
    betaPreg              = betaPreg,
    
    ## Variance components
    var_d                 = var_d,
    var_m                 = var_m,
    var_c                 = var_c,
    var_e                 = var_e,
    
    ## Genetic correlations
    rho_CMC_d             = rho_CMC_d,
    rho_CMC_m             = rho_CMC_m,
    rho_d_m               = rho_d_m,
    
    ## Causal architecture / CMC model
    nCausalPerChr         = nCausalPerChr,
    kappa                 = kappa,
    delta_CMC             = delta_CMC,
    
    ## Observed outcomes — mothers
    n_mothers             = nrow(pheno_mom),
    n_CMC_mothers         = sum(pheno_mom$CMC),
    obs_prev_CMC_mother   = mean(pheno_mom$CMC),
    
    ## Observed outcomes — children
    n_children            = nrow(pheno_kid),
    n_ASC                 = sum(pheno_kid$ASC),
    obs_prev_ASC_overall  = mean(ASC_kid),
    obs_prev_ASC_male     = mean(ASC_kid[is_male]),
    obs_prev_ASC_female   = mean(ASC_kid[is_female])
  )
  
  out_file <- file.path(OUT_DIR, sprintf("sim_run_summary_%s.csv", file_description))
  write.csv(summary_df, file = out_file, row.names = FALSE)
  cat(sprintf("  [OK] Run summary saved: %s\n", basename(out_file)))
  
  invisible(summary_df)
}