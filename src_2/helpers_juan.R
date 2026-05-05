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