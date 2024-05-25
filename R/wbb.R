wbb.Sampler <- function(y, means, covs, probs, S, a, b, c, r,
                        alphas, betas, lambdas, nus, psis, randomPrior,
                        cores, seed = 1){

  base::RNGkind("L'Ecuyer-CMRG")
  set.seed(seed) # For replication purposes

  out <- parallel::mclapply(1:S, function(s)
    wbbOneShot(y, means, covs, probs, a, b, c, r, alphas, betas, lambdas, nus,
               psis, randomPrior), mc.set.seed = TRUE, mc.cores = cores)

  parallel::mc.reset.stream()

  matrix(unlist(out), nrow = S, byrow = TRUE)
}
