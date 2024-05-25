bob_batch <- function(x_likelihood, x_means, x_covs, x_probs,
                      y, means, covs, probs, Sbatch, a, b, c, r,
                      alphas, betas, lambdas, nus, psis, cores, seed = 1){

  base::RNGkind("L'Ecuyer-CMRG")
  set.seed(seed) # For replication purposes

  out <- parallel::mclapply(1:Sbatch, function(s)
    batchOneShot(x_likelihood, x_means, x_covs, x_probs,
      y, means, covs, probs, a, b, c, r, alphas, betas, lambdas, nus, psis),
        mc.set.seed = TRUE, mc.cores = cores)

  parallel::mc.reset.stream()

  do.call(rbind, out)
}


get_indices <- function(p, K, drop.covs = TRUE){

  seqCovs <- seq_len(K * p^(2)) + (K * p)

  if (drop.covs) {
    out <- c(seq_len(p * K),                               # Means indices
             seqCovs[c(diag(p) == 1)],                     # Variances indices
             (K * (p^(2) + p) + 1):(K * (p^(2) + p) + K))  # Probs. indices
  } else {
    out <- seqCovs[c(upper.tri(matrix(NA, p, p)))]         # Covs. indices
  }

  out
}


eval_kdes <- function(draws, main.indices, covs.indices, coresKDEs){

  # Means, variances and probs
  kdes.1 <- parallel::mclapply(main.indices,
              \(j) log(FKSUM::fk_density(draws[,j], x_eval = draws[,j])$y),
                mc.cores = coresKDEs)

  # Covariances
  kdes.2 <- parallel::mclapply(covs.indices,
              \(j) log(FKSUM::fk_density(draws[,j], x_eval = draws[,j])$y),
                mc.cores = coresKDEs)

  do.call(sum, kdes.1) + (2 * do.call(sum, kdes.2))
}


KL_div <- function(x, Sbatch, y, means, covs, probs, a, b, c, r,
                   alphas, betas, lambdas, nus, psis,
                   main.indices, covs.indices, coresBO, coresKDEs, seed){

  K <- length(alphas)

  x_likelihood <- x[1]
  x_means      <- x[2:(2 + K - 1)]
  x_covs       <- x[(2 + K):(2 + 2 * K - 1)]
  x_probs      <- x[2 * K + 2]

  # Obtain approx draws for given value of x (in parallel)
  approx_draws <- bob_batch(x_likelihood, x_means, x_covs, x_probs,
                            y, means, covs, probs, Sbatch, a, b, c, r,
                            alphas, betas, lambdas, nus, psis, coresBO, seed)

  #cat("Batch Sampling: Done \n")

  theta.dim <- ncol(approx_draws)

  # Evaluate KDEs (in parallel)
  kdes <- eval_kdes(approx_draws, main.indices, covs.indices, coresKDEs)

  #cat("KDEs Evaluation: Done \n")

  # Evaluate target posterior (in parallel)
  lpost <- parallel::mclapply(1:Sbatch, \(s)
    log_post(approx_draws[s,], y, alphas, betas, lambdas, psis, nus),
      mc.cores = coresBO)

  #cat("Log-posterior Evaluation: Done \n")

  (kdes - do.call(sum, lpost)) / (Sbatch * theta.dim)
}


bob.Sampler <- function(y, means, covs, probs, S, Sbatch, a, b, c, r,
                        alphas, betas, lambdas, nus, psis,
                        lowerBound, upperBound, coresBO, coresSamp,
                        BOiterations = 120, seed = 1){

  K <- length(alphas)
  p <- ncol(y)

  main.indices <- get_indices(p, K, drop.covs = T) # Means, variances and probs
  covs.indices <- get_indices(p, K, drop.covs = F) # Covariances

  coresKDEs <- if (p > 10 && K > 2) min(coresBO, 5L) else 1L

  # Define the objective Lhat(x)
  objective <- function(x){

    loss <- tryCatch(
      {
        KL_div(x, Sbatch, y, means, covs, probs, a, b, c, r,
               alphas, betas, lambdas, nus, psis,
               main.indices, covs.indices, coresBO, coresKDEs, seed)
      }, error = function(e){
        999
      }
    )

    if (is.finite(loss)) {return(min(loss, 999))} else {return(999)}
  }

  # Bayesian Optimization
  set.seed(seed)
  x.star <- bayesOpt(lowerBound, upperBound, y, K, objective, BOiterations)
  rm(main.indices, covs.indices, coresKDEs)

  # Extract optimal values
  x.star_likelihood <- x.star[1]
  x.star_means      <- x.star[2:(2 + K - 1)]
  x.star_covs       <- x.star[(2 + K):(2 + 2 * K - 1)]
  x.star_probs      <- x.star[2 * K + 2]

  # Posterior Sampler
  out <- bob_batch(x.star_likelihood, x.star_means, x.star_covs, x.star_probs,
    y, means, covs, probs, S, a, b, c, r, alphas, betas, lambdas,
      nus, psis, coresSamp, seed)

  list("Draws" = out, "x.optim" = x.star)
}
