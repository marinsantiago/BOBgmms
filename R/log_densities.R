extractCovs <- function(theta, d, K){

  start_index <- d * K
  out <- vector("list", length = K)

  for (k in seq_len(K)) {
    start_k  <- start_index + ((k - 1) * d^(2)) + 1
    end_k    <- start_index + (k * d^(2))
    out[[k]] <- matrix(theta[start_k:end_k], nrow = d, ncol = d)
  }

  return(out)
}


log_post <- function(theta, y, alphas, betas, lambdas, psis, nus){

  d <- ncol(y)
  K <- length(alphas)

  # Extract model parameters
  means <- matrix(theta[1:(d * K)], ncol = d, byrow = TRUE)
  covs  <- extractCovs(theta, d, K)
  probs <- tail(theta, K)

  # Log-likelihood
  log.likelihood <- sum(mvnfast::dmixn(y, means, covs, w = probs, log = TRUE))

  # Log-prior
  log.prior <- sum(sapply(1:K, \(k){
    mvnfast::dmvn(means[k,], betas[[k]], covs[[k]] / lambdas[[k]], log = T) +
      LaplacesDemon::dinvwishart(covs[[k]], nus[k], psis[[k]], log = T)
  })) + LaplacesDemon::ddirichlet(probs, alphas, log = T)

  log.likelihood + log.prior
}


log_lik <- function(theta, y, d, K){

  # Extract model parameters
  means <- matrix(theta[1:(d * K)], ncol = d, byrow = TRUE)
  covs  <- extractCovs(theta, d, K)
  probs <- tail(theta, K)

  sum(mvnfast::dmixn(y, means, covs, w = probs, log = TRUE))
}
