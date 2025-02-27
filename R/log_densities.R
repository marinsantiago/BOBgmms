#' Log-posterior Density
#'
#' This function evaluates the log-posterior density of a Gaussian mixture model
#' with conjugate priors.
#' 
#' @param means A list of size \eqn{K} containing the mean vectors at which we 
#'   want to evaluate the log-posterior density. Each vector of should be of
#'   size \eqn{d}.
#' @param covs A list of of size \eqn{K} containing the covariance matrices at
#'   which we want to evaluate the log-posterior density. Each matrix should be 
#'   of size \eqn{d}-by-\eqn{d}.
#' @param probs A vector of size \eqn{d} of cluster probabilities at which we 
#'   want to evaluate the log-posterior density. The sum of all the 
#'   probabilities should equals one.
#' @param y A matrix of observations of dimension \eqn{n}-by-\eqn{d}, where each 
#'   of the \eqn{n} rows is an observation vector.
#' @param betas A list of size \eqn{K} containing the prior means of each mean
#'  vector. Each prior mean should be a vector of size \eqn{d}.
#' @param lambdas A vector of size \eqn{K} containing the \eqn{\lambda}
#'   shrinkage hyper-parameters associated with each mean vector. Each element 
#'   in the vector must be positive.
#' @param nus A vector of size \eqn{K} containing the \eqn{\nu} hyper-parameters 
#'   associated with each covariance matrix. Each element in the vector must be 
#'   larger than \eqn{d + 1}.
#' @param psis A list of size \eqn{K} containing the prior scale matrices of 
#'   each covariance matrix. Each element in the list must be a symmetric
#'   positive semidefinite matrix of size \eqn{d}-by-\eqn{d}.
#' @param alphas A vector of size \eqn{K} of Dirichlet prior hyper-parameters.
#'   All the elements in the vector must be larger than one.
#'
#' @return The evaluated log-posterior density as a scalar (numeric).
#' 
#' @author Santiago Marin
#' 
lpost <- function(means, covs, probs, y, 
                  betas, lambdas, nus, psis, alphas) {
  
  # Input validation -----------------------------------------------------------
  if (!is.numeric(y) || !is.matrix(y)) stop("y must be a numeric matrix.")
  dims_y <- dim(y)
  n <- dims_y[1]
  p <- dims_y[2]
  is.correct.init(means, covs, probs, p)
  K <- length(means)
  is.correct.hyperpriorparam(betas, lambdas, nus, psis, alphas, K, p)
  gc() # Collect garbage from input validation
  
  # log-posterior density ------------------------------------------------------
  theta <- list(means, covs, probs) |> unlist()
  log_posterior(
    theta = theta, y = y, betas = betas, lambdas = lambdas, 
    nus = nus, psis = psis, alphas = alphas
  )
}


#' Log-likelihood Function
#'
#' This function evaluates the log-likelihood of a Gaussian mixture model.
#' 
#' @param means A list of size \eqn{K} containing the mean vectors at which we 
#'   want to evaluate the log-posterior density. Each vector of should be of
#'   size \eqn{d}.
#' @param covs A list of of size \eqn{K} containing the covariance matrices at
#'   which we want to evaluate the log-posterior density. Each matrix should be 
#'   of size \eqn{d}-by-\eqn{d}.
#' @param probs A vector of size \eqn{d} of cluster probabilities at which we 
#'   want to evaluate the log-posterior density. The sum of all the 
#'   probabilities should equals one.
#' @param y A matrix of observations of dimension \eqn{n}-by-\eqn{d}, where each 
#'   of the \eqn{n} rows is an observation vector.
#'
#' @return The evaluated log-likelihood as a scalar (numeric).
#' 
#' @author Santiago Marin
#' 
loglik <- function(means, covs, probs, y) {
  # Input validation -----------------------------------------------------------
  if (!is.numeric(y) || !is.matrix(y)) stop("y must be a numeric matrix.")
  dims_y <- dim(y)
  n <- dims_y[1]
  p <- dims_y[2]
  is.correct.init(means, covs, probs, p)
  K <- length(means)
  gc() # Collect garbage from input validation
  
  # log-likelihood -------------------------------------------------------------
  theta <- list(means, covs, probs) |> unlist()
  log_likelihood(theta = theta, y = y, p = p, K = K)
}


# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------


extract_covs <- function(theta, p, K){
  start_index <- p * K
  out <- vector("list", length = K)
  for (k in seq_len(K)) {
    start_k <- start_index + ((k - 1) * p^(2)) + 1
    end_k <- start_index + (k * p^(2))
    out[[k]] <- matrix(theta[start_k:end_k], nrow = p, ncol = p)
  }
  out
}


log_likelihood <- function(theta, y, p, K){
  # Extract model parameters
  means <- matrix(theta[1:(p * K)], ncol = p, byrow = TRUE)
  covs  <- extract_covs(theta, p, K)
  probs <- utils::tail(theta, K)
  sum(mvnfast::dmixn(y, means, covs, w = probs, log = TRUE))
}


log_posterior <- function(theta, y, betas, lambdas, nus, psis, alphas) {
  # Prepare the returns
  p <- ncol(y)
  K <- length(alphas)
  # Extract model parameters
  means <- matrix(theta[1:(p * K)], ncol = p, byrow = TRUE)
  covs <- extract_covs(theta, p, K)
  probs <- utils::tail(theta, K)
  # Log-likelihood
  log.likelihood <- sum(mvnfast::dmixn(y, means, covs, w = probs, log = TRUE))
  # Log-prior
  log.prior <- sapply(
    seq_len(K), \(k) {
      mvnfast::dmvn(means[k,], betas[[k]], covs[[k]] / lambdas[[k]], log = T) +
        ldinvwishart(covs[[k]], psis[[k]], nus[k])
    }
  ) |> sum() + lddirichlet(probs, alphas)
  log.likelihood + log.prior
}
