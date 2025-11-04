#' Log-posterior optimization via expectation-maximization (EM).
#'
#' This function maximizes the weighted log-posterior density associated with a 
#' Gaussian mixture model with the EM algorithm as described in Marin et 
#' al. (2026).
#'
#' @param y A matrix of observations of dimension \eqn{n}-by-\eqn{d}, where each 
#'   of the \eqn{n} rows is an observation vector.
#' @param means.init A list of size \eqn{K} containing the initial mean 
#'   vectors, each vector of size \eqn{d}.
#' @param covs.init A list of of size \eqn{K} containing the initial covariance 
#'   matrices, each matrix of size \eqn{d}-by-\eqn{d}.
#' @param probs.init A vector of size \eqn{d} of initial probabilities. The sum
#'   of all the probabilities should equals one.
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
#' @param lik.weights A vector of size \eqn{n} of likelihood weights. All 
#'   weights must be non-negative. Default is 1 for each weight, i.e., no 
#'   weighting.
#' @param means.weights A vector of size \eqn{K} of weights associated with each
#'   mean vector. All weights must be non-negative. Default is 1 for each 
#'   weight, i.e., no weighting.
#' @param covs.weights A vector of size \eqn{K} of weights associated with each
#'   covariance matrix. All weights must be non-negative. Default is 1 for each 
#'   weight, i.e., no weighting. 
#' @param probs.weights A vector of size \eqn{K} of weights associated with each
#'   probability. All weights must be non-negative. Default is 1 for each 
#'   weight, i.e., no weighting.  
#' @param a A scalar corresponding to the \eqn{a} hyper-parameter in the 
#'   tempering profile. Constraint: \eqn{a\in [0,1)}. Default is 0.
#' @param b A scalar corresponding to the \eqn{b} hyper-parameter in the 
#'   tempering profile. Constraint: \eqn{b\in \mathbb{R}}. Default is 0.
#' @param c A scalar corresponding to the \eqn{c} hyper-parameter in the 
#'   tempering profile. Constraint: \eqn{c>0}. Default is 1.
#' @param r A scalar corresponding to the \eqn{r} hyper-parameter in the 
#'   tempering profile. Constraint: \eqn{r>0}. Default is 1.
#' @param max.iters Maximum number of EM iterations. Default is 500.
#' @param epss Convergence threshold for the EM algorithm. Default is 
#'   \code{1e-4}
#'
#' @return A list containing the EM solution for each mean vector, each 
#'   covariance matrix, and each mixture component.
#' 
#' @references
#'
#' S. Marin, B. Long,and A. H. Westveld (2026), BOB: Bayesian optimized 
#' bootstrap for approximate posterior sampling in Gaussian mixture models. 
#' \emph{Statistics and Computing}, 36, 14.
#' 
#' @author Santiago Marin
#' 
em.optim <- function(y, means.init, covs.init, probs.init,
                     betas, lambdas, nus, psis, alphas,
                     lik.weights = rep(1, nrow(y)), 
                     means.weights = rep(1, length(means.init)), 
                     covs.weights = rep(1, length(means.init)), 
                     probs.weights = rep(1, length(means.init)),
                     a = 0, b = 0, c = 1, r = 1,
                     max.iters = 500, epss = 1e-4) {
  
  # Input validation -----------------------------------------------------------
  if (!is.numeric(y) || !is.matrix(y)) stop("y must be a numeric matrix.")
  dims_y <- dim(y)
  n <- dims_y[1]
  p <- dims_y[2]
  is.correct.init(means.init, covs.init, probs.init, p)
  K <- length(means.init)
  is.correct.hyperpriorparam(betas, lambdas, nus, psis, alphas, K, p)
  is.correct.weighting(
    lik.weights, means.weights, covs.weights, probs.weights, n, K
  )
  if (!is.numeric(a) || a < 0 || a > 1) stop("Incorrect value for a.")
  if (!is.numeric(b)) stop("Incorrect value for b.")
  if (!is.numeric(c) || c <= 0) stop("Incorrect value for c.")
  if (!is.numeric(r) || r <= 0) stop("Incorrect value for r.")
  if (!is.numeric(epss) || epss <= 0) stop("Incorrect value for epss.")
  if (max.iters %% 1 != 0 || max.iters <= 0) stop("Incorrect  max.iters.")
  gc() # Collect garbage from input validation
  
  # EM optimization ------------------------------------------------------------
  em_out <- em_optim(
    y = y, means.init = means.init, covs.init = covs.init, 
    probs.init = probs.init, alphas = alphas, betas = betas, lambdas = lambdas,
    nus = nus, psis = psis, lik.weights = lik.weights, 
    means.weights = means.weights, covs.weights = covs.weights, 
    probs.weights = probs.weights, a = a, b = b, c = c, r = r, 
    max.iters = max.iters, epss = epss
  )
  gc()
  names(em_out) <- c("means", "covs", "probs")
  class(em_out) <- "EM"
  em_out
}

# ------------------------------------------------------------------------------
# Helpers 
# ------------------------------------------------------------------------------

em_optim <- function(y, means.init, covs.init, probs.init, alphas, 
                     betas, lambdas, nus, psis, lik.weights, means.weights, 
                     covs.weights, probs.weights, a = 0, b = 0, c = 1, r = 1,
                     max.iters = 500, epss = 1e-4) {
  em_optim_cpp_(
    y = y, means_init = means.init, covs_init = covs.init, 
    probs_init = probs.init, alphas = alphas, betas = betas, lambdas = lambdas,
    nus = nus, Psis = psis, lik_weights = lik.weights, 
    means_weights = means.weights, covs_weights = covs.weights,
    probs_weights = probs.weights, a = a, b = b, c = c, r = r,
    max_iters = max.iters, epss = epss
  )
}
