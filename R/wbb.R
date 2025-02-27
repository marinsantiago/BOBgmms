#' Sample from the posterior distribution of a GMM via WBB
#'
#' This function generates approximate posterior draws from Gaussian mixture
#' models using the Weighted Bayesian Bootstrap (WBB, Newton et al. 2021). The 
#' function gives the user two option of random prior weights: (i) Random prior
#' weights or (ii) fixed prior weights. Additionally, the function leverages 
#' parallel computing across multiple CPU cores with the package 
#' \code{parallel}. Please be aware that the parallelization is conducted via 
#' \emph{forking} rather than \emph{sockets}, so it only works on "unix" 
#' systems. One can verify the OS type by running the following R code:
#' \code{.Platform$OS.type}.
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
#' @param wbb.scheme A character string denoting which weighting scheme to use.
#'   The options are either (i) \code{"wbb1"} for random prior weights or (ii) 
#'   \code{"wbb2"} for fixed prior weights. Default is \code{"wbb1"}.
#' @param max.iters A positive integer, corresponding to the total number of 
#'   approximate posterior draws. Default is 20000. 
#' @param a A scalar corresponding to the \eqn{a} hyper-parameter in the 
#'   tempering profile. Note that \eqn{a\in [0,1)}. Default is 0.
#' @param b A scalar corresponding to the \eqn{b} hyper-parameter in the 
#'   tempering profile. Note that \eqn{b\in \mathbb{R}}. Default is 0.
#' @param c A scalar corresponding to the \eqn{c} hyper-parameter in the 
#'   tempering profile. Note that \eqn{c>0}. Default is 1.
#' @param r A scalar corresponding to the \eqn{r} hyper-parameter in the 
#'   tempering profile. Note that \eqn{r>0}. Default is 1.
#' @param cores The number of CPU cores to use during the sampling process.
#'   Default is \code{parallel::detectCores() - 1}.
#' @param seed The seed for random number generation. Default is
#'   \code{sample.int(.Machine$integer.max, 1)}   
#'
#' @return A matrix of size \eqn{S}-by-\eqn{\text{dim}(\boldsymbol{\theta})}
#'   of approximate posterior draws, where \eqn{S} is the total number of
#'   posterior draws and \eqn{\text{dim}(\boldsymbol{\theta})} is the total
#'   number of parameters in the model. The first \eqn{K \times d} columns are 
#'   posterior draws from the mean parameters, the next \eqn{K\times d^2} 
#'   columns are posterior draws from variance and covariance parameters, and 
#'   the last \eqn{K} columns are posterior draws from the probability 
#'   parameters.
#' 
#' @references
#' 
#' M.A. Newton, N.G. Polson, and J. Xu (2021), Weighted Bayesian bootstrap for 
#' scalable posterior distributions. \emph{Can J Statistics}, 49:421-437.
#' 
#' @author Santiago Marin
#' 
wbb.gmm <- function(y, means.init, covs.init, probs.init,
                    betas, lambdas, nus, psis, alphas,
                    wbb.scheme = "wbb1", max.iters = 20000,
                    a = 0, b = 0, c = 1, r = 1,
                    cores = parallel::detectCores() - 1,
                    seed = sample.int(.Machine$integer.max, 1)) {
  
  # Input validation -----------------------------------------------------------
  if (!is.numeric(y) || !is.matrix(y)) stop("y must be a numeric matrix.")
  dims_y <- dim(y)
  n <- dims_y[1]
  p <- dims_y[2]
  is.correct.init(means.init, covs.init, probs.init, p)
  K <- length(means.init)
  is.correct.hyperpriorparam(betas, lambdas, nus, psis, alphas, K, p)
  if (!(wbb.scheme %in% c("wbb1", "wbb2"))) {
    stop("Incorrect WBB scheme. It should be either 'wbb1' or 'wbb2'.")
  }
  if (max.iters %% 1 != 0 || max.iters <= 0) {
    stop("Incorrect value for max.iters. It should be a positive integer.")
  }
  if (!is.numeric(a) || a < 0 || a > 1) stop("Incorrect value for a.")
  if (!is.numeric(b)) stop("Incorrect value for b.")
  if (!is.numeric(c) || c <= 0) stop("Incorrect value for c.")
  if (!is.numeric(r) || r <= 0) stop("Incorrect value for r.")
  if (cores %% 1 != 0 || cores <= 0) stop("Incorrect value for cores.")
  if (seed %% 1 != 0) stop("Supplied 'seed' is not an integer.")
  gc() # Collect garbage from input validation
  
  # Pre-compute constants ------------------------------------------------------
  dims_y <- dim(y)
  n <- dims_y[1]
  p <- dims_y[2]
  K <- length(means.init)
  
  cat("Model: Gaussian mixture with conjugate priors", "\n")
  cat("Approximate posterior sampling via Weighted Bayesian Bootstrap:", "\n")
  # Set random number generator and seed
  base::RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  parallel::mc.reset.stream()
  # Iterate
  wbb_out <- pbmcapply::pbmclapply(
    seq_len(max.iters), \(iter) {
      # Random weights  --------------------------------------------------------
      if (wbb.scheme == "wbb1") {
        # Random prior weights
        lik.weights <- rexp(n, rate = 1)
        means.weights <- rexp(K, rate = 1)
        covs.weights <- rexp(K, rate = 1)
        probs.weights <- rexp(K, rate = 1)
      } else if (wbb.scheme == "wbb2") {
        # Fixed prior weights
        lik.weights <- rexp(n, rate = 1)
        means.weights <- rep(1, K)
        covs.weights <- rep(1, K)
        probs.weights <- rep(1, K)
      }
      # EM optimization --------------------------------------------------------
      em_optim(
        y = y, means.init = means.init, covs.init = covs.init,
        probs.init = probs.init, alphas = alphas, betas = betas,
        lambdas = lambdas, nus = nus, psis = psis, lik.weights = lik.weights,
        means.weights = means.weights, covs.weights = covs.weights,
        probs.weights = probs.weights, a = a, b = b, c = c, r = r, 
        max.iters = 500, epss = 1e-4
      ) |> unlist(x = _)
    },
    mc.set.seed = TRUE,
    mc.cores = cores
  ) |> unlist(x = _) |> matrix(data = _, nrow = max.iters, byrow = TRUE)
  # Re-set random generator
  base::RNGkind("default", "default", "default")
  
  class(wbb_out) <- "wbb"
  wbb_out
}
