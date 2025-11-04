#' Approximately sample from the posterior distribution of a GMM via BOB
#'
#' This function generates approximate draws for the posterior distribution of a
#' Gaussian mixture model under conjugate prior using the Bayesian optimized 
#' bootstrap (BOB, Marin et al. 2026).
#' 
#' \emph{Note}: The function leverages parallel computing across multiple CPU 
#' cores through the package \code{parallel}. Please be aware that 
#' parallelization over multiple CPU workers is conducted via \emph{forking} 
#' rather than \emph{sockets}, so it is only available on POSIX ("unix") 
#' systems. On non-POSIX platforms, the function is still operational, but the 
#' number of CPU workers will be automatically set to one. To verify your 
#' operating system (OS), simply run the following R code: 
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
#' @param lower_bound A vector of size \eqn{3K + 1} containing the lower bounds 
#'   of the search space for the Bayesian optimization procedure. All elements 
#'   in the vector should be positive.
#' @param upper_bound A vector of size \eqn{3K + 1} containing the upper bounds 
#'   of the search space for the Bayesian optimization procedure. All elements 
#'   in the vector should be positive. 
#' @param max.iters A positive integer, corresponding to the total number of 
#'   approximate posterior draws. Default is 20000. 
#' @param size.batch A positive integer, corresponding to the size of the batch
#'   of approximate posterior draws used to estimate the KL divergence between
#'   the approximate posterior and the true Bayesian posterior. Default is 1000.
#' @param bo.iters A positive integer, corresponding to the total number of 
#'   iterations in the Bayesian optimization procedure. Default is 100.
#' @param a A scalar corresponding to the \eqn{a} hyper-parameter in the 
#'   tempering profile. Constraint: \eqn{a\in [0,1)}. Default is 0.
#' @param b A scalar corresponding to the \eqn{b} hyper-parameter in the 
#'   tempering profile. Constraint: \eqn{b\in \mathbb{R}}. Default is 0.
#' @param c A scalar corresponding to the \eqn{c} hyper-parameter in the 
#'   tempering profile. Constraint: \eqn{c>0}. Default is 1. 
#' @param r A scalar corresponding to the \eqn{r} hyper-parameter in the 
#'   tempering profile. Constraint: \eqn{r>0}. Default is 1.
#' @param cores The number of CPU cores to use during the sampling process.
#'   Default is \code{parallel::detectCores() - 1}.
#' @param seed The seed for random number generation. Default is
#'   \code{sample.int(.Machine$integer.max, 1)}   
#'
#' @return A list of class "bob" containing:
#' \itemize{
#'   \item \code{post.draws}: A matrix of size 
#'   \eqn{S}-by-\eqn{\text{dim}(\boldsymbol{\theta})} of approximate posterior 
#'   draws, where \eqn{S} is the total number of posterior draws and 
#'   \eqn{\text{dim}(\boldsymbol{\theta})} is the total number of parameters in 
#'   the model. The first \eqn{K \times d} columns in the matrix are posterior 
#'   draws from the mean parameters, the next \eqn{K\times d^2} columns are 
#'   posterior draws from variance and covariance parameters, and the last 
#'   \eqn{K} columns are posterior draws from the mixture proportions.
#'   \item \code{x.optim}: The optimal \eqn{\boldsymbol{x}^{*}} values used in 
#'   the random weighting. See Marin et al. (2026) for additional details.
#' }
#' 
#' @references
#'
#' S. Marin, B. Long,and A. H. Westveld (2026), BOB: Bayesian optimized 
#' bootstrap for approximate posterior sampling in Gaussian mixture models. 
#' \emph{Statistics and Computing}, 36, 14.
#' 
#' @author Santiago Marin
#'
bob.gmm <- function(y, means.init, covs.init, probs.init, betas, lambdas, nus, 
                    psis, alphas,lower_bound, upper_bound, 
                    max.iters = 20000, size.batch = 1000, bo.iters = 100,
                    a = 0, b = 0, c = 1, r = 1,
                    cores = parallel::detectCores() - 1,
                    seed = sample.int(.Machine$integer.max, 1)) {
  
  # Input validation -----------------------------------------------------------
  if (!is.numeric(y) || !is.matrix(y)) stop("y must be a numeric matrix")
  dims_y <- dim(y)
  n <- dims_y[1]
  p <- dims_y[2]
  is.correct.init(means.init, covs.init, probs.init, p)
  K <- length(means.init)
  is.correct.hyperpriorparam(betas, lambdas, nus, psis, alphas, K, p)
  if (!is.numeric(lower_bound)) stop("'lower_bound' should be numeric")
  if (min(lower_bound) <= 0) stop("'lower_bound' should be positive")
  if (!is.numeric(upper_bound)) stop("'upper_bound' should be numeric")
  if (min(upper_bound) <= 0) stop("'upper_bound' should be positive")
  if (!is.numeric(a) || a < 0 || a > 1) stop("Incorrect value for a")
  if (!is.numeric(b)) stop("Incorrect value for b")
  if (!is.numeric(c) || c <= 0) stop("Incorrect value for c")
  if (!is.numeric(r) || r <= 0) stop("Incorrect value for r")
  if (cores %% 1 != 0 || cores <= 0) stop("Incorrect value for cores")
  if (seed %% 1 != 0) stop("'seed' should be an integer")
  nb <- c("max.iters", "size.batch", "bo.iters")
  for (i in nb) {
    if (get(i) %% 1 != 0 || get(i) <= 0) {
      stop(paste(i, " should be a positive integer"))
    }
  }
  rm(i, nb); gc() # Collect cache from input validation
  
  # Pre-compute constants ------------------------------------------------------
  K <- length(alphas)
  p <- ncol(y)
  main.idxs <- get_indices(p, K, covs.idx = FALSE) # Means, variances and probs.
  covs.idxs <- get_indices(p, K, covs.idx = TRUE)  # Covariances
  
  # Set number of parallel workers ---------------------------------------------
  if (.Platform$OS.type == "unix") {
    # Cores to evaluate the KDEs
    cores_KDEs <- if (p > 10 && K > 2) min(cores, 5L) else 1L
  } else {
    cores <- cores_KDEs <- 1L
  }

  # Define the loss function L.hat(x) ------------------------------------------
  objective <- \(x) {
    loss <- tryCatch(
      {
        KL_divergence(
          x_cov = x, 
          y = y, means.init = means.init, covs.init = covs.init,
          probs.init = probs.init, main.idxs = main.idxs, covs.idxs = covs.idxs,
          betas = betas, lambdas = lambdas, nus = nus, psis = psis, 
          alphas = alphas, batch.iters = size.batch,   
          a = a, b = b, c = c, r = r, cores.BO = cores,
          cores.kdes = cores_KDEs, seed = seed
        )
      }, error = \(e) 999
    )
    loss <- if (is.finite(loss)) min(loss, 999) else 999
    loss
  }
  
  # Bayesian optimization ------------------------------------------------------
  set.seed(seed)
  x.star <- bayes_optim(
    lower_bound = lower_bound, upper_bound = upper_bound,
    y = y, K = K, objective = objective,
    tot_iterations = bo.iters, init_samples = 30
  ) 
  rm(main.idxs, covs.idxs, cores_KDEs)
  gc() # Collect garbage 
  
  # Extract optimal x values ---------------------------------------------------
  x.star_covs <- c(x.star)

  # Posterior sampler ----------------------------------------------------------
  bob.out <- bob.batch(
    x_cov = x.star_covs,
    y = y, means.init = means.init, covs.init = covs.init,
    probs.init = probs.init, betas = betas, lambdas = lambdas, nus = nus,
    psis = psis, alphas = alphas, batch.iters = max.iters,
    a = a, b = b, c = c, r = r, cores = cores, seed = seed
  )
  
  out <- list("post.draws" = bob.out, "x.optim" = x.star)
  class(out) <- "bob"
  out
}

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

bob.batch <- function(x_cov,
                      y, means.init, covs.init, probs.init,
                      betas, lambdas, nus, psis, alphas,
                      batch.iters = 1000, a = 0, b = 0, c = 1, r = 1,
                      cores = parallel::detectCores() - 1,
                      seed = sample.int(.Machine$integer.max, 1)) {
  
  # Pre-compute constants ------------------------------------------------------
  dims_y <- dim(y)
  n <- dims_y[1]
  p <- dims_y[2]
  K <- length(means.init)
  
  # Sampler --------------------------------------------------------------------
  base::RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  parallel::mc.reset.stream()
  # Iterate
  bob_batch_out <- parallel::mclapply(
    seq_len(batch.iters), \(iter) {
      # Random data weights
      lik.weights <- rexp(n, rate = 1)
      # Random prior weights
      means.weights <- rexp(K, rate = 1)
      covs.weights <- rep(x_cov, K)
      probs.weights <- rexp(K, rate = 1)
      # EM optimization
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
  ) |> unlist(x = _) |> matrix(data = _, nrow = batch.iters, byrow = TRUE)
  # Re-set random generator
  base::RNGkind("default", "default", "default")
  bob_batch_out
}


get_indices <- function(p, K, covs.idx = FALSE) {
  pK <- p * K
  seq_covs <- seq_len(pK * p) + pK
  out <- if (!covs.idx) {
    c(
      seq_len(pK),                  # means indices
      seq_covs[c(diag(p) == 1)],    # variances indices
      seq_len(K) + pK * (p + 1)     # probs. indices
    )
  } else seq_covs[c(upper.tri(matrix(NA, p, p)))] # covs. indices
  out
}


eval_kdes <- function(post.draws, main.idxs, covs.idxs,
                      main.cores = 1, covs.cores = 1) {
  
  # Means, variances and probabilities -----------------------------------------
  kdes.1 <- parallel::mclapply(
    main.idxs, \(j) {
      FKSUM::fk_density(x = post.draws[,j], x_eval = post.draws[,j])$y |> log()
    },
    mc.cores = main.cores
  ) |> unlist(x = _) |> sum()
  
  # Covariances ----------------------------------------------------------------
  kdes.2 <- parallel::mclapply(
    covs.idxs, \(j) {
      FKSUM::fk_density(x = post.draws[,j], x_eval = post.draws[,j])$y |> log()
    },
    mc.cores = covs.cores
  ) |> unlist(x = _) |> sum()
  
  # Multiply kdes.2 by two because we are only counting the upper diag. covs. 
  kdes.1 + 2 * kdes.2
}


KL_divergence <- function(x_cov, 
                          y, means.init, covs.init, probs.init,
                          main.idxs, covs.idxs,
                          betas, lambdas, nus, psis, alphas, 
                          batch.iters = 1000, a = 0, b = 0, c = 1, r = 1,
                          cores.BO = parallel::detectCores() - 1,
                          cores.kdes = 1L,
                          seed = sample.int(.Machine$integer.max, 1)) {
  
  # Prepare the returns
  K <- length(means.init)
  p <- ncol(y)
  
  # Batch of approximate posterior draws given the current x value -------------
  approx_draws <- bob.batch(
    x_cov = x_cov, y = y, means.init = means.init, covs.init = covs.init,
    probs.init = probs.init, betas = betas, lambdas = lambdas, nus = nus,
    psis = psis, alphas = alphas, batch.iters = batch.iters, 
    a = a, b = b, c = c, r = r, cores = cores.BO, seed = seed
  )
  #cat("Batch sampling: Done \n")
  
  # Evaluate KDEs (in parallel) ------------------------------------------------
  kdes <- eval_kdes(
    post.draws = approx_draws, main.idxs = main.idxs, covs.idxs = covs.idxs, 
    main.cores = cores.kdes, covs.cores = cores.kdes
  )
  #cat("KDEs evaluation: Done \n")
  
  # Evaluate target posterior (in parallel) ------------------------------------
  lpost <- parallel::mclapply(
    seq_len(batch.iters), \(iter) {
      log_posterior(
        theta = approx_draws[iter,], y = y, betas = betas, 
        lambdas = lambdas, nus = nus, psis = psis, alphas = alphas
      )
    },
    mc.cores = cores.BO
  ) |> unlist(x = _) |> sum()
  #cat("Log-posterior evaluation: Done \n")
  
  # Empirical KL divergence ----------------------------------------------------
  theta.dim <- K * (p^(2) + p + 1)
  (kdes - lpost) / (batch.iters * theta.dim)
}
