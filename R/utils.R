
# Check if a matrix is symmetric
is.symmetric.matrix <- function (x) {
  sum(x == t(x)) == (nrow(x)^2)
}


# Check if a matrix is symmetric semidefinite positive
is.positive.semidefinite <- function (x) {
  if (!is.matrix(x)) stop("x is not a matrix.")
  if (nrow(x) != ncol(x)) stop("x is not a square matrix.")
  if (!is.symmetric.matrix(x)) stop("x is not a symmetric matrix.")
  eigs <- eigen(x, symmetric = TRUE)$values
  if (any(is.complex(eigs))) return(FALSE)
  if (all(eigs >= 0)) pd <- TRUE
  else pd <- FALSE
  pd
}


# Validate the user-supplied initial values for the EM algorithm.
is.correct.init <- function(means.init, covs.init, probs.init, p) {
  if (!is.list(means.init)) {
    stop("means.init must be a list of numeric vectors.")
  } 
  K <- length(means.init)
  if (!all(sapply(means.init, \(v) is.vector(v) & is.numeric(v)))) {
    stop("means.init must be a list of numeric vectors.")
  }
  mssg <- "does not match the dimension of the data."
  if (!all(sapply(means.init, length) == p)) {
    stop(paste("The size of the vectors in means.init", mssg))
  }
  if (!is.list(covs.init)) stop("covs.init must be a list of numeric matrices.")
  if (length(covs.init) != K) {
    stop(
      paste(
        "The number of matrices in covs.init does not",
        "match the number of vectors in means.init."
      )
    )
  }
  if (!all(sapply(covs.init, is.positive.semidefinite))) {
    paste("Not all the matrices in covs.init are semidefinite positive.")
  }
  if (!all(sapply(covs.init, ncol) == p)) {
    stop(paste("The size of the matrices in covs.init", mssg))
  }
  if (!is.numeric(probs.init) || !is.vector(probs.init)) {
    stop("probs.init must be a numeric vector.")
  }
  if (length(probs.init) != K) {
    stop("The length of probs.init does not match the number of clusetrs.")
  }
  chck.probs <- all(
    max(probs.init) <= 1, min(probs.init) >= 0, abs(sum(probs.init)-1) <= 1e-10
  )
  if (!chck.probs) stop("probs.init is not a valid vector of probabilities.")
  rm(chck.probs, K)
}


# Validate the user-supplied values for the prior hyper-parameters
is.correct.hyperpriorparam <- function(betas, lambdas, nus, Psis, 
                                       alphas, K, p) {
  
  if (!is.numeric(lambdas) || !is.vector(lambdas)) {
    stop("lambdas must be a numeric vector.")
  }
  if (!(min(lambdas) > 0)) stop("All entries in lambdas must be poitive.")
  if (length(lambdas) != K){
    stop("The length of lambdas does not match the number of clusters.")
  }
  if (!is.numeric(nus) || !is.vector(nus)) stop("nus must be a numeric vector.")
  if (!(min(nus) > (p + 1))) {
    stop("All entries in nus must be larger than p + 1.")
  }
  if (length(nus) != K){
    stop("The length of nus does not match the number of clusters.")
  }
  if (!is.list(betas) || length(betas) != K) {
    stop("betas must be a list of numeric vectors.")
  }
  if (!all(sapply(betas, \(v) is.vector(v) & is.numeric(v)))) {
    stop("betas must be a list of numeric vectors.")
  }
  mssg <- "does not match the dimension of the data."
  if (!all(sapply(betas, length) == p)) {
    stop(paste("The size of the vectors in betas", mssg))
  }
  if (!is.list(Psis) || length(Psis) != K) {
    stop("Psis must be a list of numeric matrices.")
  } 
  if (!all(sapply(Psis, is.positive.semidefinite))) {
    paste("Not all the matrices in Psis are semidefinite positive.")
  }
  if (!all(sapply(Psis, ncol) == p)) {
    stop(paste("The size of the matrices in Psis", mssg))
  }
  if (!is.numeric(alphas) || !is.vector(alphas)) {
    stop("alphas must be a numeric vector.")
  }
  if (!(min(alphas) > 1)) stop("All entries in lambdas must larger than one.")
  if (length(alphas) != K) {
    stop("The length of alphas does not match the number of clusters.")
  }
}


# Validate the user-supplied weighting scheme
is.correct.weighting <- function(lik.weights, means.weights,
                                 covs.weights, probs.weights, n, K) {
  
  if (!is.numeric(lik.weights) || !is.vector(lik.weights)) {
    stop("lik.weights must be a numeric vector.")
  }
  if (!(min(lik.weights) > 0)) stop("All lik.weights must be poitive.")
  if (length(lik.weights) != n){
    stop("The length of lik.weights does not match the sample size.")
  }
  if (!is.numeric(means.weights) || !is.vector(means.weights)) {
    stop("means.weights must be a numeric vector.")
  }
  if (!(min(means.weights) > 0)) stop("All means.weights must be poitive.")
  if (length(means.weights) != K){
    stop("The length of means.weights does not match the number of clusters.")
  }
  if (!is.numeric(covs.weights) || !is.vector(covs.weights)) {
    stop("covs.weights must be a numeric vector.")
  }
  if (!(min(covs.weights) > 0)) stop("All covs.weights must be poitive.")
  if (length(covs.weights) != K){
    stop("The length of covs.weights does not match the number of clusters.")
  }
  if (!is.numeric(probs.weights) || !is.vector(probs.weights)) {
    stop("probs.weights must be a numeric vector.")
  }
  if (!(min(probs.weights) > 0)) stop("All probs.weights must be poitive.")
  if (length(probs.weights) != K){
    stop("The length of probs.weights does not match the number of clusters.")
  }
}
