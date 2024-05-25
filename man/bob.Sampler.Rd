\name{bob.Sampler}
\alias{bob.Sampler}
\title{
Sample from the posterior distribution of GMMs using BOB
}
\description{
Sample from the joint posterior distribution of Gaussian Mixture Models using
the Bayesian Optimized Bootstrap (BOB) as in Marin et al. (2024).
}
\usage{
bob.Sampler(y, means, covs, probs, S, Sbatch, a, b, c, r,
            alphas, betas, lambdas, nus, psis, lowerBound, upperBound,
            coresBO, coresSamp, BOiterations = 120, seed = 1)
}
\arguments{
  \item{y}{Matrix of size \code{n*d} of observations, where each column is a
  feature.}
  \item{means}{List of size \code{K} containing initial mean vectors, each of size
  \code{d}.}
  \item{covs}{List of size \code{K} containing initial covariance matrices, each of
  size \code{d*d}.}
  \item{probs}{Vector of size \code{K} of initial probabilities.}
  \item{S}{Total number of posterior draws}
  \item{Sbatch}{Size of the mini-batch used to approximate the KL divergence.}
  \item{a}{Value of the tempering hyper-parameter \code{a}. It MUST be
  larger or equal than zero and smaller than one.}
  \item{b}{Value of the tempering hyper-parameter \code{b}. It MUST be real.}
  \item{c}{Value of the tempering hyper-parameter \code{c}. It MUST be larger
  than zero.}
  \item{r}{Value of the tempering hyper-parameter \code{r}. It MUST be larger
  than zero.}
  \item{alphas}{Vector of size \code{K} of Dirichlet prior hyper-parameters.
  The elements of the vector MUST be larger than one.}
  \item{betas}{List of size \code{K} containing vectors of size \code{d} of
  normal mean hyper-parameters.}
  \item{lambdas}{Vector of size \code{K} containing the lambda regularization
  parameters for each component.}
  \item{nus}{Vector of size \code{K} containing the nu regularization
  parameters for each component.}
  \item{psis}{List of size \code{K} containing matrices of size \code{d*d} of
  Inverse-wishart prior scale matrices. Each matrix MUST be symmetric positive
  semidefinite.}
  \item{lowerBound}{Vector of size \code{2*(K+1)} containing the lower bounds
  of the Bayesian Optimization problem.}
  \item{upperBound}{Vector of size \code{2*(K+1)} containing the upper bounds
  of the Bayesian Optimization problem.}
  \item{coresBO}{The number of cores to use when sampling the mini-batch of
  posterior draws to approximate the KL divergence.}
  \item{coresSamp}{The number of cores to use when sampling the total posterior
  draws after performing the Bayesian Optimization.}
  \item{BOiterations}{The number of iterations employed in the Bayesian
  Optimizer. Default is 120.}
  \item{seed}{The seed for random number generation as an integer. Default is 1.}
}
\details{
Sample from the joint posterior distribution of Gaussian Mixture Models using
the Bayesian Optimized Bootstrap (BOB) as in Marin et al. (2024). IMPORTANT:
The function builds on \code{mclapply} (i.e., a parallelized version of
\code{lapply}), which relies on forking and hence is not available on Windows
unless \code{cores = 1}.
}
\value{
\code{bob.Sampler} returns a list of size 2. The first element is a matrix of size
\code{S*theta.dim}, where \code{theta.dim} is the dimension of the parameter space.
The first \code{K*d} columns are posterior draws from the mean parameters, the next
\code{K*(d^2)} columns are posterior draws from variance and covariance parameters,
and the last \code{K} columns are posterior draws from the probability parameters.
The second elements in the list is a vector containg the optimal \code{x} values
used in the random weighting.
}
\references{
 Marin, S., Loong, B., & Westveld, A. (2024). BOB: Bayesian Optimized Bootstrap
  for Uncertainty Quantification in Gaussian Mixture Models. arXiv preprint arXiv:2311.03644.
}
\author{
Santiago Marin
}
\examples{
\dontrun{
library(mclust)

set.seed(1)

# Data
y <- mclust::banknote[,-1]
K <- 2
d <- ncol(y)
n <- nrow(y)

# Prior hyper-parameters
alphas <- rep(1.1, K)
betas <- lapply(1:K, function(k) rep(0, d))
lambdas <- rep(0.2, K)
nus <- rep(d + 1 + 0.2, K)
psis <- lapply(1:K, function(k) diag(d))

# Initial values
clust.init <- stats::kmeans(y, 2, nstart = 10)$cluster
means <- lapply(1:K, function(k) colMeans(y[clust.init == k,]))
covs <- lapply(1:K, function(k) var(y[clust.init == k,]))
probs <- sapply(1:K, function(k) sum(clust.init == k)/n)

# BO Bounds
lowerBound <- c(1.0, rep(1e-5, K), rep(1e-5, K), 1e-5)
upperBound <- c(1.5, rep(1.5, K),  rep(1.5,  K), 1.5)

bob.Sampler(y, means, covs, probs, S = 1000, Sbatch = 100,
            a = 0, b = 0, c = 1, r = 1, alphas, betas, lambdas, nus, psis,
            lowerBound, upperBound, coresBO = 1L, coresSamp = 4L)
}
}
