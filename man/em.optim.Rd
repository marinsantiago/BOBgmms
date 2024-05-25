\name{em.optim}
\alias{em.optim}
\title{
Randomly Weighted EM algorithm
}
\description{
Compute the MAP estimate of a randomly weighted posterior distribution from a
Gaussian Mixture Model via EM.
}
\usage{
em.optim(y, means, covs, probs, a, b, c, r,
         lik_weights, mean_weights, cov_weights, prob_weights,
         alphas, betas, lambdas, nus, psis, max_iters = 500, epss = 1e-4)
}
\arguments{
  \item{y}{Matrix of size \code{n*d} of observations, where each column is a
  feature.}
  \item{means}{List of size \code{K} containing initial mean vectors, each of size
  \code{d}.}
  \item{covs}{List of size \code{K} containing initial covariance matrices, each of
  size \code{d*d}.}
  \item{probs}{Vector of size \code{K} of initial probabilities.}
  \item{a}{Value of the tempering hyper-parameter \code{a}. It MUST be
  larger or equal than zero and smaller than one.}
  \item{b}{Value of the tempering hyper-parameter \code{b}. It MUST be real.}
  \item{c}{Value of the tempering hyper-parameter \code{c}. It MUST be larger
  than zero.}
  \item{r}{Value of the tempering hyper-parameter \code{r}. It MUST be larger
  than zero.}
  \item{lik_weights}{\code{n}-dimensional vector of likelihood weights.}
  \item{mean_weights}{\code{K}-dimensional vector of weights associated with
  the prior distribution of the mean vectors.}
  \item{cov_weights}{\code{K}-dimensional vector of weights associated with
  the prior distribution of the covariance matrices.}
  \item{prob_weights}{Scalar weight associated with the prior distribution of
  the component probabilities.}
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
  \item{max_iters}{Maximum number of EM iterations. Default is 500.}
  \item{epss}{Threshold for convergence. Default is 1e-4.}
}
\details{
Compute the MAP estimate of a randomly weighted posterior distribution from a
Gaussian Mixture Model using an Expectation-Maximization algorithm.
}
\value{
\code{em.optim} returns a vector of size \code{theta.dim}, where
\code{theta.dim} is the dimension of the parameter space. The first \code{K*d}
elements are mean parameters, the next \code{K*(d^2)} elements are variance and
covariance parameters, and the last \code{K} elements are probability parameters.
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

# Weights
lik_weights <- rep(1, n)
mean_weights <- cov_weights <- rep(1, 3)
prob_weights <- 1

em.optim(y, means, covs, probs, a = 0, b = 0, c = 1, r = 1,
         lik_weights, mean_weights, cov_weights, prob_weights,
         alphas, betas, lambdas, nus, psis)
}
}
