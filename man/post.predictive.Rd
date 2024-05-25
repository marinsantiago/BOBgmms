\name{post.predictive}
\alias{post.predictive}
\title{
Posterior predictive distrubtion of Gaussian Mixture Models
}
\description{
Computes the posterior predictive distrubtion of Gaussian Mixture Models.
}
\usage{
post.predictive(post.draws, d, K)
}
\arguments{
  \item{post.draws}{Matrix of size \code{S*theta.dim}, where
\code{theta.dim} is the dimension of the parameter space. The first \code{K*d}
columns are posterior draws from the mean parameters, the next \code{K*(d^2)}
columns are posterior draws from variance and covariance parameters, and the
last \code{K} columns are posterior draws from the probability parameters.}
  \item{d}{Dimension of the Gaussian Mixture.}
  \item{K}{Number of components in the Gaussian Mixture.}
}
\details{Given a matrix of posterior draws, this function computes the resulting
posterior predictive distrubtion of a Gaussian Mixture Model.}
\value{
\code{bob.Sampler} returns a matrix of size \code{S*d} of posterior predictive
draws, where each column corresponds to a feature.}
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

draws <- wbb.Sampler(y, means, covs, probs, S = 1000,
                     a = 0, b = 0, c = 1, r = 1, alphas, betas, lambdas, nus,
                     psis, randomPrior = T, cores = 4L)

post.predictive(draws, 4, 3)
}
}
