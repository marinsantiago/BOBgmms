\name{log.densities}
\alias{log_post}
\alias{log_lik}
\title{
Log-likelihood and log-posterior densities
}
\description{
Computes the log-likelihood and log-posterior densities from a
Gaussian Mixture Model (GMM) with conjugate priors.
}
\usage{
log_post(theta, y, alphas, betas, lambdas, psis, nus)
log_lik(theta, y, d, K)
}
\arguments{
  \item{theta}{Vector of size \code{theta.dim}, where
\code{theta.dim} is the dimension of the parameter space. The first \code{K*d}
elements are mean parameters, the next \code{K*(d^2)} elements are variance and
covariance parameters, and the last \code{K} elements are probability parameters.}
  \item{y}{Matrix of size \code{n*d} of observations, where each column is a
  feature.}
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
  \item{d}{Dimension of the Gaussian mixture.}
  \item{K}{Number of clusters in the Gaussian mixture.}
}
\details{
Computes the log-likelihood and log-posterior densities from a
Gaussian Mixture Model (GMM) with conjugate priors.
}
\value{
\code{log_post} and \code{log_lik} return the log-likelihood and
log-posterior densities as a scalar.
}
\references{
 Marin, S., Loong, B., & Westveld, A. (2024). BOB: Bayesian Optimized Bootstrap
  for Uncertainty Quantification in Gaussian Mixture Models. arXiv preprint arXiv:2311.03644.
}
\author{
Santiago Marin
}
