\name{distances}
\alias{distance.tv}
\alias{distance.ks}
\title{
TV and KS distances
}
\description{
Computes the Total Variation (TV) and Kolmogorov-Smirnov (KS) distances between
two posterior predictive distributions.
}
\usage{
distance.tv(target, approx, step = 0.001)
distance.ks(target, approx, step = 0.001)
}
\arguments{
  \item{target}{Matrix of size \code{S*d} of target posterior predictive
draws (usually obtained via MCMC), where each column corresponds to a feature.}
  \item{approx}{Matrix of size \code{S*d} of approximate posterior predictive
draws, where each column corresponds to a feature.}
  \item{step}{Step size used in the trapezoidal rule. Default is 0.001.}
}
\details{
Computes the Total Variation (TV) and Kolmogorov-Smirnov (KS) distances between
two posterior predictive distributions.
}
\value{
\code{distance.tv} and \code{distance.ks} return the distance between
two posterior predictive distributions as a scalar.
}
\references{
 Marin, S., Loong, B., & Westveld, A. (2024). BOB: Bayesian Optimized Bootstrap
  for Uncertainty Quantification in Gaussian Mixture Models. arXiv preprint arXiv:2311.03644.
}
\author{
Santiago Marin
}
