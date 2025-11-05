#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include "ddirichlet.h"


/**
 * Log-density function of the Dirichlet distribution                          \ 
 *                                                                             \
 * @param x Eigen::VectorXd: Point at which we want to evaluate the            \
 *  log-density function.                                                      \
 * @param alphas Eigen::VectorXd: The vector of shape parameters.              \
 *                                                                             \
 * @return double: Log-density function of the Dirichlet distribution          \
 *  evaluated at x.                                                            \
 */
// [[Rcpp::export]]
double lddirichlet(const Eigen::VectorXd & x,
                   const Eigen::VectorXd & alphas) {
  
  // Pre-compute constants
  double lmvbeta_alphas = lmvbeta(alphas);
  Eigen::VectorXd alphas_minus_one = alphas.array() - 1;
  // log-density
  double out = -lmvbeta_alphas;
  out += (x.array().log() * alphas_minus_one.array()).sum();
  return out;
}


// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------


/**
 * Log-multivariate beta function, which can be expressed in terms of the      \
 *  gamma function. For details see:                                           \
 *  https://en.wikipedia.org/wiki/Dirichlet_distribution                       \ 
 *                                                                             \
 * @param alpha Eigen::VectorXd: Point at which we want to evaluate the        \  
 *  log-multivariate beta function                                             \
 *                                                                             \
 * @return double: Log-multivariate beta function evaluated at alpha.          \
 */
double lmvbeta(const Eigen::VectorXd & alpha) {
  int p = alpha.size();
  double sum_alpha = alpha.sum();
  double log_gamma_sum = 0.0;
  for (int j = 0; j < p; ++j) {
    log_gamma_sum += std::lgamma(alpha[j]);
  }
  double log_gamma_sum_alpha = std::lgamma(sum_alpha);
  double out = log_gamma_sum - log_gamma_sum_alpha;
  return out;
}
