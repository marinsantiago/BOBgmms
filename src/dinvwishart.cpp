#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include "dinvwishart.h"


/**
 * Log-density function of the inverse-Wishart distribution                    \ 
 *                                                                             \
 * @param Sigma Eigen::MatrixXd: Symmetric positive semidefinite matrix at     \
 *  which we want to evaluate the log-density function.                        \
 * @param S Eigen::MatrixXd: Symmetric positive semidefinite scale matrix.     \
 * @param nu double: Degrees of freedom.                                       \
 *                                                                             \
 * @return double: Log-density function of the inverse-Wishart distribution    \
 *  evaluated at x.                                                            \
 */
// [[Rcpp::export]]
double ldinvwishart(const Eigen::MatrixXd & Sigma,
                    const Eigen::MatrixXd & S,
                    const double nu) {
  
  // Pre-compute constants -----------------------------------------------------
  const int p = Sigma.cols();
  const double nu_2 = nu / 2;
  const double log_det_S = log_det(S);
  const double lmvgamma_nu2 = lmvgamma(nu_2, p);
  // Cholesky decomposition of Sigma
  typedef Eigen::LLT<Eigen::MatrixXd> Chol;
  Chol chol(Sigma);
  const Chol::Traits::MatrixL & L = chol.matrixL();
  // log-determinant of Sigma
  double log_det_Sigma = 2.0;
  log_det_Sigma *=  L.toDenseMatrix().diagonal().array().log().sum();
  // S * Sigma^(-1)
  Eigen::MatrixXd Y = L.solve(S);
  Eigen::MatrixXd S_Sigmainv = L.transpose().solve(Y);
  double tr_S_Sigmainv = S_Sigmainv.trace();
  
  // log-density ---------------------------------------------------------------
  double out = (log_det_S - std::log(2) * p) * nu_2 - lmvgamma_nu2;
  out -= (nu + p + 1) * 0.5 * log_det_Sigma;
  out -= tr_S_Sigmainv * 0.5;
  return out;
}


// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------


/**
 * Log-multivariate gamma function. For details see:                           \
 *  https://en.wikipedia.org/wiki/Multivariate_gamma_function                  \ 
 *                                                                             \
 * @param x double: Point at which we want to evaluate the log-multivariate    \
 *  gamma function.                                                            \
 * @param p int: Dimension of multivariate gamma function.                     \
 *                                                                             \
 * @return double: Log-multivariate gamma function evaluated at x.             \
 */
// [[Rcpp::export]]
double lmvgamma(double x, int p) {
  double out = 0.25 * p * (p - 1) * log(M_PI);
  for (int i = 0; i < p; ++i) {
    out += std::lgamma(x - 0.5 * i);
  }
  return out;
}


// Compute the log-determinant of a non-singular matrix
/**
 * Compute the log-determinant of a non-singular matrix                        \
 *                                                                             \
 * @param x Eigen::MatrixXd: Non-singular matrix                               \
 *                                                                             \
 * @return double: Log-determinant of x.                                       \
 */
double log_det(const Eigen::MatrixXd & x) {
  Eigen::LLT<Eigen::MatrixXd> llt(x);
  double log_det = 2.0;
  log_det *= llt.matrixL().toDenseMatrix().diagonal().array().log().sum();
  return log_det;
}
