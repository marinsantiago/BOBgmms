#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include "e_step.h"


Eigen::VectorXd e_step_numerator(const Eigen::MatrixXd & y,
                                 const Eigen::VectorXd & mu,
                                 const Eigen::MatrixXd & Sigma,
                                 const double prob,
                                 const Eigen::VectorXd & lik_weights) {
  
  // Prepare the returns
  const int n = y.rows();
  const int p = y.cols();
  Eigen::VectorXd out(n);
  
  // Pre-compute constants
  const double lnSqrt2Pi = 0.5 * std::log(2 * M_PI);
  const double log_prob = std::log(prob);
  // Cholesky decomposition
  typedef Eigen::LLT<Eigen::MatrixXd> Chol;
  Chol chol(Sigma);
  const Chol::Traits::MatrixL & L = chol.matrixL();
  // Pre-compute log-determinant
  double log_det_L = L.toDenseMatrix().diagonal().array().log().sum();
 
  // Iterate over observations
  for (int i = 0; i < n; ++i) {
    // Compute quadratic form for observation i
    double quadform = (L.solve(y.row(i).transpose() - mu)).squaredNorm();
    // Evaluate log-density for observation i
    double log_dmvn = -p * lnSqrt2Pi - log_det_L - 0.5 * quadform;
    // E-step numerator for observation i
    out[i] = std::exp(log_dmvn * lik_weights[i] + log_prob * lik_weights[i]);
  }
  
  return out;
}

// [[Rcpp::export]]
Eigen::MatrixXd e_step(const Eigen::MatrixXd & y,
                       const std::vector<Eigen::VectorXd> & means,
                       const std::vector<Eigen::MatrixXd> & covs,
                       const Eigen::VectorXd probs,
                       const Eigen::VectorXd & lik_weights,
                       const double temp) {
  
  const int n = y.rows();
  const int K = means.size();
  Eigen::MatrixXd out(n, K);
  
  // Un-normalized e-step
  for (int k = 0; k < K; ++k) {
    out.col(k) = e_step_numerator(y, means[k], covs[k], probs[k], lik_weights);
  }
  // Iterate over observations
  for (int i = 0; i < n; ++i) {
    // Apply tempering
    if (std::abs(temp - 1.0) > 1e-10) {
      out.row(i) = out.row(i).array().pow(1/temp);
    }
    // Normalize
    double row_sum = out.row(i).sum();
    if (row_sum <= 0) {
      row_sum = 1e-15;
    }
    out.row(i) /= row_sum;
  }
  
  return out;
}
