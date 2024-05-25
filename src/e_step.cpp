#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include "e_step.h"


Eigen::VectorXd weighted_dens(const Eigen::MatrixXd &y,
                              const Eigen::VectorXd &mu,
                              const Eigen::MatrixXd &Sigma,
                              const double prob,
                              const Eigen::VectorXd &lik_weights) {
  
  int n = y.rows(), p = y.cols();
  Eigen::VectorXd out(n);
  
  const double lnSqrt2Pi = 0.5 * std::log(2 * M_PI);
  
  // Cholesky decomposition
  typedef Eigen::LLT<Eigen::MatrixXd> Chol;
  Chol chol(Sigma);
  const Chol::Traits::MatrixL& L = chol.matrixL();
  // Pre-compute determinant
  const double Ldet = L.determinant();
  
  // Iterate over multiple evaluations
  for(int i = 0; i < n; ++i) {
    
    // Compute the quadratic form of mvn density for observation i
    double quadform = (L.solve(y.row(i).transpose() - mu)).squaredNorm();
    // Evaluate density for observation i
    double dmvn = std::max(1e-15,
                           std::exp(-p * lnSqrt2Pi - 0.5 * quadform) / Ldet);
    // Weighted density for observation i 
    out[i] = std::pow(prob * dmvn, lik_weights[i]);
  }
  
  return out;
}


Eigen::MatrixXd e_step(const Eigen::MatrixXd &y,
                       const std::vector<Eigen::VectorXd>& mus,
                       const std::vector<Eigen::MatrixXd>& Sigmas,
                       const Eigen::VectorXd probs,
                       const Eigen::VectorXd &lik_weights,
                       const double temp) {
  
  const int n = y.rows();
  const int K = mus.size();
  Eigen::MatrixXd out(n, K);
  
  // Un-normalized e-step
  for(int k = 0; k < K; ++k) {
    out.col(k) = weighted_dens(y, mus[k], Sigmas[k], probs[k], lik_weights);
  }
  
  for(int i = 0; i < n; ++i) {
    
    // Normalize
    double rowSum = out.row(i).sum() + 1e-15;
    out.row(i) /= rowSum;
    
    // Apply tempering
    out.row(i) = out.row(i).array().pow(1/temp);
    
    // Normalize (again)
    rowSum = out.row(i).sum() + 1e-15;
    out.row(i) /= rowSum;
    
  }

  return out;
}
