#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include "dmvnorm_eigen.h"


double dmvnorm_eigen(const Eigen::VectorXd &x,
                     const Eigen::VectorXd &mu,
                     const Eigen::MatrixXd &Sigma) {
  
  int p = x.size();
  
  const double logSqrt2Pi = 0.5 * std::log(2 * M_PI);
  
  // Cholesky decomposition
  typedef Eigen::LLT<Eigen::MatrixXd> Chol;
  Chol chol(Sigma);
  const Chol::Traits::MatrixL& L = chol.matrixL();
  
  // Compute the squared quadratic form
  double quadform = (L.solve(x - mu)).squaredNorm();
  
  return std::exp(-p * logSqrt2Pi - 0.5 * quadform) / L.determinant();
}


Eigen::VectorXd dmvnorm_multievals_eigen(const Eigen::MatrixXd &x,
                                         const Eigen::VectorXd &mu,
                                         const Eigen::MatrixXd &Sigma) {
  
  int n = x.rows(), p = x.cols();
  Eigen::VectorXd out(n);
  
  const double logSqrt2Pi = 0.5 * std::log(2 * M_PI);
  
  // Cholesky decomposition
  typedef Eigen::LLT<Eigen::MatrixXd> Chol;
  Chol chol(Sigma);
  const Chol::Traits::MatrixL& L = chol.matrixL();
  
  // Pre-compute determinant
  const double Ldet = L.determinant();
  
  // Iterate over multiple evaluations
  double quadform = 0;
  for(int i = 0; i < n; ++i) {
    
    // Compute the squared quadratic form
    quadform = (L.solve(x.row(i).transpose() - mu)).squaredNorm();
    
    // Evaluate density
    out[i] = std::exp(-p * logSqrt2Pi - 0.5 * quadform) / Ldet;
  }
  
  return out;
}
