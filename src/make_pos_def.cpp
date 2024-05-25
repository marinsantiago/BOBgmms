#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include "make_pos_def.h"


Eigen::MatrixXd make_positive_definite(const Eigen::MatrixXd& x) {
  
  const int p = x.cols();
  
  // Compute the smallest eigenvalue
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(x);
  const double smallestEigenval = eigenSolver.eigenvalues().minCoeff();
  
  // Compute the adjusting quantity
  const double adjst  = std::max(0.0, -smallestEigenval) + 1e-10;
  Eigen::MatrixXd out = x + Eigen::MatrixXd::Identity(p, p) * adjst;
  
  return out;
}


Eigen::MatrixXd make_symmetric(const Eigen::MatrixXd& x) {
  
  Eigen::MatrixXd out = 0.5 * (x + x.transpose());
  
  return out;
}
