#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include "make_pos_def.h"
#include "weighted_moments.h"


Eigen::VectorXd weighted_mean(const Eigen::MatrixXd& y,
                              const Eigen::VectorXd& uq,
                              const double n_weights) {
  
  
  const int p = y.cols();
  
  Eigen::VectorXd out;
  
  // Copy 'y' to avoid modifying the original
  Eigen::MatrixXd y_uq = y;
  
  // Multiply the cols of 'y' element-wise with vector 'weights'
  y_uq.array().colwise() *= uq.array();
  
  // Weighted sample mean vector 
  out = y_uq.colwise().sum() / (n_weights + 1e-15);
  
  // Check all the values are real (smaller than abs(1e+14))
  const double thold = 1e+10;
  double chck        = out.cwiseAbs().maxCoeff();
  if (chck < thold) {
    return out;
  } else {
    return Eigen::VectorXd::Zero(p);
  }
  
}


Eigen::MatrixXd weighted_ssr(const Eigen::MatrixXd& y,
                             const Eigen::VectorXd& ybar_w,
                             const Eigen::VectorXd& uq) {
  
  const double p = y.cols();
  const double n = y.rows();
  
  Eigen::MatrixXd ssr = Eigen::MatrixXd::Zero(p, p);
  
  for(int i = 0; i < n; ++i) {
    
    // Subtract weighted mean from each observation
    Eigen::VectorXd diff = (y.row(i).transpose() - ybar_w) + 
      Eigen::VectorXd::Ones(p) * 1e-10;
    
    // Compute SSR
    ssr += (diff * diff.transpose()) * std::max(1e-15, uq[i]);
  }
  
  // Check if all values are finite (non-NaN, non-infinity)
  if (ssr.allFinite()) {
    return make_symmetric(ssr);
  } else {
    return Eigen::MatrixXd::Zero(p, p);
  }
  
}
