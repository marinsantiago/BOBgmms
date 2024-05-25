#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include "distances.h"


// [[Rcpp::export]]
double tv_distance(const Eigen::MatrixXd &target,
                   const Eigen::MatrixXd &approx,
                   const double step) {
  
  const int p = target.cols();
  
  // Set up environments
  Rcpp::Environment stats_pkg = Rcpp::Environment::namespace_env("stats");
  Rcpp::Function ecdf         = stats_pkg["ecdf"];
  
  // Initialize returns
  double out = 0.0;
  
  for(int j = 0; j < p; ++j) {
    
    Eigen::VectorXd target_variable = target.col(j);
    Eigen::VectorXd approx_variable = approx.col(j);
  
    // Define range for the TV distance
    double minVal         = target_variable.minCoeff() - 0.5;
    double maxVal         = target_variable.maxCoeff() + 0.5;
    int size              = static_cast<int>((maxVal - minVal) / step) + 1;
    Eigen::VectorXd range = Eigen::VectorXd::LinSpaced(size, minVal, maxVal);
    
    // Compute and evaluate ECDFs
    Rcpp::Function Ftarget        = ecdf(target_variable);
    Rcpp::Function Fapprox        = ecdf(approx_variable);
    Eigen::VectorXd Ftarget_evals = Rcpp::as<Eigen::VectorXd>(Ftarget(range));
    Eigen::VectorXd Fapprox_evals = Rcpp::as<Eigen::VectorXd>(Fapprox(range));
    
    // Trapezoidal rule
    double diff = (Ftarget_evals - Fapprox_evals).array().abs().sum();
    double tv   = diff * (step / 2);
    
    out += tv;
    
  }
  
  // Normalize
  return out / p;
  
}


// [[Rcpp::export]]
double ks_distance(const Eigen::MatrixXd &target,
                   const Eigen::MatrixXd &approx,
                   const double step) {
  
  const int p = target.cols();
  
  // Set up environments
  Rcpp::Environment stats_pkg = Rcpp::Environment::namespace_env("stats");
  Rcpp::Function ecdf         = stats_pkg["ecdf"];
  
  // Initialize returns
  double out = 0.0;
  
  for(int j = 0; j < p; ++j) {
    
    Eigen::VectorXd target_variable = target.col(j);
    Eigen::VectorXd approx_variable = approx.col(j);
    
    // Define range for the TV distance
    double minVal         = target_variable.minCoeff() - 0.5;
    double maxVal         = target_variable.maxCoeff() + 0.5;
    int size              = static_cast<int>((maxVal - minVal) / step) + 1;
    Eigen::VectorXd range = Eigen::VectorXd::LinSpaced(size, minVal, maxVal);
    
    // Compute and evaluate ECDFs
    Rcpp::Function Ftarget        = ecdf(target_variable);
    Rcpp::Function Fapprox        = ecdf(approx_variable);
    Eigen::VectorXd Ftarget_evals = Rcpp::as<Eigen::VectorXd>(Ftarget(range));
    Eigen::VectorXd Fapprox_evals = Rcpp::as<Eigen::VectorXd>(Fapprox(range));
    
    // Compute supremum
    Eigen::VectorXd diff = (Ftarget_evals - Fapprox_evals).array().abs();
    double ks            = diff.maxCoeff();
    
    out += ks;
    
  }
  
  // Normalize
  return out / p;
  
}
