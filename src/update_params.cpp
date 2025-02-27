#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include "update_params.h"


double temperature(const int t,
                   const double a,
                   const double b,
                   const double c,
                   const double r) {
  
  double tau  = (t + (c * r)) / r;
  double temp = 1 + std::pow(a, tau) + b * (std::sin(tau) / tau);
  return temp;
} 


Rcpp::List weighted_moments(const Eigen::MatrixXd & y,
                            const Eigen::VectorXd & uq,
                            const double n_tilde) {
  
  // Prepare the returns
  const int p = y.cols();
  Rcpp::List out(2);
  
  // Weighted mean
  Eigen::VectorXd w_mean = (y.array().colwise() * uq.array()).colwise().sum();
  w_mean /= (n_tilde + 1e-10);
  out[0] = w_mean;
  
  // Subtract the weighted mean from each observation
  Eigen::MatrixXd y_centered = y.array().rowwise() - w_mean.transpose().array();
  
  // Weighted SSR matrix
  Eigen::MatrixXd y_centered_w;
  Eigen::VectorXd sqrt_uq = uq.array().abs().sqrt(); // abs. for num. stability
  y_centered_w = y_centered.array().colwise() * sqrt_uq.array();
  Eigen::MatrixXd w_ssr = y_centered_w.transpose() * y_centered_w;
  out[1] = w_ssr;
  
  return out;
}

// [[Rcpp::export]]
Eigen::VectorXd get_beta_bar(const Eigen::VectorXd & beta_vec,
                             const Eigen::VectorXd & w_mean,
                             const double lambda_tilde,
                             const double lambda_bar,
                             const double n_tilde) {
  
  // Prepare the returns
  const int p = beta_vec.size();
  Eigen::VectorXd out;
  
  // Updated beta
  out = (lambda_tilde/lambda_bar) * beta_vec + (n_tilde/lambda_bar) * w_mean;
  out = out.array().isNaN().select(beta_vec, out);

  return out;
} 


// [[Rcpp::export]]
Eigen::MatrixXd get_Psi_bar(const Eigen::MatrixXd & Psi_tilde,
                            const Eigen::MatrixXd & w_ssr,
                            const Eigen::VectorXd & beta_vec,
                            const Eigen::VectorXd & w_mean,
                            const double lambda_tilde,
                            const double lambda_bar,
                            const double n_tilde) {
  
  // Prepare the returns
  const int p = beta_vec.size();
  Eigen::MatrixXd out = Eigen::MatrixXd::Identity(p, p) * 1e-10;
  
  // Update psi
  Eigen::VectorXd w_beta = w_mean - beta_vec;
  out += (lambda_tilde * n_tilde) / lambda_bar * (w_beta * w_beta.transpose());
  out += Psi_tilde + w_ssr;
  if (out.array().isNaN().sum() > 0) {
    return Psi_tilde;
  }

  return out;
}
