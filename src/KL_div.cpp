#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include "em_optim.h"
#include "KL_div.h"


// [[Rcpp::export]]
Eigen::VectorXd batchOneShot(double x_likelihood,
                             Eigen::VectorXd x_mus,
                             Eigen::VectorXd x_sigmas,
                             double x_probs,
                             const Eigen::MatrixXd &y,
                             std::vector<Eigen::VectorXd>& mus,
                             std::vector<Eigen::MatrixXd>& Sigmas,
                             Eigen::VectorXd probs,
                             const double a,
                             const double b,
                             const double c,
                             const double r,
                             const Eigen::VectorXd &alphas,
                             const std::vector<Eigen::VectorXd>& betas,
                             const Eigen::VectorXd &lambdas,
                             const Eigen::VectorXd &nus,
                             const std::vector<Eigen::MatrixXd>& Psis) {

  const int n = y.rows();

  Rcpp::Environment stats_pkg = Rcpp::Environment::namespace_env("stats");
  Rcpp::Function rexpp        = stats_pkg["rexp"];

  // Step 1: Sample weights associated with the likelihood
  Eigen::VectorXd u           = Rcpp::as<Eigen::VectorXd>(rexpp(n, 1));
  Eigen::VectorXd u_a         = u.array().pow(x_likelihood);
  Eigen::VectorXd lik_weights = (u_a * n) / u_a.sum();

  // Step 2: Sample weights associated with the prior
  double prob_weights          = x_probs;
  Eigen::VectorXd cov_weights  = x_sigmas;
  Eigen::VectorXd mean_weights = x_mus;

  // Step 3: Find the posterior mode
  int max_iters = 500;     // Max iterations EM algorithm
  double epss   = 1e-4;    // Stopping rule EM algorithm

  // EM algorithm (tempered)
  Eigen::VectorXd Theta_out = em_optim(y,
                                       mus,
                                       Sigmas,
                                       probs,
                                       a,
                                       b,
                                       c,
                                       r,
                                       lik_weights,
                                       mean_weights,
                                       cov_weights,
                                       prob_weights,
                                       alphas,
                                       betas,
                                       lambdas,
                                       nus,
                                       Psis,
                                       max_iters,
                                       epss);

  return Theta_out;
}
