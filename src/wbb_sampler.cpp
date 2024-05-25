#include <RcppEigen.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include "helpers.h"
#include "em_optim.h"
#include "wbb_sampler.h"


Eigen::MatrixXd wbb(const Eigen::MatrixXd &y,
                    std::vector<Eigen::VectorXd>& mus,
                    std::vector<Eigen::MatrixXd>& Sigmas,
                    Eigen::VectorXd probs,
                    const int S,
                    const double a,
                    const double b,
                    const double c,
                    const double r,
                    const Eigen::VectorXd &alphas,
                    const std::vector<Eigen::VectorXd>& betas,
                    const Eigen::VectorXd &lambdas,
                    const Eigen::VectorXd &nus,
                    const std::vector<Eigen::MatrixXd>& Psis,
                    const bool randPrior_w) {

  const int n = y.rows();
  const int p = y.cols();
  const int K = mus.size();

  const int theta_dim = (pow(p, 2) + p + 1) * K;
  Eigen::MatrixXd Theta_out(S, theta_dim);

  for(int s = 0; s < S; ++s) {

    // Step 1: Sample weights associated with the likelihood
    Rcpp::NumericVector u       = Rcpp::rexp(n, 1);
    Eigen::VectorXd lik_weights = convertEigenVec(u);

    // Step 2: Sample weights associated with the prior
    double prob_weights;
    Eigen::VectorXd cov_weights;
    Eigen::VectorXd mean_weights;

    if (randPrior_w) {

      // Random prior weights
      Rcpp::NumericVector u_sigma = Rcpp::rexp(K, 1);
      Rcpp::NumericVector u_mu    = Rcpp::rexp(K, 1);
      prob_weights                = Rcpp::rexp(1, 1)[0];
      mean_weights                = convertEigenVec(u_mu);
      cov_weights                 = convertEigenVec(u_sigma);

    } else {

      // Fixed prior weights
      prob_weights = 1;
      cov_weights  = Eigen::VectorXd::Ones(K);
      mean_weights = Eigen::VectorXd::Ones(K);

    }

    // Step 3: Find the posterior mode

    int max_iters = 500;     // Max iterations EM algorithm
    double epss   = 1e-4;    // Stopping rule EM algorithm

    // EM algorithm (tempered)
    Theta_out.row(s) = em_optim(y,
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
  }

  return Theta_out;
}


// [[Rcpp::export]]
Eigen::VectorXd wbbOneShot(const Eigen::MatrixXd &y,
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
                           const std::vector<Eigen::MatrixXd>& Psis,
                           const bool randPrior_w) {

  const int n = y.rows();
  const int K = mus.size();

  // Step 1: Sample weights associated with the likelihood
  Rcpp::NumericVector u       = Rcpp::rexp(n, 1);
  Eigen::VectorXd lik_weights = convertEigenVec(u);

  // Step 2: Sample weights associated with the prior
  double prob_weights;
  Eigen::VectorXd cov_weights;
  Eigen::VectorXd mean_weights;

  if (randPrior_w) {

    // Random prior weights
    Rcpp::NumericVector u_sigma = Rcpp::rexp(K, 1);
    Rcpp::NumericVector u_mu    = Rcpp::rexp(K, 1);
    prob_weights                = Rcpp::rexp(1, 1)[0];
    mean_weights                = convertEigenVec(u_mu);
    cov_weights                 = convertEigenVec(u_sigma);

  } else {

    // Fixed prior weights
    prob_weights = 1;
    cov_weights  = Eigen::VectorXd::Ones(K);
    mean_weights = Eigen::VectorXd::Ones(K);

  }

  // Step 3: Find the posterior mode

  int max_iters = 500;     // Max iterations EM algorithm
  double epss   = 1e-4;    // Stopping rule EM algorithm

  // EM algorithm (tempered)
  Eigen::VectorXd theta_out = em_optim(y,
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

  return theta_out;
}
