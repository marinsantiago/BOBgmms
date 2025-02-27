#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include "e_step.h"
#include "update_params.h"
#include "utils.h"
#include "em_optim.h"


/**
 * This function maximizes the weighted log-posterior density associated with  \
 * a Gaussian Mixture Model via the EM algorithm, as in Marin et al. (2025+).  \ 
 *                                                                             \
 * @param y Eigen::MatrixXd: A matrix of observations of dimension             \
 *  \eqn{n}-by-\eqn{d}, where each  of the \eqn{n} rows is an                  \
 *  observation vector.                                                        \
 * @param means_init std::vector<Eigen::VectorXd>: A list of size \eqn{K}      \
 *  containing the initial mean vectors, each vector of size \eqn{d}.          \
 * @param covs_init std::vector<Eigen::MatrixXd>: A list of of size \eqn{K}    \
 *  containing the initial covariance matrices, each matrix                    \
 *  of size \eqn{d}-by-\eqn{d}.                                                \
 * @param probs_init Eigen::VectorXd: A vector of size \eqn{d} of initial      \
 *  probabilities. The sum of all the probabilities should equals one.         \
 * @param alphas Eigen::VectorXd: A vector of size \eqn{K} of Dirichlet prior  \
 *  hyper-parameters. All the elements in the vector must be larger than one.  \
 * @param betas std::vector<Eigen::VectorXd>: A list of size \eqn{K}           \
 *  containing the prior means of each mean vector. Each prior mean should     \
 *  be a vector of size \eqn{d}.                                               \
 * @param lambdas Eigen::VectorXd: A vector of size \eqn{K} containing the     \
 *  \eqn{\lambda} shrinkage hyper-parameters associated with each mean vector. \
 *  Each element in the vector must be positive.                               \
 * @param nus Eigen::VectorXd: A vector of size \eqn{K} containing the         \
 *  \eqn{\nu} hyper-parameters associated with each covariance matrix. Each    \
 *  element in the vector must be larger than \eqn{d + 1}.                     \
 * @param Psis std::vector<Eigen::MatrixXd>: A list of size \eqn{K} containing \ 
 *  the prior scale matrices of each covariance matrix. Each element in the    \
 *  list must be a symmetric positive semidefinite matrix of size              \
 *  \eqn{d}-by-\eqn{d}.                                                        \
 * @param lik_weights Eigen::VectorXd: A vector of size \eqn{n} of likelihood  \
 *  weights. All weights must be non-negative.                                 \
 * @param means_weights: A vector of size \eqn{K} of weights associated with   \
 *  each mean vector. All weights must be non-negative.                        \
 * @param covs_weights Eigen::VectorXd: A vector of size \eqn{K} of weights    \
 *  associated with each covariance matrix. All weights must be non-negative.  \
 * @params probs_weights Eigen::VectorXd: A vector of size \eqn{K} of weights  \
 *  associated with each mixture proportion. All weights must be non-negative. \
 * @param a double: A scalar corresponding to the \eqn{a} hyper-parameter in   \
 *  the tempering profile. Note that \eqn{a\in [0,1)}.                         \
 * @param b double: A scalar corresponding to the \eqn{b} hyper-parameter in   \
 *  the tempering profile. Note that \eqn{b\in \mathbb{R}}.                    \
 * @param c double: A scalar corresponding to the \eqn{c} hyper-parameter in   \
 *  the tempering profile. Note that \eqn{c>0}.                                \  
 * @param r double: A scalar corresponding to the \eqn{r} hyper-parameter in   \
 *  the tempering profile. Note that \eqn{r>0}.                                \
 * @param max_iters int: Maximum number of EM iterations.                      \
 * @param epss double: Convergence threshold for the EM algorithm.             \
 *                                                                             \
 * @return Rcpp::List: A list containing the EM solution for each mean vector, \
 *  each covariance matrix, and each mixture proportion in the mixture model.  \
 */
// [[Rcpp::export]]
Rcpp::List em_optim_cpp_(const Eigen::MatrixXd & y,
                         std::vector<Eigen::VectorXd> & means_init,
                         std::vector<Eigen::MatrixXd> & covs_init,
                         Eigen::VectorXd probs_init,
                         const Eigen::VectorXd & alphas,
                         const std::vector<Eigen::VectorXd> & betas,
                         const Eigen::VectorXd & lambdas,
                         const Eigen::VectorXd & nus,
                         const std::vector<Eigen::MatrixXd> & Psis,
                         const Eigen::VectorXd & lik_weights,
                         const Eigen::VectorXd & means_weights,
                         const Eigen::VectorXd & covs_weights,
                         const Eigen::VectorXd & probs_weights,
                         const double a,
                         const double b,
                         const double c,
                         const double r,
                         int max_iters,    // 500
                         double epss) {    // 1e-4
  
  // Prepare the returns
  const int n = y.rows();
  const int p = y.cols();
  const int K = means_init.size();
  Rcpp::List out(3); // means, covs, and probs.
  
  // Initialize updated version of the parameters
  std::vector<Eigen::VectorXd> & means_new = means_init;
  std::vector<Eigen::MatrixXd> & covs_new = covs_init;
  Eigen::VectorXd probs_new = probs_init;
  
  // Pre-compute constants
  Eigen::VectorXd lambdas_tilde;
  lambdas_tilde = (lambdas.array() * means_weights.array()).array() + 1e-10;
  Eigen::VectorXd nus_tilde = nus.array() + p + 2;
  nus_tilde = nus_tilde.array() * covs_weights.array();
  nus_tilde = nus_tilde.array() - 2 - p;
  Eigen::VectorXd alphas_tilde = alphas.array() - 1;
  alphas_tilde = alphas_tilde.array() * probs_weights.array();
  alphas_tilde = alphas_tilde.array() + 1;
  std::vector<Eigen::MatrixXd> Psis_tilde;
  Psis_tilde.resize(K);
  for (int k = 0; k < K; ++k) {
    Eigen::MatrixXd psi_tilde_k = Eigen::MatrixXd::Identity(p, p) * 1e-10;
    psi_tilde_k = (Psis[k] * covs_weights[k]) + psi_tilde_k;
    Psis_tilde[k] = psi_tilde_k;
  }
  
  // EM iterations -------------------------------------------------------------
  Eigen::VectorXd alphas_bar(K);
  for (int t = 0; t < max_iters; ++t) {
    
    // Expectation step --------------------------------------------------------
    // Loop starts at 0 (add 1 to the current temperature)
    double temp = temperature(t + 1, a, b, c, r);
    // Weighted expected value of latent variables
    Eigen::MatrixXd Q = e_step(
      y, means_init, covs_init, probs_init, lik_weights, temp
    );
    
    // Maximization step -------------------------------------------------------
    for (int k = 0; k < K; ++k) {
      Eigen::VectorXd uq_k = lik_weights.array() * Q.col(k).array();
      double n_tilde_k = uq_k.sum() + 1e-10;
      // Weighted mean vector and SSR matrix
      Rcpp::List w_moments_k = weighted_moments(y, uq_k, n_tilde_k);
      Eigen::VectorXd ybar_tilde_k = w_moments_k[0];
      Eigen::MatrixXd ssr_tilde_k = w_moments_k[1];
      // Posterior hyper-parameters
      double lambda_bar_k = lambdas_tilde[k] + n_tilde_k;
      double nu_bar_k = nus_tilde[k] + n_tilde_k;
      Eigen::VectorXd beta_bar_k = get_beta_bar(
        betas[k], ybar_tilde_k, lambdas_tilde[k], lambda_bar_k, n_tilde_k
      );
      Eigen::MatrixXd Psi_bar_k = get_Psi_bar(
        Psis_tilde[k], ssr_tilde_k, betas[k], ybar_tilde_k, 
        lambdas_tilde[k], lambda_bar_k, n_tilde_k
      );
      // Parameter updates
      Eigen::MatrixXd cov_new_k;
      cov_new_k = Psi_bar_k / std::max(nu_bar_k + p + 1, 1e-10);
      cov_new_k = (cov_new_k + cov_new_k.transpose()) / 2.0;  // Make symmetric
      cov_new_k += Eigen::MatrixXd::Identity(p, p) * 1e-10;
      covs_new[k] = cov_new_k;
      means_new[k] = beta_bar_k;
      alphas_bar[k] = alphas_tilde[k] + n_tilde_k; //Q.col(k).sum() + 1e-10;
    }
    probs_new = alphas_bar.array() - 1;
    probs_new = probs_new / (alphas_bar.sum() - K + 1e-10);
    
    // Check for EM convergence ------------------------------------------------
    double err = chck_convergence(
      means_init, covs_init, probs_init, means_new, covs_new, probs_new
    );
    if (err < epss) {
      break;
    }
    // Otherwise, keep updating the model parameters
    means_init = means_new;
    covs_init = covs_new;
    probs_init = probs_new;
  }
  
  out[0] = means_new;
  out[1] = covs_new;
  out[2] = probs_init;
  return out;
}
