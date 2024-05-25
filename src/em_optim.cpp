#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include "e_step.h"
#include "temperature.h"
#include "make_pos_def.h"
#include "weighted_moments.h"
#include "helpers.h"
#include "em_optim.h"


Eigen::VectorXd get_betaBar(const Eigen::VectorXd& beta_vec,
                            const Eigen::VectorXd& ybar_w,
                            const double lambdaTilde,
                            const double lambdaBar,
                            const double n_w) {

  Eigen::VectorXd out;

  // Output
  Eigen::VectorXd betaTilde = (lambdaTilde / lambdaBar) * beta_vec;
  out                       = betaTilde + (n_w / lambdaBar) * ybar_w;

  // Check all the values are real (smaller than abs(1e+14))
  const double thold = 1e+10;
  double chck        = out.cwiseAbs().maxCoeff();
  if (chck < thold) {
    return out;
  } else {
    return betaTilde;
  }

}


Eigen::MatrixXd get_psiBar(const Eigen::MatrixXd& psiTilde,
                           const Eigen::MatrixXd& SRRTilde,
                           const Eigen::VectorXd& beta_vec,
                           const Eigen::VectorXd& ybar_w,
                           const double lambdaTilde,
                           const double lambdaBar,
                           const double n_w) {

  const int p = psiTilde.cols();

  Eigen::MatrixXd out;

  double adj_diff           = (lambdaTilde * n_w) / (lambdaBar + 1e-15);
  Eigen::VectorXd diff_mean = (ybar_w - beta_vec) +
    Eigen::VectorXd::Ones(p) * 1e-10;;

  // Output
  out = psiTilde + SRRTilde + (diff_mean * diff_mean.transpose()) * adj_diff;

  // Check if all values are finite (non-NaN, non-infinity)
  if (out.allFinite()) {
    return make_symmetric(out);
  } else {
    return psiTilde;
  }

}


Eigen::MatrixXd update_sigma(const Eigen::MatrixXd& Psi_bar,
                             const double nu_bar) {

  const int p = Psi_bar.cols();

  double sigma_adj    = nu_bar + p + 1 + 1e-15;
  Eigen::MatrixXd out = make_symmetric(Psi_bar / sigma_adj) +
    Eigen::MatrixXd::Identity(p, p) * 1e-15;

  return make_positive_definite(out);
}


Eigen::VectorXd update_probs(const Eigen::VectorXd& alpha_bar) {

  const int K = alpha_bar.size();
  Eigen::VectorXd out(K);

  double sum_alpha = alpha_bar.sum();
  for(int k = 0; k < K; ++k) {

    double pi_new_k = (alpha_bar[k] - 1) / (sum_alpha - K + 1e-25);

    // Make sure the probabilities are in (0,1)
    double max_prob = 1 - 1e-15;
    double min_prob = 1e-15;
    pi_new_k        = std::min(max_prob, pi_new_k);
    pi_new_k        = std::max(min_prob, pi_new_k);
    out[k]          = pi_new_k;
  }

  return out / out.sum();
}


// [[Rcpp::export]]
Eigen::VectorXd em_optim(const Eigen::MatrixXd &y,
                         std::vector<Eigen::VectorXd>& mus,
                         std::vector<Eigen::MatrixXd>& Sigmas,
                         Eigen::VectorXd probs,
                         const double a,
                         const double b,
                         const double c,
                         const double r,
                         const Eigen::VectorXd &lik_weights,
                         const Eigen::VectorXd &mean_weights,
                         const Eigen::VectorXd &cov_weights,
                         const double prob_weights,
                         const Eigen::VectorXd &alphas,
                         const std::vector<Eigen::VectorXd>& betas,
                         const Eigen::VectorXd &lambdas,
                         const Eigen::VectorXd &nus,
                         const std::vector<Eigen::MatrixXd>& Psis,
                         int max_iters,    // 500
                         double epss) {    // 1e-4

  const int n = y.rows();
  const int p = y.cols();
  const int K = mus.size();

  // Copy the initial values of the parameters to avoid modifying the originals
  std::vector<Eigen::VectorXd>& mus_current    = mus;
  std::vector<Eigen::MatrixXd>& Sigmas_current = Sigmas;
  Eigen::VectorXd probs_current                = probs;

  // Initialize updated version of the parameters
  std::vector<Eigen::VectorXd>& mus_new        = mus;
  std::vector<Eigen::MatrixXd>& Sigmas_new     = Sigmas;
  Eigen::VectorXd probs_new                    = probs;

  Eigen::VectorXd theta_new;

  for(int t = 0; t < max_iters; ++t) {

    // E-step (tempered)
    double temp            = temperature(t + 1, // Loop starts at 0 (add 1)
                                         a, b, c, r);
    Eigen::MatrixXd Q_temp = e_step(y, mus, Sigmas, probs, lik_weights, temp);

    // M-steps
    Eigen::VectorXd alphaBar(K);

    for(int k = 0; k < K; ++k) {

      // Define useful intermediate quantities
      Eigen::VectorXd uq_k       = lik_weights.cwiseProduct(Q_temp.col(k));
      double n_k                 = std::max(1e-15, uq_k.sum());
      double lambTilde_k         = std::max(1e-15,
                                            lambdas[k] * mean_weights[k]);
      double nnTild              = (cov_weights[k] * (nus[k] + p + 2)) - 2 - p;
      double nuTilde_k           = std::max(1e-15, nnTild);
      Eigen::MatrixXd PsiTilde_k = Psis[k] * cov_weights[k];

      // Weighted mean vector and SSR matrix
      Eigen::VectorXd ybarTilde_k = weighted_mean(y, uq_k, n_k);
      Eigen::MatrixXd STilde_k    = weighted_ssr(y, ybarTilde_k, uq_k);

      // Posterior hyper-parameters
      double lambBar_k          = lambTilde_k + n_k + 1e-15;
      double nuBar_k            = nuTilde_k + n_k + 1e-15;
      alphaBar[k]               = prob_weights * (alphas[k] - 1) + 1 + n_k;
      Eigen::VectorXd betaBar_k = get_betaBar(betas[k],
                                              ybarTilde_k,
                                              lambTilde_k,
                                              lambBar_k,
                                              n_k);

      Eigen::MatrixXd PsiBar_k  = get_psiBar(PsiTilde_k,
                                             STilde_k,
                                             betas[k],
                                             ybarTilde_k,
                                             lambTilde_k,
                                             lambBar_k,
                                             n_k);

      // Parameter updates
      mus_new[k]    = betaBar_k;
      Sigmas_new[k] = update_sigma(PsiBar_k, nuBar_k);

    }
    probs_new = update_probs(alphaBar);

    // Collection of all current parameter values
    Eigen::VectorXd theta_current = collect_params(vec_vec(mus_current),
                                                   vec_mat(Sigmas_current),
                                                   probs_current);

    // Collection of all updated parameter values
    theta_new = collect_params(vec_vec(mus_new),
                               vec_mat(Sigmas_new),
                               probs_new);

    // Check for convergence
    double er = ((theta_current - theta_new).squaredNorm()) / theta_new.size();

    if (er < epss) {
      break;
    }

    // Otherwise, keep updating the parameters
    mus_current    = mus_new;
    probs_current  = probs_new;
    Sigmas_current = Sigmas_new;
  }

  return theta_new;
}


Eigen::MatrixXd recover_Z(const Eigen::MatrixXd &y,
                          const Eigen::VectorXd& theta,
                          std::vector<Eigen::VectorXd>& betas,
                          std::vector<Eigen::MatrixXd>& Psis,
                          const bool recover_Q) {

  const int n = y.rows();
  const int p = y.cols();
  const int K = betas.size();

  Eigen::MatrixXd Z(n, K);

  // Recover parameters
  std::vector<Eigen::VectorXd> mus    = recov_means(theta,
                                                    betas,
                                                    p,
                                                    K);

  std::vector<Eigen::MatrixXd> sigmas = recov_covs(theta,
                                                   Psis,
                                                   p,
                                                   K);

  Eigen::VectorXd probs               = recov_probs(theta,
                                                    p,
                                                    K);

  // Fixed Likelihood weights
  Eigen::VectorXd lik_weights         = Eigen::VectorXd::Ones(n);

  // No tempering
  const double temp = 1.0;

  // Recover Q matrix
  Eigen::MatrixXd Q = e_step(y, mus, sigmas, probs, lik_weights, temp);
  if (recover_Q) {
    return Q;
  }

  // Otherwise recover Z
  for(int i = 0; i < n; ++i) {

    Eigen::VectorXd z_row = Eigen::VectorXd::Zero(K);
    // Find the index of the maximum element in each row of Q
    Eigen::Index maxIndex;
    double maxValue       = Q.row(i).maxCoeff(&maxIndex);
    // Set the element at maxIndex to one
    z_row(maxIndex)       = 1.0;
    Z.row(i)              = z_row;

  }

  return Z;
}
