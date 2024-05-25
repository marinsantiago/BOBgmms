#ifndef EM_OPTIM_H
#define EM_OPTIM_H

#include <RcppEigen.h>

Eigen::VectorXd get_betaBar(const Eigen::VectorXd& beta_vec,
                            const Eigen::VectorXd& ybar_w,
                            const double lambdaTilde,
                            const double lambdaBar,
                            const double n_w);

Eigen::MatrixXd get_psiBar(const Eigen::MatrixXd& psiTilde,
                           const Eigen::MatrixXd& SRRTilde,
                           const Eigen::VectorXd& beta_vec,
                           const Eigen::VectorXd& ybar_w,
                           const double lambdaTilde,
                           const double lambdaBar,
                           const double n_w);

Eigen::MatrixXd update_sigma(const Eigen::MatrixXd& Psi_bar,
                             const double nu_bar); 

Eigen::VectorXd update_probs(const Eigen::VectorXd& alpha_bar);

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
                         double epss);

Eigen::MatrixXd recover_Z(const Eigen::MatrixXd &y,
                          const Eigen::VectorXd& theta,
                          std::vector<Eigen::VectorXd>& betas,
                          std::vector<Eigen::MatrixXd>& Psis,
                          const bool recover_Q);

#endif  // EM_OPTIM_H
