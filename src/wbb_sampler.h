#ifndef WBB_SAMPLER_H
#define WBB_SAMPLER_H

#include <RcppEigen.h>
#include <Rcpp.h>

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
                    const bool randPrior_w);

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
                           const bool randPrior_w);

#endif  // WBB_SAMPLER_H
