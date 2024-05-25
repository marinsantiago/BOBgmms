#ifndef KL_DIV_H
#define KL_DIV_H

#include <RcppEigen.h>

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
                             const std::vector<Eigen::MatrixXd>& Psis);

#endif  // KL_DIV_H
