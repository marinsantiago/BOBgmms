#ifndef E_STEP_H
#define E_STEP_H

#include <RcppEigen.h>

Eigen::VectorXd weighted_dens(const Eigen::MatrixXd &y,
                              const Eigen::VectorXd &mu,
                              const Eigen::MatrixXd &Sigma,
                              const double prob,
                              const Eigen::VectorXd &lik_weights);

Eigen::MatrixXd e_step(const Eigen::MatrixXd &y,
                       const std::vector<Eigen::VectorXd>& mus,
                       const std::vector<Eigen::MatrixXd>& Sigmas,
                       const Eigen::VectorXd probs,
                       const Eigen::VectorXd &lik_weights,
                       const double temp);

#endif  // E_STEP_H
