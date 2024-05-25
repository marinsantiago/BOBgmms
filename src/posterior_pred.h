#ifndef POSTERIOR_PRED_H
#define POSTERIOR_PRED_H

#include <RcppEigen.h>

Eigen::MatrixXd post_pred(const Eigen::MatrixXd& post_draws,
                          std::vector<Eigen::VectorXd> betas,
                          std::vector<Eigen::MatrixXd>& Psis);

#endif  // POSTERIOR_PRED_H
