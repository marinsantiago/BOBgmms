#ifndef WEIGHTED_MOMENTS_H
#define WEIGHTED_MOMENTS_H

#include <RcppEigen.h>

Eigen::VectorXd weighted_mean(const Eigen::MatrixXd& y,
                              const Eigen::VectorXd& uq,
                              const double n_weights); // Function declaration

Eigen::MatrixXd weighted_ssr(const Eigen::MatrixXd& y,
                             const Eigen::VectorXd& ybar_w,
                             const Eigen::VectorXd& uq); // Function declaration

#endif  // WEIGHTED_MOMENTS_H
