#ifndef UPDATE_PARAMS_H
#define UPDATE_PARAMS_H

#include <RcppEigen.h>


double temperature(const int t,
                   const double a,
                   const double b,
                   const double c,
                   const double r);

Rcpp::List weighted_moments(const Eigen::MatrixXd & y,
                            const Eigen::VectorXd & uq,
                            const double n_tilde);

Eigen::VectorXd get_beta_bar(const Eigen::VectorXd & beta_vec,
                             const Eigen::VectorXd & w_mean,
                             const double lambda_tilde,
                             const double lambda_bar,
                             const double n_tilde);

Eigen::MatrixXd get_Psi_bar(const Eigen::MatrixXd & Psi_tilde,
                            const Eigen::MatrixXd & w_ssr,
                            const Eigen::VectorXd & beta_vec,
                            const Eigen::VectorXd & w_mean,
                            const double lambda_tilde,
                            const double lambda_bar,
                            const double n_tilde);

#endif  // UPDATE_PARAMS_H
