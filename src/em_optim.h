#ifndef EM_OPTIM_H
#define EM_OPTIM_H

#include <RcppEigen.h>


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
                         int max_iters,      
                         double epss);

#endif  // EM_OPTIM_H
