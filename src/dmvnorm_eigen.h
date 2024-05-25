#ifndef DMVNORM_EIGEN_H
#define DMVNORM_EIGEN_H

#include <RcppEigen.h>

double dmvnorm_eigen(const Eigen::VectorXd &x,
                     const Eigen::VectorXd &mu,
                     const Eigen::MatrixXd &Sigma); 

Eigen::VectorXd dmvnorm_multievals_eigen(const Eigen::MatrixXd &x,
                                         const Eigen::VectorXd &mu,
                                         const Eigen::MatrixXd &Sigma);

#endif  // DMVNORM_EIGEN_H
