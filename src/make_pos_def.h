#ifndef MAKE_POS_DEF_H
#define MAKE_POS_DEF_H

#include <RcppEigen.h>

Eigen::MatrixXd make_positive_definite(const Eigen::MatrixXd& x); 

Eigen::MatrixXd make_symmetric(const Eigen::MatrixXd& x);

#endif  // MAKE_POS_DEF_H
