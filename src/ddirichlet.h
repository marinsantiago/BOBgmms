#ifndef DDIRICHLET_H
#define DDIRICHLET_H

#include <RcppEigen.h>


double lddirichlet(const Eigen::VectorXd & x,
                   const Eigen::VectorXd & alphas);

double lmvbeta(const Eigen::VectorXd & alpha);

#endif  // DDIRICHLET_H
