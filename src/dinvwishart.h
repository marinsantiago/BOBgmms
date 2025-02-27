#ifndef DINVWISHART_H
#define DINVWISHART_H

#include <RcppEigen.h>


double ldinvwishart(const Eigen::MatrixXd & Sigma,
                    const Eigen::MatrixXd & S,
                    const double nu);

double lmvgamma(double x, int p);

double log_det(const Eigen::MatrixXd & x);

#endif  // DINVWISHART_H
