#ifndef DISTANCES_H
#define DISTANCES_H

#include <RcppEigen.h>

double tv_distance(const Eigen::MatrixXd &target,
                   const Eigen::MatrixXd &approx,
                   const double step);

double ks_distance(const Eigen::MatrixXd &target,
                   const Eigen::MatrixXd &approx,
                   const double step);

#endif  // DISTANCES_H
