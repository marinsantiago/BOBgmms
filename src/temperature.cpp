#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::depends(RcppProgress)]]

// #include <progress.hpp>
// #include <progress_bar.hpp>

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include "temperature.h"


double temperature(const double t,
                   const double a,
                   const double b,
                   const double c,
                   const double r) {

  double tau  = (t + (c * r)) / r;
  double temp = 1 + std::pow(a, tau) + b * (std::sin(tau) / tau);

  return temp;
}
