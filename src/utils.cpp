#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include "utils.h"

// Function that takes a list of vectors and "unlist" it
Eigen::VectorXd unlist_vecs(const std::vector<Eigen::VectorXd> & list_of_vecs) {
  // Compute the total size of the resulting vector
  int totalSize = 0;
  for (const auto& vec : list_of_vecs) {
    totalSize += vec.size();
  }
  // Create a single vector by concatenating all the vectors
  Eigen::VectorXd out(totalSize);
  int idx = 0;
  for (const auto & vec : list_of_vecs) {
    Eigen::Map<Eigen::VectorXd>(out.data() + idx, vec.size()) = vec;
    idx += vec.size();
  }
  return out;
}


// Function that takes a list of matrices and "unlist" it
Eigen::VectorXd unlist_mats(const std::vector<Eigen::MatrixXd> & list_of_mats) {
  // Compute the total size of the resulting vector
  int total_size = 0;
  for (const auto & matrix : list_of_mats) {
    total_size += matrix.size();
  }
  // Create a single vector by concatenating the matrices
  Eigen::VectorXd out(total_size);
  int idx = 0;
  for (const auto& matrix : list_of_mats) {
    Eigen::Map<Eigen::VectorXd>(out.data() + idx, matrix.size()) =
      Eigen::Map<const Eigen::VectorXd>(matrix.data(), matrix.size());
    idx += matrix.size();
  }
  return out;
}


// Collect all model parameters in a single vector.
// v1: unlisted means
// v2: unlisted covs
// v3: vector of probabilities
Eigen::VectorXd collect_params(const Eigen::VectorXd & v1,
                               const Eigen::VectorXd & v2,
                               const Eigen::VectorXd & v3) {
  
  Eigen::VectorXd out;
  out.resize(v1.size() + v2.size() + v3.size());
  // Copy the data directly without dynamic memory allocation
  std::memcpy(out.data(), v1.data(), v1.size() * sizeof(double));
  std::memcpy(out.data() + v1.size(), v2.data(), v2.size() * sizeof(double));
  std::memcpy(out.data() + v1.size() + 
    v2.size(), v3.data(), v3.size() * sizeof(double));
  return out;
}


// Euclidean distance between theta_old and theta_new
double chck_convergence(std::vector<Eigen::VectorXd> & means_old,
                        std::vector<Eigen::MatrixXd> & covs_old,
                        Eigen::VectorXd probs_old,
                        std::vector<Eigen::VectorXd> & means_new,
                        std::vector<Eigen::MatrixXd> & covs_new,
                        Eigen::VectorXd probs_new) {
  
  // theta parameters old and new
  Eigen::VectorXd theta_old = collect_params(
    unlist_vecs(means_old), unlist_mats(covs_old), probs_old
  );
  Eigen::VectorXd theta_new = collect_params(
    unlist_vecs(means_new), unlist_mats(covs_new), probs_new
  );
  const int theta_dim = theta_old.size();
  double err = ((theta_old - theta_new).squaredNorm()) / theta_new.size();
  return err;
}
