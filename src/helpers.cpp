#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include "helpers.h"


// vec-vec concatenation:
// Function that takes a vector of vectors and returns the concatenated vector.
Eigen::VectorXd vec_vec(const std::vector<Eigen::VectorXd>& vec_of_vec) {

  // Compute the total size of the resulting vector
  int totalSize = 0;
  for (const auto& vec : vec_of_vec) {
    totalSize += vec.size();
  }

  // Create a single vector by concatenating all the vectors
  Eigen::VectorXd out(totalSize);
  int Indexx = 0;
  for (const auto& vec : vec_of_vec) {
    Eigen::Map<Eigen::VectorXd>(out.data() + Indexx, vec.size()) = vec;
    Indexx += vec.size();
  }

  return out;
}


// vec-mat concatenation:
// Function that takes a vector of matrices and returns the concatenated vector.
Eigen::VectorXd vec_mat(const std::vector<Eigen::MatrixXd>& vec_of_mats) {
  // Compute the total size of the resulting vector
  int totalSize = 0;
  for (const auto& matrix : vec_of_mats) {
    totalSize += matrix.size();
  }

  // Create a single vector by concatenating the matrices
  Eigen::VectorXd out(totalSize);
  int currentIndex = 0;
  for (const auto& matrix : vec_of_mats) {
    Eigen::Map<Eigen::VectorXd>(out.data() + currentIndex, matrix.size()) =
      Eigen::Map<const Eigen::VectorXd>(matrix.data(), matrix.size());
    currentIndex += matrix.size();
  }

  return out;
}


Eigen::VectorXd collect_params(const Eigen::VectorXd& v1,
                               const Eigen::VectorXd& v2,
                               const Eigen::VectorXd& v3) {

  Eigen::VectorXd out;
  out.resize(v1.size() + v2.size() + v3.size());

  // Copy the data directly without dynamic memory allocation
  std::memcpy(out.data(), v1.data(), v1.size() * sizeof(double));
  std::memcpy(out.data() + v1.size(), v2.data(), v2.size() * sizeof(double));
  std::memcpy(out.data() + v1.size() +
    v2.size(), v3.data(), v3.size() * sizeof(double));

  return out;
}


// Function that converts Rcpp::NumericVector to Eigen::VectorXd
Eigen::VectorXd convertEigenVec(Rcpp::NumericVector& v1) {

  // Use Rcpp::as to convert NumericVector to Eigen::Map
  Eigen::Map<Eigen::VectorXd> out(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(v1));

  // Return the Eigen::VectorXd
  return out;
}


// Function that converts Eigen::VectorXd into Eigen::MatrixXd
Eigen::MatrixXd vecToMat(const Eigen::VectorXd& v1,
                         int nrows,
                         int ncols) {

  return Eigen::Map<const Eigen::MatrixXd>(v1.data(), nrows, ncols);
}


Eigen::VectorXd recov_probs(const Eigen::VectorXd& theta,
                            int p,
                            int K) {

  return theta.tail(K);

}


std::vector<Eigen::MatrixXd>& recov_covs(const Eigen::VectorXd& theta,
                                         std::vector<Eigen::MatrixXd>& Psis,
                                         int p,
                                         int K) {

  std::vector<Eigen::MatrixXd>& out = Psis;

  for(int k = 0; k < K; ++k) {

    const int init = K * p + (((k + 1) - 1) * std::pow(p, 2));

    Eigen::VectorXd sigma_vec = theta.segment(init, p * p);
    out[k] = vecToMat(sigma_vec, p, p);

  }

  return out;
}


std::vector<Eigen::VectorXd>& recov_means(const Eigen::VectorXd& theta,
                                          std::vector<Eigen::VectorXd>& betas,
                                          int p,
                                          int K) {

  std::vector<Eigen::VectorXd>& out = betas;

  for(int k = 0; k < K; ++k) {

    const int init = p * (k);
    out[k]         = theta.segment(init, p);

  }

  return out;
}


// Concatenate a vector of vectors into a single long vector
Eigen::VectorXd concVecs(const std::vector<Eigen::VectorXd>& vecOfVecs) {
  // Calculate the total size needed for the concatenated vector
  int totalSize = 0;
  for (const auto& vector : vecOfVecs) {
    totalSize += vector.size();
  }

  // Create the concatenated vector
  Eigen::VectorXd concatenatedVector(totalSize);

  // Copy the data from individual vectors to the concatenated vector
  int currentIndex = 0;
  for (const auto& vector : vecOfVecs) {
    concatenatedVector.segment(currentIndex, vector.size()) = vector;
    currentIndex += vector.size();
  }

  return concatenatedVector;
}

// Each entry of the vector to be repeated K times sequentially
Eigen::VectorXd repeatEntries(const Eigen::VectorXd& originalVector,
                              int K) {
  // Calculate the size of the repeated vector
  int repeatedSize = originalVector.size() * K;

  // Create a vector to store the repeated values
  Eigen::VectorXd repeatedVector(repeatedSize);

  // Populate the repeated vector
  for (int i = 0; i < originalVector.size(); ++i) {
    for (int j = 0; j < K; ++j) {
      repeatedVector[i * K + j] = originalVector[i];
    }
  }

  return repeatedVector;
}
