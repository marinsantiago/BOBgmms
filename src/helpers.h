#ifndef HELPERS_H
#define HELPERS_H

#include <RcppEigen.h>

Eigen::VectorXd vec_vec(const std::vector<Eigen::VectorXd>& vec_of_vec);

Eigen::VectorXd vec_mat(const std::vector<Eigen::MatrixXd>& vec_of_mats);

Eigen::VectorXd collect_params(const Eigen::VectorXd &mus_vec,
                               const Eigen::VectorXd &sigmas_vec,
                               const Eigen::VectorXd &probs_vec); 

Eigen::VectorXd convertEigenVec(Rcpp::NumericVector& v1);

Eigen::MatrixXd vecToMat(const Eigen::VectorXd& v1,
                         int nrows,
                         int ncols);

Eigen::VectorXd recov_probs(const Eigen::VectorXd& theta,
                            int p,
                            int K);

std::vector<Eigen::MatrixXd>& recov_covs(const Eigen::VectorXd& theta,
                                         std::vector<Eigen::MatrixXd>& Psis,
                                         int p,
                                         int K);

std::vector<Eigen::VectorXd>& recov_means(const Eigen::VectorXd& theta,
                                          std::vector<Eigen::VectorXd>& betas,
                                          int p,
                                          int K);

Eigen::VectorXd concVecs(const std::vector<Eigen::VectorXd>& vecOfVecs);

Eigen::VectorXd repeatEntries(const Eigen::VectorXd& originalVector,
                              int K);

#endif  // HELPERS_H
