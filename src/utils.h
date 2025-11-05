#ifndef UTILS_H
#define UTILS_H

#include <RcppEigen.h>


Eigen::VectorXd unlist_vecs(const std::vector<Eigen::VectorXd> & list_of_vecs);

Eigen::VectorXd unlist_mats(const std::vector<Eigen::MatrixXd> & list_of_mats);

Eigen::VectorXd collect_params(const Eigen::VectorXd & v1,
                               const Eigen::VectorXd & v2,
                               const Eigen::VectorXd & v3);

double chck_convergence(std::vector<Eigen::VectorXd> & means_old,
                        std::vector<Eigen::MatrixXd> & covs_old,
                        Eigen::VectorXd probs_old,
                        std::vector<Eigen::VectorXd> & means_new,
                        std::vector<Eigen::MatrixXd> & covs_new,
                        Eigen::VectorXd probs_new);

#endif  // UTILS_H
