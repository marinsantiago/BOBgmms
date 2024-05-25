#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include "helpers.h"
#include "posterior_pred.h"


// [[Rcpp::export]]
Eigen::MatrixXd post_pred(const Eigen::MatrixXd& post_draws,
                          std::vector<Eigen::VectorXd> betas,
                          std::vector<Eigen::MatrixXd>& Psis) {
  
 const int S = post_draws.rows();
 const int K = betas.size();
 const int p = betas[0].size();
 
 
 // Set up environments
 Rcpp::Environment base_pkg = Rcpp::Environment::namespace_env("base");
 Rcpp::Environment MASS_pkg = Rcpp::Environment::namespace_env("MASS");
 Rcpp::Function sample      = base_pkg["sample"];
 Rcpp::Function mvrnorm     = MASS_pkg["mvrnorm"];
 
 // Prepare output
 Eigen::MatrixXd out(S, p);
 
 for(int s = 0; s < S; ++s) {
   
   // Extract model parameters
   Eigen::VectorXd theta               = post_draws.row(s);
   std::vector<Eigen::VectorXd> mus    = recov_means(theta, betas, p, K);
   std::vector<Eigen::MatrixXd> sigmas = recov_covs(theta, Psis, p, K);
   Eigen::VectorXd probs               = recov_probs(theta, p, K);
   
   // Sample latent cluster variables
   int z = Rcpp::as<int>(sample(Eigen::VectorXd::LinSpaced(K, 1, K),
                                1, 
                                Rcpp::Named("prob", probs)));
   
   // Adjust z, because C++ starts counting at zero
   z -= 1;
   
   // Sample from the corresponding multivariate normal
   out.row(s) = Rcpp::as<Eigen::VectorXd>(mvrnorm(1, mus[z], sigmas[z]));
   
 }
 
 return out;  
}
