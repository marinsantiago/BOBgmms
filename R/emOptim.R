em.optim <- function(y, means, covs, probs, a, b, c, r,
                    lik_weights, mean_weights, cov_weights, prob_weights,
                    alphas, betas, lambdas, nus, psis,
                    max_iters = 500, epss = 1e-4){

  em_optim(y, means, covs, probs, a, b, c, r, lik_weights, mean_weights,
           cov_weights, prob_weights,  alphas, betas, lambdas, nus, psis,
           max_iters, epss)
}
