lhsSample <- function(samples, dim, lowerBound, upperBound){

  x <- lhs::randomLHS(samples, dim)
  t(apply(x, 1, \(r) qunif(r, lowerBound, upperBound)))
}


nextEvals <- function(evals, outcomes, lBound, uBound){

  # Fit a GP model
  gp <- DiceKriging::km(~1, design = evals, response = outcomes,
                        covtype = "matern5_2", multistart = 1,
                        control = list(trace = FALSE, pop.size = 20))

  # Optimize the acquisition function
  opt.acq <- DiceOptim::max_qEI(gp, npoints = 2L, lower = lBound,
              upper = uBound, crit = "exact", minimization = TRUE,
                optimcontrol = list(nStarts = 1, methods = "BFGS"))

  as.matrix(opt.acq$par)
}


lhsAugment <- function(evals, lBound, uBound){

  std.evals <- t(apply(evals, 1, \(r) punif(r, lBound, uBound)))

  out <- tryCatch(
    {
      tail(lhs::augmentLHS(std.evals, 2), 2, keepnums = FALSE)
    }, error = function(e){
      matrix(runif(ncol(evals) * 2), nrow = 2)
    }
  )

  t(apply(out, 1, \(r) qunif(r, lBound, uBound)))
}


bayesOpt <- function(lowerBound, upperBound, y, K, objective,
                     tot_iterations, init_samples = 20 + (5 * K)){

  dim_bayesOpt <- 2 * (K + 1)
  totalEvals   <- init_samples + 3 + (tot_iterations * 2)

  # Initial Random Latin Hyper-cube design
  x_evals <- matrix(1, totalEvals, dim_bayesOpt)
  y_evals <- rep(999, totalEvals)

  cat("Latin Hypercube Initialization: \n")

  # The first row corresponds to the vector of ones
  # The second and third rows are the boundaries of the search space
  x_evals[2,] <- lowerBound
  x_evals[3,] <- upperBound
  x_evals[4:(init_samples + 3),] <- lhsSample(init_samples, dim_bayesOpt,
                                              lowerBound, upperBound)

  # Evaluate the loss at the initial points
  y_evals[1:(init_samples + 3)] <- pbapply::pbapply(
    x_evals[1:(init_samples + 3),], 1, objective
  )

  # Extract optimal x and loss
  best_index <- which.min(y_evals[1:(init_samples + 3)])
  best_loss  <- y_evals[best_index]
  best_x     <- x_evals[best_index,]

  cat("Best Initial Loss:", best_loss, "\n")
  cat("Best Initial Input:", best_x, "\n")

  # Bayesian Optimizer
  bo.loop <- seq(from = init_samples + 4, to = totalEvals, by = 2)
  bo.iter <- 0
  for (i in bo.loop) {

    bo.iter <- bo.iter + 1

    # Get the candidate point for evaluation
    candidates_x <- tryCatch(
      {
        suppressWarnings(nextEvals(x_evals[seq_len(i - 1),],
                                   y_evals[seq_len(i - 1) ],
                                   lowerBound, upperBound))
      }, error = function(e){
        lhsAugment(x_evals[seq_len(i - 1),], lowerBound, upperBound)
      }
    )

    #cat("BO: Done \n")

    # Evaluate the candidate
    candidates_loss <- apply(candidates_x, 1, objective)
    best_cand_index <- which.min(candidates_loss)
    best_cand_loss  <- candidates_loss[best_cand_index]
    best_cand_input <- c(candidates_x[best_cand_index,])

    # Update design matrix and outcomes vector
    x_evals[i:(i + 1),] <- candidates_x
    y_evals[i:(i + 1)]  <- candidates_loss

    if (best_cand_loss < best_loss) {
      best_loss <- best_cand_loss      # Update best loss
      best_x    <- best_cand_input     # Update best input
    }

    cat("---------------------------------------\n")
    cat("BO iteration:", bo.iter,              "\n")
    cat("Current Loss:", best_cand_loss,       "\n")
    cat("Best Loss:", best_loss,               "\n")
    cat("---------------------------------------\n")

  }

  cat("---------------------------------------\n")
  cat("BO Summary:",                         "\n")
  cat("Best Loss:", best_loss,               "\n")
  cat("Best input:", best_x,                 "\n")
  cat("---------------------------------------\n")

  return(best_x)
}
