bayes_optim <- function(lower_bound, upper_bound, y, K, objective,
                        tot_iterations, init_samples = 30) {
  
  # Prepare output 
  total_evals <- init_samples + 3 + tot_iterations
  
  # Latin hyper-cube initialization --------------------------------------------
  cat("Latin hyper-cube initialization: \n")
  x_evals <- matrix(1, nrow = total_evals, ncol = 1)
  y_evals <- rep(999, total_evals)
  # The first row of x corresponds to the vector of ones (i.e., wbb)
  # The second and third rows are the boundaries of the search space
  x_evals[2,] <- lower_bound
  x_evals[3,] <- upper_bound
  x_evals[4:(init_samples + 3),] <- lhs_sample(
    samples = init_samples, lower_bound = lower_bound, 
    upper_bound = upper_bound
  )
  # Evaluate the loss at the initial points
  y_evals[1:(init_samples + 3)] <- pbapply::pblapply(
    x_evals[1:(init_samples + 3),], objective
  ) |> unlist()
  # Extract optimal x and loss
  best_index <- which.min(y_evals[1:(init_samples + 3)])
  best_loss <- y_evals[best_index]
  best_x <- x_evals[best_index,]
  cat("Best initial loss:", best_loss, "\n")
  cat("Best initial input:", best_x, "\n")
  
  # Bayesian optimizer ---------------------------------------------------------
  bo.loop <- (init_samples + 4):total_evals
  bo.iter <- 0
  for (i in bo.loop) {
    bo.iter <- bo.iter + 1
    # Get the candidate point for evaluation
    candidate_x <- tryCatch(
      {
        next_evals(
          evals = as.matrix(x_evals[seq_len(i - 1),]),
          outcomes = y_evals[seq_len(i - 1)], 
          lower_bound = lower_bound, upper_bound = upper_bound
        ) |> suppressWarnings()
      }, error = \(e) {
        lhs_augment(
          evals = as.matrix(x_evals[seq_len(i - 1),]), 
          lower_bound = lower_bound, upper_bound = upper_bound
        )
      }
    ) |> c()
    # Make sure that "candidate_x" is not in "x_evals"
    if (sum(abs(x_evals[seq_len(i - 1),] - candidate_x) < 1e-10) > 0) {
      candidate_x <- lhs_augment(
        evals = as.matrix(x_evals[seq_len(i - 1),]), 
        lower_bound = lower_bound, upper_bound = upper_bound
      )
    }
    # Make sure that "candidate_x" is a real number
    if (!is.finite(candidate_x)) {
      candidate_x <- lhs_augment(
        evals = as.matrix(x_evals[seq_len(i - 1),]), 
        lower_bound = lower_bound, upper_bound = upper_bound
      )
    }
    # Evaluate the loss at the candidate points
    candidates_loss <- objective(candidate_x)
    # Update design matrix and outcomes vector
    x_evals[i,] <- candidate_x
    y_evals[i] <- candidates_loss
    # Update best current points
    if (candidates_loss < best_loss) {
      best_loss <- candidates_loss        # Update best loss
      best_x <- candidate_x               # Update best input
    }
    cat("--------------------------------------------\n")
    cat("BO iteration:", bo.iter, "\n")
    cat("Current Loss:", candidates_loss, "\n")
    cat("Best Loss:", best_loss, "\n")
    cat("--------------------------------------------\n")
  }
  cat("--------------------------------------------\n")
  cat("BO Summary:", "\n")
  cat("Best Loss:", best_loss, "\n")
  cat("Best input:", best_x, "\n")
  cat("--------------------------------------------\n")
  
  best_x
}


# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------


# Generate a latin hyper-cube design
lhs_sample <- function(samples, lower_bound, upper_bound) {
  # The dimension of the search space is 1
  lhs::randomLHS(samples, 1) |> qunif(p = _, lower_bound, upper_bound)
}


# Augment the existing latin hyper-cube design
lhs_augment <- function(evals, lower_bound, upper_bound) {
  # Make sure evals are in (0, 1)
  stdz.evals <- punif(evals, lower_bound, upper_bound)
  # Augment the evaluations
  out <- utils::tail(lhs::augmentLHS(stdz.evals, 1), 1, keepnums = FALSE)
  # Return the evaluations to the original scale
  qunif(out, lower_bound, upper_bound)
}


# Next points for evaluation
next_evals <- function(evals, outcomes, lower_bound, upper_bound) {
  # Fit the GP model -----------------------------------------------------------
  gp <- DiceKriging::km(
    formula = ~ 1, design = evals, response = outcomes, covtype = "matern5_2",
    multistart = 1, control = list(trace = FALSE, pop.size = 20)
  )
  
  # Define the acquisition function --------------------------------------------
  EI_fun <- \(x) DiceOptim::EI(x, gp) # Expected Improvement
  
  # Optimize the acquisition function ------------------------------------------
  EI.x <- seq(lower_bound, upper_bound, by = 0.02)
  EI.y <- sapply(EI.x, \(r) EI_fun(r))
  EI.x[which.max(EI.y)]
}
