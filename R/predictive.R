post.predictive <- function(post.draws, d, K){

  list.vectors  <- lapply(1:K, function(m) rep(0, d))
  list.matrices <- lapply(1:K, function(m) diag(d))
  out           <- post_pred(post.draws, list.vectors, list.matrices)

  out
}
