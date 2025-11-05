## Analytical expression of the Expected Improvement criterion
## 
## Computes the Expected Improvement at current location. The current minimum
## of the observations can be replaced by an arbitrary value (plugin), which
## is usefull in particular in noisy frameworks.
## 
## 
## @param x a vector representing the input for which one wishes to calculate
## EI,
## @param model an object of class \code{\link[DiceKriging]{km}},
## @param plugin optional scalar: if provided, it replaces the minimum of the
## current observations,
## @param type "SK" or "UK" (by default), depending whether uncertainty
## related to trend estimation has to be taken into account,
## @param minimization logical specifying if EI is used in minimiziation or in
## maximization,
## @param envir an optional environment specifying where to assign
## intermediate values for future gradient calculations. Default is NULL.
## @param proxy an optional Boolean, if TRUE EI is replaced by the kriging mean (to minimize)
## @return The expected improvement, defined as \deqn{EI(x) := E[( min(Y(X)) -
## Y(x))^{+} | Y(X)=y(X)],} where X is the current design of experiments and Y
## is the random process assumed to have generated the objective function y.
## If a plugin is specified, it replaces \deqn{min(Y(X))} in the previous
## formula.
## @author David Ginsbourger 
## 
## Olivier Roustant 
## 
## Victor Picheny 
##
## @references
## 
## D.R. Jones, M. Schonlau, and W.J. Welch (1998), Efficient global
## optimization of expensive black-box functions, \emph{Journal of Global
## Optimization}, 13, 455-492.
## 
## T.J. Santner, B.J. Williams, and W.J. Notz (2003), \emph{The design and
## analysis of computer experiments}, Springer.
## 
## M. Schonlau (1997), \emph{Computer experiments and global optimization},
## Ph.D. thesis, University of Waterloo.
##
EI <- function (x, model, plugin = NULL, type = "UK", 
                minimization = TRUE, envir = NULL, proxy = FALSE) {
  
  if (is.null(plugin)){ 
    if (minimization) {
      plugin <- min(model@y)
    } else {
      plugin <- -max(model@y)
    }
  }
  m <- plugin
  
  ##############################################################################
  # Convert x in proper format(s)
  if (is.data.frame(x)) {
    d <- length(x)
    if (d != model@d){ stop("x does not have the right size") }
    newdata.num <- as.numeric(x)
    newdata <- data.frame(t(newdata.num))
  } else {
    if (is.null(dim(x))) {
      d <- length(x)
      if (d != model@d){ stop("x does not have the right size") }
      newdata <- data.frame(t(as.numeric(x)))
    } else {
      d <- ncol(x)
      if (d != model@d){ stop("x does not have the right size") }
      newdata <- data.frame(x)
    }
  }
  colnames(newdata) <- colnames(model@X)
  
  ##############################################################################
  if (!is.null(envir)) {
    predx <- predict(
      object=model, newdata=newdata, type=type, checkNames = FALSE
    )
  } else {
    predx <- predict(
      object=model, newdata=newdata, type=type, 
      checkNames = FALSE, light.return=TRUE
    )
  }
  
  kriging.mean <- predx$mean
  if(!minimization) {
    kriging.mean <- -kriging.mean
  }
  kriging.sd   <- predx$sd
  
  if (proxy) {
    xcr <- xcr.prob <- xcr.dens <- NULL
    res <- m - kriging.mean
  } else {
    xcr <- (m - kriging.mean) / kriging.sd
    
    xcr.prob <- pnorm(xcr)
    xcr.dens <- dnorm(xcr)	        
    res <- (m - kriging.mean) * xcr.prob + kriging.sd * xcr.dens
  }
  
  ##############################################################################
  if (!is.null(envir)) 
  { assign("xcr", xcr, envir=envir)
    assign("xcr.prob", xcr.prob, envir=envir)
    assign("xcr.dens", xcr.dens, envir=envir)
    assign("kriging.sd", kriging.sd, envir=envir)
    assign("c", predx$c, envir=envir)
    assign("Tinv.c", predx$Tinv.c, envir=envir)
  }
  return(res)
}