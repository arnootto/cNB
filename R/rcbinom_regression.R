#' Generate Negative Binomial Regression Data
#'
#' This function generates synthetic data for a negative binomial regression model.
#'
#' @param n Number of observations to generate. Default is 1.
#' @param off Offset term. Default is 0.
#' @param alpha Dispersion parameter. Default is 0.5.
#' @param xv Vector of regression coefficients. Default is c(1, 0.75, -1.5).
#'
#' @return A data frame containing the generated negative binomial response variable and the predictor variables.
#'
#' @import MASS
#' 
#' @export 
 rnbinom_regression <- function(n = 1, off = 0,
                               alpha = 0.5,
                               xv = c(1, 0.75, -1.5)) {
  p <- length(xv) - 1
  X <- cbind(1, matrix(rbinom(n * p, 1, 0.5), ncol = p))
  
  xb <- X %*% xv
  r <- 1 / alpha
  exb <- exp(xb + off)
  nby <- rnegbin(n, mu = exb, theta = r)
  out <- data.frame(cbind(nby, X[, -1, drop = FALSE]))
  names(out) <- c("nby", paste("x", 1:p, sep = ""))
  return(out)
}

#' Generate Contaminated Negative Binomial Regression Data
#'
#' This function generates synthetic data for a contaminated negative binomial regression model.
#'
#' @param n Number of observations to generate. Default is 50000.
#' @param off Offset term. Default is 0.
#' @param alpha Dispersion parameter. Default is 0.5.
#' @param delta Proportion of good data points. Default is 0.95.
#' @param eta Inflation parameter for bad data points. Default is 2.
#' @param xv Vector of regression coefficients. Default is c(1, 0.75, -1.5).
#'
#' @return A data frame containing the generated contaminated negative binomial response variable and the predictor variables.
#'
#' @import MASS
#' 
#' @export 
rcnbinom_regression <- function(n = 50000, off = 0,
                                alpha = 0.5,
                                delta = 0.95,
                                eta = 2,
                                xv = c(1, 0.75, -1.5)) {
  p <- length(xv) - 1
  bernoulli <- rbinom(n, 1, delta)
  cnby_good <- rnbinom_regression(n = sum(bernoulli == 1), off = off, alpha = alpha, xv = xv)
  if (sum(bernoulli == 0) != 0) {
    cnby_bad <- rnbinom_regression(n = sum(bernoulli == 0), off = off, alpha = alpha * eta, xv = xv)
    cnby <- rbind(cnby_good, cnby_bad)
  } else {
    cnby <- cnby_good
  }
  
  names(cnby) <- c("cnby", paste("x", 1:p, sep = ""))
  return(cnby)
}