#' contaminated Negative Binomial density
#'
#' Density function of contaminated negative binomial distribution
#'
#' @param x vector of (non negative integer) quantities
#' @param mu mean parameter. mu > 0
#' @param alpha dispersion parameter. alpha > 0.
#' @param delta proportion of mild outliers (bad points). 0 < delta < 1.
#' @param eta inflation parameter. eta > 1.
#' @param log loglical; if True probabilities p are given as log(p). Default is False.
#' 
#' @details The contaminated negative binomial distribution with mu=\eqn{\mu}, alpha=\eqn{\alpha}, delta=\eqn{\delta}, and eta=\eqn{\eta} has density
#' \deqn{
#' f(x;\mu,\alpha,\delta,\eta)=(1-\delta)\binom{x-1/\alpha-1}{x}(\frac{\alpha\mu}{1+\alpha\mu})^x(\frac{1}{1+\alpha\mu})^{1/\alpha}+\delta\binom{x-1/(\eta\alpha)-1}{x}(\frac{\eta\alpha\mu}{1+\eta\alpha\mu})^x(\frac{1}{1+\eta\alpha\mu})^{1/(\eta\alpha)}
#' }
#' for, \eqn{x=0,1,\dots}, \eqn{\mu>0, \alpha>0, 0<\delta<1} and \eqn{\eta>1}. 
#'
#' @return An list of elements:
#'    \item{d}{The density of the contaminated negative binomial distribution }
#'    
#'
#' @import stats
#'
#' @export 
dcnbinom<-function(x, mu, alpha, delta, eta, log=F){
  log_term1 = log(1-delta) + dnbinom(x, size = 1/alpha, mu = mu, log = T)
  log_term2 = log(delta) + dnbinom(x, size = 1/(alpha*eta), mu = mu, log = T)
  d = (exp(log_term1) + exp(log_term2))
  if (log == F) {
    return(d)
  } else {
    return(log(d))
  }
}