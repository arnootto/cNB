#' Maximum likelihood estimation of contaminated negative binomial regression
#'
#' Maximum likelihood estimation of the contaminated negative binomial regression model via an Expectation-Maximization algorithm.
#'
#' @param formula an object of class 'formula': a symbolic description of the model to be fitted. 
#' @param data a mandatory data frame containing the variables in the model.
#' @param start vector of initial values. If NULL, then values produced by glm.nb are used as initial values with delta=0.05 and eta=1.1.
#' @param method optimization method to be used. Default is "BFGS". Other options are "Nelder-Mead".
#' @param reltol relative convergence tolerance in each optimization process. Defaults to 1e-15.
#' @param maxit the number of inner iterations in each optimization process. Defaults to 10000.
#' @param em.tol  the EM convergence tolerance. Defaults to 1e-10.
#' @param em.maxit  the number of EM iterations. Defaults to 1000.
#'
#' @details The \code{ml.cnb} function fits the contaminated negative binomial regression model (see Otto et. al (2024)). 
#'
#' @return An list of elements:
#'    \item{results}{A data frame with parameter estimates and standard errors. }
#'    \item{alpha}{Maximum likelihood estimate of alpha.}
#'    \item{delta}{Maximum likelihood estimate of delta.}
#'    \item{beta}{Maximum likelihood estimates of the regression coefficients.}
#'    \item{mu}{Minimum likelihood estimates of the mean paramter.}
#'    \item{X}{The design matrix used in the model.}
#'    \item{y}{The response variable.}
#'    \item{loglike}{The log-likelihood value at convergence.}
#'    \item{AIC}{Akaike Information Criterion (AIC) for the fitted model.}
#'    \item{BIC}{Bayesian Information Criterion (BIC) for the fitted model.}
#'    \item{LRpvalue}{p-value of Likelihood Ratio test.}
#'
#' @import MASS
#' 
#' 
#'
#' @export 

ml.cnb <- function(formula, data, start = NULL, method = "BFGS", reltol=1e-15, maxit=10000, em.tol=1e-5, em.maxit=1000) {
  mf <- model.frame(formula, data)
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  nb2X <- model.matrix(formula, data = data)  #design matrix
  nb2.reg.ml <- function(b.hat, X, y,delta, v) {
    alpha.hat <- exp(b.hat[1]) #alpha
    eta.hat <- exp(b.hat[2])+1 #degree of contamination
    xb.hat <- X %*% b.hat[-c(1:2)] #coefficients
    mu.hat <- exp(xb.hat) #link
    ll <- sum(ll(y,mu=mu.hat, alpha = alpha.hat, eta=eta.hat,v=v))
    #print(c(alpha.hat,delta,eta.hat,ll))
    return(ll)
  }
  if (is.null(start)) {#initial param
    nb <-MASS::glm.nb(formula = formula, data = data, control = glm.control(maxit=25))
    start <- matrix(NA,(ncol(mf)+2),1)
    start[1] <- 1/nb$theta #alpha
    start[2] <- 2 # eta
    start[3:(ncol(mf)+2)] <- nb$coefficients
    xb.start <- nb2X %*% start[-c(1:2)] #coefficients
    mu.start <- exp(xb.start) #link
    nbloglik <- logLik(nb)[1] #loglikelihood of nb
    nbterms <- length(coefficients(nb)) #number of parameters of nb
  }
  delta=0.05  
  loglik    <- NULL
  cond=F
  i=1
  loglik[1]<-1000000000
  while(cond==F){
    i=i+1
    v_i = delta * dnbinom(y, size = 1/(start[1]*start[2]), mu = mu.start, log = F)/dcnbinom(y,mu= mu.start, alpha=start[1], delta=delta, eta=start[2], log=F)
    delta=min(0.5,mean(v_i))
    start[1] <- log(start[1]) #alpha
    start[2] <- log(start[2]-1) #eta
    fit <- optim(par = start,
                 fn = nb2.reg.ml,
                 X = nb2X,
                 y = y,
                 method = method,
                 control = list(
                   fnscale = -1,
                   maxit = maxit,
                   reltol = reltol),
                 hessian = F,
                 delta = delta,
                 v = v_i
    )
    beta.hat <- fit$par
    start[1] <- exp(beta.hat[1])
    start[2] <- exp(beta.hat[2])+1
    start[3:(ncol(mf)+2)] <- beta.hat[-c(1:2)]
    beta <- beta.hat[(3:length(beta.hat))]
    mu.start <- exp(nb2X %*% beta)
    lc=sum(dcnbinom(y, mu=mu.start, alpha = start[1], delta=delta, eta=start[2], log = T))
    loglik[i]=lc
    #print(i)
    if(abs(loglik[i]-loglik[i-1])<em.tol||i>=em.maxit){
      cond=T
      
    }
  }
  
  beta.hat <- fit$par
  beta.hat[1] <- exp(beta.hat[1])
  beta.hat[2] <- exp(beta.hat[2])+1
  
  par <-c(beta.hat,delta)
  dcnbinomHess<-function(X, par, y,log=F){
    alpha <- par[1] #alpha
    eta <- par[2]#degree of contamination
    xb.hat <- X %*% par[-c(1:2,length(par))] #coefficients
    mu <- exp(xb.hat) #link
    delta <- par[length(par)]
    
    log_term1 = log(1-delta) + dnbinom(y, size = 1/alpha, mu = mu, log = T)
    log_term2 = log(delta) + dnbinom(y, size = 1/(alpha*eta), mu = mu, log = T)
    d = (exp(log_term1) + exp(log_term2))
    if (log == F) {
      return(sum(d))
    } else {
      return(sum(log(d)))
    }
  }
  Hess=optimHess(par=par,fn=dcnbinomHess,log=T, X=nb2X, y=y,control = list(
    fnscale = -1,
    maxit = maxit,
    reltol = reltol))
  se.beta.hat <- sqrt(diag(solve(-Hess)))
  
  results <- data.frame(Estimate = cbind(c(beta.hat,delta),se.beta.hat))
  rownames(results) <- c("alpha", "eta", colnames(nb2X),"delta")
  colnames(results) <- c("estimates", "se")
  alpha <- beta.hat[1];  eta <- beta.hat[2]
  beta <- beta.hat[(3:length(beta.hat))]; delta <- delta
  mu <- exp(nb2X %*% beta)
  results <- results[c(3:(nrow(results)-1), 1,2,nrow(results)),,drop=F]
  AIC <- -2*lc+nrow(results)*2
  BIC <- -2*lc+nrow(results)*log(nrow(data))
  LR <- 1-pchisq(-2*(nbloglik-lc), nrow(results)-nbterms)
  return(list(
    results = results,
    alpha =alpha,
    delta = delta,
    eta = eta,
    beta = beta,
    mu = mu,
    X = nb2X,
    y = y,
    loglike = lc,
    AIC = AIC,
    BIC = BIC,
    LRpvalue =LR
  ))
}



ll<-function(x, mu, alpha, delta, eta, v){
  (1-v)*(lgamma(x+1/alpha)-lgamma(x+1)-lgamma(1/alpha)+x*log(alpha*mu)-x*log(1+alpha*mu)-1/alpha*log(1+alpha*mu))+(v)*(lgamma(x+1/(eta*alpha))-lgamma(x+1)-lgamma(1/(eta*alpha))+x*log(eta*alpha*mu)-x*log(1+eta*alpha*mu)-1/(eta*alpha)*log(1+eta*alpha*mu))
}

