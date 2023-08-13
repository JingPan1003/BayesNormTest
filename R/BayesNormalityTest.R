#############################################################################
## Descriptions for 'BayesNormTestLocal' function:
#  This is a Bayesian hypothesis test of normality with local prior of gamma.
#############################################################################
## Arguments:
  # data:the data set for testing
  # m0: the mean of normal distribution for prior of mu
  # s0: the standard deviation of normal distribution for prior of mu
  # a0: the shape parameter of gamma distribution for prior of sigma
  # b0: the rate parameter of gamma distribution for prior of sigma
  # m0l: the mean of normal distribution for local prior of gamma
  # s0l: the standard deviation of normal distribution for local prior of gamma
  # initN: numeric vector, initial values for the parameters to be optimized 
  #         for normal model
  # init2L: numeric vector, initial values for the parameters to be optimized 
  #         for two piece model
###############################################################################

rm(list = ls())
# Load the twopiece R package
#library(devtools)
#install_github("FJRubio67/twopiece")
library(twopiece)
library(numDeriv)


BayesNormTestLocal <- function(data, m0, s0, a0, b0, m0l,s0l,initN,init2L){
  # Normal model
  lpriorn <- function(par){
    mu = par[1]; sigma = exp(par[2])
    out <- dnorm(par[1], m0, s0, log = TRUE) +
      dgamma(sigma, shape = a0, scale = 1/b0, log = TRUE) + par[2]
    return(out)
  }
  
  
  # Two piece normal model: Local prior
  lprior2L <- function(par){
    mu = par[1]; sigma = exp(par[2]);
    out <- dnorm(par[1], m0, s0, log = TRUE) +
      dgamma(sigma, shape = a0, scale = 1/b0, log = TRUE) + par[2] +
      dnorm(par[3], m0l, s0l, log = TRUE)
    return(out)
  }
  
  
  # Log posterior function: normal model
  log_postn <- function(par){
    mu = par[1]; sigma = exp(par[2])
    ll <- sum(dnorm(data, mean = mu, sd = sigma, log = TRUE))
    lprior <- lpriorn(par)
    return(-ll-lprior)
  }
  
  # Log posterior function: twopiece normal model + local prior
  log_post2L <- function(par){
    mu = par[1]; sigma = exp(par[2]); skew <- tanh(par[3])
    ll <- sum( dtp3(data, mu, sigma, skew, param = "eps", FUN = dnorm, log = TRUE))
    lprior <- lprior2L(par)
    return(-ll-lprior)
  }
  # get the optimisation
  OPTPN <- nlminb(start = initN, objective = log_postn, control = list(iter.max = 1000),
                  hessian = TRUE)
  OPTP2L <- nlminb(start = init2L, objective = log_post2L, control = list(iter.max = 1000),
                   hessian = TRUE)
  
  MAPN <- c(OPTPN$par[1], exp(OPTPN$par[2]))
  MAP2L <- c(OPTP2L$par[1], exp(OPTP2L$par[2]), tanh(OPTP2L$par[3]))
  
  # Compute hessian matrix
  hessn <- hessian(func = log_postn, x = OPTPN$par)
  hess2L <- hessian(func = log_post2L, x = OPTP2L$par)
  
  # Compute the logarithm of Laplace Approximation
  Lap2L <- 0.5*3*log(2*pi) - 0.5*log(det(hess2L)) - OPTP2L$objective
  LapN <- 0.5*3*log(2*pi) - 0.5*log(det(hessn)) - OPTPN$objective
  
  # Compute the posterior probability of Model 1
  probH1 <- 1/( 1 + exp(LapN-Lap2L) )
  # Compute the posterior probability of Model 0
  probH0 <-  exp(LapN-Lap2L)/( 1 + exp(LapN-Lap2L) )
  
  out <- list( probNorm = probH0, probTPN = probH1)
  
  return(out)
  
}

#############################################################################
## Descriptions for 'BayesNormTestNonLocal' function:
#  This is a Bayesian hypothesis test of normality with non-local prior of gamma.
#############################################################################
## Arguments:
# data:the data set applied for testing
# m0: the mean of normal distribution for prior of mu
# s0: the standard deviation of normal distribution for prior of mu
# a0: the shape parameter of gamma distribution for prior of sigma
# b0: the rate parameter of gamma distribution for prior of sigma
# m0nl: the mean of normal distribution for non-local prior of gamma
# s0nl: the standard deviation of normal distribution for non-local prior of gamma
# initN: numeric vector, initial values for the parameters to be optimized 
#         for normal model
# init2NL: numeric vector, initial values for the parameters to be optimized 
#         for two piece model
###############################################################################


BayesNormTestNonLocal <- function(data, m0, s0, a0, b0, m0nl,s0nl,initN,init2NL){
  # Normal model
  lpriorn <- function(par){
    mu = par[1]; sigma = exp(par[2])
    out <- dnorm(par[1], m0, s0, log = TRUE) + 
      dgamma(sigma, shape = a0, scale = 1/b0, log = TRUE) + par[2]
    return(out)
  }
  
  
  # Two piece normal model: Non-local prior
  lprior2NL <- function(par){
    mu = par[1]; sigma = exp(par[2]); skew <- tanh(par[3])
    out <- dnorm(par[1], m0, s0, log = TRUE) + 
      dgamma(sigma, shape = a0, scale = 1/b0, log = TRUE) + par[2] + 
      dnorm(par[3], m0nl, s0nl, log = TRUE)  + 2*log(abs(par[3])) - 2*log(s0nl)
    return(out)
  }
  
  
  
  # Log posterior function: normal model
  log_postn <- function(par){
    mu = par[1]; sigma = exp(par[2])
    ll <- sum(dnorm(data, mean = mu, sd = sigma, log = TRUE))
    lprior <- lpriorn(par)
    return(-ll-lprior)
  }
  
  # Log posterior function: twopiece normal model + nonlocal prior
  log_post2NL <- function(par){
    mu = par[1]; sigma = exp(par[2]); skew <- tanh(par[3])
    ll <- sum( dtp3(data, mu, sigma, skew, param = "eps", FUN = dnorm, log = TRUE))
    lprior <- lprior2NL(par)
    return(-ll-lprior)
  }
  # get the optimisation  
  OPTPN <- nlminb(start = initN, objective = log_postn, control = list(iter.max = 1000), 
                  hessian = TRUE)
  OPTP2NL <- nlminb(start = init2NL, objective = log_post2NL, control = list(iter.max = 1000), 
                    hessian = TRUE)
  
  MAPN <- c(OPTPN$par[1], exp(OPTPN$par[2]))
  MAP2NL <- c(OPTP2NL$par[1], exp(OPTP2NL$par[2]), tanh(OPTP2NL$par[3]))
  
  # Compute hessian matrix  
  hessn <- hessian(func = log_postn, x = OPTPN$par)
  hess2NL <- hessian(func = log_post2NL, x = OPTP2NL$par)
  
  # Compute the logarithm of Laplace Approximation  
  Lap2NL <- 0.5*3*log(2*pi) - 0.5*log(det(hess2NL)) - OPTP2NL$objective
  LapN <- 0.5*3*log(2*pi) - 0.5*log(det(hessn)) - OPTPN$objective
  
  # Compute the posterior probability of Model 1  
  probH1 <- 1/( 1 + exp(LapN-Lap2NL) )
  # Compute the posterior probability of Model 0  
  probH0 <-  exp(LapN-Lap2NL)/( 1 + exp(LapN-Lap2NL) )
  
  out <- list( NLprobNorm = probH0, NLprobTPN = probH1)
  
  return(out)
  
}

