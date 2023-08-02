#############################################################################
## Descriptions for 'BayesNormTestLocal' function:
#  This is a Bayesian hypothesis test of normality with local prior of gamma.
#############################################################################
## Arguments:
  # data:the data set used for testing
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

  OPTPN <- nlminb(start = initN, objective = log_postn, control = list(iter.max = 1000),
                  hessian = TRUE)
  OPTP2L <- nlminb(start = init2L, objective = log_post2L, control = list(iter.max = 1000),
                   hessian = TRUE)

  MAPN <- c(OPTPN$par[1], exp(OPTPN$par[2]))
  MAP2L <- c(OPTP2L$par[1], exp(OPTP2L$par[2]), tanh(OPTP2L$par[3]))


  hessn <- hessian(func = log_postn, x = OPTPN$par)
  hess2L <- hessian(func = log_post2L, x = OPTP2L$par)

  Lap2L <- 0.5*3*log(2*pi) - 0.5*log(det(hess2L)) - OPTP2L$objective
  LapN <- 0.5*3*log(2*pi) - 0.5*log(det(hessn)) - OPTPN$objective

  probH1 <- 1/( 1 + exp(LapN-Lap2L) )

  probH0 <-  exp(LapN-Lap2L)/( 1 + exp(LapN-Lap2L) )

  out <- list( probNorm = probH0, probTPN = probH1)

  return(out)

}

#############################################################################
## Descriptions for 'BayesNormTestLocal' function:
#  This is a Bayesian hypothesis test of normality with non-local prior of gamma.
#############################################################################
## Arguments:
# data:the dataa set applied for testing
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

  OPTPN <- nlminb(start = initN, objective = log_postn, control = list(iter.max = 1000),
                  hessian = TRUE)
  OPTP2NL <- nlminb(start = init2NL, objective = log_post2NL, control = list(iter.max = 1000),
                    hessian = TRUE)

  MAPN <- c(OPTPN$par[1], exp(OPTPN$par[2]))
  MAP2NL <- c(OPTP2NL$par[1], exp(OPTP2NL$par[2]), tanh(OPTP2NL$par[3]))


  hessn <- hessian(func = log_postn, x = OPTPN$par)
  hess2NL <- hessian(func = log_post2NL, x = OPTP2NL$par)

  Lap2NL <- 0.5*3*log(2*pi) - 0.5*log(det(hess2NL)) - OPTP2NL$objective
  LapN <- 0.5*3*log(2*pi) - 0.5*log(det(hessn)) - OPTPN$objective

  probH1 <- 1/( 1 + exp(LapN-Lap2NL) )

  probH0 <-  exp(LapN-Lap2NL)/( 1 + exp(LapN-Lap2NL) )

  out <- list( probNorm = probH0, probTPN = probH1)

  return(out)

}

