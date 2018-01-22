#---------------
# CMP functions
#---------------

# Unnormalised CMP density function
logDensCMP1 <- function(y, gamma, nu, useC = FALSE) {
  ## 'y'
  stopifnot(round(y)==y)
  stopifnot(identical(length(y), 1L))
  stopifnot(!is.na(y))
  stopifnot(y >= 0L)
  ## 'gamma'
  stopifnot(is.double(gamma))
  stopifnot(identical(length(gamma), 1L))
  stopifnot(!is.na(gamma))
  stopifnot(gamma >= 0L)
  ## 'nu'
  stopifnot(is.double(nu))
  stopifnot(identical(length(nu), 1L))
  stopifnot(!is.na(nu))
  stopifnot(nu >= 0L)
  gamma <- ifelse(gamma == 0, exp(-100), gamma)
  if (useC) {
    .Call(logDensCMP1_R, y, gamma, nu)
  }
  else {
    nu * (y * log(gamma) - lgamma(y + 1))
  }
}

# Function to generate underdispersed data
rcmpUnder <- function(mu, nu, max){
  ## 'mu'
  stopifnot(is.double(mu))
  stopifnot(identical(length(mu), 1L))
  stopifnot(!is.na(mu))
  stopifnot(mu >= 0)
  ## 'nu'
  stopifnot(is.double(nu))
  stopifnot(identical(length(nu), 1L))
  stopifnot(!is.na(nu))
  stopifnot(nu > 0)
  ## 'max'
  stopifnot(is.double(max))
  stopifnot(identical(length(max), 1L))
  stopifnot(!is.na(max))
  stopifnot(max >= 0L)
  #
  fl <- floor(mu)
  logm <- lgamma(fl + 1)
  lmu <- ifelse(mu==0, exp(-100), log(mu))
  for(i in seq_len(max)){
    ynew <- rpois(n = 1, lambda = mu)
    logy <- lgamma(ynew + 1)
    log_a <- (nu - 1) * (lmu * (ynew - fl) - logy + logm)
    u <- log(runif(n = 1))
    if(u < log_a)
      return(ynew)
  }
  return(-Inf)
}

# Function to generate overdispersed data
rcmpOver <- function(mu, nu, max){
  ## 'mu'
  stopifnot(is.double(mu))
  stopifnot(identical(length(mu), 1L))
  stopifnot(!is.na(mu))
  stopifnot(mu >= 0L)
  ## 'nu'
  stopifnot(is.double(nu))
  stopifnot(identical(length(nu), 1L))
  stopifnot(!is.na(nu))
  stopifnot(nu >= 0L)
  ## 'max'
  stopifnot(is.double(max))
  stopifnot(identical(length(max), 1L))
  stopifnot(!is.na(max))
  stopifnot(max >= 0L)
  #
  if( mu == 0){
    lmu <- exp(-100)
    p <- 0.99999
  } else {
    lmu <- log(mu)
    p <- 2 * nu / (2 * mu * nu + 1 + nu)
  }
  fl <- floor(mu / ((1 - p)^(1 / nu)))
  logfl <- lgamma(fl + 1)
  for(i in seq_len(max)){
    ynew <- rgeom(n = 1, prob = p)
    logy <- lgamma(ynew + 1)
    log_a <- (ynew - fl) * (nu * lmu - log(1 - p)) + nu * (logfl - logy)
    u <- log(runif(n = 1))
    if(u < log_a)
      return(ynew)
  }
  return(-Inf)
}

# Random generation function for CMP data
rcmp1 <- function(mu, nu, max){
  ## 'mu'
  stopifnot(is.double(mu))
  stopifnot(identical(length(mu), 1L))
  stopifnot(!is.na(mu))
  stopifnot(mu >= 0L)
  ## 'nu'
  stopifnot(is.double(nu))
  stopifnot(identical(length(nu), 1L))
  stopifnot(!is.na(nu))
  stopifnot(nu >= 0L)
  ## 'max'
  stopifnot(is.double(max))
  stopifnot(identical(length(max), 1L))
  stopifnot(!is.na(max))
  stopifnot(max >= 0L)
  #
  if( nu < 1){
    rcmpOver(mu = mu , nu = nu, max = max)
  } else {
    rcmpUnder(mu = mu , nu = nu, max = max)
  }
}
