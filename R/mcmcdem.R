library(LaplacesDemon) # to use half t distribution
#---------------
# CMP functions
#---------------

# Unnormalised CMP density function
logDensCMP1 <- function(y, gamma, nu, useC = FALSE) {
  ## 'y'
  stopifnot(is.integer(y))
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

#---------------------
# Exchange algorithms
#---------------------

# For model type: y_i ~ CMP(mu_i, nu_i)
# The function applies to one observation at time
a.exch <- function(y, lambda.curr, omega.curr, lambda.cand, omega.cand, q, mu,
                   sigma, eta, tau, max = 10^4){
  # CMP parameters
  mu.curr <- exp(lambda.curr)*q
  nu.curr <- exp(omega.curr)
  #
  mu.cand <- exp(lambda.cand)*q
  nu.cand <- exp(omega.cand)
  #
  # Generate the auxiliary data
  y.new <- rcmp1(mu.cand, nu.cand, max = max)
  #
  # Log-likelihood ratio
  llik.ratio <- sum(logDensCMP1(y,mu.cand,nu.cand),
                    - logDensCMP1(y,mu.curr,nu.curr),
                    logDensCMP1(y.new,mu.curr,nu.curr),
                    - logDensCMP1(y.new,mu.cand,nu.cand))
  # Prior ratio
  prior.ratio <- sum(dnorm(lambda.cand, mu, sigma, log=TRUE),
                     - dnorm(lambda.curr, mu, sigma, log=TRUE),
                     dnorm(omega.cand, eta, tau, log=TRUE),
                     - dnorm(omega.curr, eta, tau, log=TRUE))
  # Posterior
  post <- llik.ratio + prior.ratio
  post > log(runif(1))
}

# For model type: y_i ~ CMP(mu_i, nu)
# The function applies to all the data y
a.exch1 <- function(y, lambda.curr, omega.curr, lambda.cand, omega.cand, q, mu,
                    sigma, eta, tau, max = 10^4){
  n<- length(y)
  # CMP parameters
  #
  mu.curr <- exp(lambda.curr)*q
  nu.curr <- exp(omega.curr)
  #
  mu.cand <- exp(lambda.cand)*q
  nu.cand <- exp(omega.cand)
  #
  # Generate the auxiliary data
  y.new <- c()
  for(k in 1:n){
    y.new[k] <- rcmp1(mu.cand[k], nu.cand, max = max)
  }
  #
  # Log-likelihood ratio
  llik1 <- c()
  for(a in 1:n){
    llik1[a] <- sum(logDensCMP1(y[a],mu.cand[a],nu.cand),
                    - logDensCMP1(y[a],mu.curr[a],nu.curr),
                    logDensCMP1(y.new[a],mu.curr[a],nu.curr),
                    - logDensCMP1(y.new[a],mu.cand[a],nu.cand))

  }
  llik.ratio <- sum(llik1)
  # Prior ratio
  prior.ratio <- sum(dnorm(omega.cand, eta, tau, log=TRUE),
                     - dnorm(omega.curr, eta, tau, log=TRUE))
  # Posterior
  post <- llik.ratio + prior.ratio
  post > log(runif(1))
}

#------
# MCMC
#------
#---------------------------------
# First simple model, sd ~ half t
#---------------------------------
# Function for variance M-H algorithm
lpost.var <- function(var.cand, ldata.curr, mean.curr, var.curr, scale.var, nu.var){
  lprior.cand <- dhalft(var.cand, scale = scale.var, nu = nu.var, log=T)
  llike.cand <- sum(dnorm(ldata.curr, mean = mean.curr, sd = var.cand, log = T))
  lprior.curr <- dhalft(var.curr, scale = scale.var, nu = nu.var, log=T)
  llike.curr <- sum(dnorm(ldata.curr, mean = mean.curr, sd = var.curr, log = T))
  lpost <- lprior.cand + llike.cand - lprior.curr - llike.curr
  return(lpost>log(runif(1)))
}

mcmc.dem.simple <- function(y, iter, expo, gamma0 = rep(1,n), nu0 = rep(1,n),
                            lgamma0 = log(gamma0), lnu0 = log(nu0)){

  gamma.mat <- matrix(NA, iter, n)
  nu.mat <- matrix(NA, iter, n)


  # Initialise

  gamma.mat[1,] <-gamma0
  nu.mat[1,] <- nu0

  # Acceptance

  gamma.nu.acc <- 0
  gamma.acc <- 0
  nu.acc <- 0

  # nu and gamma are updated with the exchange algorithm
  # for mu and eta Gibbs sampling
  # sigma and tau MH?
  for(i in 2:iter){
    if(i%%(0.1*iter)==0) print(i)
    gamma.curr <- gamma.mat[i-1,]
    nu.curr <- nu.mat[i-1,]
    lgamma.curr <- log(gamma.curr)
    lnu.curr <- log(nu.curr)


    # generate candidates

    # candidate for c(gamma, nu):

    gamma.cand <- exp(rnorm(n, lgamma.curr, 0.1))
    nu.cand <- exp(rnorm(n, lnu.curr, 0.1))

    # Exchange algorithm for gamma and nu
    par.post<- c()
    for(k in 1:n){
      par.post[k] <- a.exch(y[k], gamma.curr[k], nu.curr[k],
                            gamma.cand[k], nu.cand[k], expo[k])
    }

    gamma.mat[i,] <- ifelse(par.post, gamma.cand, gamma.curr)
    nu.mat[i,] <- ifelse(par.post, nu.cand, nu.curr)
    gamma.nu.acc <- ifelse(par.post, gamma.nu.acc+1, gamma.nu.acc)
    lgamma.curr <- ifelse(par.post, log(gamma.cand), log(gamma.curr))
    lnu.curr <- ifelse(par.post, log(nu.cand), log(nu.curr))

    # candidate for c(gamma, nu):

    gamma.cand <- exp(rnorm(n, lgamma.curr, 0.1))
    nu.cand <- exp(rnorm(n, lnu.curr, 0.1))

    # Exchange algorithm for gamma:
    par.post.g<- c()
    for(k in 1:n){
      par.post.g[k] <- a.exch(y[k], gamma.curr[k], nu.curr[k],
                              gamma.cand[k], nu.curr[k], expo[k])
    }

    gamma.mat[i,] <- ifelse(par.post.g, gamma.cand, gamma.curr)
    gamma.acc <- ifelse(par.post.g, gamma.acc+1, gamma.acc)

    # Exchange algorithm for nu:
    par.post.n <- c()
    for(k in 1:n){
      par.post.n[k] <- a.exch(y[k], gamma.curr[k], nu.curr[k],
                              gamma.curr[k], nu.cand[k], expo[k])
    }

    nu.mat[i,] <- ifelse(par.post.n, nu.cand, nu.curr)
    nu.acc <- ifelse(par.post.n, nu.acc+1, nu.acc)


  }

  return(list (parameters = list(gamma = gamma.mat, nu = nu.mat),
               acceptance = list(gamma.nu = gamma.nu.acc/iter,
                                 gamma = gamma.acc/iter,
                                 nu = nu.acc/iter))
  )
}

#-------------------------------------
# Model assuming IG prior on Variance
# (i.e. Gibbs sampling for sigma and tau as conjugate model)
#-------------------------------------

mcmc.dem.IG <- function(y, iter, expo, lambda0 = rep(0,length(y)), omega0 = rep(0,length(y)),
                        mu0 = 0, sigma0 = 0.1, eta0 = 0, tau0 = 0.5,
                        alpha0 = 0, beta0 = 0.1, delta0 = 0, xi0 = 0.1,
                        sh.sigma0 = length(y), sh.tau0 = length(y),
                        sc.sigma0 = 0.1, sc.tau0 = 0.1,
                        mean.lambda.cand = 0, mean.omega.cand = 0,
                        var.lambda.cand = 1, var.omega.cand = 1,
                        model.upd){

  n <- length(y)
  # Space for storing values
  lambda.mat <- matrix(NA, iter, n)
  omega.mat <- matrix(NA, iter, n)
  mu.vec <- c()
  sigma.vec <- c()
  eta.vec <- c()
  tau.vec <- c()

  # Initialise

  lambda.mat[1,] <- lambda0
  omega.mat[1,] <- omega0
  mu.vec[1] <- mu0
  sigma.vec[1] <- sigma0
  eta.vec[1] <- eta0
  tau.vec[1] <- tau0

  # Acceptance

  acceptance <- 0

  # Begin loop
  for(i in 2:iter){
    if(i%%(0.1*iter)==0) print(i)
    lambda.curr <- lambda.mat[(i-1),]
    omega.curr <- omega.mat[(i-1),]
    #
    # Parameters for Gibbs sampling

    # Update mu
    # mean
    num.alpha <- alpha0 / beta0 + sum(lambda.mat[(i-1),])/sigma.vec[(i-1)]
    den.alpha <- (1 / beta0 + n / sigma.vec[(i-1)])
    alpha <-  num.alpha / den.alpha
    # sd
    beta <- sqrt(1 / (1 / beta0 + n / sigma.vec[(i-1)]))
    #
    mu.vec[i] <- rnorm(1, alpha, beta)

    # Update sigma
    # shape
    sh.sigma <- sh.sigma0 + n / 2
    # scale
    sc.sigma <- sc.sigma0 + 0.5 * sum((lambda.mat[(i-1),] - mu.vec[i])^2)
    #
    sigma.vec[i] <- 1 / rgamma(1, shape = sh.sigma, scale = sc.sigma)

    # Update eta
    # mean
    num.delta <- delta0 / xi0 + sum(omega.mat[(i-1),])/tau.vec[(i-1)]
    den.delta <- 1 / xi0 + n / tau.vec[(i-1)]
    delta <-  num.delta / den.delta
    # sd
    xi <- sqrt(1 / (1 / xi0 + n / tau.vec[(i-1)]))
    #
    eta.vec[i] <- rnorm(1, delta, xi)

    # Update tau
    # shape
    sh.tau <- sh.tau0 + n / 2
    # scale
    sc.tau <- sc.tau0 + 0.5 * sum((omega.mat[(i-1),] - eta.vec[i])^2)
    #
    tau.vec[i] <- 1 / rgamma(1, shape = sh.tau, scale = sc.tau)

    mu.curr <- mu.vec[i]
    sigma.curr <- sigma.vec[i]
    eta.curr <- eta.vec[i]
    tau.curr <- tau.vec[i]

    # candidate for lambda and omega:
    # by default standard normal or if model.upd == "RW" it is a random walk step
    mean.lambda.cand <- ifelse(model.upd=="RW", lambda.mat[(i-1),], mean.lambda.cand)
    mean.omega.cand <- ifelse(model.upd=="RW", omega.mat[(i-1),], mean.omega.cand)
    #
    lambda.cand <- rnorm(n, mean.lambda.cand, var.lambda.cand)
    omega.cand <- rnorm(n, mean.omega.cand, var.omega.cand)


    # Exchange algorithm for lambda and omega
    par.post <- c()
    for(k in 1:n){
      par.post[k] <- a.exch(y[k], lambda.curr[k], omega.curr[k],
                            lambda.cand[k], omega.cand[k],
                            expo[k], mu.curr, sigma.curr,
                            eta.curr, tau.curr)
    }

    lambda.mat[i,] <- ifelse(par.post, lambda.cand, lambda.curr)
    omega.mat[i,] <- ifelse(par.post, omega.cand, omega.curr)
    acceptance <- ifelse(par.post, acceptance + 1, acceptance)

  }

  return(list (parameters = list(gamma = exp(lambda.mat), nu = exp(omega.mat),
                                 mu = mu.vec, sigma = sigma.vec,
                                 eta = eta.vec, tau = tau.vec),
               acceptance = acceptance/iter)
  )
}


#-----------------------------------
# Model assuming variances as known
#-----------------------------------
mcmc.dem.KV <- function(y, iter, expo, lambda0 = rep(0,length(y)), omega0 = rep(0,length(y)),
                        mu0 = 0, sigma, eta0 = 0, tau,
                        alpha0 = 0, beta0 = 0.1, delta0 = 0, xi0 = 0.1,
                        mean.lambda.cand = 0, mean.omega.cand = 0,
                        var.lambda.cand = 1, var.omega.cand = 1,
                        model.upd){

  n <- length(y)
  # Space for storing values
  lambda.mat <- matrix(NA, iter, n)
  omega.mat <- matrix(NA, iter, n)
  mu.vec <- c()
  eta.vec <- c()

  # Initialise

  lambda.mat[1,] <- lambda0
  omega.mat[1,] <- omega0
  mu.vec[1] <- mu0
  eta.vec[1] <- eta0


  # Acceptance

  acceptance <- 0

  # Begin loop
  for(i in 2:iter){
    if(i%%(0.1*iter)==0) print(i)
    lambda.curr <- lambda.mat[(i-1),]
    omega.curr <- omega.mat[(i-1),]
    #
    # Parameters for Gibbs sampling

    # Update mu
    # mean
    num.alpha <- alpha0 / beta0 + sum(lambda.mat[(i-1),])/sigma
    den.alpha <- (1 / beta0 + n / sigma)
    alpha <-  num.alpha / den.alpha
    # sd
    beta <- sqrt(1 / (1 / beta0 + n / sigma))
    #
    mu.vec[i] <- rnorm(1, alpha, beta)

    # Update eta
    # mean
    num.delta <- delta0 / xi0 + sum(omega.mat[(i-1),])/tau
    den.delta <- 1 / xi0 + n / tau
    delta <-  num.delta / den.delta
    # sd
    xi <- sqrt(1 / (1 / xi0 + n / tau))
    #
    eta.vec[i] <- rnorm(1, delta, xi)

    mu.curr <- mu.vec[i]
    eta.curr <- eta.vec[i]

    # candidate for lambda and omega:
    # by default standard normal or if model.upd == "RW" it is a random walk step
    # not using mu and eta values from the Gibbs sampler in the proposal
    mean.lambda.cand <- ifelse(model.upd=="RW", lambda.mat[(i-1),], mean.lambda.cand)
    mean.omega.cand <- ifelse(model.upd=="RW", omega.mat[(i-1),], mean.omega.cand)
    #
    lambda.cand <- rnorm(n, mean.lambda.cand, var.lambda.cand)
    omega.cand <- rnorm(n, mean.omega.cand, var.omega.cand)


    # Exchange algorithm for lambda and omega
    par.post <- c()
    for(k in 1:n){
      par.post[k] <- a.exch(y[k], lambda.curr[k], omega.curr[k],
                            lambda.cand[k], omega.cand[k],
                            expo[k], mu.curr, sigma,
                            eta.curr, tau)
    }

    lambda.mat[i,] <- ifelse(par.post, lambda.cand, lambda.curr)
    omega.mat[i,] <- ifelse(par.post, omega.cand, omega.curr)
    acceptance <- ifelse(par.post, acceptance + 1, acceptance)

  }

  return(list (parameters = list(gamma = exp(lambda.mat), nu = exp(omega.mat),
                                 mu = mu.vec, eta = eta.vec),
               acceptance = acceptance/iter)
  )
}

#-------------------------------------------------------
# Only estimate lambda and omega, everything else known
# For model type: y_i ~ CMP(mu_i, nu_i)
#-------------------------------------------------------
mcmc.dem.lo <- function(y, iter, expo, lambda0 = rep(0,length(y)),
                        omega0 = rep(0,length(y)),
                        mu, sigma, eta, tau,
                        mean.lambda.cand = 0, mean.omega.cand = 0,
                        var.lambda.cand = 1, var.omega.cand = 1,
                        model.upd){

  n <- length(y)
  # Space for storing values
  lambda.mat <- matrix(NA, iter, n)
  omega.mat <- matrix(NA, iter, n)

  # Initialise

  lambda.mat[1,] <- lambda0
  omega.mat[1,] <- omega0


  # Acceptance

  acceptance <- 0

  # Begin loop
  for(i in 2:iter){
    if(i%%(0.1*iter)==0) print(i)
    lambda.curr <- lambda.mat[(i-1),]
    omega.curr <- omega.mat[(i-1),]
    #

    # candidate for lambda and omega:
    # by default standard normal or if model.upd == "RW" it is a random walk step
    mean.lambda.cand <- ifelse(model.upd=="RW", lambda.mat[(i-1),], mean.lambda.cand)
    mean.omega.cand <- ifelse(model.upd=="RW", omega.mat[(i-1),], mean.omega.cand)
    #
    lambda.cand <- rnorm(n, mean.lambda.cand, var.lambda.cand)
    omega.cand <- rnorm(n, mean.omega.cand, var.omega.cand)


    # Exchange algorithm for lambda and omega
    par.post <- c()
    for(k in 1:n){
      par.post[k] <- a.exch(y[k], lambda.curr[k], omega.curr[k],
                            lambda.cand[k], omega.cand[k],
                            expo[k], mu, sigma,
                            eta, tau)
    }

    lambda.mat[i,] <- ifelse(par.post, lambda.cand, lambda.curr)
    omega.mat[i,] <- ifelse(par.post, omega.cand, omega.curr)
    acceptance <- ifelse(par.post, acceptance + 1, acceptance)

  }

  return(list (parameters = list(gamma = exp(lambda.mat), nu = exp(omega.mat)),
               acceptance = acceptance/iter)
  )
}

#-------------------------------------
# Only estimate lambda and omega
# For model type: y_i ~ CMP(mu_i, nu)
#-------------------------------------

mcmc.dem.lo1 <- function(y, iter, expo, lambda0 = rep(0,length(y)), omega0 = 0,
                         mu, sigma, eta, tau,
                         mean.lambda.cand = 0, mean.omega.cand = 0,
                         var.lambda.cand = 1, var.omega.cand = 1,
                         model.upd){

  n <- length(y)
  # Space for storing values
  lambda.mat <- matrix(NA, iter, n)
  omega.vec <- c()

  # Initialise

  lambda.mat[1,] <- lambda0
  omega.vec[1] <- omega0


  # Acceptance

  acceptance.lambda <- 0
  acceptance.omega <- 0
  # Begin loop
  for(i in 2:iter){
    if(i%%(0.1*iter)==0) print(i)
    lambda.curr <- lambda.mat[(i-1),]
    omega.curr <- omega.vec[(i-1)]
    #

    # candidate for lambda and omega:
    # by default standard normal or if model.upd == "RW" it is a random walk step
    mean.lambda.cand <- ifelse(model.upd=="RW", lambda.mat[(i-1),], mean.lambda.cand)
    mean.omega.cand <- ifelse(model.upd=="RW", omega.vec[(i-1)], mean.omega.cand)
    #
    lambda.cand <- rnorm(n, mean.lambda.cand, var.lambda.cand)
    omega.cand <- rnorm(1, mean.omega.cand, var.omega.cand)


    # Exchange algorithm for lambda and omega
    par.post.lambda <- c()
    for(k in 1:n){
      par.post.lambda[k] <- a.exch(y[k], lambda.curr[k], omega.curr,
                                   lambda.cand[k], omega.cand,
                                   expo[k], mu, sigma,
                                   eta, tau)
    }


    lambda.mat[i,] <- ifelse(par.post.lambda, lambda.cand, lambda.curr)
    acceptance.lambda <- ifelse(par.post.lambda, acceptance.lambda + 1,acceptance.lambda)

    par.post.omega <- a.exch1(y, lambda.curr, omega.curr,
                              lambda.curr, omega.cand,
                              expo, mu, sigma,
                              eta, tau)
    if(par.post.omega){
      omega.vec[i] <- omega.cand
      acceptance.omega <- acceptance.omega + 1
    } else {

      omega.vec[i] <- omega.curr
      acceptance.omega <- acceptance.omega
    }
  }

  return(list (parameters = list(gamma = exp(lambda.mat), nu = exp(omega.vec)),
               acceptance = list(gamma=acceptance.lambda/iter,
                                 nu=acceptance.omega/iter))
  )
}

#---------------------------------------
# Only estimate omega
# For model type: y_i ~ CMP(mu_i, nu_i)
#---------------------------------------

mcmc.dem.om <- function(y, iter, expo, lambda, omega0 = rep(0,length(y)),
                        mu, sigma, eta, tau,
                        mean.omega.cand = 0,
                        var.omega.cand = 1,
                        model.upd){

  n <- length(y)
  # Space for storing values
  omega.mat <- matrix(NA, iter, n)

  # Initialise
  omega.mat[1,] <- omega0


  # Acceptance

  acceptance <- 0

  # Begin loop
  for(i in 2:iter){
    if(i%%(0.1*iter)==0) print(i)

    omega.curr <- omega.mat[(i-1),]
    #

    # candidate for lambda and omega:
    # by default standard normal or if model.upd == "RW" it is a random walk step
    mean.omega.cand <- ifelse(model.upd=="RW", omega.mat[(i-1),], mean.omega.cand)
    #
    omega.cand <- rnorm(n, mean.omega.cand, var.omega.cand)


    # Exchange algorithm for lambda and omega
    par.post <- c()
    for(k in 1:n){
      par.post[k] <- a.exch(y[k], lambda[k], omega.curr[k],
                            lambda[k], omega.cand[k],
                            expo[k], mu, sigma, eta, tau)
    }

    omega.mat[i,] <- ifelse(par.post, omega.cand, omega.curr)
    acceptance <- ifelse(par.post, acceptance + 1, acceptance)

  }

  return(list (parameters = list(nu = exp(omega.mat)),
               acceptance = acceptance/iter)
  )
}


#-------------------------------------
# Only estimate unique omega
# For model type: y_i ~ CMP(mu_i, nu)
#-------------------------------------

mcmc.dem.om1 <- function(y, iter, expo, lambda, omega0 = 0,
                         mu, sigma, eta, tau,
                         mean.omega.cand = 0,
                         var.omega.cand = 1,
                         model.upd){

  n <- length(y)
  # Space for storing values
  omega.vec <- c()

  # Initialise
  omega.vec[1] <- omega0


  # Acceptance

  acceptance <- 0

  # Begin loop
  for(i in 2:iter){
    if(i%%(0.1*iter)==0) print(i)

    omega.curr <- omega.vec[i-1]
    #

    # candidate for lambda and omega:
    # by default standard normal or if model.upd == "RW" it is a random walk step
    mean.omega.cand <- ifelse(model.upd=="RW", omega.vec[(i-1)], mean.omega.cand)
    #
    omega.cand <- rnorm(1, eta, var.omega.cand)


    # Exchange algorithm for lambda and omega
    par.post <- a.exch1(y, lambda, omega.curr,
                        lambda, omega.cand,
                        expo, mu, sigma, eta, tau)

    if(par.post){
      omega.vec[i] <- omega.cand
      acceptance <- acceptance + 1
    } else {
      omega.vec[i] <- omega.curr
      acceptance <-acceptance
    }

  }

  return(list (parameters = list(nu = exp(omega.vec)),
               acceptance = acceptance/iter)
  )
}


#------------------------------------------------------------------
# Variance known, unique omega, separate updating lambda and omega
# For model type: y_i ~ CMP(mu_i, nu)
#------------------------------------------------------------------

mcmc.dem.kvslom1 <- function(y, iter, expo, lambda0 = rep(0,length(y)), omega0 = 0,
                             mu0 = 0, sigma0 = 0.1, eta0 = 0, tau0 = 0.1,
                             alpha0 = 0, beta0 = 0.1, delta0 = 0, xi0 = 0.1,
                             mean.lambda.cand = 0, mean.omega.cand = 0,
                             var.lambda.cand = 1, var.omega.cand = 1,
                             model.upd){

  n <- length(y)
  # Space for storing values
  lambda.mat <- matrix(NA, iter, n)
  omega.vec <- c()
  mu.vec <- c()
  eta.vec <- c()

  # Initialise

  lambda.mat[1,] <- lambda0
  omega.vec[1] <- omega0
  mu.vec[1] <- mu0
  eta.vec[1] <- eta0


  # Acceptance

  acceptance.lambda <- 0
  acceptance.omega <- 0

  # Begin loop
  for(i in 2:iter){
    if(i%%(0.1*iter)==0) print(i)
    lambda.curr <- lambda.mat[(i-1),]
    omega.curr <- omega.vec[(i-1)]
    #
    # Parameters for Gibbs sampling

    # Update mu
    # mean
    num.alpha <- alpha0 / beta0 + sum(lambda.mat[(i-1),])/sigma
    den.alpha <- (1 / beta0 + n / sigma)
    alpha <-  num.alpha / den.alpha
    # sd
    beta <- sqrt(1 / (1 / beta0 + n / sigma))
    #
    mu.vec[i] <- rnorm(1, alpha, beta)

    # Update eta
    # mean
    num.delta <- delta0 / xi0 + omega.vec[(i-1)]/tau
    den.delta <- 1 / xi0 + 1 / tau
    delta <-  num.delta / den.delta
    # sd
    xi <- sqrt(1 / (1 / xi0 + 1 / tau))
    #
    eta.vec[i] <- rnorm(1, delta, xi)

    mu.curr <- mu.vec[i]
    eta.curr <- eta.vec[i]

    # candidate for lambda and omega:
    # by default standard normal or if model.upd == "RW" it is a random walk step
    mean.lambda.cand <- ifelse(model.upd=="RW", lambda.mat[(i-1),], mean.lambda.cand)
    mean.omega.cand <- ifelse(model.upd=="RW", omega.vec[(i-1)], mean.omega.cand)
    #
    lambda.cand <- rnorm(n, mean.lambda.cand, var.lambda.cand)
    omega.cand <- rnorm(1, mean.omega.cand, var.omega.cand)


    # Exchange algorithm for lambda
    par.post.lambda <- c()
    for(k in 1:n){
      par.post.lambda[k] <- a.exch(y[k], lambda.curr[k], omega.curr,
                                   lambda.cand[k], omega.curr,
                                   expo[k], mu.curr, sigma,
                                   eta.curr, tau)
    }

    lambda.mat[i,] <- ifelse(par.post.lambda, lambda.cand, lambda.curr)
    acceptance.lambda <- ifelse(par.post.lambda, acceptance.lambda + 1, acceptance.lambda)


    # Exchange algorithm for omega
    par.post.omega <- a.exch1(y, lambda.curr, omega.curr,
                              lambda.curr, omega.cand,
                              expo, mu.curr, sigma,
                              eta.curr, tau)

    if(par.post.omega){
      omega.vec[i] <- omega.cand
      acceptance.omega <- acceptance.omega + 1
    } else {
      omega.vec[i] <- omega.curr
      acceptance.omega <- acceptance.omega}


  }


  return(list (parameters = list(gamma = exp(lambda.mat), nu = exp(omega.vec),
                                 mu = mu.vec, eta = eta.vec),
               acceptance = list(gamma= acceptance.lambda/iter,
                                 nu = acceptance.omega/iter))
  )
}


#----------------------------------------------------
# Variance known, separate updating lambda and omega
# For model type: y_i ~ CMP(mu_i, nu_i)
#----------------------------------------------------

mcmc.dem.kvslom <- function(y, iter, expo, lambda0 = rep(0,length(y)), omega0 = rep(0,length(y)),
                            mu0 = 0, sigma, eta0 = 0, tau0 = 0,
                            alpha0 = 0, beta0 = 0.1, delta0 = 0, xi0 = 0.1,
                            mean.lambda.cand = 0, mean.omega.cand = 0,
                            var.lambda.cand = 1, var.omega.cand = 1,
                            model.upd){

  n <- length(y)
  # Space for storing values
  lambda.mat <- matrix(NA, iter, n)
  omega.mat <- matrix(NA, iter, n)
  mu.vec <- c()
  eta.vec <- c()

  # Initialise

  lambda.mat[1,] <- lambda0
  omega.mat[1,] <- omega0
  mu.vec[1] <- mu0
  eta.vec[1] <- eta0


  # Acceptance

  acceptance.lambda <- 0
  acceptance.omega <- 0

  # Begin loop
  for(i in 2:iter){
    if(i%%(0.1*iter)==0) print(i)
    lambda.curr <- lambda.mat[(i-1),]
    omega.curr <- omega.mat[(i-1),]
    #
    # Parameters for Gibbs sampling

    # Update mu
    # mean
    num.alpha <- alpha0 / beta0 + sum(lambda.mat[(i-1),])/sigma
    den.alpha <- (1 / beta0 + n / sigma)
    alpha <-  num.alpha / den.alpha
    # sd
    beta <- sqrt(1 / (1 / beta0 + n / sigma))
    #
    mu.vec[i] <- rnorm(1, alpha, beta)

    # Update eta
    # mean
    num.delta <- delta0 / xi0 + omega.mat[(i-1),]/tau
    den.delta <- 1 / xi0 + 1 / tau
    delta <-  num.delta / den.delta
    # sd
    xi <- sqrt(1 / (1 / xi0 + 1 / tau))
    #
    eta.vec[i] <- rnorm(1, delta, xi)

    mu.curr <- mu.vec[i]
    eta.curr <- eta.vec[i]

    # candidate for lambda and omega:
    # by default standard normal or if model.upd == "RW" it is a random walk step
    mean.lambda.cand <- ifelse(model.upd=="RW", lambda.mat[(i-1),], mean.lambda.cand)
    mean.omega.cand <- ifelse(model.upd=="RW", omega.mat[(i-1),], mean.omega.cand)
    #
    lambda.cand <- rnorm(n, mean.lambda.cand, var.lambda.cand)
    omega.cand <- rnorm(n, mean.omega.cand, var.omega.cand)


    # Exchange algorithm for lambda
    par.post.lambda <- c()
    for(k in 1:n){
      par.post.lambda[k] <- a.exch(y[k], lambda.curr[k], omega.curr[k],
                                   lambda.cand[k], omega.curr[k],
                                   expo[k], mu.curr, sigma,
                                   eta.curr, tau)
    }

    lambda.mat[i,] <- ifelse(par.post, lambda.cand, lambda.curr)
    acceptance.lambda <- ifelse(par.post.lambda, acceptance.lambda + 1, acceptance.lambda)


    # Exchange algorithm for omega

    par.post.omega <- c()
    for(k in 1:n){
      par.post.omega[k] <- a.exch(y[k], lambda.curr[k], omega.curr[k],
                                  lambda.curr[k], omega.cand[k],
                                  expo[k], mu.curr, sigma,
                                  eta.curr, tau)
    }

    omega.mat[i,] <- ifelse(par.post, omega.cand, omega.curr)
    acceptance.omega <- ifelse(par.post.omega, acceptance.omega + 1, acceptance.omega)

  }


  return(list (parameters = list(gamma = exp(lambda.mat), nu = exp(omega.mat),
                                 mu = mu.vec, eta = eta.vec),
               acceptance = list(gamma=acceptance.lambda/iter,
                                 nu = acceptance.omega/iter))
  )
}

#-------------------------------------------------------------------
# Variance and second level known, double updating lambda and omega
# For model type: y_i ~ CMP(mu_i, nu)
#-------------------------------------------------------------------

mcmc.dem.double <- function(y, iter, expo, lambda0 = rep(0,length(y)), omega0 = rep(0,length(y)),
                            mu0 = 0, sigma, eta0 = 0, tau0 = 0,
                            alpha0 = 0, beta0 = 0.1, delta0 = 0, xi0 = 0.1,
                            mean.lambda.cand = 0, mean.omega.cand = 0,
                            var.lambda.cand = 10, var.omega.cand = 10,
                            model.upd){
  n <- length(y)
  # Space for storing values
  lambda.mat <- matrix(NA, iter, n)
  omega.mat <- matrix(NA, iter, n)
  mu.vec <- c()
  eta.vec <- c()

  # Initialise

  lambda.mat[1,] <- lambda0
  omega.mat[1,] <- omega0
  mu.vec[1] <- mu0
  eta.vec[1] <- eta0


  # Acceptance

  acceptance <- 0
  acceptance.lambda <- 0
  acceptance.omega <- 0

  # Begin loop
  for(i in 2:iter){
    if(i%%(0.1*iter)==0) print(i)
    lambda.curr <- lambda.mat[(i-1),]
    omega.curr <- omega.mat[(i-1),]
    #
    # Parameters for Gibbs sampling

    # Update mu
    # mean
    num.alpha <- alpha0 / beta0 + sum(lambda.mat[(i-1),])/sigma
    den.alpha <- (1 / beta0 + n / sigma)
    alpha <-  num.alpha / den.alpha
    # sd
    beta <- sqrt(1 / (1 / beta0 + n / sigma))
    #
    mu.vec[i] <- rnorm(1, alpha, beta)

    # Update eta
    # mean
    num.delta <- delta0 / xi0 + omega.mat[(i-1),]/tau
    den.delta <- 1 / xi0 + 1 / tau
    delta <-  num.delta / den.delta
    # sd
    xi <- sqrt(1 / (1 / xi0 + 1 / tau))
    #
    eta.vec[i] <- rnorm(1, delta, xi)

    mu.curr <- mu.vec[i]
    eta.curr <- eta.vec[i]

    # candidate for lambda and omega:

    # by default standard normal or if model.upd == "RW" it is a random walk step
    mean.lambda.cand <- ifelse(model.upd=="RW", lambda.mat[(i-1),], mu.curr)
    mean.omega.cand <- ifelse(model.upd=="RW", omega.mat[(i-1),], eta.curr)
    #
    lambda.cand <- rnorm(n, mean.lambda.cand, var.lambda.cand)
    omega.cand <- rnorm(n, mean.omega.cand, var.omega.cand)


    # Exchange algorithm for lambda and omega
    par.post <- c()
    for(k in 1:n){
      par.post[k] <- a.exch(y[k], lambda.curr[k], omega.curr[k],
                            lambda.cand[k], omega.cand[k],
                            expo[k], mu.curr, sigma,
                            eta.curr, tau)
    }

    lambda.mat[i,] <- ifelse(par.post, lambda.cand, lambda.curr)
    omega.mat[i,] <- ifelse(par.post, omega.cand, omega.curr)
    acceptance <- ifelse(par.post, acceptance + 1, acceptance)

    # Exchange algorithm for lambda

    par.post.lambda <- a.exch1(y, lambda.curr, omega.curr, lambda.cand, omega.curr, expo,
                               mu.curr, sigma, eta.curr, tau)


    lambda.mat[i,] <- ifelse(par.post.lambda, lambda.cand, lambda.curr)
    acceptance.lambda <- ifelse(par.post.lambda, acceptance.lambda + 1, acceptance.lambda)

    # Exchange algorithm for omega

    par.post.omega <- a.exch1(y, lambda.curr, omega.curr,
                              lambda.curr, omega.cand, expo,
                              mu.curr, sigma, eta.curr, tau)

    omega.mat[i,] <- ifelse(par.post.omega, omega.cand, omega.curr)
    acceptance.omega <- ifelse(par.post.omega, acceptance.omega + 1, acceptance.omega)

  }


  return(list (parameters = list(gamma = exp(lambda.mat), nu = exp(omega.mat),
                                 mu = mu.vec, eta = eta.vec),
               acceptance = list(couple = acceptance/iter,
                                 lambda = acceptance.lambda/iter,
                                 omega = acceptance.omega/iter)))
}

library(devtools)
devtools::load_all()
