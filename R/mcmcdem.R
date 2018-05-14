#--------------------------------------------
# MCMC For model type: y_i = CMP(mu_i, nu_i)
#--------------------------------------------
library(LaplacesDemon) # to use half t distribution
#---------------------------------
# First simple model, sd = half t
#---------------------------------
# Function for variance M-H algorithm
lpost.sd <- function(sd.cand, ldata.curr, mean.curr, sd.curr, scale.sd, nu.sd){
  lprior.cand <- dhalft(sd.cand, scale = scale.sd, nu = nu.sd, log = T)
  llike.cand <- sum(dnorm(ldata.curr, mean = mean.curr, sd = sd.cand, log = T))
  lprior.curr <- dhalft(sd.curr, scale = scale.sd, nu = nu.sd, log = T)
  llike.curr <- sum(dnorm(ldata.curr, mean = mean.curr, sd = sd.curr, log = T))
  lpost <- lprior.cand + llike.cand - lprior.curr - llike.curr
  return(lpost>log(runif(1)))
}

mcmc.dem.simple <- function(y, iter, expo, gamma0 = rep(1,length(y)), nu0 = rep(1,length(y)),
                            lgamma0 = log(gamma0), lnu0 = log(nu0),
                            mu = 0, sigma = 10, eta = 0, tau =10){
  n <-length(y)
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
    lgamma.curr<- log(gamma.curr)
    lnu.curr <- log(nu.curr)


    # generate candidates

    # candidate for c(gamma, nu):

    lgamma.cand <- rnorm(n, lgamma.curr, 0.1)
    lnu.cand <- rnorm(n, lnu.curr, 0.1)

    # Exchange algorithm for gamma and nu
    par.post<- c()
    for(k in 1:n){
      par.post[k] <- a.exch(y[k], lgamma.curr[k], lnu.curr[k],
                            lgamma.cand[k], lnu.cand[k], expo[k], mu,
                            sigma, eta, tau)$test
    }

    gamma.mat[i,] <- ifelse(par.post, exp(lgamma.cand), exp(lgamma.curr))
    nu.mat[i,] <- ifelse(par.post, exp(lnu.cand), exp(lnu.curr))
    gamma.nu.acc <- ifelse(par.post, gamma.nu.acc+1, gamma.nu.acc)
    lgamma.curr <- ifelse(par.post, lgamma.cand, lgamma.curr)
    lnu.curr <- ifelse(par.post, lnu.cand, lnu.curr)

    # candidate for c(gamma, nu):

    lgamma.cand2 <- rnorm(n, lgamma.curr, 0.1)
    lnu.cand2 <- rnorm(n, lnu.curr, 0.1)

    # Exchange algorithm for gamma:
    par.post.g<- c()
    for(k in 1:n){
      par.post.g[k] <- a.exch(y[k], lgamma.curr[k], lnu.curr[k],
                              lgamma.cand2[k], lnu.curr[k], expo[k],
                              mu, sigma, eta, tau)$test
    }

    gamma.mat[i,] <- ifelse(par.post.g, exp(lgamma.cand2), exp(lgamma.curr))
    gamma.acc <- ifelse(par.post.g, gamma.acc+1, gamma.acc)

    # Exchange algorithm for nu:
    par.post.n <- c()
    for(k in 1:n){
      par.post.n[k] <- a.exch(y[k], lgamma.curr[k], lnu.curr[k],
                              lgamma.curr[k], lnu.cand2[k], expo[k],
                              mu, sigma, eta, tau)$test
    }

    nu.mat[i,] <- ifelse(par.post.n, exp(lnu.cand2), exp(lnu.curr))
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
                        sh.sigma0 = 0.001, sh.tau0 = 0.001,
                        sc.sigma0 = 0.001, sc.tau0 = 0.001,
                        mean.lambda.cand = 0, mean.omega.cand = 0,
                        sd.lambda.cand = 1, sd.omega.cand = 1,
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

  lambda.mat[1,]<- lambda0
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
    num.alpha <- alpha0 / beta0^2 + sum(lambda.mat[(i-1),])/(sigma.vec[(i-1)])^2
    den.alpha <- (1 / beta0^2 + n / (sigma.vec[(i-1)])^2)
    alpha <-  num.alpha / den.alpha
    # sd
    beta <- sqrt(1 / (1 / beta0^2 + n / (sigma.vec[(i-1)])^2))
    #
    mu.vec[i] <- rnorm(1, alpha, beta)

    # Update sigma
    # shape
    sh.sigma <- sh.sigma0 + n / 2
    # scale
    sc.sigma <- sc.sigma0 + 0.5 * sum((lambda.mat[(i-1),] - mu.vec[i])^2)
    #
    sigma.vec[i] <- sqrt(1 / rgamma(1, shape = sh.sigma, scale = sc.sigma))

    # Update eta
    # mean
    num.delta <- delta0 / xi0^2 + sum(omega.mat[(i-1),])/(tau.vec[(i-1)])^2
    den.delta <- 1 / xi0^2 + n / (tau.vec[(i-1)])^2
    delta <-  num.delta / den.delta
    # sd
    xi <- sqrt(1 / (1 / xi0^2 + n / (tau.vec[(i-1)])^2 ))
    #
    eta.vec[i] <- rnorm(1, delta, xi)

    # Update tau
    # shape
    sh.tau <- sh.tau0 + n / 2
    # scale
    sc.tau <- sc.tau0 + 0.5 * sum((omega.mat[(i-1),] - eta.vec[i])^2)
    #
    tau.vec[i] <- sqrt(1 / rgamma(1, shape = sh.tau, scale = sc.tau))

    mu.curr <- mu.vec[i]
    sigma.curr <- sigma.vec[i]
    eta.curr <- eta.vec[i]
    tau.curr <- tau.vec[i]

    # candidate for lambda and omega:
    # by default standard normal or if model.upd == "RW" it is a random walk step
    # or with hyperparameters if Gibbs
    if(model.upd=="RW"){
      mean.lambda.cand <- lambda.mat[(i-1),]
      mean.omega.cand <- omega.mat[(i-1),]
      sd.lambda.cand <- sd.lambda.cand
      sd.omega.cand <- sd.omega.cand
    } else if(model.upd=="Gibbs"){
      mean.lambda.cand <- mu.curr
      mean.omega.cand <- eta.curr
      sd.lambda.cand <- sigma.curr
      sd.omega.cand <- tau.curr
    } else{
      mean.lambda.cand <- mean.lambda.cand
      mean.omega.cand <- mean.omega.cand
      sd.lambda.cand <- sd.lambda.cand
      sd.omega.cand <- sd.omega.cand
    }


    #
    lambda.cand<- rnorm(n, mean.lambda.cand, sd.lambda.cand)
    omega.cand <- rnorm(n, mean.omega.cand, sd.omega.cand)


    # Exchange algorithm for lambda and omega
    par.post <- c()
    for(k in 1:n){
      par.post[k] <- a.exch(y[k], lambda.curr[k], omega.curr[k],
                            lambda.cand[k], omega.cand[k],
                            expo[k], mu.curr, sigma.curr,
                            eta.curr, tau.curr)$test
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
                        sd.lambda.cand = 1, sd.omega.cand = 1,
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
    num.alpha <- alpha0 / beta0^2 + sum(lambda.mat[(i-1),])/sigma^2
    den.alpha <- (1 / beta0^2 + n / sigma^2)
    alpha <-  num.alpha / den.alpha
    # sd
    beta <- sqrt(1 / (1 / beta0^2 + n / sigma^2))
    #
    mu.vec[i] <- rnorm(1, alpha, beta)

    # Update eta
    # mean
    num.delta <- delta0 / xi0^2 + sum(omega.mat[(i-1),])/tau^2
    den.delta <- 1 / xi0^2 + n / tau^2
    delta <-  num.delta / den.delta
    # sd
    xi <- sqrt(1 / (1 / xi0^2 + n / tau^2))
    #
    eta.vec[i] <- rnorm(1, delta, xi)

    mu.curr <- mu.vec[i]
    eta.curr <- eta.vec[i]

    # candidate for lambda and omega:
    # by default standard normal or if model.upd == "RW" it is a random walk step
    # or with hyperparameters if Gibbs
    if(model.upd=="RW"){
      mean.lambda.cand <- lambda.mat[(i-1),]
      mean.omega.cand <- omega.mat[(i-1),]
      sd.lambda.cand <- sd.lambda.cand
      sd.omega.cand <- sd.omega.cand
    } else if(model.upd=="Gibbs"){
      mean.lambda.cand <- mu.curr
      mean.omega.cand <- eta.curr
      sd.lambda.cand <- sigma
      sd.omega.cand <- tau
    } else{
      mean.lambda.cand <- mean.lambda.cand
      mean.omega.cand <- mean.omega.cand
      sd.lambda.cand <- sd.lambda.cand
      sd.omega.cand <- sd.omega.cand
    }
    #
    lambda.cand <- rnorm(n, mean.lambda.cand, sd.lambda.cand)
    omega.cand <- rnorm(n, mean.omega.cand, sd.omega.cand)


    # Exchange algorithm for lambda and omega
    par.post <- c()
    for(k in 1:n){
      par.post[k] <- a.exch(y[k], lambda.curr[k], omega.curr[k],
                            lambda.cand[k], omega.cand[k],
                            expo[k], mu.curr, sigma,
                            eta.curr, tau)$test
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
#-------------------------------------------------------
mcmc.dem.lo <- function(y, iter, expo, lambda0 = rep(0,length(y)),
                        omega0 = rep(0,length(y)),
                        mu, sigma, eta, tau,
                        mean.lambda.cand = 0, mean.omega.cand = 0,
                        sd.lambda.cand = 1, sd.omega.cand = 1,
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
    # or with hyperparameters if Gibbs
    if(model.upd=="RW"){
      mean.lambda.cand <- lambda.mat[(i-1),]
      mean.omega.cand <- omega.mat[(i-1),]
      sd.lambda.cand <- sd.lambda.cand
      sd.omega.cand <- sd.omega.cand
    } else if(model.upd=="Gibbs"){
      mean.lambda.cand <- mu
      mean.omega.cand <- eta
      sd.lambda.cand <- sigma
      sd.omega.cand <- tau
    } else{
      mean.lambda.cand <- mean.lambda.cand
      mean.omega.cand <- mean.omega.cand
      sd.lambda.cand <- sd.lambda.cand
      sd.omega.cand <- sd.omega.cand
    }
    #
    lambda.cand <- rnorm(n, mean.lambda.cand, sd.lambda.cand)
    omega.cand <- rnorm(n, mean.omega.cand, sd.omega.cand)


    # Exchange algorithm for lambda and omega
    par.post <- c()
    for(k in 1:n){
      par.post[k] <- a.exch(y[k], lambda.curr[k], omega.curr[k],
                            lambda.cand[k], omega.cand[k],
                            expo[k], mu, sigma,
                            eta, tau)$test
    }

    lambda.mat[i,] <- ifelse(par.post, lambda.cand, lambda.curr)
    omega.mat[i,] <- ifelse(par.post, omega.cand, omega.curr)
    acceptance <- ifelse(par.post, acceptance + 1, acceptance)

  }

  return(list (parameters = list(gamma = exp(lambda.mat), nu = exp(omega.mat)),
               acceptance = acceptance/iter)
  )
}


#---------------------------------------
# Only estimate omega
#---------------------------------------
mcmc.dem.om <- function(y, iter, expo, lambda, omega0 = rep(0,length(y)),
                        mu, sigma, eta, tau,
                        mean.omega.cand = 0,
                        sd.omega.cand = 1,
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
    # or with hyperparameters if Gibbs
    if(model.upd=="RW"){
      mean.omega.cand <- omega.mat[(i-1),]
      sd.omega.cand <- sd.omega.cand
    } else if(model.upd=="Gibbs"){
      mean.omega.cand <- eta
      sd.omega.cand <- tau
    } else{
      mean.omega.cand <- mean.omega.cand
      sd.omega.cand <- sd.omega.cand
    }
    #
    omega.cand <- rnorm(n, mean.omega.cand, sd.omega.cand)


    # Exchange algorithm for lambda and omega
    par.post <- c()
    for(k in 1:n){
      par.post[k] <- a.exch(y[k], lambda[k], omega.curr[k],
                            lambda[k], omega.cand[k],
                            expo[k], mu, sigma, eta, tau)$test
    }

    omega.mat[i,] <- ifelse(par.post, omega.cand, omega.curr)
    acceptance <- ifelse(par.post, acceptance + 1, acceptance)

  }

  return(list (parameters = list(nu = exp(omega.mat)),
               acceptance = acceptance/iter)
  )
}


#----------------------------------------------------
# Variance known, separate updating lambda and omega
#----------------------------------------------------
mcmc.dem.kvslom <- function(y, iter, expo, lambda0 = rep(0,length(y)), omega0 = rep(0,length(y)),
                            mu0 = 0, sigma, eta0 = 0, tau0 = 0,
                            alpha0 = 0, beta0 = 10, delta0 = 0, xi0 = 10,
                            mean.lambda.cand = 0, mean.omega.cand = 0,
                            sd.lambda.cand = 1, sd.omega.cand = 1,
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
    num.alpha <- alpha0 / beta0^2 + sum(lambda.mat[(i-1),])/sigma^2
    den.alpha <- (1 / beta0^2 + n / sigma^2)
    alpha <-  num.alpha / den.alpha
    # sd
    beta <- sqrt(1 / (1 / beta0^2 + n / sigma^2))
    #
    mu.vec[i] <- rnorm(1, alpha, beta)

    # Update eta
    # mean
    num.delta <- delta0 / xi0^2 + sum(omega.mat[(i-1),])/tau^2
    den.delta <- 1 / xi0^2 + 1 / tau^2
    delta <-  num.delta / den.delta
    # sd
    xi <- sqrt(1 / (1 / xi0^2 + 1 / tau^2))
    #
    eta.vec[i] <- rnorm(1, delta, xi)

    mu.curr <- mu.vec[i]
    eta.curr <- eta.vec[i]

    # candidate for lambda and omega:
    # by default standard normal or if model.upd == "RW" it is a random walk step
    # or with hyperparameters if Gibbs
    if(model.upd=="RW"){
      mean.lambda.cand <- lambda.mat[(i-1),]
      mean.omega.cand <- omega.mat[(i-1),]
      sd.lambda.cand <- sd.lambda.cand
      sd.omega.cand <- sd.omega.cand
    } else if (model.upd=="Gibbs") {
      mean.lambda.cand <- mu.curr
      mean.omega.cand <- eta.curr
      sd.lambda.cand <- sigma
      sd.omega.cand <- tau
    } else{
      mean.lambda.cand <- mean.lambda.cand
      mean.omega.cand <- mean.omega.cand
      sd.lambda.cand <- sd.lambda.cand
      sd.omega.cand <- sd.omega.cand
    }
  #
  lambda.cand <- rnorm(n, mean.lambda.cand, sd.lambda.cand)
  omega.cand <- rnorm(n, mean.omega.cand, sd.omega.cand)


    # Exchange algorithm for lambda
    par.post.lambda <- c()
    for(k in 1:n){
      par.post.lambda[k] <- a.exch(y[k], lambda.curr[k], omega.curr[k],
                                   lambda.cand[k], omega.curr[k],
                                   expo[k], mu.curr, sigma,
                                   eta.curr, tau)$test
    }

    lambda.mat[i,] <- ifelse(par.post.lambda, lambda.cand, lambda.curr)
    acceptance.lambda <- ifelse(par.post.lambda, acceptance.lambda + 1, acceptance.lambda)


    # Exchange algorithm for omega

    par.post.omega <- c()
    for(k in 1:n){
      par.post.omega[k] <- a.exch(y[k], lambda.curr[k], omega.curr[k],
                                  lambda.curr[k], omega.cand[k],
                                  expo[k], mu.curr, sigma,
                                  eta.curr, tau)$test
    }

    omega.mat[i,] <- ifelse(par.post.omega, omega.cand, omega.curr)
    acceptance.omega <- ifelse(par.post.omega, acceptance.omega + 1, acceptance.omega)

  }


  return(list (parameters = list(gamma = exp(lambda.mat), nu = exp(omega.mat),
                                 mu = mu.vec, eta = eta.vec),
               acceptance = list(gamma=acceptance.lambda/iter,
                                 nu = acceptance.omega/iter)))
}
