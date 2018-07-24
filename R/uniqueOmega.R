# For model type: y_i = CMP(mu_i, nu)

#-------------------------------------
# Only estimate lambda and omega
#-------------------------------------
mcmc.dem.lo.un <- function(y, iter, expo, lambda0 = rep(0,length(y)), omega0 = 0,
                         mu, sigma, eta, tau,
                         mean.lambda.cand = 0, mean.omega.cand = 0,
                         sd.lambda.cand = 1, sd.omega.cand = 1,
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
    # or with hyperparameters if Gibbs
    if(model.upd=="RW"){
      mean.lambda.cand <- lambda.mat[(i-1),]
      mean.omega.cand <- omega.vec[i-1]
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
    omega.cand <- rnorm(1, mean.omega.cand, sd.omega.cand)


    # Exchange algorithm for lambda and omega
    par.post.lambda <- c()
    for(k in 1:n){
      par.post.lambda[k] <- a.exch(y[k], lambda.curr[k], omega.curr,
                                   lambda.cand[k], omega.curr,
                                   expo[k], mu, sigma,
                                   eta, tau)$test
    }


    lambda.mat[i,] <- ifelse(par.post.lambda, lambda.cand, lambda.curr)
    acceptance.lambda <- ifelse(par.post.lambda, acceptance.lambda + 1,acceptance.lambda)

    par.post.omega <- a.exch1(y, lambda.curr, omega.curr,
                              lambda.curr, omega.cand,
                              expo, mu, sigma,
                              eta, tau)$test
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


#-------------------------------------
# Model assuming IG prior on Variance
# (i.e. Gibbs sampling for sigma and tau as conjugate model)
#-------------------------------------

mcmc.dem.IG.un <- function(y, iter, expo, lambda0 = rep(0,length(y)), omega0 = 0,
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
  omega.vec <- c()
  mu.vec <- c()
  sigma.vec <- c()
  eta.vec <- c()
  tau.vec <- c()
  
  # Initialise
  
  lambda.mat[1,]<- lambda0
  omega.vec[1] <- omega0
  mu.vec[1] <- mu0
  sigma.vec[1] <- sigma0
  eta.vec[1] <- eta0
  tau.vec[1] <- tau0
  
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
    num.delta <- delta0 / xi0^2 + sum(omega.vec[(i-1)])/(tau.vec[(i-1)])^2
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
    sc.tau <- sc.tau0 + 0.5 * sum((omega.vec[(i-1)] - eta.vec[i])^2)
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
      mean.omega.cand <- omega.vec[(i-1)]
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
    omega.cand <- rnorm(1, mean.omega.cand, sd.omega.cand)
    
    
    # Exchange algorithm for lambda and omega
    par.post.lambda <- c()
    for(k in 1:n){
      par.post.lambda[k] <- a.exch(y[k], lambda.curr[k], omega.curr,
                                   lambda.cand[k], omega.curr,
                                   expo[k], mu.curr, sigma.curr,
                                   eta.curr, tau.curr)$test
    }
    
    lambda.mat[i,] <- ifelse(par.post.lambda, lambda.cand, lambda.curr)
    acceptance.lambda <- ifelse(par.post.lambda, acceptance.lambda + 1, acceptance.lambda)
    lambda.curr <- lambda.mat[i,]
    
    par.post.omega <- c()
    par.post.omega <- a.exch1(y, lambda.curr, omega.curr,
                              lambda.curr, omega.cand,
                              expo, mu.curr, sigma.curr,
                              eta.curr, tau.curr)$test
    if(par.post.omega){
      omega.vec[i] <- omega.cand
      acceptance.omega <- acceptance.omega + 1
    } else {
      omega.vec[i] <- omega.curr
      acceptance.omega <- acceptance.omega
    }
    
  }
  
  return(list (parameters = list(gamma = exp(lambda.mat), nu = exp(omega.vec),
                                 mu = mu.vec, sigma = sigma.vec,
                                 eta = eta.vec, tau = tau.vec),
               acceptance = list(gamma=acceptance.lambda/iter,
                                 nu=acceptance.omega/iter))
  )
}

# Model assuming uninformative prior on Variance
# (i.e. Gibbs sampling for sigma and tau as conjugate model)
#-------------------------------------

mcmc.dem.UnInf.un <- function(y, iter, expo, lambda0 = rep(0,length(y)), omega0 = 0,
                           mu0, sigma0, eta0, tau0,
                           alpha0, beta0, delta0, xi0,
                            mean.lambda.cand, mean.omega.cand,
                           sd.lambda.cand, sd.omega.cand,
                           model.upd){
  
  n <- length(y)
  # Space for storing values
  lambda.mat <- matrix(NA, iter, n)
  omega.vec <- c()
  mu.vec <- c()
  sigma.vec <- c()
  eta.vec <- c()
  tau.vec <- c()
  
  # Initialise
  
  lambda.mat[1,]<- lambda0
  omega.vec[1] <- omega0
  mu.vec[1] <- mu0
  sigma.vec[1] <- sigma0
  eta.vec[1] <- eta0
  tau.vec[1] <- tau0
  
  # Acceptance
  
  acceptance.lambda <- 0
  acceptance.omega <- 0
  
  # Begin loop
  for(i in 2:iter){
    if(i%%(0.1*iter)==0) print(i)
    lambda.curr <- lambda.mat[(i-1),]
    omega.curr <- omega.vec[(i-1)]
    sigma.curr <- sigma.vec[(i-1)]
    tau.curr <- tau.vec[(i-1)]
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
    # generate new sigma 
    sigma.cand <- runif(1, sigma.curr - 0.5, sigma.curr + 0.5)
    if(sigma.cand > 0){
      lpost.sigma.curr <- dunif(sigma.curr, ) + sum(dnorm(lambda.mat[(i-1),], mu.vec[i], sigma.curr, log = TRUE)) - log(sigma.curr)
      lpost.sigma.cand <- sum(dnorm(lambda.mat[(i-1),], mu.vec[i], sigma.cand, log = TRUE)) - log(sigma.cand)
      diff.sigma <- lpost.sigma.cand - lpost.sigma.curr
      sigma.vec[i] <- ifelse(diff.sigma > runif(1), sigma.cand, sigma.curr)
    } else {
      sigma.vec[i] <- sigma.curr
    }
      
    # posterior and metropolis
    
    # Update eta
    # mean
    num.delta <- delta0 / xi0^2 + sum(omega.vec[(i-1)])/(tau.vec[(i-1)])^2
    den.delta <- 1 / xi0^2 + n / (tau.vec[(i-1)])^2
    delta <-  num.delta / den.delta
    # sd
    xi <- sqrt(1 / (1 / xi0^2 + n / (tau.vec[(i-1)])^2 ))
    #
    eta.vec[i] <- rnorm(1, delta, xi)
    
    # Update tau
    # generate new tau 
    tau.cand <- runif(1, sigma.curr - 0.5, sigma.curr + 0.5)
    # posterior and metropolis
    if(tau.cand > 0){
      lpost.tau.curr <- sum(dnorm(omega.vec[(i-1)], eta.vec[i], tau.curr, log = TRUE)) - log(tau.curr)
      lpost.tau.cand <- sum(dnorm(omega.vec[(i-1)], eta.vec[i], tau.cand, log = TRUE)) - log(tau.cand)
      diff.tau <- lpost.tau.cand - lpost.tau.curr
      tau.vec[i] <- ifelse(diff.tau > runif(1), tau.cand, tau.curr)
    } else {
      tau.vec[i] <- tau.curr
    }
    
    mu.curr <- mu.vec[i]
    sigma.curr <- sigma.vec[i]
    eta.curr <- eta.vec[i]
    tau.curr <- tau.vec[i]
    
    # candidate for lambda and omega:
    # by default standard normal or if model.upd == "RW" it is a random walk step
    # or with hyperparameters if Gibbs
    if(model.upd=="RW"){
      mean.lambda.cand <- lambda.mat[(i-1),]
      mean.omega.cand <- omega.vec[(i-1)]
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
    omega.cand <- rnorm(1, mean.omega.cand, sd.omega.cand)
    
    
    # Exchange algorithm for lambda and omega
    par.post.lambda <- c()
    for(k in 1:n){
      par.post.lambda[k] <- a.exch(y[k], lambda.curr[k], omega.curr,
                                   lambda.cand[k], omega.curr,
                                   expo[k], mu.curr, sigma.curr,
                                   eta.curr, tau.curr)$test
    }
    
    lambda.mat[i,] <- ifelse(par.post.lambda, lambda.cand, lambda.curr)
    acceptance.lambda <- ifelse(par.post.lambda, acceptance.lambda + 1, acceptance.lambda)
    lambda.curr <- lambda.mat[i,]
    
    par.post.omega <- c()
    par.post.omega <- a.exch1(y, lambda.curr, omega.curr,
                              lambda.curr, omega.cand,
                              expo, mu.curr, sigma.curr,
                              eta.curr, tau.curr)$test
    if(par.post.omega){
      omega.vec[i] <- omega.cand
      acceptance.omega <- acceptance.omega + 1
    } else {
      omega.vec[i] <- omega.curr
      acceptance.omega <- acceptance.omega
    }
    
  }
  
  return(list (parameters = list(gamma = exp(lambda.mat), nu = exp(omega.vec),
                                 mu = mu.vec, sigma = sigma.vec,
                                 eta = eta.vec, tau = tau.vec),
               acceptance = list(gamma=acceptance.lambda/iter,
                                 nu=acceptance.omega/iter))
  )
}
#----------------------------
# Model assuming uninformative prior on Variance and nu
# (i.e. Gibbs sampling for sigma and tau as conjugate model)
#-------------------------------------

mcmc.dem.UnInf.un.bis <- function(y, iter, expo, lambda0 = rep(0,length(y)), omega0 = 0,
                              mu0 = 0, sigma0 = 0.1, eta0 = 0, tau0 = 0.5,
                              alpha0 = 0, beta0 = 0.1, delta0 = 0, xi0 = 0.1,
                              mean.lambda.cand = 0, mean.omega.cand = 0,
                              sd.lambda.cand = 1, sd.omega.cand = 1,
                              model.upd){
  
  n <- length(y)
  # Space for storing values
  lambda.mat <- matrix(NA, iter, n)
  omega.vec <- c()
  mu.vec <- c()
  sigma.vec <- c()
  eta.vec <- c()
  tau.vec <- c()
  
  # Initialise
  
  lambda.mat[1,]<- lambda0
  omega.vec[1] <- omega0
  mu.vec[1] <- mu0
  sigma.vec[1] <- sigma0
  eta.vec[1] <- eta0
  tau.vec[1] <- tau0
  
  # Acceptance
  
  acceptance.lambda <- 0
  acceptance.omega <- 0
  
  # Begin loop
  for(i in 2:iter){
    if(i%%(0.1*iter)==0) print(i)
    lambda.curr <- lambda.mat[(i-1),]
    omega.curr <- omega.vec[(i-1)]
    sigma.curr <- sigma.vec[(i-1)]
    tau.curr <- tau.vec[(i-1)]
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
    # generate new sigma 
    sigma.cand <- runif(1, sigma.curr - 0.5, sigma.curr + 0.5)
    if(sigma.cand > 0){
      lpost.sigma.curr <- dunif(sigma.curr, ) + sum(dnorm(lambda.mat[(i-1),], mu.vec[i], sigma.curr, log = TRUE)) - log(sigma.curr)
      lpost.sigma.cand <- sum(dnorm(lambda.mat[(i-1),], mu.vec[i], sigma.cand, log = TRUE)) - log(sigma.cand)
      diff.sigma <- lpost.sigma.cand - lpost.sigma.curr
      sigma.vec[i] <- ifelse(diff.sigma > runif(1), sigma.cand, sigma.curr)
    } else {
      sigma.vec[i] <- sigma.curr
    }
    
    # posterior and metropolis
    
    # Update eta
    # mean
    num.delta <- delta0 / xi0^2 + sum(omega.vec[(i-1)])/(tau.vec[(i-1)])^2
    den.delta <- 1 / xi0^2 + n / (tau.vec[(i-1)])^2
    delta <-  num.delta / den.delta
    # sd
    xi <- sqrt(1 / (1 / xi0^2 + n / (tau.vec[(i-1)])^2 ))
    #
    eta.vec[i] <- rnorm(1, delta, xi)
    
    # Update tau
    # generate new tau 
    tau.cand <- runif(1, sigma.curr - 0.5, sigma.curr + 0.5)
    # posterior and metropolis
    if(tau.cand > 0){
      lpost.tau.curr <- sum(dnorm(omega.vec[(i-1)], eta.vec[i], tau.curr, log = TRUE)) - log(tau.curr)
      lpost.tau.cand <- sum(dnorm(omega.vec[(i-1)], eta.vec[i], tau.cand, log = TRUE)) - log(tau.cand)
      diff.tau <- lpost.tau.cand - lpost.tau.curr
      tau.vec[i] <- ifelse(diff.tau > runif(1), tau.cand, tau.curr)
    } else {
      tau.vec[i] <- tau.curr
    }
    
    mu.curr <- mu.vec[i]
    sigma.curr <- sigma.vec[i]
    eta.curr <- eta.vec[i]
    tau.curr <- tau.vec[i]
    
    # candidate for lambda and omega:
    # by default standard normal or if model.upd == "RW" it is a random walk step
    # or with hyperparameters if Gibbs
    if(model.upd=="RW"){
      mean.lambda.cand <- lambda.mat[(i-1),]
      mean.omega.cand <- omega.vec[(i-1)]
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
    omega.cand <- rnorm(1, mean.omega.cand, sd.omega.cand)
    
    
    # Exchange algorithm for lambda and omega
    par.post.lambda <- c()
    for(k in 1:n){
      par.post.lambda[k] <- a.exch(y[k], lambda.curr[k], omega.curr,
                                   lambda.cand[k], omega.curr,
                                   expo[k], mu.curr, sigma.curr,
                                   eta.curr, tau.curr)$test
    }
    
    lambda.mat[i,] <- ifelse(par.post.lambda, lambda.cand, lambda.curr)
    acceptance.lambda <- ifelse(par.post.lambda, acceptance.lambda + 1, acceptance.lambda)
    lambda.curr <- lambda.mat[i,]
    
    par.post.omega <- c()
    par.post.omega <- a.exch.uninf.nu(y, lambda.curr, omega.curr,
                                      lambda.curr, omega.cand,
                                      expo, mu.curr, sigma.curr,
                                      eta.curr, tau.curr)$test
    if(par.post.omega){
      omega.vec[i] <- omega.cand
      acceptance.omega <- acceptance.omega + 1
    } else {
      omega.vec[i] <- omega.curr
      acceptance.omega <- acceptance.omega
    }
    
  }
  
  return(list (parameters = list(gamma = exp(lambda.mat), nu = exp(omega.vec),
                                 mu = mu.vec, sigma = sigma.vec,
                                 eta = eta.vec, tau = tau.vec),
               acceptance = list(gamma=acceptance.lambda/iter,
                                 nu=acceptance.omega/iter))
  )
}


#-------------------------------------
# Only estimate unique omega
#-------------------------------------
mcmc.dem.om.un <- function(y, iter, expo, lambda, omega0 = 0,
                         mu, sigma, eta, tau,
                         mean.omega.cand = 0,
                         sd.omega.cand = 1,
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
    # or with hyperparameters if Gibbs
    if(model.upd=="RW"){
      mean.omega.cand <- omega.vec[i-1]
      sd.omega.cand <- sd.omega.cand
    } else if(model.upd=="Gibbs"){
      mean.omega.cand <- eta
      sd.omega.cand <- tau
    } else{
      mean.omega.cand <- mean.omega.cand
      sd.omega.cand <- sd.omega.cand
    }

    #
    omega.cand <- rnorm(1, eta, sd.omega.cand)


    # Exchange algorithm for lambda and omega
    par.post <- a.exch1(y, lambda, omega.curr,
                        lambda, omega.cand,
                        expo, mu, sigma, eta, tau)$test

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
#------------------------------------------------------------------
mcmc.dem.kvslom.un <- function(y, iter, expo, lambda0 = rep(0,length(y)), omega0 = 0,
                             mu0 = 0, sigma, eta0 = 0, tau,
                             alpha0 = 0, beta0 = 10, delta0 = 0, xi0 = 10,
                             mean.lambda.cand = 0, mean.omega.cand = 0,
                             sd.lambda.cand = 1, sd.omega.cand = 1,
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
    num.alpha <- alpha0 / beta0^2 + sum(lambda.mat[(i-1),])/sigma^2
    den.alpha <- (1 / beta0^2 + n / sigma^2)
    alpha <-  num.alpha / den.alpha
    # sd
    beta <- sqrt(1 / (1 / beta0^2 + n / sigma^2))
    #
    mu.vec[i] <- rnorm(1, alpha, beta)

    # Update eta
    # mean
    num.delta <- delta0 / xi0^2 + sum(omega.vec[(i-1)])/tau^2
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
      mean.omega.cand <- omega.vec[i-1]
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
    omega.cand <- rnorm(1, mean.omega.cand, sd.omega.cand)


    # Exchange algorithm for lambda
    par.post.lambda <- c()
    for(k in 1:n){
      par.post.lambda[k] <- a.exch(y[k], lambda.curr[k], omega.curr,
                                   lambda.cand[k], omega.curr,
                                   expo[k], mu.curr, sigma,
                                   eta.curr, tau)$test
    }

    lambda.mat[i,] <- ifelse(par.post.lambda, lambda.cand, lambda.curr)
    acceptance.lambda <- ifelse(par.post.lambda, acceptance.lambda + 1, acceptance.lambda)


    # Exchange algorithm for omega
    par.post.omega <- a.exch1(y, lambda.curr, omega.curr,
                              lambda.curr, omega.cand,
                              expo, mu.curr, sigma,
                              eta.curr, tau)$test

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

#------------------------------------------------------------------
# Variance known, unique omega, separate updating lambda and omega
#------------------------------------------------------------------
mcmc.dem.kvslom.un.bis <- function(y, iter, expo, lambda0 = rep(0,length(y)), omega0 = 0,
                               mu0 = 0, sigma, eta0 = 0, tau,
                               alpha0 = 0, beta0 = 10, delta0 = 0, xi0 = 10,
                               mean.lambda.cand = 0, mean.omega.cand = 0,
                               sd.lambda.cand = 1, sd.omega.cand = 1,
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
    num.alpha <- alpha0 / beta0^2 + sum(lambda.mat[(i-1),])/sigma^2
    den.alpha <- (1 / beta0^2 + n / sigma^2)
    alpha <-  num.alpha / den.alpha
    # sd
    beta <- sqrt(1 / (1 / beta0^2 + n / sigma^2))
    #
    mu.vec[i] <- rnorm(1, alpha, beta)

    # Update eta
    # mean
    num.delta <- delta0 / xi0^2 + sum(omega.vec[(i-1)])/tau^2
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
      mean.omega.cand <- omega.vec[i-1]
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
    omega.cand <- rnorm(1, mean.omega.cand, sd.omega.cand)


    # Exchange algorithm for lambda
    par.post.lambda <- c()
    for(k in 1:n){
      par.post.lambda[k] <- a.exch(y[k], lambda.curr[k], omega.curr,
                                   lambda.cand[k], omega.curr,
                                   expo[k], mu.curr, sigma,
                                   eta.curr, tau)$test
    }

    lambda.mat[i,] <- ifelse(par.post.lambda, lambda.cand, lambda.curr)
    acceptance.lambda <- ifelse(par.post.lambda, acceptance.lambda + 1, acceptance.lambda)


    # Exchange algorithm for omega
    par.post.omega <- a.exch1.bis(y, lambda.curr, omega.curr,
                              lambda.curr, omega.cand,
                              expo, mu.curr, sigma,
                              eta.curr, tau)$test

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


#-------------------------------------------------------------------
# Variance and second level known, double updating lambda and omega
#-------------------------------------------------------------------
mcmc.dem.double <- function(y, iter, expo, lambda0 = rep(0,length(y)), omega0 = rep(0,length(y)),
                            mu0 = 0, sigma, eta0 = 0, tau,
                            alpha0 = 0, beta0 = 1, delta0 = 0, xi0 = 1,
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
      mean.omega.cand <- omega.vec[i-1]
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

    # Exchange algorithm for lambda

    par.post.lambda <- c()
    for(k in 1:n){
      par.post.lambda[k] <- a.exch1(y[k], lambda.curr[k], omega.curr[k],
                                    lambda.cand[k], omega.curr[k], expo[k],
                                    mu.curr, sigma, eta.curr, tau)$test
      }

    lambda.mat[i,] <- ifelse(par.post.lambda, lambda.cand, lambda.curr)
    acceptance.lambda <- ifelse(par.post.lambda, acceptance.lambda + 1, acceptance.lambda)

    # Exchange algorithm for omega

    par.post.omega <- c()
    for(k in 1:n){
      par.post.omega[k] <- a.exch1(y[k], lambda.curr[k], omega.curr[k],
                                    lambda.curr[k], omega.cand[k], expo[k],
                                    mu.curr, sigma, eta.curr, tau)$test
    }
    omega.mat[i,] <- ifelse(par.post.omega, omega.cand, omega.curr)
    acceptance.omega <- ifelse(par.post.omega, acceptance.omega + 1, acceptance.omega)

  }


  return(list (parameters = list(gamma = exp(lambda.mat), nu = exp(omega.mat),
                                 mu = mu.vec, eta = eta.vec),
               acceptance = list(couple = acceptance/iter,
                                 lambda = acceptance.lambda/iter,
                                 omega = acceptance.omega/iter)))
}

