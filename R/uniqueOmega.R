# For model type: y_i = CMP(mu_i, nu)

#-------------------------------------
# Only estimate lambda and omega
#-------------------------------------

mcmc.dem.lo.un <- function(y, iter, expo, lambda0 = rep(0,length(y)), omega0 = 0,
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
    # or with hyperparameters if Gibbs
    if(model.upd=="RW"){
      mean.lambda.cand <- lambda.mat[(i-1),]
      mean.omega.cand <- omega.vec[i-1]
      var.lambda.cand <- var.lambda.cand
      var.omega.cand <- var.omega.cand
    } else {
      if(model.upd=="Gibbs"){
        mean.lambda.cand <- mu.curr
        mean.omega.cand <- eta.curr
        var.lambda.cand <- sigma.curr
        var.omega.cand <- tau.curr
      } else{
        mean.lambda.cand <- mean.lambda.cand
        mean.omega.cand <- mean.omega.cand
        var.lambda.cand <- var.lambda.cand
        var.omega.cand <- var.omega.cand
      }
    }
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

#-------------------------------------
# Only estimate unique omega
#-------------------------------------

mcmc.dem.om.un <- function(y, iter, expo, lambda, omega0 = 0,
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
    # or with hyperparameters if Gibbs
    if(model.upd=="RW"){
      mean.omega.cand <- omega.vec[i-1]
      var.omega.cand <- var.omega.cand
    } else {
      if(model.upd=="Gibbs"){
        mean.omega.cand <- eta
        var.omega.cand <- tau
      } else{
        mean.omega.cand <- mean.omega.cand
        var.omega.cand <- var.omega.cand
      }
    }

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
#------------------------------------------------------------------

mcmc.dem.kvslom.un <- function(y, iter, expo, lambda0 = rep(0,length(y)), omega0 = 0,
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
    # or with hyperparameters if Gibbs
    if(model.upd=="RW"){
      mean.lambda.cand <- lambda.mat[(i-1),]
      mean.omega.cand <- omega.vec[i-1]
      var.lambda.cand <- var.lambda.cand
      var.omega.cand <- var.omega.cand
    } else {
      if(model.upd=="Gibbs"){
        mean.lambda.cand <- mu.curr
        mean.omega.cand <- eta.curr
        var.lambda.cand <- sigma.curr
        var.omega.cand <- tau.curr
      } else{
        mean.lambda.cand <- mean.lambda.cand
        mean.omega.cand <- mean.omega.cand
        var.lambda.cand <- var.lambda.cand
        var.omega.cand <- var.omega.cand
      }
    }
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



#-------------------------------------------------------------------
# Variance and second level known, double updating lambda and omega
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
    # or with hyperparameters if Gibbs
    if(model.upd=="RW"){
      mean.lambda.cand <- lambda.mat[(i-1),]
      mean.omega.cand <- omega.vec[i-1]
      var.lambda.cand <- var.lambda.cand
      var.omega.cand <- var.omega.cand
    } else {
      if(model.upd=="Gibbs"){
        mean.lambda.cand <- mu.curr
        mean.omega.cand <- eta.curr
        var.lambda.cand <- sigma.curr
        var.omega.cand <- tau.curr
      } else{
        mean.lambda.cand <- mean.lambda.cand
        mean.omega.cand <- mean.omega.cand
        var.lambda.cand <- var.lambda.cand
        var.omega.cand <- var.omega.cand
      }
    }


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

