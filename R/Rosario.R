#--------------------------------------------
# MCMC For model type: y_i = CMP(mu_i, nu_i)
# Uniform prior on nu_i
#--------------------------------------------
library(LaplacesDemon) # to use half t distribution
#---------------------------------
# First simple model, sd = half t
#---------------------------------

mcmc.dem.simple.unif <- function(y, iter, expo, gamma0 = rep(1,length(y)), nu0 = rep(1,length(y)),
                                 lgamma0 = log(gamma0), lnu0 = log(nu0),
                                 mu = 0, sigma = 10,
                                 lnumin = -9, lnumax = 3){
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
    lnu.cand <- runif(n, lnumin, lnumax)
    
    # Exchange algorithm for gamma and nu
    par.post<- c()
    for(k in 1:n){
      par.post[k] <- a.exch.unif(y[k], lgamma.curr[k], lnu.curr[k],
                                 lgamma.cand[k], lnu.cand[k], expo[k], mu,
                                 sigma)$test
    }
    
    gamma.mat[i,] <- ifelse(par.post, exp(lgamma.cand), exp(lgamma.curr))
    nu.mat[i,] <- ifelse(par.post, exp(lnu.cand), exp(lnu.curr))
    gamma.nu.acc <- ifelse(par.post, gamma.nu.acc+1, gamma.nu.acc)
    lgamma.curr <- ifelse(par.post, lgamma.cand, lgamma.curr)
    lnu.curr <- ifelse(par.post, lnu.cand, lnu.curr)
    
    # candidate for c(gamma, nu):
    
    lgamma.cand2 <- rnorm(n, lgamma.curr, 0.1)
    lnu.cand2 <- runif(n, lnumin, lnumax)
    
    # Exchange algorithm for gamma:
    par.post.g<- c()
    for(k in 1:n){
      par.post.g[k] <- a.exch.unif(y[k], lgamma.curr[k], lnu.curr[k],
                                   lgamma.cand2[k], lnu.curr[k], expo[k],
                                   mu, sigma)$test
    }
    
    gamma.mat[i,] <- ifelse(par.post.g, exp(lgamma.cand2), exp(lgamma.curr))
    gamma.acc <- ifelse(par.post.g, gamma.acc+1, gamma.acc)
    
    # Exchange algorithm for nu:
    par.post.n <- c()
    for(k in 1:n){
      par.post.n[k] <- a.exch.unif(y[k], lgamma.curr[k], lnu.curr[k],
                                   lgamma.curr[k], lnu.cand2[k], expo[k],
                                   mu, sigma)$test
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

mcmc.dem.IG.unif <- function(y, iter, expo, lambda0 = rep(0,length(y)), omega0 = rep(0,length(y)),
                             mu0 = 0, sigma0 = 0.1,
                             alpha0 = 0, beta0 = 0.1,
                             sh.sigma0 = 0.001, sc.sigma0 = 0.001,
                             mean.lambda.cand = 0, sd.lambda.cand = 1,
                             model.upd, ommin = 0.0001, ommax = 5){
  
  n <- length(y)
  # Space for storing values
  lambda.mat <- matrix(NA, iter, n)
  omega.mat <- matrix(NA, iter, n)
  mu.vec <- c()
  sigma.vec <- c()
  
  # Initialise
  
  lambda.mat[1,]<- lambda0
  omega.mat[1,] <- omega0
  mu.vec[1] <- mu0
  sigma.vec[1] <- sigma0
  
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
    
    mu.curr <- mu.vec[i]
    sigma.curr <- sigma.vec[i]
    
    # candidate for lambda and omega:
    # by default standard normal or if model.upd == "RW" it is a random walk step
    # or with hyperparameters if Gibbs
    if(model.upd=="RW"){
      mean.lambda.cand <- lambda.mat[(i-1),]
      sd.lambda.cand <- sd.lambda.cand
      omega.cand <- log(runif(n, ommin, ommax))
    } else{
      mean.lambda.cand <- mean.lambda.cand
      sd.lambda.cand <- sd.lambda.cand
      omega.cand <- log(runif(n, ommin, ommax))
    }
    
    
    #
    lambda.cand<- rnorm(n, mean.lambda.cand, sd.lambda.cand)
    
    
    
    # Exchange algorithm for lambda and omega
    par.post <- c()
    for(k in 1:n){
      par.post[k] <- a.exch.unif(y[k], lambda.curr[k], omega.curr[k],
                                 lambda.cand[k], omega.cand[k],
                                 expo[k], mu.curr, sigma.curr, ommin, ommax)$test
    }
    
    lambda.mat[i,] <- ifelse(par.post, lambda.cand, lambda.curr)
    omega.mat[i,] <- ifelse(par.post, omega.cand, omega.curr)
    acceptance <- ifelse(par.post, acceptance + 1, acceptance)
    
  }
  
  return(list (parameters = list(gamma = exp(lambda.mat), nu = exp(omega.mat),
                                 mu = mu.vec, sigma = sigma.vec),
               acceptance = acceptance/iter)
  )
}

#-------------------------------------
# Model assuming IG prior on Variance separate updating
# (i.e. Gibbs sampling for sigma and tau as conjugate model)
#-------------------------------------
#---------------------------------------------
#
mcmc.dem.IG.unif1 <- function(y, iter, expo, lambda0 = 0, nu0 = 0,
                           mu0 = 0, sigma0 = 0.1,
                           alpha0 = 0, beta0 = 0.1,
                           sh.sigma0, sc.sigma0, 
                           in.sd.lambda = 1, in.sd.nu = 1){
  
  n <- length(y)
  # Space for storing values
  lambda.vec <- c()
  nu.vec <- c()
  mu.vec <- c()
  sigma.vec <- c()
  
  # Initialise
  
  lambda.vec[1]<- lambda0
  nu.vec[1] <- nu0
  mu.vec[1] <- mu0
  sigma.vec[1] <- sigma0
  
  # Acceptance
  
  acceptance <- 0
  
  # Begin loop
  for(i in 2:iter){
    if(i%%(0.1*iter)==0) print(i)
    lambda.curr <- lambda.vec[(i-1)]
    nu.curr <- nu.vec[(i-1)]
    #
    # Parameters for Gibbs sampling
    
    # Update mu
    # mean
    num.alpha <- alpha0 / beta0^2 + sum(lambda.vec[(i-1)])/(sigma.vec[(i-1)])^2
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
    sc.sigma <- sc.sigma0 + 0.5 * sum((lambda.vec[(i-1)] - mu.vec[i])^2)
    #
    sigma.vec[i] <- sqrt(1 / rgamma(1, shape = sh.sigma, scale = sc.sigma))
    
    
    mu.curr <- mu.vec[i]
    sigma.curr <- sigma.vec[i]
    
    # candidate for lambda and nu:
    # by default standard normal or if model.upd == "RW" it is a random walk step
    # or with hyperparameters if Gibbs
    
    if(i<500){
           q.var.lambda <- c(in.sd.lambda,0)
             q.var.nu <- c(0,in.sd.nu)
           }else{
             q.var.lambda <- c(var(lambda.vec[(i-201):(i-1)]),cov(nu.vec[(i-201):(i-1)],lambda.vec[(i-201):(i-1)]))
             q.var.nu <- c(cov(nu.vec[(i-201):(i-1),],lambda.vec[(i-201):(i-1)]),var(nu.vec[(i-201):(i-1)])) 
           }
           q.mean <- c(lambda.curr,nu.curr)
         q.var.matrix <- cbind(q.var.lambda,q.var.nu)
           #Proposal
           
           Q <- rmvnorm(1,q.mean, q.var.matrix)
         lambda.cand <- Q[1]
           nu.cand <- Q[2]
           if(nu.cand < 0 | nu.cand > 20){
             lambda.vec[i] <- lambda.curr
             nu.vec[i] <- nu.curr
           acceptance <- acceptance 
             } else {
             #ACCEPTANCE STEP
             
             # Exchange algorithm for gamma and nu
             
             par.post <- a.exch1.unif2(y, lambda.curr, nu.curr,
                                       lambda.cand, nu.cand, expo, mu,
                                       sigma, ommin, ommax)$test
             
             
             lambda.vec[i] <- ifelse(par.post, lambda.cand), lambda.curr))
             nu.vec[i] <- ifelse(par.post, nu.cand, nu.curr)
             acceptance <- ifelse(par.post, acceptance+1, acceptance)
           }
    
  }
  
  return(list (parameters = list(gamma = exp(lambda.vec), nu = nu.vec,
                                 mu = mu.vec, sigma = sigma.vec),
               acceptance = list(acceptance = acceptance/iter))
  )
}


#------------------------------------
# No omega only nu
#------------------------------------



a.exch.unif2 <- function(y, lambda.curr, nu.curr, lambda.cand, nu.cand, q, mu,
                         sigma, ommin=0.001, ommax=10, max = 10^4){
  # CMP parameters
  mu.curr <- exp(lambda.curr)*q
  #
  mu.cand <- exp(lambda.cand)*q
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
                     -log(nu.cand),
                     log(nu.curr)
                     #dunif(nu.cand, ommin, ommax, log=TRUE),
                     #- dunif(nu.curr, ommin, ommax, log = TRUE)
                     # dgamma(nu.cand, shape = 0.15, rate = .075, log=TRUE),
                     # - dgamma(nu.curr, shape = .15, rate = .075, log=TRUE)
                     )
  # Posterior
  post <- llik.ratio + prior.ratio
  test <- post > log(runif(1))
  return(list(postdens = post, test = test))

  }



library(mvtnorm)

mcmc.dem.simple.unif2 <- function(y, iter, expo, gamma0 = rep(1,length(y)), nu0 = rep(1,length(y)),
                                  lgamma0 = log(gamma0),
                                  mu = 0, sigma = 10,
                                  ommin = 0.0001, ommax = 10,
                                  in.sd.lgamma, in.sd.nu,
                                  model.upd){
  library(mvtnorm)
  n <-length(y)
  gamma.mat <- matrix(NA, iter, n)
  nu.mat <- matrix(NA, iter, n)
  
  
  # Initialise
  
  gamma.mat[1,] <-gamma0
  nu.mat[1,] <- nu0
  
  # Acceptance
  
  gamma.nu.acc <- 0
  
  # nu and gamma are updated with the exchange algorithm
  # for mu and eta Gibbs sampling
  # sigma and tau MH?
  for(i in 2:iter){
    
    if(i%%(0.1*iter)==0) print(i)
    par.post<- matrix(NA, ncol=n)
    for(obs in 1:n){
    gamma.curr <- gamma.mat[i-1,obs]
    nu.curr <- nu.mat[i-1,obs]
    lgamma.curr<- log(gamma.curr)
    
    # generate candidates
    # Sample from a bivariate normal distribution
    if(i<500){
    q.var.lgamma <- c(in.sd.lgamma,0)
    q.var.nu <- c(0,in.sd.nu)
    }else{
      q.var.lgamma <- c(var(gamma.mat[(i-201):(i-1),obs]),cov(nu.mat[(i-201):(i-1),obs],gamma.mat[(i-201):(i-1),obs]))
      q.var.nu <- c(cov(nu.mat[(i-201):(i-1),obs],gamma.mat[(i-201):(i-1),obs]),var(nu.mat[(i-201):(i-1),obs])) 
      }
    q.mean <- c(lgamma.curr,nu.curr)
    q.var.matrix <- cbind(q.var.lgamma,q.var.nu)
    #Proposal
    
    Q <- rmvnorm(1,q.mean, q.var.matrix)
    lgamma.cand <- Q[1]
    nu.cand <- Q[2]
    if(nu.cand < 0 | nu.cand > 40){
      gamma.mat[i,obs] <- exp(lgamma.curr)
      nu.mat[i,obs] <- nu.curr
      gamma.nu.acc <- gamma.nu.acc
    } else {
    #ACCEPTANCE STEP
    
    # Exchange algorithm for gamma and nu
  
      par.post[obs] <- a.exch.unif2(y[obs], lgamma.curr, nu.curr,
                                  lgamma.cand, nu.cand, expo[obs], mu,
                                  sigma)$test

    
    gamma.mat[i,obs] <- ifelse(par.post[obs], exp(lgamma.cand), exp(lgamma.curr))
    nu.mat[i,obs] <- ifelse(par.post[obs], nu.cand, nu.curr)
gamma.nu.acc <- ifelse(par.post[obs], gamma.nu.acc+1, gamma.nu.acc)
}

  }
  }
  return(list (parameters = list(gamma = gamma.mat, nu = nu.mat),
               acceptance = list(gamma.nu = gamma.nu.acc/(iter*n)))
         )
}




a.exch1.unif2 <- function(y, lambda.curr, nu.curr, lambda.cand, nu.cand, q, mu,
                          sigma, ommin, ommax, max = 10^4){
  n<- length(y)
  # CMP parameters
  #
  mu.curr <- exp(lambda.curr)*q
  #
  mu.cand <- exp(lambda.cand)*q
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
  prior.ratio <- sum(dnorm(lambda.cand,mu,sigma, log = TRUE),
                     -dnorm(lambda.curr,mu,sigma, log = TRUE),
                     #dunif(nu.cand, ommin, ommax),
                     #- dunif(nu.curr, ommin, ommax)
                     -log(nu.cand),
                     log(nu.curr))
# Posterior
post <- llik.ratio + prior.ratio
test <- post > log(runif(1))
return(list(postdens = post, test = test))
}

a.exch1.comp <- function(y, lambda.curr, nu.curr, lambda.cand, nu.cand, q, mu,
                          sigma, ommin, ommax, max = 10^4){
  n<- length(y)
  # CMP parameters
  #
  mu.curr <- exp(lambda.curr)*q
  gamma.curr <- mu.curr^nu.curr

  #
  mu.cand <- exp(lambda.cand)*q
  gamma.cand <- mu.cand^nu.cand
  #
  # Generate the auxiliary data
  
  #
  # Log-likelihood ratio
  lik.cand <- c()
  lik.curr <- c()
  for(i in 1:n){
  lik.cand[i] <- log(dcom(x = y[i], lambda = gamma.cand[i], nu = nu.cand))
  lik.curr[i] <- log(dcom(x = y[i], lambda = gamma.curr[i], nu = nu.curr))
  }
  llik.ratio <- sum(lik.cand) - sum(lik.curr)
  # Prior ratio
  prior.ratio <- sum(dnorm(lambda.cand,mu,sigma, log = TRUE),
                     -dnorm(lambda.curr,mu,sigma, log = TRUE),
                     #dunif(nu.cand, ommin, ommax),
                     #- dunif(nu.curr, ommin, ommax)
                     -log(nu.cand),
                     log(nu.curr))
  # Posterior
  post <- llik.ratio + prior.ratio
  test <- post > log(runif(1))
  return(list(postdens = post, test = test))
}



mcmc.dem.simple.unif12 <- function(y, iter, expo, gamma0 = 1, nu0 = 1,
                                  lgamma0 = log(gamma0),
                                  mu = 0, sigma = 10,
                                  ommin = 0.0001, ommax = 10,
                                  in.sd.lgamma, in.sd.nu,
                                  model.upd){
  library(mvtnorm)
  n <-length(y)
  gamma.mat <- matrix(NA, iter, 1)
  nu.mat <- matrix(NA, iter, 1)
  
  
  # Initialise
  
  gamma.mat[1,] <- gamma0
  nu.mat[1,] <- nu0
  
  # Acceptance
  
  gamma.nu.acc <- 0
  
  # nu and gamma are updated with the exchange algorithm
  # for mu and eta Gibbs sampling
  # sigma and tau MH?
  for(i in 2:iter){
    
    if(i%%(0.1*iter)==0) print(i)
    #par.post<- matrix(NA, ncol=n)
    
      gamma.curr <- gamma.mat[i-1,]
      nu.curr <- nu.mat[i-1,]
      lgamma.curr<- log(gamma.curr)
      
      # generate candidates
      # Sample from a bivariate normal distribution
      if(i<500){
        q.var.lgamma <- c(in.sd.lgamma,0)
        q.var.nu <- c(0,in.sd.nu)
      }else{
        q.var.lgamma <- c(var(gamma.mat[(i-201):(i-1),]),cov(nu.mat[(i-201):(i-1),],gamma.mat[(i-201):(i-1),]))
        q.var.nu <- c(cov(nu.mat[(i-201):(i-1),],gamma.mat[(i-201):(i-1),]),var(nu.mat[(i-201):(i-1),])) 
      }
      q.mean <- c(lgamma.curr,nu.curr)
      q.var.matrix <- cbind(q.var.lgamma,q.var.nu)
      #Proposal
      Q <- rmvnorm(1,q.mean, q.var.matrix)
      lgamma.cand <- Q[1]
      nu.cand <- Q[2]
      if(nu.cand < 0 | nu.cand > 20){
        gamma.mat[i,] <- exp(lgamma.curr)
        nu.mat[i,] <- nu.curr
        gamma.nu.acc <- gamma.nu.acc
      } else {
        #ACCEPTANCE STEP
        
        # Exchange algorithm for gamma and nu
        
        par.post <- a.exch1.unif2(y, lgamma.curr, nu.curr,
                                      lgamma.cand, nu.cand, expo, mu,
                                      sigma, ommin, ommax)$test
        
        
        gamma.mat[i,] <- ifelse(par.post, exp(lgamma.cand), exp(lgamma.curr))
        nu.mat[i,] <- ifelse(par.post, nu.cand, nu.curr)
        gamma.nu.acc <- ifelse(par.post, gamma.nu.acc+1, gamma.nu.acc)
      }
      
    }
  return(list (parameters = list(gamma = gamma.mat, nu = nu.mat),
               acceptance = list(gamma.nu = gamma.nu.acc/iter))
  )
}


mcmc.dem.simple.comp <- function(y, iter, expo, gamma0 = 1, nu0 = 1,
                                   lgamma0 = log(gamma0),
                                   mu = 0, sigma = 10,
                                   ommin = 0.0001, ommax = 10,
                                   in.sd.lgamma, in.sd.nu,
                                   model.upd){
  library(mvtnorm)
  n <-length(y)
  gamma.mat <- matrix(NA, iter, 1)
  nu.mat <- matrix(NA, iter, 1)
  
  
  # Initialise
  
  gamma.mat[1,] <- gamma0
  nu.mat[1,] <- nu0
  
  # Acceptance
  
  gamma.nu.acc <- 0
  
  # nu and gamma are updated with the exchange algorithm
  # for mu and eta Gibbs sampling
  # sigma and tau MH?
  for(i in 2:iter){
    
    if(i%%(0.1*iter)==0) print(i)
    #par.post<- matrix(NA, ncol=n)
    
    gamma.curr <- gamma.mat[i-1,]
    nu.curr <- nu.mat[i-1,]
    lgamma.curr<- log(gamma.curr)
    
    # generate candidates
    # Sample from a bivariate normal distribution
    if(i<500){
      q.var.lgamma <- c(in.sd.lgamma,0)
      q.var.nu <- c(0,in.sd.nu)
    } else {
      q.var.lgamma <- c(var(gamma.mat[(i-201):(i-1),]),cov(nu.mat[(i-201):(i-1),],gamma.mat[(i-201):(i-1),]))
      q.var.nu <- c(cov(nu.mat[(i-201):(i-1),],gamma.mat[(i-201):(i-1),]),var(nu.mat[(i-201):(i-1),])) 
    }
    q.mean <- c(lgamma.curr,nu.curr)
    q.var.matrix <- cbind(q.var.lgamma,q.var.nu)
    #Proposal
    const = (2.38^2)/2 #al denominatore c'è la dimensione della distribuzione multivariata
    Q <- rmvnorm(1,q.mean, const*q.var.matrix)
    lgamma.cand <- Q[1]
    nu.cand <- Q[2]
    
    if(nu.cand < 0 | nu.cand > 20){
      gamma.mat[i,] <- exp(lgamma.curr)
      nu.mat[i,] <- nu.curr
      gamma.nu.acc <- gamma.nu.acc
    } else {
      #ACCEPTANCE STEP
      
      # Exchange algorithm for gamma and nu
      
      par.post <- a.exch1.comp(y, lgamma.curr, nu.curr,
                                lgamma.cand, nu.cand, expo, mu,
                                sigma, ommin, ommax)$test
      
      
      gamma.mat[i,] <- ifelse(par.post, exp(lgamma.cand), exp(lgamma.curr))
      nu.mat[i,] <- ifelse(par.post, nu.cand, nu.curr)
      gamma.nu.acc <- ifelse(par.post, gamma.nu.acc+1, gamma.nu.acc)
    }
    
  }
  return(list (parameters = list(gamma = gamma.mat, nu = nu.mat),
               acceptance = list(gamma.nu = gamma.nu.acc/iter))
  )
}

#-------------------------------------
# Model assuming IG prior on Variance
# (i.e. Gibbs sampling for sigma and tau as conjugate model)
#-------------------------------------

a.exch.new <- function(y, lambda.curr, omega.curr, lambda.cand, omega.cand, q, mu,
                         sigma, ommin=0.001, ommax=10, max = 10^4){
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
                     - dnorm(omega.curr, eta, tau, log=TRUE)
                     #-log(nu.cand),
                     #log(nu.curr)
                     #dunif(nu.cand, ommin, ommax, log=TRUE),
                     #- dunif(nu.curr, ommin, ommax, log = TRUE)
                     # dgamma(nu.cand, shape = 0.15, rate = .075, log=TRUE),
                     # - dgamma(nu.curr, shape = .15, rate = .075, log=TRUE)
  )
  # Posterior
  post <- llik.ratio + prior.ratio
  test <- post > log(runif(1))
  return(list(postdens = post, test = test))
  
}




mcmc.dem.IG.new <- function(y, iter, expo, lambda0 = rep(0,length(y)), omega0 = rep(0,length(y)),
                        mu0 = 0, sigma0 = 1, eta0 = 0, tau0 = 0.5,
                        alpha0 = 0, beta0 = 10, delta0 = 0, xi0 = 10,
                        sh.sigma0 = 3.5, sh.tau0 = 3.5,
                        sc.sigma0 = 5, sc.tau0 = 5,
                        mean.lambda.cand = 0, mean.omega.cand = 0,
                        sd.lambda.cand = 1, sd.omega.cand = 1,
                        model.upd, in.sd.lambda, in.sd.omega){
  library(mvtnorm)
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
    par.post<- matrix(NA, ncol=n)
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
    
    for(obs in 1:n){
      
      lambda.curr.un <- lambda.mat[i-1,obs]
      omega.curr.un <- omega.mat[i-1,obs]
      gamma.curr.un <- exp(lambda.curr)
      
      # generate candidates
      # Sample from a bivariate normal distribution
      if(i<500){
        q.var.lambda <- c(in.sd.lambda,0)
        q.var.omega <- c(0,in.sd.omega)
      }else{
        q.var.lambda <- c(var(lambda.mat[(i-201):(i-1),obs]),
                          cov(omega.mat[(i-201):(i-1),obs],lambda.mat[(i-201):(i-1),obs]))
        q.var.omega <- c(cov(omega.mat[(i-201):(i-1),obs],lambda.mat[(i-201):(i-1),obs]),
                      var(omega.mat[(i-201):(i-1),obs])) 
      }
      q.mean <- c(lambda.curr.un,omega.curr.un)
      q.var.matrix <- cbind(q.var.lambda,q.var.omega)
      #Proposal
      const <- (2.38^2)/2
      Q <- rmvnorm(1,q.mean, const * q.var.matrix)
      lambda.cand.un <- Q[1]
      omega.cand.un <- Q[2]
      if(omega.cand.un < 0 | omega.cand.un > log(20)){
        lambda.mat[i,obs] <- lambda.curr.un
        omega.mat[i,obs] <- omega.curr.un
        acceptance <- acceptance
      } else {
        #ACCEPTANCE STEP
        
        # Exchange algorithm for gamma and nu
        
        par.post[obs] <- a.exch.new(y[obs], lambda.curr.un, omega.curr.un,
                                      lambda.cand.un, omega.cand.un, expo[obs], mu.curr,
                                      sigma.curr, eta.curr, tau.curr)$test
        
        
        lambda.mat[i,obs] <- ifelse(par.post[obs], lambda.cand.un, lambda.curr.un)
        omega.mat[i,obs] <- ifelse(par.post[obs], omega.cand.un, omega.curr.un)
        acceptance <- ifelse(par.post[obs], acceptance+1, acceptance)
      }
      
    }
    
    
    
    # by default standard normal or if model.upd == "RW" it is a random walk step
    # or with hyperparameters if Gibbs
  #   if(model.upd=="RW"){
  #     mean.lambda.cand <- lambda.mat[(i-1),]
  #     mean.omega.cand <- omega.mat[(i-1),]
  #     sd.lambda.cand <- sd.lambda.cand
  #     sd.omega.cand <- sd.omega.cand
  #   } else if(model.upd=="Gibbs"){
  #     mean.lambda.cand <- mu.curr
  #     mean.omega.cand <- eta.curr
  #     sd.lambda.cand <- sigma.curr
  #     sd.omega.cand <- tau.curr
  #   } else{
  #     mean.lambda.cand <- mean.lambda.cand
  #     mean.omega.cand <- mean.omega.cand
  #     sd.lambda.cand <- sd.lambda.cand
  #     sd.omega.cand <- sd.omega.cand
  #   }
  #   
  #   
  #   #
  #   lambda.cand<- rnorm(n, mean.lambda.cand, sd.lambda.cand)
  #   omega.cand <- rnorm(n, mean.omega.cand, sd.omega.cand)
  #   
  #   
  #   # Exchange algorithm for lambda and omega
  #   par.post <- c()
  #   for(k in 1:n){
  #     par.post[k] <- a.exch(y[k], lambda.curr[k], omega.curr[k],
  #                           lambda.cand[k], omega.cand[k],
  #                           expo[k], mu.curr, sigma.curr,
  #                           eta.curr, tau.curr)$test
  #   }
  #   
  #   lambda.mat[i,] <- ifelse(par.post, lambda.cand, lambda.curr)
  #   omega.mat[i,] <- ifelse(par.post, omega.cand, omega.curr)
  #   acceptance <- ifelse(par.post, acceptance + 1, acceptance)
  #   
  # }
  # 
  
    }
  return(list (parameters = list(gamma = exp(lambda.mat), nu = exp(omega.mat),
                                 mu = mu.vec, sigma = sigma.vec,
                                 eta = eta.vec, tau = tau.vec),
               acceptance = acceptance/(iter*n))
  )
}
