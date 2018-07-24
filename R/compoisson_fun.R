#-----------------------
# Metropolis-Hastings
#-----------------------

a.exch <- function(y, lambda.curr, omega.curr, lambda.cand, omega.cand, q, mu,
                   sigma, eta, tau, max = 10^4){
  n<- length(y)
  # CMP parameters
  #
  mu.curr <- exp(lambda.curr)*q
  nu.curr <- exp(omega.curr)
  gamma.curr <- mu.curr^nu.curr
  #
  mu.cand <- exp(lambda.cand)*q
  nu.cand <- exp(omega.cand)
  gamma.cand <- mu.cand^nu.cand
  #
  # Log-likelihood ratio
  
  lik.cand <- log(dcom(x = y, lambda = gamma.cand, nu = nu.cand))
  lik.curr <- log(dcom(x = y, lambda = gamma.curr, nu = nu.curr))
  
  llik.ratio <- sum(lik.cand) - sum(lik.curr)
  # Prior ratio
  prior.ratio <- sum(dnorm(lambda.cand, mu, sigma, log=TRUE),
                     - dnorm(lambda.curr, mu, sigma, log=TRUE),
                     dnorm(omega.cand, eta, tau, log=TRUE),
                     - dnorm(omega.curr, eta, tau, log=TRUE))
  # Posterior
  post <- llik.ratio + prior.ratio
  test <- post > log(runif(1))
  return(list(postdens = post, test = test))
}


a.exch1 <- function(y, lambda.curr, omega.curr, lambda.cand, omega.cand, q, mu,
                    sigma, eta, tau, max = 10^4){
  n<- length(y)
  # CMP parameters
  #
  mu.curr <- exp(lambda.curr)*q
  nu.curr <- exp(omega.curr)
  gamma.curr <- mu.curr^nu.curr
  
  #
  mu.cand <- exp(lambda.cand)*q
  nu.cand <- exp(omega.cand)
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
  prior.ratio <- sum(sum(dnorm(lambda.cand,mu,sigma, log = TRUE)),
                     -sum(dnorm(lambda.curr,mu,sigma, log = TRUE)),
                     dnorm(omega.cand, eta, tau, log = TRUE),
                     - dnorm(omega.curr, eta, tau, log = TRUE))
  # Posterior
  post <- llik.ratio + prior.ratio
  test <- post > log(runif(1))
  return(list(postdens = post, test = test))
}


a.exch.comp <- function(y, lambda.curr, nu.curr, lambda.cand, nu.cand, q, mu,
                         sigma, ommin=0.001, ommax=10, max = 10^4){
  
  n<- length(y)
  # CMP parameters
  #
  mu.curr <- exp(lambda.curr)*q
  gamma.curr <- mu.curr^nu.curr
  
  #
  mu.cand <- exp(lambda.cand)*q
  gamma.cand <- mu.cand^nu.cand
  #
  # Log-likelihood ratio
  
    lik.cand <- log(dcom(x = y, lambda = gamma.cand, nu = nu.cand))
    lik.curr <- log(dcom(x = y, lambda = gamma.curr, nu = nu.curr))
  
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


a.exch1.comp <- function(y, lambda.curr, nu.curr, lambda.cand, nu.cand, q, mu,
                         sigma, ommin, ommax, nu.prior){
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



#------------------------------------
# MCMC
#-----------------------------------

mcmc.dem.simple.comp <- function(y, iter, expo, gamma0 = rep(1,length(y)), nu0 = rep(1,length(y)),
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
      if(i < 500){
        q.var.lgamma <- c(in.sd.lgamma,0)
        q.var.nu <- c(0,in.sd.nu)
       }else{
         q.var.lgamma <- c(var(gamma.mat[(i-201):(i-1),obs]),
                           cov(nu.mat[(i-201):(i-1),obs],gamma.mat[(i-201):(i-1),obs]))
         q.var.nu <- c(cov(nu.mat[(i-201):(i-1),obs],gamma.mat[(i-201):(i-1),obs]),
                       var(nu.mat[(i-201):(i-1),obs])) 
       }
      q.mean <- c(lgamma.curr,nu.curr)
      q.var.matrix <- cbind(q.var.lgamma,q.var.nu)
      #Proposal
      
      # const <- (2.38^2)/2
      Q <- rmvnorm(1, q.mean, q.var.matrix)
      lgamma.cand <- Q[1]
      nu.cand <- Q[2]
      if(nu.cand < 0 | nu.cand > 40){
        gamma.mat[i,obs] <- exp(lgamma.curr)
        nu.mat[i,obs] <- nu.curr
        gamma.nu.acc <- gamma.nu.acc
      } else {
        #ACCEPTANCE STEP
        
        # Exchange algorithm for gamma and nu
        
        par.post[obs] <- a.exch.comp(y[obs], lgamma.curr, nu.curr,
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




mcmc.dem.simple.comp1 <- function(y, iter, expo, gamma0 = 1, nu0 = 1,
                                 lgamma0 = log(gamma0),
                                 mu = 0, sigma = 10,
                                 ommin = 0.001, ommax = 10,
                                 in.sd.lgamma, in.sd.nu){
  library(mvtnorm)
  library(compoisson)
  n <-length(y)
  gamma.mat <- matrix(NA, iter, 1)
  nu.mat <- matrix(NA, iter, 1)
  VarCov <- array(NA, dim = c(2,2,iter))
  
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
    VarCov[,,i] <- q.var.matrix
    #Proposal
    const <- (2.38^2)/2 #al denominatore c'è la dimensione della distribuzione multivariata
    Q <- rmvnorm(1,q.mean, const*q.var.matrix)
    lgamma.cand <- Q[1]
    nu.cand <- Q[2]
    
    if(nu.cand < 0 | nu.cand > 10){
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
               acceptance = list(gamma.nu = gamma.nu.acc/iter),
               varcov = VarCov))
}


mcmc.dem.IG.comp1 <- function(y, iter, expo, lambda0 = 0, nu0 = 0,
                              mu0 = 0, sigma0 = 0.1,
                              alpha0 = 0, beta0 = 10,
                              sh.sigma0 = 0.001, 
                              sc.sigma0 = 0.001, 
                              in.sd.lambda = 1, in.sd.nu = 1,
                              ommin = 0.0001, ommax = 20){
  
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
      q.var.nu <- c(cov(nu.vec[(i-201):(i-1)],lambda.vec[(i-201):(i-1)]),var(nu.vec[(i-201):(i-1)])) 
    }
    q.mean <- c(lambda.curr,nu.curr)
    q.var.matrix <- cbind(q.var.lambda,q.var.nu)
    
    # if(!is.finite(q.var.matrix)){
    #   print("problems VarCov matrix")
    # }
    
    #q.var.matrix <-ifelse(is.finite(q.var.matrix), q.var.matrix, matrix(c(in.sd.lambda,0,0,in.sd.nu), 2,2))

    #Proposal
    const <- (2.38^2)/2
    Q <- rmvnorm(1,q.mean, const * q.var.matrix)
    lambda.cand <- Q[1]
    nu.cand <- Q[2]
    if(nu.cand < 0 | nu.cand > 20){
      lambda.vec[i] <- lambda.curr
      nu.vec[i] <- nu.curr
      acceptance <- acceptance 
    } else {
      #ACCEPTANCE STEP
      
      # Exchange algorithm for gamma and nu
      
      par.post <- a.exch1(y, lambda.curr, nu.curr,
                                lambda.cand, nu.cand, expo, mu,
                                sigma, ommin, ommax)$test
      
      
      lambda.vec[i] <- ifelse(par.post, lambda.cand, lambda.curr)
nu.vec[i] <- ifelse(par.post, nu.cand, nu.curr)
acceptance <- ifelse(par.post, acceptance+1, acceptance)
    }
    
  }
  
  return(list (parameters = list(gamma = exp(lambda.vec), nu = nu.vec,
                                 mu = mu.vec, sigma = sigma.vec),
               acceptance = list(acceptance = acceptance/iter))
  )
}

