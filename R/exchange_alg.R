#---------------------
# Exchange algorithms
#---------------------

# For model type: y_i = CMP(mu_i, nu_i)
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
  test <- post > log(runif(1))
  return(list(postdens = post, test = test))
}

# For model type: y_i = CMP(mu_i, nu)
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
  prior.ratio <- sum(sum(dnorm(lambda.cand,mu,sigma, log = TRUE)),
                     -sum(dnorm(lambda.curr,mu,sigma, log = TRUE)),
                     dnorm(omega.cand, eta, tau, log = TRUE),
                     - dnorm(omega.curr, eta, tau, log = TRUE))
  # Posterior
  post <- llik.ratio + prior.ratio
  test <- post > log(runif(1))
  return(list(postdens = post, test = test))
}

#-------------------------------------------
# For model type: y_i = CMP(mu_i, nu)
# The function applies to all the data y
#-------------------------------------------

a.exch1.bis <- function(y, lambda.curr, omega.curr, lambda.cand, omega.cand, q, mu,
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
  y.new <- sapply(mu.cand, rcmp1, nu = nu.cand, max = max)
  #
  # Log-likelihood ratio
  dens1 <- mapply(logDensCMP1, y = y, gamma = mu.cand ,nu = nu.cand)
  dens2 <- mapply(logDensCMP1, y = y, gamma = mu.curr ,nu = nu.curr)
  dens3 <- mapply(logDensCMP1, y = y.new, gamma = mu.curr ,nu = nu.curr)
  dens4 <-mapply(logDensCMP1, y = y.new, gamma = mu.cand ,nu = nu.cand)

  llik.ratio <- sum(cbind(dens1, - dens2, dens3, -dens4))
  # Prior ratio
  prior.ratio <- sum(dnorm(omega.cand, eta, tau, log = TRUE),
                     - dnorm(omega.curr, eta, tau, log = TRUE))
  # Posterior
  post <- llik.ratio + prior.ratio
  test <- post > log(runif(1))
  return(list(postdens = post, test = test))
}

#----------------------------------------
# Exchange alg uninformative prior on nu
#----------------------------------------

a.exch.uninf.nu <- function(y, lambda.curr, omega.curr, lambda.cand, omega.cand, q, mu,
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
  prior.ratio <- 0
  
  # Posterior
  post <- llik.ratio + prior.ratio
  test <- post > log(runif(1))
  return(list(postdens = post, test = test))
}

#--------------------------------------------
# For model type: y_i = CMP(mu_i, nu)
# The function applies to all the data y
a.exch1.bis <- function(y, lambda.curr, omega.curr, lambda.cand, omega.cand, q, mu,
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
  y.new <- sapply(mu.cand, rcmp1, nu = nu.cand, max = max)
  #
  # Log-likelihood ratio
  dens1 <- mapply(logDensCMP1, y = y, gamma = mu.cand ,nu = nu.cand)
  dens2 <- mapply(logDensCMP1, y = y, gamma = mu.curr ,nu = nu.curr)
  dens3 <- mapply(logDensCMP1, y = y.new, gamma = mu.curr ,nu = nu.curr)
  dens4 <-mapply(logDensCMP1, y = y.new, gamma = mu.cand ,nu = nu.cand)
  
  llik.ratio <- sum(cbind(dens1, - dens2, dens3, -dens4))
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
