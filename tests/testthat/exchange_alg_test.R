# Test exchange algorithm functions

n.test <- 5
test.identity <- FALSE
test.extended <- TRUE

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


test_that("Posterior density is right", {
  a.exch <- demCMP:::a.exch
  for (seed in seq_len(n.test)) {
    set.seed(seed)
    y <- rpois(n = 1, lambda = 10)
    lambda.curr <- runif(n = 1)
    lambda.cand <- runif(n = 1)
    omega.curr <- runif(n = 1)
    omega.curr <- runif(n = 1)
    q <- rpois(n = 1, lambda = 10)
    mu <- 0
    sigma <- 1
    eta <- 0
    tau <- 1
    max = 10^4
    ans.obtained <- a.exch(y, lambda.curr, omega.curr, lambda.cand, omega.cand, q, mu,
                     sigma, eta, tau, max = 10^4)$postdens
    ans.expected <- mydcmp(y = y, gamma = gamma, nu = nu, log = TRUE)
    if (test.identity)
      expect_identical(ans.obtained, ans.expected)
    else
      expect_equal(ans.obtained, ans.expected)
  }
})
