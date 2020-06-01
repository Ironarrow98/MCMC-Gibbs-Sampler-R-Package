#--- normal-normal-mixture (NNM) model -------------------------------------------
#
# this file contains tests of elementary functions to set up the MCMC
# sampler for the posterior distribution of
#
#

# Contents of file `test-nnm.R`

context("nnm-gibbs-sampler")

test_that("", {

  result <- array(NA, dim = c(6, 15))

  #--- conditional updates for theta_i ----------------------------

  # simulate data
  p <- 2
  N <- 50
  K <- 4

  rho <- runif(K)
  rho <- rho / sum(rho) # 1 * 4
  z <- replicate(n = N,expr = rcategorical(matrix(rho)))

  mu <- matrix(rnorm(K * p),K)
  Sigma <- replicate(K, crossprod(matrix(rnorm(p ^ 2), p, p)))
  theta <- rmNorm(n = N, mu = mu[z, ], Sigma = Sigma[, , z])

  V <- replicate(N, crossprod(matrix(rnorm(p ^ 2), p, p)))
  y <- rmNorm(n = N, mu = theta, Sigma = V)


  ntest <- 15

  Theta <- replicate(ntest, {
    list(theta = rmNorm(n = N, mu = mu[z, ], Sigma = Sigma[, , z]))
  }, simplify = FALSE)


  # 1. check that (simplified) conditional posterior = (unsimplified) loglik
  # assume flat prior

  # unsimplified
  llu <- sapply(1:ntest, function(ii) {
    theta1 <- Theta[[ii]]$theta
    nnm_loglik_26(mu = mu, Sigma = Sigma, rho = rho, y = y, V = V, theta = theta1, z = z)
  })

  # simplified
  lls <- sapply(1:ntest, function(ii) {
    theta1 <- Theta[[ii]]$theta
    sum(dmNorm(theta1, mu = mu[z, ], Sigma = Sigma[, , z], log = TRUE) +
          dmNorm(y, mu = theta1, Sigma = V, log = TRUE))
  })


  result[1, ] <- llu - lls




  #--- conditional updates for mu_k ----------------------------

  # simulate data
  p <- 2
  N <- 50
  K <- 4

  rho <- runif(K)
  rho <- rho / sum(rho) # 1 * 4
  z <- replicate(n = N, expr = rcategorical(matrix(rho)))

  mu <- matrix(rnorm(K * p), K)
  Sigma <- replicate(K, crossprod(matrix(rnorm(p ^ 2), p, p)))
  theta <- rmNorm(n = N, mu = mu[z, ], Sigma = Sigma[, , z])

  V <- replicate(N, crossprod(matrix(rnorm(p ^ 2), p, p)))
  y <- rmNorm(n = N, mu = theta, Sigma = V)


  ntest <- 15

  Mu <- replicate(ntest, {
    list(mu = matrix(rnorm(K * p), K))
  }, simplify = FALSE)


  # unsimplified
  llu <- sapply(1:ntest, function(ii) {
    mu1 <- Mu[[ii]]$mu
    nnm_loglik_26(mu = mu1, Sigma = Sigma, rho = rho, y = y, V = V, theta = theta, z = z)
  })


  # simplified
  lls <- sapply(1:ntest, function(ii) {
    mu1 <- Mu[[ii]]$mu
    N_k <- sapply(X = 1:K, FUN = function(kk){sum(z == kk)}, simplify = "array")
    sum(dmNorm(x = mu1,
               mu = t(sapply(X = 1:K, FUN = function(kk){colMeans(subset(theta, z == kk))}, simplify = "array")),
               Sigma = sweep(Sigma, 3, N_k, "/"),log = TRUE))
  })


  result[2, ] <- llu - lls




  #--- conditional updates for sigma_k ----------------------------

  # simulate data
  p <- 2
  N <- 50
  K <- 4

  rho <- runif(K)
  rho <- rho / sum(rho) # 1 * 4
  z <- replicate(n = N, expr = rcategorical(matrix(rho)))

  mu <- matrix(rnorm(K * p),K)
  Sigma <- replicate(K, crossprod(matrix(rnorm(p ^ 2), p, p)))
  theta <- rmNorm(n = N, mu = mu[z, ], Sigma = Sigma[, , z])

  V <- replicate(N, crossprod(matrix(rnorm(p ^ 2), p, p)))
  y <- rmNorm(n = N, mu = theta, Sigma = V)


  ntest <- 15

  Sigma <- replicate(ntest, {
    list(Sigma = replicate(K, crossprod(matrix(rnorm(p^2), p, p))))
  }, simplify = FALSE)

  # unsimplified
  llu <- sapply(1:ntest, function(ii) {
    Sigma <- Sigma[[ii]]$Sigma
    nnm_loglik(mu = mu, Sigma = Sigma, rho = rho, y = y, V = V, theta = theta, z = z)
  })


  # simplified
  lls <- sapply(1:ntest, function(ii) {
    Sigma <- Sigma[[ii]]$Sigma
    omega <- var(y)
    vk <- p + 2
    N_k <- sapply(X = 1:K, simplify = "array", FUN = function(i){sum(z==i)})
    Psi <- sapply(X = 1:K, simplify = "array", FUN = function(kk){crossprod(sweep(subset(theta,z==kk), 2, mu[kk,],'-'))})
    Psi <- sapply(X = 1:K, simplify = "array", FUN = function(kk){Psi[,,kk]+omega})
    det(Sigma[,,1]) + sum(diag(solve(Sigma[,,1]) %*% omega))
    prior <- sapply(X = 1:K, FUN = function(kk){(vk+p+1)*log(det(Sigma[,,kk]))+sum(diag(solve(Sigma[,,kk])%*%omega))})
    sum(diwish(X = Sigma, Psi = Psi, nu = N_k+vk, log = TRUE)) + (1/2) * sum(prior)
  })


  result[3, ] <- llu - lls




  #--- conditional updates for rho ----------------------------

  # simulate data
  p <- 2
  N <- 50
  K <- 4

  rho <- runif(K)
  rho <- rho / sum(rho) # 1 * 4
  z <- replicate(n = N,expr = rcategorical(matrix(rho)))

  mu <- matrix(rnorm(K * p),K)
  Sigma <- replicate(K, crossprod(matrix(rnorm(p ^ 2), p, p)))
  theta <- rmNorm(n = N, mu = mu[z, ], Sigma = Sigma[, , z])

  V <- replicate(N, crossprod(matrix(rnorm(p ^ 2), p, p)))
  y <- rmNorm(n = N, mu = theta, Sigma = V)


  ntest <- 15

  Rho <- replicate(ntest, {
    list(rho = rdirichlet(1,rep(N / K, K)))
  }, simplify = FALSE)


  # unsimplified
  llu <- sapply(1:ntest, function(ii) {
    rho1 <- Rho[[ii]]$rho
    nnm_loglik_26(mu = mu, Sigma = Sigma, rho = rho1, y = y, V = V, theta = theta, z = z)

  })


  # simplified
  lls <- sapply(1:ntest, function(ii) {
    rho1 <- Rho[[ii]]$rho
    N_k <- sapply(X = 1:K, simplify = "array", FUN = function(i){sum(z == i)})
    sum(ddirichlet(x = rho1, alpha = N_k + 1, log = TRUE))
  })

  lls1 <- sapply(1:ntest, function(ii) {
    rho1 <- Rho[[ii]]$rho
    sum(log(rho1[z]))
  })


  result[4, ] <- llu - lls
  result[5, ] <- llu - lls1




  #--- conditional updates for z ----------------------------

  # simulate data
  p <- 2
  N <- 50
  K <- 4

  rho <- runif(K)
  rho <- rho / sum(rho) # 1 * 4
  z <- replicate(n = N,expr = rcategorical(matrix(rho)))

  mu <- matrix(rnorm(K * p),K)
  Sigma <- replicate(K, crossprod(matrix(rnorm(p ^ 2), p, p)))
  theta <- rmNorm(n = N, mu = mu[z, ], Sigma = Sigma[, , z])

  V <- replicate(N, crossprod(matrix(rnorm(p ^ 2), p, p)))
  y <- rmNorm(n = N, mu = theta, Sigma = V)


  ntest <- 15

  Z <- replicate(ntest, {
    list(z = replicate(n = N,expr = rcategorical(matrix(rho))))
  }, simplify = FALSE)


  # unsimplified
  llu <- sapply(1:ntest, function(ii) {
    z1 <- Z[[ii]]$z
    nnm_loglik_26(mu = mu, Sigma = Sigma, rho = rho, y = y, V = V, theta = theta, z = z1)
  })


  # simplified
  lls <- sapply(1:ntest, function(ii) {
    z1 <- Z[[ii]]$z
    sum(log(rho[z1]) + dmNorm(x = theta, mu = mu[z1, ], Sigma = Sigma[, , z1],log = TRUE))
  })


  result[6, ] <- llu - lls

  for(ii in 1:(ntest - 1)) {
    expect_equal(result[, ii], result[, ii+1])
  }


})







