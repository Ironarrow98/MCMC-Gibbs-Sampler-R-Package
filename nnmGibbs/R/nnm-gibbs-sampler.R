
#--- conditional distrubution for theta_i ----------------------------

#' Conditional Distrubution of theta given other parameters
#'
#' @param theta A `N x p` matrix of observation means.
#' @param mu A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Sigma A `p x p x K` matrix of group variances.
#' @param z A length-`N` vector of integers between 1 and `K` indicating the group membership of each column of `theta`.
#' @param V A `p x p x N` matrix of variances corresponding to each observation.
#' @param y An `N x p` matrix of observations, each of which is a column.
#'
#' @return The conditional log-PDF of theta given other parameters
#' @export
theta.loglik <- function(theta, mu, Sigma, z, V, y) {
  result <- sum(dmNorm(theta, mu = mu[z, ], Sigma = Sigma[, , z], log = TRUE) + dmNorm(y, mu = theta, Sigma = V, log = TRUE))
  return(result)
}




#--- conditional updates for theta_i ----------------------------

#' Conditional Draw of theta given other parameters
#'
#' @param mu_curr A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Sigma_curr A `p x p x K` matrix of group variances.
#' @param Z_curr A length-`N` vector of integers between 1 and `K` indicating the group membership of each column of `theta`.
#' @param Y An `N x p` matrix of observations, each of which is a column.
#' @param V A `p x p x N` matrix of variances corresponding to each observation.
#'
#' @return The conditional draw of theta given other parameters
#' @export
update_theta <- function(mu_curr, Sigma_curr, Z_curr, Y, V) {
  N <- nrow(Y)
  result <- rRxNorm(n = N, x = Y, V = V, lambda = mu_curr[Z_curr, ], Sigma = Sigma_curr[, , Z_curr])
  return(result)
}

# update_theta(mu_curr = mu,Sigma_curr = sigma, Z_curr = z, Y = y, V = V)




#--- conditional distrubution for mu_k ----------------------------

#' Conditional Distrubution of mu given other parameters
#'
#' @param theta A `N x p` matrix of observation means.
#' @param mu A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Sigma A `p x p x K` matrix of group variances.
#' @param z A length-`N` vector of integers between 1 and `K` indicating the group membership of each column of `theta`.
#'
#' @return The conditional log-PDF of mu given other parameters
#' @export
mu.loglik <- function(theta, mu, Sigma, z) {
  K <- nrow(mu)
  N_k <- sapply(X = 1:K, FUN = function(kk){sum(z == kk)}, simplify = "array")
  result <- sum(dmNorm(x = mu,
                       mu = t(sapply(X = 1:K, FUN = function(kk){colMeans(subset(theta, z == kk))}, simplify = "array")),
                       Sigma = sweep(Sigma, 3, N_k, "/"),log = TRUE))
  return(result)
}




#--- conditional updates for mu_k ----------------------------

#' Conditional Draw of mu given other parameters
#'
#' @param theta_curr A `N x p` matrix of observation means.
#' @param Sigma_curr A `p x p x K` matrix of group variances.
#' @param Z_curr A length-`N` vector of integers between 1 and `K` indicating the group membership of each column of `theta`.
#' @param rho_curr A length-`K` probability vector of group membership.
#' @param Y An `N x p` matrix of observations, each of which is a column.
#'
#' @return The conditional draw of mu given other parameters
#' @export
update_mu <- function(theta_curr, Sigma_curr, Z_curr, rho_curr, Y) {
  K <- length(rho_curr) # number of groups
  N <- nrow(Y)
  N_k <- sapply(X = 1:K, FUN = function(kk){sum(Z_curr==kk)},simplify = "array")

  mu = t(sapply(X = 1:K, FUN = function(kk){colMeans(subset(theta_curr,Z_curr==kk))},simplify = "array")) # ????

  result <- rmNorm(n = K, mu = mu, Sigma = sweep(Sigma_curr, 3, N_k, "/")) # ??? colmean or rowmean?
  return(result)
}

# update_mu(theta_curr=theta, Sigma_curr=sigma, Z_curr=z, rho_curr=rho, Y=y)




#--- conditional distrubution for sigma_k ----------------------------

#' Conditional Distrubution of sigma given other parameters
#'
#' @param theta A `N x p` matrix of observation means.
#' @param mu A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Sigma A `p x p x K` matrix of group variances.
#' @param z A length-`N` vector of integers between 1 and `K` indicating the group membership of each column of `theta`.
#' @param y An `N x p` matrix of observations, each of which is a column.
#'
#' @return The conditional log-PDF of sigma given other parameters
#' @export
Sigma.loglik <- function(theta, mu, Sigma, z, y) {
  K <- nrow(mu) # number of groups
  p <- ncol(mu) #?????
  N <- length(y)
  omega <- var(y)
  vk <- p + 2
  N_k <- sapply(X = 1:K, simplify = "array", FUN = function(i){sum(z==i)})
  Psi <- sapply(X = 1:K, simplify = "array", FUN = function(kk){crossprod(sweep(subset(theta,z==kk), 2, mu[kk,],'-'))})
  Psi <- sapply(X = 1:K, simplify = "array", FUN = function(kk){Psi[,,kk]+omega})
  det(Sigma[,,1]) + sum(diag(solve(Sigma[,,1]) %*% omega))
  prior <- sapply(X = 1:K, FUN = function(kk){(vk+p+1)*log(det(Sigma[,,kk]))+sum(diag(solve(Sigma[,,kk])%*%omega))})
  result <- sum(diwish(X = Sigma, Psi = Psi, nu = N_k+vk, log = TRUE)) + (1/2) * sum(prior)
  return(result)
}




#--- conditional updates for sigma_k ----------------------------

#' Conditional Draw of sigma given other parameters
#'
#' @param theta_curr A `N x p` matrix of observation means.
#' @param mu_curr A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Z_curr A length-`N` vector of integers between 1 and `K` indicating the group membership of each column of `theta`.
#' @param Y An `N x p` matrix of observations, each of which is a column.
#'
#' @return The conditional draw of sigma given other parameters
#' @export
update_Sigma <- function(theta_curr, mu_curr, Z_curr, Y) {
  K <- nrow(mu_curr) # number of groups
  p <- ncol(mu_curr) #?????
  N <- length(Y)

  omega <- var(Y)
  vk <- p + 2

  N_k <- sapply(X = 1:K, simplify = "array",
                FUN = function(i){sum(Z_curr==i)})

  Psi <- sapply(X = 1:K, simplify = "array",
                FUN = function(kk){crossprod(sweep(subset(theta_curr,Z_curr==kk), 2, mu_curr[kk,],'-'))})

  Psi <- sapply(X = 1:K, simplify = "array",
                FUN = function(kk){Psi[,,kk]+omega})

  result <- riwish(n = K, Psi = Psi, nu = N_k+vk)
  return(result)
}

# update_Sigma(theta_curr=theta, mu_curr=mu, Z_curr=z, Y=y)




#--- conditional distrubution for rho_k ----------------------------

#' Conditional Distrubution of rho given other parameters
#'
#' @param rho A length-`K` probability vector of group membership.
#' @param mu A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param z A length-`N` vector of integers between 1 and `K` indicating the group membership of each column of `theta`.
#'
#' @return The conditional log-PDF of rho given other parameters
#' @export
rho.loglik <- function(rho, mu, z) {
  K <- nrow(mu) # number of groups
  N_k <- sapply(X = 1:K, simplify = "array", FUN = function(i){sum(z == i)})
  result <- sum(ddirichlet(x = rho, alpha = N_k + 1, log = TRUE))
  return(result)
}




#--- conditional updates for rho ----------------------------

#' Conditional Draw of rho given other parameters
#'
#' @param mu_curr A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Z_curr A length-`N` vector of integers between 1 and `K` indicating the group membership of each column of `theta`.
#'
#' @return The conditional draw of rho given other parameters
#' @export
update_rho <- function(mu_curr, Z_curr) {
  K <- nrow(mu_curr) # number of groups
  N <- nrow(Z_curr)
  N_k <- sapply(X = 1:K, simplify = "array", FUN = function(i){sum(Z_curr==i)})
  result <- rdirichlet(1, N_k+1)
  return(result)
}

# update_rho(mu_curr = mu, Z_curr = z)




#--- conditional distrubution for z ----------------------------

#' Conditional Distrubution of z given other parameters
#'
#' @param theta A `N x p` matrix of observation means.
#' @param mu A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Sigma A `p x p x K` matrix of group variances.
#' @param z A length-`N` vector of integers between 1 and `K` indicating the group membership of each column of `theta`.
#' @param rho A length-`K` probability vector of group membership.
#'
#' @return The conditional log-PDF of z given other parameters
#' @export
z.loglik <- function(theta, mu, Sigma, z, rho) {
  result <- sum(log(rho[z]) + dmNorm(x = theta, mu = mu[z, ], Sigma = Sigma[, , z], log = TRUE))
  return(result)
}




#--- conditional updates for z ----------------------------

#' Conditional Draw of z given other parameters
#'
#' @param theta_curr A `N x p` matrix of observation means.
#' @param mu_curr A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Sigma_curr A `p x p x K` matrix of group variances.
#' @param rho_curr A length-`K` probability vector of group membership.
#' @param Y An `N x p` matrix of observations, each of which is a column.
#'
#' @return The conditional draw of z and lambda given other parameters
#' @export
update_Z <- function(theta_curr, mu_curr, Sigma_curr, rho_curr, Y) {
  K <- nrow(mu_curr) # number of groups
  N <- nrow(Y)
  k_ik <- sapply(X=1:K, simplify = 'array', FUN = function(kk){dmNorm(x = theta_curr, mu = mu_curr[kk,],Sigma = Sigma_curr[,,kk],log = TRUE)})
  k_ik <- sweep(k_ik,2,log(rho_curr),'+')
  lambda <- exp(k_ik)

  z <- sapply(X = 1:N,simplify = 'array', FUN = function(ii){rcategorical(matrix(lambda[ii,]))})
  N_k <- sapply(X = 1:K, simplify = "array",
                FUN = function(i){sum(z==i)})
  while (any(N_k==0)) {
    z <- sapply(X = 1:N,simplify = 'array', FUN = function(ii){rcategorical(matrix(lambda[ii,]))})
    N_k <- sapply(X = 1:K, simplify = "array",
                  FUN = function(i){sum(z==i)})
  }

  result <- list(Z = z, lambda = lambda)
  return(result)
}

# update_Z(theta_curr = theta, mu_curr = mu, Sigma_curr = sigma, rho_curr = rho, Y = y)




#--- gibbs sampler for all parameters ---------------------------------------

#' Formal Gibbs Sampler for all five parameters of normal-normal-mixture (NNM) model for clustering
#'
#' @param nsamples An integer represents the number of samples the Gibbs Sampler will perform
#' @param burn An integer represents the burn-in sampler of the Gibbs Sampler
#' @param K An integer represents the number of mixture components.
#' @param Y An `N x p` matrix of observations, each of which is a column.
#' @param V A `p x p x N` matrix of variances corresponding to each observation.
#' @param theta_init A `N x p` matrix of observation means.
#' @param mu_init A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Sigma_init A `p x p x K` matrix of group variances.
#' @param Z_init A length-`N` vector of integers between 1 and `K` indicating the group membership of each column of `theta`.
#' @param rho_init A length-`K` probability vector of group membership.
#' @param return_Z An True/False parameter indicates whether to store z during every sample process or not
#' @param return_theta An True/False parameter indicates whether to store theta during every sample process or not
#'
#' @return A list contains all updated parameters
#' @export
gibbs_sampler <- function(nsamples, burn, K, Y, V, theta_init, mu_init, Sigma_init, Z_init, rho_init, return_Z = FALSE, return_theta = FALSE) {
  K <- nrow(mu_init)
  p <- ncol(mu_init)
  N <- length(Y)

  # allocate memory for output
  mu_out <- array(NA, dim = c(K, p, nsamples))
  Sigma_out <- array(NA, dim = c(p, p, K, nsamples))
  rho_out <- matrix(NA, K, nsamples)
  lambda_out <- 0
  if(return_Z) {
    # may not want to do this at each step (or at all) for large problems,
    # since this is an `N x nsamples` matrix which might take up a ton of memory
    Z_out <- matrix(NA, N, nsamples)
  }
  if(return_theta) {
    theta_out <- array(NA, dim = c(N, p, nsamples))
  }

  # initialize sampler
  theta_curr <- theta_init
  mu_curr <- mu_init
  Sigma_curr <- Sigma_init
  Z_curr <- Z_init
  rho_curr <- rho_init

  # gibbs begin
  for(ii in -burn:nsamples) {
    # updates
    theta_curr <- update_theta(mu_curr, Sigma_curr, Z_curr, Y, V)
    mu_curr <- update_mu(theta_curr, Sigma_curr, Z_curr, rho_curr, Y)
    Sigma_curr <- update_Sigma(theta_curr, mu_curr, Z_curr, Y)
    Z_list <- update_Z(theta_curr, mu_curr, Sigma_curr, rho_curr, Y)
    Z_curr <- Z_list$Z
    rho_curr <- update_rho(mu_curr, Z_curr)
    # store
    if(ii > 0) {
      # storage
      mu_out[,,ii] <- mu_curr
      Sigma_out[,,,ii] <- Sigma_curr
      rho_out[,ii] <- rho_curr
      lambda_out <- lambda_out + Z_list$lambda
      if(return_Z) {
        Z_out[,ii] <- Z_curr
      }
      if(return_theta) {
        theta_out[,,ii] <- theta_curr
      }
    }
  }

  result <- list(rho = rho_out, mu = mu_out, Sigma = Sigma_out, theta = NA, z = NA, lamdba = lambda_out)
  if(return_Z) {
    result$z <- Z_out
  }
  if(return_theta) {
    result$theta <- theta_out
  }
  return(result)
}

