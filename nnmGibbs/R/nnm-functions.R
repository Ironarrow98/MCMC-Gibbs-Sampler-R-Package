

#' Prior data loglikelihood for the normal-normal-mixture (NNM) model.
#'
#' @param K An integer represents the number of mixture components.
#' @param p An integer represents is the dimension of each observation.
#' @param vk An integer used in Inverse-Wishart distribution which set to `p + 2`.
#' @param omega An `p x p` symmetric positive-definite matrix used in Inverse-Wishart distribution.
#' @param Sigma A `p x p x K` matrix of group variances.
#'
#' @return The Prior data loglikelihood evaluated at the inputs (scalar).
#' @export
nnm_rest <- function(K, p, vk, omega, Sigma) {
  result <- 0
  for(kk in 1:K) {
    result <- result + mniw::diwish(Sigma[, , kk], omega, vk, log = TRUE)
  }
  return(result)
}




#' Complete data loglikelihood for the normal-normal-mixture (NNM) model.
#'
#' @param mu A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Sigma A `p x p x K` matrix of group variances.
#' @param rho A length-`K` probability vector of group membership.
#' @param y An `N x p` matrix of observations, each of which is a column.
#' @param V A `p x p x N` matrix of variances corresponding to each observation.
#' @param theta A `N x p` matrix of observation means.
#' @param z A length-`N` vector of integers between 1 and `K` indicating the group membership of each column of `theta`.
#'
#' @return The complete data loglikelihood evaluated at the inputs (scalar).
#' @export
nnm_loglik <- function(mu, Sigma, rho, y, V, theta, z) {
  K <- nrow(mu) # number of groups
  ll <- 0
  for(kk in 1:K) {
    idk <- z == kk # T/F membership in group k
    Nk <- sum(idk) # number of observations in group k
    if(Nk > 0) {
      # group membership contribution
      ll <- ll + Nk * log(rho[kk])
      # mean contribution
      ll <- ll + sum(dmNorm(theta[idk,], log = TRUE, mu = mu[kk,], Sigma = Sigma[,,kk]))
    }
  }
  # observation contribution
  result <- ll + sum(dmNorm(y, mu = theta, Sigma = V, log = TRUE))
  return(result)
}




#' Complete data loglikelihood with prior for the normal-normal-mixture (NNM) model.
#'
#' @param mu A `K x p` matrix of group means, where `K` is the number of mixture components and `p` is the dimension of each observation.
#' @param Sigma A `p x p x K` matrix of group variances.
#' @param rho A length-`K` probability vector of group membership.
#' @param y An `N x p` matrix of observations, each of which is a column.
#' @param V A `p x p x N` matrix of variances corresponding to each observation.
#' @param theta A `N x p` matrix of observation means.
#' @param z A length-`N` vector of integers between 1 and `K` indicating the group membership of each column of `theta`.
#'
#' @return The complete data loglikelihood with prior evaluated at the inputs (scalar).
#' @export
nnm_loglik_26 <- function(mu, Sigma, rho, y, V, theta, z) {
  K <- nrow(mu) # number of groups
  p <- ncol(mu)
  vk <- p + 2
  omega <- var(y)
  result <- nnm_loglik(mu, Sigma, rho, y, V, theta, z) - nnm_rest(K, p, vk, omega, Sigma)
  return(result)
}




#' Initial cluster allocation.
#'
#' Initializes the clusters using the kmeans++ algorithm of Arthur & Vassilvitskii (2007).
#'
#' @param y An `N x p` matrix of observations, each of which is a column.
#' @param K Number of clusters (positive integer).
#' @param max_iter The maximum number of steps in the [stats::kmeans()] algorithm.
#' @param nstart The number of random starts in the [stats::kmeans()] algorithm.
#'
#' @return A vector of length `N` consisting of integers between `1` and `K` specifying an initial clustering of the observations.
#'
#' @references Arthur, D., Vassilvitskii, S. "k-means++: the advantages of careful seeding" *Proceedings of the 18th Annual ACM-SIAM Symposium on Discrete Algorithms.* (2007): 1027–1035. <http://ilpubs.stanford.edu:8090/778/1/2006-13.pdf>.
#' @export
init_z <- function(y, K, max_iter = 10, nstart = 10) {
  # init++
  N <- nrow(y)
  p <- ncol(y)
  x <- t(y) # easier to use columns
  centers <- matrix(NA, p, K) # centers
  icenters <- rep(NA, K-1) # indices of centers
  minD <- rep(Inf, N) # current minimum distance vector
  # initialize
  icenters[1] <- sample(N,1)
  centers[,1] <- x[,icenters[1]]
  for(ii in 1:(K-1)) {
    D <- colSums((x - centers[,ii])^2) # new distance
    minD <- pmin(minD, D)
    icenters[ii+1] <- sample(N, 1, prob = minD)
    centers[,ii+1] <- x[,icenters[ii+1]]
  }
  centers <- t(centers)
  colnames(centers) <- rownames(x)
  # run kmeans with these centers
  km <- kmeans(x = y, centers = centers, nstart = nstart, iter.max = max_iter)
  km$cluster
}




#' Sample from a categorical distribution.
#'
#' Performs one draw from a categorical distribution (i.e., a multinomial distribution with size `n = 1`) for each of multiple probability vectors.
#'
#' @param prob An `n_cat x n_prob` matrix of probability vectors, each of which is a column of `prob`.  The entries of `prob` must be nonnegative, but are internally normalized such that each column sums to one.
#' @return A vector of length `n_prob` of integers between 1 and `n_cat` indicating the category draw for each column of `prob`.
#' @export
rcategorical <- function(prob) {
  if(any(prob < 0)) stop("prob must contain only nonnegative entries.")
  cprob <- apply(prob, 2, cumsum)
  u <- runif(ncol(prob)) * cprob[nrow(prob),]
  apply(sweep(cprob, 2, u) >= 0, 2, which.max)
}




#' Density of the Dirichlet distribution.
#'
#' @param x Observation vector of nonnegative entries which sum to one.
#' @param alpha Weight parameter: a vector of nonnegative entries the same length as `x`.
#' @param log Logical; whether to evaluate the density on the log scale.
#'
#' @return The density or log-density evaluated at the inputs (scalar).
#' @export
ddirichlet <- function(x, alpha, log = FALSE) {
  ld <- sum(lgamma(alpha)) - lgamma(sum(alpha))
  ld <- ld + sum((alpha-1) * log(x))
  if(!log) ld <- exp(ld)
  ld
}




#' Random sampling from the Dirichlet distribution.
#'
#' @param n Number of random draws.
#' @param alpha Weight parameter: a vector of nonnegative entries.
#' @return A matrix of size `n x length(alpha)` of which each row is a random draw.
#' @export
rdirichlet <- function(n, alpha) {
  K <- length(alpha) # number of categories
  X <- matrix(rgamma(n*K, shape = alpha), K, n)
  drop(t(sweep(X, 2, colSums(X), "/")))
}



