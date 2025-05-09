# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   https://r-pkgs.org
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

hello <- function() {
  print("Hello, world!")
}
# penalize regression
# example
set.seed(123)
n <- 100
p <- 10
X <- matrix(rnorm(n * p), n, p)
beta_true <- c(3, -2, rep(0, p - 2))
y <- X %*% beta_true + rnorm(n)

lambda <- 0.1

# lasso + GD(subgradient)
lasso_gd <- function(X, y, lambda, step_size = 0.001, max_iter = 1000) {
  n <- nrow(X)
  p <- ncol(X)
  beta <- rep(0, p)

  for (iter in 1:max_iter) {
    grad <- -t(X) %*% (y - X %*% beta) / n
    beta <- beta - step_size * (grad + lambda * sign(beta))  # Subgradient
  }
  return(beta)
}

# lasso + CDA
soft_threshold <- function(z, gamma) {
  sign(z) * pmax(abs(z) - gamma, 0)
}

lasso_cda <- function(X, y, lambda, max_iter = 1000, tol = 1e-6) {
  n <- nrow(X)
  p <- ncol(X)
  beta <- rep(0, p)

  for (iter in 1:max_iter) {
    beta_old <- beta
    for (j in 1:p) {
      r_j <- y - X[, -j] %*% beta[-j]
      zj <- sum(X[, j]^2)
      pj <- sum(X[, j] * r_j)
      beta[j] <- soft_threshold(pj / zj, lambda / zj)
    }
    if (sum(abs(beta - beta_old)) < tol) break
  }
  return(beta)
}


# lasso + PGD
lasso_pgd <- function(X, y, lambda, step_size = 0.01, max_iter = 1000) {
  n <- nrow(X)
  p <- ncol(X)
  beta <- rep(0, p)

  for (iter in 1:max_iter) {
    grad <- -t(X) %*% (y - X %*% beta) / n
    beta <- soft_threshold(beta - step_size * grad, lambda * step_size)
  }
  return(beta)
}
