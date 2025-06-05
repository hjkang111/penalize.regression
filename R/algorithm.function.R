
#' Coordinate Descent Algorithm (CDA) for Penalized Regression
#'
#' This function implements the Coordinate Descent Algorithm (CDA) for solving
#' penalized regression problems, supporting LASSO, Ridge, Elastic Net, SCAD,
#' and MCP penalties.
#'
#' @param X A numeric matrix of predictors (n x p).
#' @param y A numeric response vector of length n.
#' @param method A character string specifying the penalty method: "lasso",
#' "ridge", "elasticnet", "scad", or "mcp".
#' @param lambda A non-negative regularization parameter.
#' @param learning_rate A learning rate parameter (currently not used in CDA,
#' but reserved for consistency).
#' @param max_iter An integer specifying the maximum number of iterations.
#' @param alpha Elastic Net mixing parameter (0 <= alpha <= 1). Used only when
#' \code{method = "elasticnet"}.
#' @param gamma A tuning parameter for SCAD and MCP (typically >= 3).
#'
#' @return A numeric vector of estimated regression coefficients of length p.
#' @export
perform_CDA <- function(X, y, method, lambda, learning_rate = 0.01, max_iter = 1000, alpha = 0.5, gamma = 3.7) {
  n <- nrow(X)
  p <- ncol(X)
  beta <- rep(0, p)
  Xy <- t(X) %*% y
  XX <- colSums(X^2)

  soft_threshold <- function(z, t) {
    sign(z) * pmax(0, abs(z) - t)
  }

  for (iter in 1:max_iter) {
    for (j in 1:p) {
      r_j <- y - X %*% beta + X[, j] * beta[j]
      rho_j <- sum(X[, j] * r_j)

      if (method == "lasso") {
        beta[j] <- soft_threshold(rho_j / XX[j], lambda / XX[j])

      } else if (method == "ridge") {
        beta[j] <- rho_j / (XX[j] + 2 * lambda)

      } else if (method == "elasticnet") {
        z <- rho_j / XX[j]
        beta[j] <- soft_threshold(z, lambda * alpha / XX[j]) / (1 + lambda * (1 - alpha) / XX[j])

      } else if (method == "scad") {
        z <- rho_j / XX[j]
        if (abs(z) <= lambda) {
          beta[j] <- soft_threshold(z, lambda / XX[j])
        } else if (abs(z) <= gamma * lambda) {
          beta[j] <- soft_threshold(z, gamma * lambda / (gamma - 1) / XX[j])
        } else {
          beta[j] <- z
        }

      } else if (method == "mcp") {
        z <- rho_j / XX[j]
        abj <- abs(z)

        if (abj <= gamma * lambda) {
          beta[j] <- soft_threshold(z, lambda / XX[j]) / (1 - 1 / gamma)
        } else {
          beta[j] <- z
        }

      } else {
        stop("Unsupported method in CDA.")
      }
    }
  }

  return(beta)
}

#' Fast Iterative Shrinkage-Thresholding Algorithm (FISTA) for Penalized Regression
#'
#' This function implements the FISTA algorithm for solving penalized regression problems, supporting LASSO, Ridge, Elastic Net, SCAD, and MCP penalties.
#'
#' @param X A numeric matrix of predictors (n x p).
#' @param y A numeric response vector of length n.
#' @param method A character string specifying the penalty method: "lasso", "ridge", "elasticnet", "scad", or "mcp".
#' @param lambda A non-negative regularization parameter.
#' @param learning_rate A numeric learning rate controlling the gradient step size (default: 1e-3).
#' @param max_iter An integer specifying the maximum number of iterations.
#' @param alpha Elastic Net mixing parameter (0 <= alpha <= 1). Used only when \code{method = "elasticnet"}.
#' @param gamma A tuning parameter for SCAD and MCP (typically >= 3).
#'
#' @return A numeric vector of estimated regression coefficients of length p.
#' @export
perform_FISTA <- function(X, y, method, lambda, learning_rate = 1e-3, max_iter = 1000, alpha = 0.5, gamma = 3.7) {
  n <- nrow(X)
  p <- ncol(X)

  beta <- rep(0, p)
  beta_old <- beta
  t <- 1

  soft_threshold <- function(z, t) {
    sign(z) * pmax(0, abs(z) - t)
  }

  grad <- function(beta) {
    -t(X) %*% (y - X %*% beta) / n
  }

  penalty_grad <- function(beta_j) {
    if (method == "lasso") {
      return(lambda * sign(beta_j))
    } else if (method == "ridge") {
      return(2 * lambda * beta_j)
    } else if (method == "elasticnet") {
      return(lambda * (alpha * sign(beta_j) + 2 * (1 - alpha) * beta_j))
    } else if (method == "scad") {
      abj <- abs(beta_j)
      if (abj <= lambda) {
        return(lambda * sign(beta_j))
      } else if (abj <= gamma * lambda) {
        return(((gamma * lambda - abj) / (gamma - 1)) * sign(beta_j))
      } else {
        return(0)
      }
    } else if (method == "mcp") {
      abj <- abs(beta_j)
      if (abj <= gamma * lambda) {
        return(lambda * (1 - abj / (gamma * lambda)) * sign(beta_j))
      } else {
        return(0)
      }
    } else {
      stop("Unsupported method.")
    }
  }

  for (k in 1:max_iter) {
    z <- beta + ((t - 1) / (t + 2)) * (beta - beta_old)
    grad_z <- grad(z)

    beta_new <- numeric(p)

    for (j in 1:p) {
      if (method == "ridge") {
        beta_new[j] <- z[j] - learning_rate * (grad_z[j] + 2 * lambda * z[j])

      } else if (method %in% c("lasso", "elasticnet")) {
        thresh <- learning_rate * lambda * alpha

        beta_new[j] <- soft_threshold(z[j] - learning_rate * grad_z[j], thresh)

      } else if (method %in% c("scad", "mcp")) {
        pen_grad_val <- penalty_grad(z[j])
        thresh <- learning_rate * abs(pen_grad_val)

        if (thresh < 1e-12) thresh <- 0

        beta_new[j] <- soft_threshold(z[j] - learning_rate * grad_z[j], thresh)

      } else {
        stop("Unsupported method in FISTA.")
      }
    }

    beta_old <- beta
    beta <- beta_new
    t <- t + 1

    if (sqrt(sum((beta - beta_old)^2)) < 1e-6) break
  }

  return(beta)
}

#' Local Linear Approximation (LLA) Algorithm for Penalized Regression
#'
#' This function implements the Local Linear Approximation (LLA) algorithm for solving penalized regression problems, supporting LASSO, Ridge, Elastic Net, SCAD, and MCP penalties.
#'
#' @param X A numeric matrix of predictors (n x p).
#' @param y A numeric response vector of length n.
#' @param method A character string specifying the penalty method: "lasso", "ridge", "elasticnet", "scad", or "mcp".
#' @param lambda A non-negative regularization parameter.
#' @param learning_rate A numeric learning rate controlling the update step size (default: 0.01).
#' @param max_iter An integer specifying the maximum number of iterations.
#' @param alpha Elastic Net mixing parameter (0 <= alpha <= 1). Used only when \code{method = "elasticnet"}.
#' @param gamma A tuning parameter for SCAD and MCP (typically >= 3).
#'
#' @return A numeric vector of estimated regression coefficients of length p.
#' @export
perform_LLA <- function(X, y, method, lambda, learning_rate = 0.01, max_iter = 100, alpha = 0.5, gamma = 3.7) {
  n <- nrow(X)
  p <- ncol(X)
  beta <- rep(0, p)
  tol <- 1e-4

  for (iter in 1:max_iter) {
    weights <- rep(1, p)

    for (j in 1:p) {
      bj <- beta[j]

      if (method == "lasso") {
        weights[j] <- 1
      } else if (method == "ridge") {
        weights[j] <- 2 * abs(bj)
      } else if (method == "elasticnet") {
        weights[j] <- alpha + 2 * (1 - alpha) * abs(bj)
      } else if (method == "scad") {
        abj <- abs(bj)
        if (abj <= lambda) {
          weights[j] <- 1
        } else if (abj <= gamma * lambda) {
          weights[j] <- (gamma * lambda - abj) / ((gamma - 1) * lambda)
        } else {
          weights[j] <- 0
        }
      } else if (method == "mcp") {
        abj <- abs(bj)
        if (abj <= gamma * lambda) {
          weights[j] <- 1 - abj / (gamma * lambda)
        } else {
          weights[j] <- 0
        }
      } else {
        stop("Unsupported method in LLA.")
      }
    }

    for (j in 1:p) {
      r_j <- y - X %*% beta + X[, j] * beta[j]
      rho_j <- sum(X[, j] * r_j)
      XX_j <- sum(X[, j]^2)
      beta[j] <- sign(rho_j) * max(0, abs(rho_j) - lambda * weights[j]) / XX_j
    }
  }

  return(beta)
}
