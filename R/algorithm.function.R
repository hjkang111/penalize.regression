'''# algorithm functions

# soft threshold function
soft_threshold <- function(z, lambda) {
  return(sign(z) * pmax(0, abs(z) - lambda))
}

mcp_proximal <- function(x, lambda, gamma) {
  abs_x <- abs(x)
  prox <- numeric(length(x))
  for (i in seq_along(x)) {
    if (abs_x[i] <= lambda * gamma) {
      prox[i] <- sign(x[i]) * pmax(abs_x[i] - lambda, 0) / (1 - 1/gamma)
    } else {
      prox[i] <- x[i]
    }
  }
  return(prox)
}

scad_proximal <- function(x, lambda, gamma) {
  abs_x <- abs(x)
  prox <- numeric(length(x))
  for (i in seq_along(x)) {
    if (abs_x[i] <= 2 * lambda) {
      prox[i] <- soft_threshold(x[i], lambda)
    } else if (abs_x[i] <= lambda * gamma) {
      prox[i] <- ((gamma - 1) * x[i] - sign(x[i]) * gamma * lambda) / (gamma - 2)
    } else {
      prox[i] <- x[i]
    }
  }
  return(prox)
}

elasticnet_proximal <- function(x, lambda, alpha, learning_rate) {
  thresh <- learning_rate * lambda * alpha
  scaled_x <- x / (1 + 2 * learning_rate * lambda * (1 - alpha))
  return(soft_threshold(scaled_x, thresh))
}


# perform CDA
perform_CDA <- function(X, y,method, lambda, max_iter = 1000, alpha = 1, gamma = 3.7) {
  n <- nrow(X)
  p <- ncol(X)
  beta <- rep(0, p)
  r <- y - X %*% beta  # 초기 잔차

  for (iter in 1:max_iter) {
    for (j in 1:p) {
      tmp <- r + X[, j] * beta[j]  # 잔차 복원
      rho <- sum(X[, j] * tmp) / n

      if (method == "lasso") {
        beta_j_new <- soft_threshold(rho, lambda)
      } else if (method == "ridge") {
        beta_j_new <- rho / (1 + lambda)
      } else if (method == "elasticnet") {
        beta_j_new <- soft_threshold(rho, lambda * alpha) / (1 + lambda * (1 - alpha))
      } else if (method == "scad") {
        beta_j_new <- scad_proximal(rho, lambda, gamma)
      } else if (method == "mcp") {
        beta_j_new <- mcp_proximal(rho, lambda, gamma)
      } else {
        stop("Unknown penalty type.")
      }

      r <- tmp - X[, j] * beta_j_new  # 잔차 갱신
      beta[j] <- beta_j_new
    }
  }

  return(matrix(beta, ncol = 1, dimnames = list(colnames(X), NULL)))
}




# perform PGD
perform_PGD <- function(X, y, method, lambda, learning_rate, max_iter, alpha = 0.5, gamma = 3.7) {
  # 현재는 Lasso에만 맞춰진 형태이므로, method별로 처리할 수 있게 만들어야 함
  n <- nrow(X)
  p <- ncol(X)
  beta <- rep(0, p)

  for (iter in 1:max_iter) {
    pred <- X %*% beta
    grad <- (2/n) * t(X) %*% (pred - y)

    if (method == "lasso") {
      beta <- soft_threshold(beta - learning_rate * grad, learning_rate * lambda)
    } else if (method == "ridge") {
      beta <- beta - learning_rate * (grad + 2 * lambda * beta)
    } else if (method == "elasticnet") {
      ridge_grad <- 2 * lambda * (1 - alpha) * beta
      lasso_grad <- lambda * alpha * sign(beta)
      beta <- beta - learning_rate * (grad + ridge_grad + lasso_grad)
    } else {
      stop("PGD does not support non-convex penalties like SCAD or MCP.")
    }
  }

  return(beta)
}
######## fista 구현하기 위한 함수들


penalty_proximal <- function(x, method, lambda, learning_rate, alpha=0.5, gamma=3.7) {
  if (method == "lasso") {
    return(soft_threshold(x, learning_rate * lambda))
  } else if (method == "ridge") {
    # Ridge는 proximal 연산 대신 단순 감쇠(scaling)
    return(x / (1 + 2 * learning_rate * lambda))
  } else if (method == "elasticnet") {
    return(elasticnet_proximal(x, lambda, alpha, learning_rate))
  } else if (method == "mcp") {
    return(mcp_proximal(x, learning_rate * lambda, gamma))
  } else if (method == "scad") {
    return(scad_proximal(x, learning_rate * lambda, gamma))
  } else {
    stop("Unsupported penalty in proximal operator.")
  }
}

perform_FISTA <- function(X, y, method, lambda, learning_rate, max_iter, alpha=0.5, gamma=3.7) {
  n <- nrow(X)
  p <- ncol(X)
  beta <- rep(0, p)
  z <- beta
  t_k <- 1

  for (iter in 1:max_iter) {
    pred <- X %*% z
    grad <- (2/n) * t(X) %*% (pred - y)

    beta_new <- penalty_proximal(z - learning_rate * grad, method, lambda, learning_rate, alpha, gamma)
    t_k_new <- (1 + sqrt(1 + 4 * t_k^2)) / 2
    z <- beta_new + ((t_k - 1) / t_k_new) * (beta_new - beta)

    beta <- beta_new
    t_k <- t_k_new
  }

  return(beta)
}
################################################

# Newton
grad_penalty <- function(beta, lambda) {
  return(2 * lambda * beta)
}

perform_Newton <- function(X, y, method, lambda, max_iter, ...) {

  n <- nrow(X)
  p <- ncol(X)
  beta <- rep(0, p)

  for (iter in 1:max_iter) {
    pred <- X %*% beta
    grad_loss <- (2 / n) * t(X) %*% (pred - y)
    grad_pen <- 2 * lambda * beta
    grad <- grad_loss + grad_pen

    hess <- (2 / n) * t(X) %*% X + diag(rep(2 * lambda, p))  # Ridge Hessian

    beta <- beta - solve(hess) %*% grad
  }

  return(beta)
}




##################################################

perform_GD <- function(X, y, method, lambda, learning_rate, max_iter, alpha = 0.5, gamma = 3.7) {
  n <- nrow(X)
  p <- ncol(X)
  beta <- rep(0, p)

  for (iter in 1:max_iter) {
    pred <- X %*% beta
    residual <- pred - y
    grad_loss <- (2 / n) * t(X) %*% residual
    grad_penalty <- penalty_gradient(beta, method, lambda, alpha, gamma)
    gradient <- grad_loss + grad_penalty

    beta <- beta - learning_rate * gradient
  }
  return(beta)
}

penalty_gradient <- function(beta, method, lambda, alpha = 0.5, gamma = 3.7) {
  if (method == "lasso") {
    return(lambda * sign(beta))
  } else if (method == "ridge") {
    return(2 * lambda * beta)
  } else if (method == "elasticnet") {
    return(lambda * (alpha * sign(beta) + 2 * (1 - alpha) * beta))
  } else if (method == "scad") {
    grad <- numeric(length(beta))
    for (i in 1:length(beta)) {
      b <- beta[i]
      if (abs(b) <= lambda) {
        grad[i] <- lambda * sign(b)
      } else if (abs(b) <= gamma * lambda) {
        grad[i] <- ((gamma * lambda - abs(b)) / (gamma - 1)) * sign(b)
      } else {
        grad[i] <- 0
      }
    }
    return(grad)
  } else if (method == "mcp") {
    grad <- numeric(length(beta))
    for (i in 1:length(beta)) {
      b <- beta[i]
      if (abs(b) <= gamma * lambda) {
        grad[i] <- lambda * (1 - abs(b) / (gamma * lambda)) * sign(b)
      } else {
        grad[i] <- 0
      }
    }
    return(grad)
  } else {
    stop("Penalty gradient not implemented for this method")
  }
}

# gd X
# cda O
# FISTA X
# Non-convex algorithms'''

rm(list=ls())

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

        # MCP 업데이트 처리
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
    -t(X) %*% (y - X %*% beta)
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
        pen_grad <- lambda * (if (method == "lasso") sign(z[j])
                              else alpha * sign(z[j]) + 2 * (1 - alpha) * z[j])
        beta_new[j] <- soft_threshold(z[j] - learning_rate * grad_z[j], learning_rate * lambda * alpha)
      } else if (method %in% c("scad", "mcp")) {
        beta_new[j] <- soft_threshold(z[j] - learning_rate * grad_z[j], learning_rate * abs(penalty_grad(z[j])))
      } else {
        stop("Unsupported method in FISTA.")
      }
    }

    beta_old <- beta
    beta <- beta_new
    t <- t + 1
  }

  return(beta)
}

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
