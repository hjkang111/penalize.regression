# algorithm functions

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
perform_CDA <- function(X, y, lambda, penalty, max_iter = 1000, alpha = 1, gamma = 3.7) {
  penalty <- tolower(penalty)  # 소문자로 통일
  n <- nrow(X)
  p <- ncol(X)
  beta <- rep(0, p)
  r <- y - X %*% beta  # 초기 잔차

  for (iter in 1:max_iter) {
    for (j in 1:p) {
      tmp <- r + X[, j] * beta[j]  # 잔차 복원
      rho <- sum(X[, j] * tmp) / n

      if (penalty == "lasso") {
        beta_j_new <- soft_threshold(rho, lambda)
      } else if (penalty == "ridge") {
        beta_j_new <- rho / (1 + lambda)
      } else if (penalty == "elasticnet") {
        beta_j_new <- soft_threshold(rho, lambda * alpha) / (1 + lambda * (1 - alpha))
      } else if (penalty == "scad") {
        beta_j_new <- scad_proximal(rho, lambda, gamma)
      } else if (penalty == "mcp") {
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

