# algorithm functions

# soft threshold function
soft_threshold <- function(z, lambda) {
  return(sign(z) * pmax(0, abs(z) - lambda))
}

# CGD
perform_CDA <- function(X, y, method, lambda, learning_rate, max_iter) {
  n <- nrow(X)
  p <- ncol(X)
  beta <- rep(0, p)

  for (iter in 1:max_iter) {
    for (j in 1:p) {
      pred <- X %*% beta
      residual <- pred - y
      grad_j <- (2/n) * sum(X[, j] * residual) + grad_penalty(beta[j], method, lambda)

      beta[j] <- beta[j] - learning_rate * grad_j
    }
  }

  return(beta)
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



# FISTA
perform_FISTA <- function(X, y, lambda, learning_rate, max_iter) {
  n <- nrow(X)
  p <- ncol(X)
  beta <- rep(0, p)
  z <- beta
  t_k <- 1

  for (iter in 1:max_iter) {
    pred <- X %*% z
    grad <- (2/n) * t(X) %*% (pred - y)

    beta_new <- soft_threshold(z - learning_rate * grad, learning_rate * lambda)
    t_k_new <- (1 + sqrt(1 + 4 * t_k^2)) / 2
    z <- beta_new + ((t_k - 1) / t_k_new) * (beta_new - beta)

    beta <- beta_new
    t_k <- t_k_new
  }

  return(beta)
}

# Newton
perform_Newton <- function(X, y, lambda, max_iter) {
  n <- nrow(X)
  p <- ncol(X)
  beta <- rep(0, p)

  for (iter in 1:max_iter) {
    pred <- X %*% beta
    grad <- (2/n) * t(X) %*% (pred - y)
    hess <- (2/n) * t(X) %*% X + diag(rep(1e-6, p))  # 정칙화 추가

    beta <- beta - solve(hess) %*% grad
  }

  return(as.vector(beta))
}


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

