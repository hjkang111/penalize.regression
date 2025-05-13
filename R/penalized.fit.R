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


## 함수 구조 개요
#penalized_fit <- function(X, y, type = "lasso", method = "cda", lambda = 0.1, ...)
  # type: "lasso", "ridge", "elastic", "scad", "mcp", "adaptive"
  # method: "cda", "gd", "pgd"


soft_threshold <- function(z, gamma) {
  sign(z) * pmax(abs(z) - gamma, 0)
}

scad_derivative <- function(beta, lambda, a = 3.7) {
  deriv <- numeric(length(beta))
  for (j in seq_along(beta)) {
    b <- abs(beta[j])
    if (b <= lambda) {
      deriv[j] <- lambda * sign(beta[j])
    } else if (b <= a * lambda) {
      deriv[j] <- ((a * lambda - b) / (a - 1)) * sign(beta[j])
    } else {
      deriv[j] <- 0
    }
  }
  return(deriv)
}

mcp_derivative <- function(beta, lambda, gamma = 3) {
  deriv <- numeric(length(beta))
  for (j in seq_along(beta)) {
    b <- abs(beta[j])
    if (b <= gamma * lambda) {
      deriv[j] <- (lambda - b / gamma) * sign(beta[j])
    } else {
      deriv[j] <- 0
    }
  }
  return(deriv)
}

# 가중치 벡터를 따로?
adaptive_lasso_derivative <- function(beta, weights) {
  return(weights * sign(beta))
}

get_penalty <- function(beta, type, lambda, alpha = 0.5) {
  if (type == "ridge") return(lambda * beta)
  if (type == "lasso") return(lambda * sign(beta))
  if (type == "elastic") return(lambda * (alpha * sign(beta) + (1 - alpha) * beta))
  if (type == "scad") return(scad_derivative(beta, lambda))
  if (type == "mcp") return(mcp_derivative(beta, lambda))
  if (type == "adlasso") return(lambda * sign(beta) / (abs(beta) + 1e-6))
  stop("지원하지 않는 type입니다.")
}

penalized_fit <- function(X, y, type = "lasso", method = "cda", lambda = 0.1,
                          step_size = 0.01, max_iter = 1000, tol = 1e-6, alpha = 0.5) {
  n <- nrow(X)
  p <- ncol(X)
  beta <- rep(0, p)


  if (method == "gd") {
    for (iter in 1:max_iter) {
      grad <- -t(X) %*% (y - X %*% beta) / n
      penalty_grad <- get_penalty(beta, type, lambda, alpha)
      beta <- beta - step_size * (grad + penalty_grad)
    }
  } else if (method == "pgd") {
    for (iter in 1:max_iter) {
      grad <- -t(X) %*% (y - X %*% beta) / n
      z <- beta - step_size * grad
      if (type == "lasso") {
        beta <- soft_threshold(z, lambda * step_size)
      } else if (type == "elastic") {
        beta <- soft_threshold(z, lambda * alpha * step_size) / (1 + lambda * (1 - alpha) * step_size)
      } else {
        stop("PGD는 현재 lasso/elastic만 지원합니다.")
      }
    }
  } else if (method == "cda") {
    for (iter in 1:max_iter) {
      beta_old <- beta
      for (j in 1:p) {
        r_j <- y - X[, -j] %*% beta[-j]
        zj <- sum(X[, j]^2)
        pj <- sum(X[, j] * r_j)

        if (type == "lasso") {
          beta[j] <- soft_threshold(pj / zj, lambda / zj)
        } else if (type == "ridge") {
          beta[j] <- pj / (zj + lambda)
        } else if (type == "elastic") {
          beta[j] <- soft_threshold(pj / zj, lambda * alpha / zj) / (1 + lambda * (1 - alpha) / zj)
        } else {
          stop("CDA는 현재 ridge/lasso/elastic만 지원합니다.")
        }
      }
      if (sum(abs(beta - beta_old)) < tol) break
    }
  } else {
    stop("지원하지 않는 알고리즘 method입니다.")
  }

  return(list(coefficients = beta,
              type = type,
              method = method))
  return(beta)
}

# main function example
penalized_regression <- function(X, y, method = c("lasso", "ridge", "scad"),
                                 algorithm = c("GD", "CDA"), lambda = 1,
                                 learning_rate = 0.01, max_iter = 1000) {
  # method와 algorithm 선택
  method <- match.arg(method)
  algorithm <- match.arg(algorithm)

  # 알고리즘과 메소드에 맞는 서브 함수 호출
  if (algorithm == "GD") {
    return(perform_GD(X, y, method, lambda, learning_rate, max_iter))
  } else if (algorithm == "CDA") {
    return(perform_CDA(X, y, method, lambda, max_iter))
  }
}
