# penalize regression

# check algorithm
check_algorithm_compatibility <- function(method, algorithm) {
  # SCAD, MCP는 PGD, FISTA, NEWTON과 함께 쓰면 수렴 문제 가능
  if (method %in% c("scad", "mcp") && algorithm %in% c("PGD", "FISTA")) {
    warning(sprintf("Algorithm '%s' may have convergence issues with '%s' penalty due to non-convexity. Use with caution.", algorithm, method))
  }

  if (method %in% c("scad", "mcp") && algorithm == "NEWTON") {
    warning(sprintf("Algorithm '%s' may not converge reliably with non-convex '%s' penalty. Consider alternative algorithms.", algorithm, method))
  }

  # Ridge + PGD, FISTA는 비효율적일 수 있음
  if (method == "ridge" && algorithm %in% c("PGD", "FISTA")) {
    warning(sprintf("Algorithm '%s' is typically inefficient for convex '%s' penalty. Consider GD or NEWTON instead.", algorithm, method))
  }

  # GD + SCAD, MCP: 가능하지만 느리고 불안정할 수 있음
  if (method %in% c("scad", "mcp") && algorithm == "GD") {
    warning(sprintf("Algorithm '%s' may converge slowly or unstably with non-convex '%s' penalty. Use carefully.", algorithm, method))
  }

  # NEWTON + Lasso, Elastic Net은 경고 불필요 (볼록 문제)
  # 별도 처리 안 함 (필요 시 추가 가능)
}


penalized_regression <- function(X, y,
                                 method = c("lasso", "ridge", "scad", "mcp", "elasticnet"),
                                 algorithm = c("GD", "CDA", "PGD", "FISTA", "NEWTON"),
                                 lambda = 1, learning_rate = 0.01, max_iter = 1000, alpha = 0.5, gamma = 3.7) {
  # method와 algorithm 선택
  method <- match.arg(method)
  algorithm <- match.arg(algorithm)

  # warning
  check_algorithm_compatibility(method, algorithm)

  # 알고리즘 함수 호출 (가정: 각각 함수는 X, y, method, lambda ... 인자를 받음)
  if (algorithm == "GD") {
    return(perform_GD(X, y, method, lambda, learning_rate, max_iter))
  } else if (algorithm == "CDA") {
    return(perform_CDA(X, y, method, lambda, max_iter))
  } else if (algorithm == "PGD") {
    return(perform_PGD(X, y, method, lambda, learning_rate, max_iter, alpha, gamma))
  } else if (algorithm == "FISTA") {
    return(perform_FISTA(X, y, method, lambda, learning_rate, max_iter, alpha, gamma))
  } else if (algorithm == "NEWTON") {
    return(perform_NEWTON(X, y, method, lambda, learning_rate, max_iter, alpha, gamma))
  } else {
    stop("Unknown algorithm selected.")
  }
}


# example
## Data 2 -> Generating data
set.seed(924)

# predictor X
X1 <- rnorm(100, mean = 50, sd = 10)  # 평균 50, 표준편차 10
X2 <- rnorm(100, mean = 30, sd = 5)   # 평균 30, 표준편차 5
X3 <- rnorm(100, mean = 100, sd = 20)  # 평균 100, 표준편차 20
X4 <- rnorm(100, mean = 101, sd = 3)  # 평균 101, 표준편차 3
X5 <- rnorm(100, mean = 57, sd = 100)  # 평균 57, 표준편차 100
X6 <- rnorm(100, mean = 88, sd = 20)  # 평균 88, 표준편차 58
X7 <- rnorm(100, mean = 34, sd = 9)  # 평균 34, 표준편차 9
X8 <- rnorm(100, mean = 83, sd = 23)  # 평균 83, 표준편차 23


# Response Y
eps = rnorm(100, mean = 0, sd = 10)
Y <- 5 + 3*X1 + 2*X2 - 0.5*X3 + X4 + 1.3*X5 - 0.4*X6 + 2*X7 - X8 + eps


Data2 <- data.frame(X1 = X1, X2 = X2, X3 = X3,X4=X4, X5=X5,X6=X6,X7=X7,X8=X8, Y = Y)
head(Data2)
X <- as.matrix(Data2[, 1:8])
y <- Data2[[9]]

penalized_regression(as.matrix(Data2[, 1:8]), Data2[[9]], method="ridge", algorithm = "GD")

str(Data2[, 1:8])
str(Data2[[9]])
