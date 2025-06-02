# penalize regression

# check algorithm
check_algorithm_compatibility <- function(method, algorithm) {
  method <- tolower(method)
  algorithm <- tolower(algorithm)

  if (method == "ridge" && algorithm == "cda") {
    warning("CDA may not be the most appropriate choice for Ridge penalty.")
  }
  if (method %in% c("scad", "mcp") && algorithm == "fista") {
    stop("FISTA is not recommended for non-convex penalties like SCAD or MCP.")
  }
  if (!method %in% c("scad", "mcp") && algorithm == "lla") {
    warning("LLA is primarily designed for SCAD or MCP penalties and may not be optimal for Ridge, Lasso, or Elastic Net.")
  }
}


penalized_regression <- function(X, y,
                                 method = c("lasso", "ridge", "scad", "mcp", "elasticnet"),
                                 algorithm = c("cda", "fista", "lla"),
                                 lambda = 1, learning_rate = 0.01,
                                 max_iter = 1000, alpha = 0.5, gamma = 3.7) {
  method <- match.arg(method)
  algorithm <- match.arg(algorithm)

  check_algorithm_compatibility(method, algorithm)

  if (algorithm == "cda") {
    return(perform_CDA(X, y, method, lambda, learning_rate, max_iter, alpha, gamma))
  } else if (algorithm == "fista") {
    return(perform_FISTA(X, y, method, lambda, learning_rate, max_iter, alpha, gamma))
  } else if (algorithm == "lla") {
    return(perform_LLA(X, y, method, lambda, learning_rate, max_iter, gamma))
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

X <- scale(X)
y <- scale(y)

#ridge
penalized_regression(X, y, method="ridge", algorithm = "cda")
penalized_regression(X, y, method="ridge", algorithm = "fista") #
penalized_regression(X, y, method="ridge", algorithm = "lla")

#lasso
penalized_regression(X, y, method="lasso", algorithm = "cda")
penalized_regression(X, y, method="lasso", algorithm = "fista") #
penalized_regression(X, y, method="lasso", algorithm = "lla")

#elastic net
penalized_regression(X, y, method="elasticnet", algorithm = "cda")
penalized_regression(X, y, method="elasticnet", algorithm = "fista") #
penalized_regression(X, y, method="elasticnet", algorithm = "lla")
# gd, cda,pgd랑 값이 같고, fista값이 같음

#mcp
beta_hat <- penalized_regression(X, y, method="mcp", algorithm = "cda")
beta_hat <- penalized_regression(X, y, method="mcp", algorithm = "fista") #
beta_hat <- penalized_regression(X, y, method="mcp", algorithm = "lla") #

# scad
beta_hat <- penalized_regression(X, y, method="scad", algorithm = "cda")

y_pred <- X %*% beta_hat

# MSE 계산
mean((y - y_pred)^2)


beta_hat <- penalized_regression(X, y, method="scad", algorithm = "fista")
beta_hat <- penalized_regression(X, y, method="scad", algorithm = "lla")


lambda_seq <- seq(0.1, 2, length.out = 50)
plot_LASSO(X,y,lambda_seq)
length(X)
length(y)
