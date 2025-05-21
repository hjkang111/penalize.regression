# penalize regression

# check algorithm
check_algorithm_compatibility <- function(method, algorithm) {
  method <- tolower(method)
  algorithm <- tolower(algorithm)

  # Incompatible combinations – stop with error
  if (method == "ridge" && algorithm == "cda") {
    stop("CDA is not applicable to Ridge (L2) penalty because it assumes an L1 structure.")
  }
  if (method == "mcp" && algorithm %in% c("pgd", "fista")) {
    stop("PGD and FISTA are not suitable for non-convex MCP penalty due to lack of convergence guarantee.")
  }
  if (method == "scad" && algorithm %in% c("pgd", "fista")) {
    stop("PGD and FISTA are not suitable for non-convex SCAD penalty due to lack of convergence guarantee.")
  }
  if (method %in% c("lasso","elasticnet","mcp", "scad") && algorithm == "newton") {
    stop("Newton method currently supports only Ridge penalty")
  }

  # Inefficient or unstable combinations – give a warning
  if (method == "ridge" && algorithm %in% c("pgd", "fista")) {
    warning("PGD or FISTA may be inefficient for convex Ridge penalty. Consider using GD or NEWTON instead.")
  }
  if (method %in% c("mcp", "scad") && algorithm == "gd") {
    warning("GD may be slow or unstable with non-convex penalties such as SCAD or MCP.")
  }
  if (method %in% c("mcp", "scad") && algorithm == "newton") {
    warning("NEWTON may not converge reliably for non-convex penalties such as SCAD or MCP.")
  }
}



penalized_regression <- function(X, y,
                                 method = c("lasso", "ridge", "scad", "mcp", "elasticnet"),
                                 algorithm = c("GD", "CDA", "PGD", "FISTA", "NEWTON"),
                                 lambda = 1, learning_rate = 0.01, max_iter = 1000, alpha = 0.5, gamma = 3.7) {
  method <- match.arg(method)
  algorithm <- match.arg(algorithm)

  check_algorithm_compatibility(method, algorithm)

  if (algorithm == "GD") {
    return(perform_GD(X, y, method, lambda, learning_rate, max_iter))
  } else if (algorithm == "CDA") {
    return(perform_CDA(X, y, method, lambda, learning_rate, max_iter))
  } else if (algorithm == "PGD") {
    return(perform_PGD(X, y, method, lambda, learning_rate, max_iter, alpha, gamma))
  } else if (algorithm == "FISTA") {
    # 여기서 method, alpha, gamma 모두 넘겨줘야 함
    return(perform_FISTA(X, y, method, lambda, learning_rate, max_iter, alpha, gamma))
  } else if (algorithm == "NEWTON") {
    return(perform_Newton(X, y, method, lambda, max_iter))
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
penalized_regression(X, y, method="ridge", algorithm = "GD")
penalized_regression(X, y, method="ridge", algorithm = "CDA")
penalized_regression(X, y, method="ridge", algorithm = "PGD")
penalized_regression(X, y, method="ridge", algorithm = "FISTA")
penalized_regression(X, y, method="ridge", algorithm = "NEWTON")

#lasso
penalized_regression(X, y, method="lasso", algorithm = "GD")
penalized_regression(X, y, method="lasso", algorithm = "CDA") #
penalized_regression(X, y, method="lasso", algorithm = "PGD")
penalized_regression(X, y, method="lasso", algorithm = "FISTA")
penalized_regression(X, y, method="lasso", algorithm = "NEWTON")
# gd, cda랑 값이 같고, pgd랑 fista값이 같음.. 왜 다르지

#elastic net
penalized_regression(X, y, method="elasticnet", algorithm = "GD")
penalized_regression(X, y, method="elasticnet", algorithm = "CDA") #
penalized_regression(X, y, method="elasticnet", algorithm = "PGD")
penalized_regression(X, y, method="elasticnet", algorithm = "FISTA")
penalized_regression(X, y, method="elasticnet", algorithm = "NEWTON")
# gd, cda,pgd랑 값이 같고, fista값이 같음

#mcp
penalized_regression(X, y, method="mcp", algorithm = "GD")
penalized_regression(X, y, method="mcp", algorithm = "CDA") #
penalized_regression(X, y, method="mcp", algorithm = "PGD")
penalized_regression(X, y, method="mcp", algorithm = "FISTA")
penalized_regression(X, y, method="mcp", algorithm = "NEWTON")

# scad
penalized_regression(X, y, method="scad", algorithm = "GD")
penalized_regression(X, y, method="scad", algorithm = "CDA") #
penalized_regression(X, y, method="scad", algorithm = "PGD")
penalized_regression(X, y, method="scad", algorithm = "FISTA")
penalized_regression(X, y, method="scad", algorithm = "NEWTON")


# 1 cda 알고리즘 수정
# 2 알고리즘마다 반환값 왜 다른지
# method랑 algorithm 잘 맞는지
# 옵션값 잘 돌아가는지
