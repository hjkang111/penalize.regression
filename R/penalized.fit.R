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
