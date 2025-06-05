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


#' Penalized Regression Wrapper
#'
#' A unified function to fit penalized regression models using various penalty
#' methods and optimization algorithms. This wrapper handles Ridge, Lasso, Elastic Net,
#' SCAD, and MCP penalties, with support for Coordinate Descent (CDA), FISTA,
#' and Local Linear Approximation (LLA) algorithms.
#'
#' @param X A numeric matrix of predictors (n x p).
#' @param y A numeric response vector of length n.
#' @param method A character string specifying the penalty type:
#'   one of "lasso", "ridge", "elasticnet", "scad", or "mcp".
#' @param algorithm A character string specifying the optimization algorithm:
#'   one of "cda" (Coordinate Descent Algorithm), "fista" (Fast Iterative Shrinkage-Thresholding Algorithm),
#'   or "lla" (Local Linear Approximation).
#' @param lambda A non-negative regularization parameter controlling the overall penalty strength.
#' @param learning_rate A positive numeric learning rate for gradient-based algorithms (default: 0.01).
#' @param max_iter An integer specifying the maximum number of iterations (default: 1000).
#' @param alpha Elastic Net mixing parameter (between 0 and 1, default: 0.5). Only used if \code{method = "elasticnet"}.
#' @param gamma A tuning parameter for SCAD and MCP (typically >= 3.7, default: 3.7).
#'
#' @return A numeric vector of estimated regression coefficients of length p.
#'
#' @details
#' The function automatically checks for compatibility between the selected
#' penalty method and optimization algorithm, issuing warnings or errors if
#' inappropriate combinations are chosen.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' X <- matrix(rnorm(100*10), 100, 10)
#' y <- rnorm(100)
#' beta <- penalized_regression(X, y, method = "lasso", algorithm = "cda", lambda = 0.1)
#' }
#'
#' @export
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
