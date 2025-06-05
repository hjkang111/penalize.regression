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

install.packages("roxygen2")
library(roxygen2)
devtools::load_all()
??perform_CDA

library(roxygen2)
roxygen2::roxygenise()
detach("package:penalize.regression", unload=TRUE)

search()
install.packages(c("MASS", "rlang", "cli", "ellipsis", "glue", "R6", "fansi", "pillar"))

getOption("repos")
# named vector 형태로 https://cran.rstudio.com/ 가 보여야 정상입니다.
remotes::install_local("C:/Users/naomi/OneDrive/Desktop/penalize.regression")
remotes::install_local("C:/Users/naomi/OneDrive/Desktop/penalize.regression", repos = "https://cran.rstudio.com/")



# 패키지 소스 빌드
devtools::build("C:/Users/naomi/OneDrive/Desktop/penalize.regression")

# 빌드된 tar.gz 파일 이름을 확인 후 (예: penalize.regression_0.1.0.tar.gz)
install.packages("C:/Users/naomi/OneDrive/Desktop/penalize.regression_0.1.0.tar.gz", repos = NULL, type = "source")
getOption("repos")
# CRAN 값이 https://cran.rstudio.com/ 으로 제대로 설정되어 있어야 합니다.
install.packages(c("MASS"))

install.packages("C:/Users/naomi/OneDrive/Desktop/penalize.regression_0.1.0.tar.gz", repos = NULL, type = "source", dependencies = TRUE)
