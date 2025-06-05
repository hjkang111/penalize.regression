####################매소드 기준 시각화 연습
plot_LASSO <- function(X, y, lambda_seq, max_iter = 1000, learning_rate = 0.01) {
  coefs <- matrix(0, nrow = length(lambda_seq), ncol = ncol(X))

  for (i in 1:length(lambda_seq)) {
    coefs[i, ] <- penalized_regression(X, y, method = "lasso", algorithm = "cda", lambda = lambda_seq[i], max_iter = max_iter, learning_rate = learning_rate)
  }

  matplot(lambda_seq, t(coefs), type = "l", lty = 1, col = 1:ncol(X),
          xlab = expression(lambda), ylab = "Coefficients", main = "Lasso Path")
}

plot_Ridge <- function(X, y, lambda_seq, max_iter = 1000, learning_rate = 0.01) {
  coefs <- matrix(0, nrow = length(lambda_seq), ncol = ncol(X))

  for (i in 1:length(lambda_seq)) {
    coefs[i, ] <- penalized_regression(X, y, method = "ridge", algorithm = "cda", lambda = lambda_seq[i], max_iter = max_iter, learning_rate = learning_rate)
  }

  matplot(lambda_seq, t(coefs), type = "l", lty = 1, col = 1:ncol(X),
          xlab = expression(lambda), ylab = "Coefficients", main = "Ridge Path")
}

plot_SCAD <- function(X, y, lambda_seq, max_iter = 1000, learning_rate = 0.01) {
  coefs <- matrix(0, nrow = length(lambda_seq), ncol = ncol(X))

  for (i in 1:length(lambda_seq)) {
    coefs[i, ] <- penalized_regression(X, y, method = "scad", algorithm = "cda", lambda = lambda_seq[i], max_iter = max_iter, learning_rate = learning_rate)
  }

  matplot(lambda_seq, t(coefs), type = "l", lty = 1, col = 1:ncol(X),
          xlab = expression(lambda), ylab = "Coefficients", main = "SCAD Path")
}

plot_MCP <- function(X, y, lambda_seq, max_iter = 1000, learning_rate = 0.01) {
  coefs <- matrix(0, nrow = length(lambda_seq), ncol = ncol(X))

  for (i in 1:length(lambda_seq)) {
    coefs[i, ] <- penalized_regression(X, y, method = "mcp", algorithm = "cda", lambda = lambda_seq[i], max_iter = max_iter, learning_rate = learning_rate)
  }

  matplot(lambda_seq, t(coefs), type = "l", lty = 1, col = 1:ncol(X),
          xlab = expression(lambda), ylab = "Coefficients", main = "MCP Path")
}
