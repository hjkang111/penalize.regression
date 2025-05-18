# regularization path plot
plot_reg_path <- function(lambda_seq, beta_mat, method) {
  matplot(lambda_seq, t(beta_mat), type = "l", lty = 1, col = 1:nrow(beta_mat),
          xlab = expression(lambda), ylab = "Coefficients",
          main = paste("Regularization Path -", method))
  legend("topright", legend = paste0("Beta", 1:nrow(beta_mat)), col = 1:nrow(beta_mat), lty = 1)
}

# grid search + cross validation
tune_lambda <- function(X, y, method, algorithm, lambda_grid, folds = 5) {
  n <- nrow(X)
  fold_ids <- sample(rep(1:folds, length.out = n))
  cv_errors <- numeric(length(lambda_grid))

  for (i in seq_along(lambda_grid)) {
    lambda <- lambda_grid[i]
    fold_mse <- numeric(folds)

    for (k in 1:folds) {
      train_idx <- which(fold_ids != k)
      val_idx <- which(fold_ids == k)

      beta <- penalized_regression(X[train_idx, ], y[train_idx], method, algorithm, lambda)
      pred <- X[val_idx, ] %*% beta
      fold_mse[k] <- mean((y[val_idx] - pred)^2)
    }

    cv_errors[i] <- mean(fold_mse)
  }

  best_lambda <- lambda_grid[which.min(cv_errors)]
  list(best_lambda = best_lambda, cv_errors = cv_errors)
}

# residual plot
plot_residuals <- function(X, y, beta) {
  pred <- X %*% beta
  residuals <- y - pred
  plot(pred, residuals, xlab = "Predicted values", ylab = "Residuals", main = "Residual Plot")
  abline(h = 0, col = "red", lty = 2)
}

# variable importance
plot_var_importance <- function(beta) {
  barplot(abs(beta), main = "Variable Importance (abs coefficients)",
          ylab = "Absolute Coefficient", col = "steelblue")
}
