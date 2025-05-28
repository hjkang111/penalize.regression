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





compare_algorithms <- function(beta_cda, beta_lla, beta_fista,
                               method = "lasso",
                               feature_names = NULL,
                               title = NULL,
                               digits = 4) {
  if (is.null(feature_names)) {
    feature_names <- paste0("V", seq_along(beta_cda))
  }

  df <- data.frame(
    Feature = feature_names,
    CDA = round(beta_cda, digits),
    LLA = round(beta_lla, digits),
    FISTA = round(beta_fista, digits)
  )

  library(reshape2)
  library(ggplot2)

  df_long <- melt(df, id.vars = "Feature", variable.name = "Algorithm", value.name = "Coefficient")

  # 시각화: 계수 비교 바그래프
  ggplot(df_long, aes(x = Feature, y = Coefficient, fill = Algorithm)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(title = paste0("Coefficient Comparison: ", toupper(method),
                        if (!is.null(title)) paste0(" (", title, ")") else ""),
         y = "Estimated Coefficient") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # 수치 비교
  l2_diff_cda_lla <- sqrt(sum((beta_cda - beta_lla)^2))
  l2_diff_cda_fista <- sqrt(sum((beta_cda - beta_fista)^2))
  l2_diff_lla_fista <- sqrt(sum((beta_lla - beta_fista)^2))

  cat("\n[ L2 Distance between coefficient vectors ]\n")
  cat(sprintf("CDA vs LLA   : %.4f\n", l2_diff_cda_lla))
  cat(sprintf("CDA vs FISTA : %.4f\n", l2_diff_cda_fista))
  cat(sprintf("LLA vs FISTA : %.4f\n", l2_diff_lla_fista))

  return(invisible(df))
}


# 예시: 세 알고리즘으로 추정한 회귀계수 벡터 -- 나중에 수정
beta_cda <- penalized_regression(X,y,'lasso','cda')
beta_lla <- penalized_regression(X,y,'lasso','lla')
beta_fista <- penalized_regression(X,y,'lasso','fista')

compare_algorithms(beta_cda, beta_lla,
                   method = "lasso",
                   feature_names = c("X1", "X2", "X3", "X4","X5","X6","X7","X8"),
                   title = "lambda = 0.1")

