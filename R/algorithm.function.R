#다른 파일에 있는 penalty_function 불러와서 이용할 것임
#library(penalize.regression)
#library(penalize.regression)


# gradient descent
perform_GD <- function(X, y, method, lambda, learning_rate, max_iter) {
  n <- nrow(X)
  p <- ncol(X)

  # 초기 회귀 계수 설정
  beta <- rep(0, p)

  # 반복문을 통해 Gradient Descent 최적화
  for (iter in 1:max_iter) {
    # 예측값 계산
    pred <- X %*% beta
    residual <- pred - y

    # Gradient 계산
    gradient <- (2/n) * t(X) %*% residual + penalty_function(beta, method, lambda)

    beta <- beta - learning_rate * gradient
  }

  return(beta)
}


# coordinate descent
perform_CDA <- function(X, y, method, lambda, max_iter) {
  n <- nrow(X)
  p <- ncol(X)

  beta <- rep(0, p)

  # 반복문을 통해 Coordinate Descent 최적화
  for (iter in 1:max_iter) {
    for (j in 1:p) {
      # 예측값 계산
      pred <- X %*% beta
      residual <- y - pred + X[, j] * beta[j]

      # 새로운 계수 계산 (Soft Thresholding)
      beta[j] <- sign(t(X[, j]) %*% residual) * max(0, abs(t(X[, j]) %*% residual) - lambda)
    }
  }

  return(beta)
}
