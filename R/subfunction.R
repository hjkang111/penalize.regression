# penalty function

penalty_function <- function(beta, method, lambda, alpha = 0.5, gamma = 3.7) {
  # method: 'lasso', 'ridge', 'elasticnet', 'scad', 'mcp'
  # lambda: 패널티 강도 (정규화 계수)
  # alpha: Elastic Net에서 Lasso와 Ridge의 비율 (0은 Ridge, 1은 Lasso)
  # gamma: SCAD와 MCP에서 사용하는 매개변수 (SCAD에서 3.7이 표준)

  if (method == "lasso") {                     # lasso
    return(lambda * sum(abs(beta)))
  }

  else if (method == "ridge") {                  # ridge
    return(lambda * sum(beta^2))
  }

  else if (method == "elasticnet") {        # elastic net
    # alpha 0: Ridge, 1: Lasso
    return(lambda * (alpha * sum(abs(beta)) + (1 - alpha) * sum(beta^2)))
  }

  else if (method == "scad") {             # scad
    penalty <- 0
    for (i in 1:length(beta)) {
      if (abs(beta[i]) <= lambda) {
        penalty <- penalty + lambda * abs(beta[i])  # Lasso처럼 패널티
      } else if (abs(beta[i]) > lambda & abs(beta[i]) <= gamma * lambda) {
        # SCAD 중간 범위: 점차적으로 패널티를 줄임
        penalty <- penalty + lambda * (abs(beta[i]) - 0.5 * (gamma * lambda - abs(beta[i])))
      } else {
        # SCAD 큰 값: 패널티를 거의 주지 않음
        penalty <- penalty + 0.5 * gamma * lambda^2
      }
    }
    return(penalty)
  }


  else if (method == "mcp") {         # mcp
    penalty <- 0
    for (i in 1:length(beta)) {
      if (abs(beta[i]) <= lambda) {
        penalty <- penalty + lambda * abs(beta[i])  # Lasso처럼 패널티
      } else if (abs(beta[i]) > lambda & abs(beta[i]) <= gamma * lambda) {
        # MCP 중간 범위: 부드럽게 패널티를 감소
        penalty <- penalty + lambda * (abs(beta[i]) - (abs(beta[i])^2 / (2 * lambda)))
      } else {
        # MCP 큰 값: 거의 패널티가 없음
        penalty <- penalty + 0.5 * gamma * lambda^2
      }
    }
    return(penalty)
  } else {
    stop("Unknown method. Please choose one of 'lasso', 'ridge', 'elasticnet', 'scad', or 'mcp'.")
  }
}
