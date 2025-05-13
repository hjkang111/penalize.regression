
## Data 1 -> Boston data
install.packages("MASS")
library(MASS)
data(Boston)

Boston <- MASS::Boston
usethis::use_data(Boston, overwrite = TRUE)
# 패키지 사용할 때 data(Boston)

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
usethis::use_data(Data2, overwrite = TRUE)
