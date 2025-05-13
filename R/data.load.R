# function 2
# Data
install.packages("MASS")
library(MASS)
data(Boston)

Boston <- MASS::Boston
usethis::use_data(Boston, overwrite = TRUE)

