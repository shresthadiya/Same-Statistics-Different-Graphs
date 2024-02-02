
# Load necessary libraries
library(ggplot2)
library(MASS)
library(e1071)

# Choose a dataset from Anscombe's Quartet
data <- anscombe[, c("x1", "y1")]
data

# Perform Gram-Schmidt orthogonalization
U <- matrix(0, nrow = nrow(data), ncol = ncol(data))
U[, 1] <- data[, 1] / sqrt(sum(data[, 1]^2))

for (i in 2:ncol(data)) {
  U[, i] <- data[, i]
  for (j in 1:(i - 1)) {
    U[, i] <- U[, i] - (sum(U[, j] * U[, i]) * U[, j])
  }
  U[, i] <- U[, i] / sqrt(sum(U[, i]^2))
}

# Display the corresponding orthogonal matrix
cat("Corresponding orthogonal matrix is\n")
print(U)

# Perform the transformation
cat("After transformation, the matrix is\n")
cov_matrix <- cov(data)
Y <- sqrt(10) * U %*% chol(cov_matrix) + matrix(rep(colMeans(data), each = nrow(data)), nrow = nrow(data), ncol = ncol(data), byrow = TRUE)
print(Y)

# Scatter plot with Cook's D values
ggplot(data, aes(x = x1, y = y1)) +
  geom_point(aes(size = cooksd)) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Linear Regression with Cook's D Values")

# Scatter plot with skewness differential
ggplot(data, aes(x = fitted(model), y = residuals(model))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "Scatter Plot with Skewness Differential")

# Scatter plot with ks differential
ggplot(data.frame(value = c(X,Y), sample = rep(c("X", "Y"), each = length(X)))) +
  geom_point(aes(x = value, y = 0, color = sample)) +
  labs(x = "Variable", y = "", title = "Scatter Plot with Kolmogorov-Smirnov Differential") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Fit a linear regression model
model <- lm(y1 ~ x1, data = data)

# Cook's D values
cooksd <- cooks.distance(model)

# KS values
kstest <- ks.test(X,Y)
kstest

# Skewness of the residuals
residuals_skewness <- e1071::skewness(residuals(model))

# Print results
cat("Cook's D values:\n", cooksd, "\n")
cat("Skewness of residuals:\n", residuals_skewness, "\n")