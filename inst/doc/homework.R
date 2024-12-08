## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----cas----------------------------------------------------------------------
# Clear memory
rm(list = ls())

# Function to generate random samples from Beta(3, 3)
generate_samples <- function(n) {
  return(rbeta(n, 3, 3))
}

# Function to estimate the CDF using Monte Carlo simulation
monte_carlo_estimate <- function(x, n_samples = 10000) {
  samples <- generate_samples(n_samples)
  return(mean(samples <= x))
}

# Function to compare Monte Carlo estimates with pbeta
compare_estimates <- function(x_values) {
  monte_carlo_results <- sapply(x_values, monte_carlo_estimate)
  pbeta_results <- pbeta(x_values, 3, 3)
  
  results <- data.frame(x = x_values, 
                        Monte_Carlo_Estimate = monte_carlo_results, 
                        pbeta_Value = pbeta_results)
  
  return(results)
}

# Values to estimate
x_values <- seq(0.1, 0.9, by = 0.1)

# Compare estimates
results <- compare_estimates(x_values)
print(results)


## ----j------------------------------------------------------------------------
# Clear memory
rm(list = ls())

# Function to compute the CDF of the Rayleigh distribution
calculate_cdf <- function(x, sigma = 2) {
  1 - exp(-((x)^2) / (2 * (sigma^2)))
}

# Monte Carlo simulation without antithetic variables
mc_simulation <- function(x, r = 10000) {
  y <- runif(r)
  g <- (y * x^2) * exp(-(y * x^2) / (2 * (2^2)))
  mean(g) / (2^2)
}

# Monte Carlo simulation with antithetic variables
mc_simulation_antithetic <- function(x, r = 10000) {
  y <- runif(r / 2)
  v <- 1 - y
  y <- c(y, v)
  g <- (y * x^2) * exp(-(y * x^2) / (2 * (2^2)))
  mean(g) / (2^2)
}

# Main code execution
x <- seq(0.1, 2.5, length = 10)
set.seed(12345)

Phi <- calculate_cdf(x)  # Results from Rayleigh(2) 
MC1 <- sapply(x, mc_simulation, r = 10000)
MC2 <- sapply(x, mc_simulation_antithetic, r = 10000)

print(round(rbind(x, Phi, MC1, MC2), 5)) # Comparison of estimates

# Variance reduction assessment
m <- 1000
MC1_samples <- numeric(m)
MC2_samples <- numeric(m)
x_val <- 1.65

for (i in 1:m) {
  MC1_samples[i] <- mc_simulation(x_val, r = 1000)
  MC2_samples[i] <- mc_simulation_antithetic(x_val, r = 1000)
}

print(sd(MC1_samples))
print(sd(MC2_samples))
print((var(MC1_samples) - var(MC2_samples)) / var(MC1_samples))


## ----yu-----------------------------------------------------------------------
# Clear memory
rm(list = ls())

# Define the probability density functions
g_function <- function(x) {
  (1 / sqrt(2 * pi)) * x^2 * exp(-x^2 / 2) * (x > 1)
}

f1_function <- function(x) {
  dnorm(x, 2, sqrt(1.1))
}

f2_function <- function(x) {
  dgamma(x, 3, 2)
}

# Plotting functions
plot_density_functions <- function(x) {
  g <- g_function(x)
  f1 <- f1_function(x)
  f2 <- f2_function(x)

  # Figure 1
  plot(x, g, type = "l", main = "", ylab = "", ylim = c(0, 0.4), lwd = 2)
  lines(x, f1, lty = 2, lwd = 2, col = "red")
  lines(x, f2, lty = 3, lwd = 2, col = "blue")
  legend("topright", legend = c("g", "f1", "f2"), lty = 1:3, col = c("black", "red", "blue"), lwd = 2, inset = 0.02)

  # Figure 2
  plot(x, g, type = "l", main = "", ylab = "", ylim = c(0, 3), lwd = 2)
  lines(x, g / f1, lty = 2, lwd = 2, col = "red")
  lines(x, g / f2, lty = 3, lwd = 2, col = "blue")
  legend("topright", legend = c("g", "f1", "f2"), lty = 1:3, col = c("black", "red", "blue"), lwd = 2, inset = 0.02)
}

# Standard error calculation function
calculate_standard_errors <- function(m) {
  se <- numeric(2)
  
  set.seed(123)
  x <- rnorm(m, 2, sqrt(1.1))
  fg <- g_function(x) / f1_function(x)
  se[1] <- sd(fg)

  x <- rgamma(m, 3, 2)
  fg <- g_function(x) / f2_function(x)
  se[2] <- sd(fg)

  return(round(se, 5))
}

# Monte Carlo estimate function
monte_carlo_estimate <- function(m) {
  theta_hat <- numeric(2)
  
  set.seed(123)
  x <- rnorm(m, 2, sqrt(1.1))
  fg <- g_function(x) / f1_function(x)
  theta_hat[1] <- mean(fg)

  x <- rgamma(m, 3, 2)
  fg <- g_function(x) / f2_function(x)
  theta_hat[2] <- mean(fg)
  
  return(theta_hat)
}

# Main execution
x <- seq(1, 10, 0.1)
plot_density_functions(x)  # Plot the density functions

m <- 10000
se_results <- calculate_standard_errors(m)  # Calculate standard errors
print(se_results)  # Print the corresponding standard errors

theta_hat_results <- monte_carlo_estimate(m)  # Monte Carlo estimate
phi <- (1 / sqrt(2 * pi)) * exp(-1 / 2) + 1 - pnorm(1)  # Real value of the integral
print(round(c(theta_hat_results, phi), 5))  # Print the results



## ----lka----------------------------------------------------------------------
# Clear memory
rm(list = ls())

# Load necessary library
library(ggplot2)

# Function to time the sorting algorithm
time_sorting <- function(n) {
  data <- sample(1:n)  # Generate a random permutation
  start_time <- Sys.time()
  sort(data)  # Sort using built-in quicksort
  end_time <- Sys.time()
  return(as.numeric(end_time - start_time))
}

# Function to perform simulations and calculate average times
average_sorting_time <- function(n, num_simulations) {
  times <- numeric(num_simulations)
  for (j in 1:num_simulations) {
    times[j] <- time_sorting(n)
  }
  return(mean(times))
}

# Function to run the Monte Carlo experiment
run_experiment <- function(n_values, num_simulations) {
  average_times <- numeric(length(n_values))
  t_n <- n_values * log(n_values)
  
  for (i in seq_along(n_values)) {
    average_times[i] <- average_sorting_time(n_values[i], num_simulations)
  }
  
  results <- data.frame(n = n_values, a_n = average_times, t_n = t_n)
  
  # Perform linear regression
  model <- lm(a_n ~ t_n, data = results)
  return(list(results = results, model = model))
}

# Parameters
n_values <- c(10^4, 2 * 10^4, 4 * 10^4, 6 * 10^4, 8 * 10^4)
num_simulations <- 100

# Run the experiment
experiment_results <- run_experiment(n_values, num_simulations)

# Print summary of the regression model
summary(experiment_results$model)

# Plot scatter plot with regression line
ggplot(experiment_results$results, aes(x = t_n, y = a_n)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "Computation Time vs n log(n)",
       x = "t_n = n log(n)",
       y = "Average Computation Time (seconds)") +
  theme_minimal()

## ----ex3----------------------------------------------------------------------
# Clear memory
rm(list = ls())

# Function to generate skewness coefficients
calculate_skewness <- function(n, m) {
  xbar <- numeric(n)
  m3 <- numeric(n)
  m2 <- numeric(n)
  v <- numeric(n)

  for (i in 1:n) {
    x <- rnorm(m, mean = 0, sd = sqrt(6 / n))
    xbar[i] <- mean(x)
    m3[i] <- mean((x - xbar[i])^3)
    m2[i] <- mean((x - xbar[i])^2)
    v[i] <- m3[i] / (m2[i]^1.5)
  }
  
  return(v)
}

# Function to compute quantiles and standard errors
compute_quantiles_and_errors <- function(v, q, n) {
  u <- sort(v)
  gq <- numeric(length(q))
  lsq <- numeric(length(q))
  f <- numeric(length(q))
  ssd <- numeric(length(q))

  for (j in seq_along(q)) {
    gq[j] <- u[n * q[j]]
    lsq[j] <- qnorm(q[j], mean = 0, sd = sqrt(6 / n))
    f[j] <- (1 / sqrt(2 * pi * (6 / n))) * exp(-gq[j]^2 / (2 * (6 / n)))
    ssd[j] <- sqrt(q[j] * (1 - q[j]) / (n * f[j]^2))
  }
  
  return(data.frame(gq, lsq, ssd))
}

table_print <- function(table){
  knitr::kable(table)
}

# Main execution
n <- 1000
m <- 1000
set.seed(12345)

# Calculate skewness coefficients
v <- calculate_skewness(n, m)

# Define quantiles to compute
quantiles <- c(0.025, 0.05, 0.95, 0.975)

# Compute quantiles and standard errors
comparison <- compute_quantiles_and_errors(v, quantiles, n)

# Display results in a table
table_print(comparison)

## ----lk-----------------------------------------------------------------------
# Clear memory
rm(list = ls())

# Load necessary library
library(MASS)

# Function to calculate power for Pearson correlation
calculate_pearson_power <- function(n, m, correlations) {
  power.p <- numeric(length(correlations))
  for (i in seq_along(correlations)) {
    p.p <- replicate(m, {
      x <- mvrnorm(n, c(0, 1), matrix(c(1, correlations[i], correlations[i], 1), 2, 2))
      cortest <- cor.test(x[, 1], x[, 2], alternative = "two.sided", method = "pearson")
      cortest$p.value
    })
    power.p[i] <- mean(p.p <= 0.05)
  }
  return(power.p)
}

# Function to calculate power for Spearman's rank correlation
calculate_spearman_power <- function(n, m, correlations) {
  power.s <- numeric(length(correlations))
  for (i in seq_along(correlations)) {
    p.s <- replicate(m, {
      x <- mvrnorm(n, c(0, 1), matrix(c(1, correlations[i], correlations[i], 1), 2, 2))
      cortest <- cor.test(x[, 1], x[, 2], alternative = "two.sided", method = "spearman")
      cortest$p.value
    })
    power.s[i] <- mean(p.s <= 0.05)
  }
  return(power.s)
}

# Function to calculate power for Kendall's coefficient
calculate_kendall_power <- function(n, m, correlations) {
  power.k <- numeric(length(correlations))
  for (i in seq_along(correlations)) {
    p.k <- replicate(m, {
      x <- mvrnorm(n, c(0, 1), matrix(c(1, correlations[i], correlations[i], 1), 2, 2))
      cortest <- cor.test(x[, 1], x[, 2], alternative = "two.sided", method = "kendall")
      cortest$p.value
    })
    power.k[i] <- mean(p.k <= 0.05)
  }
  return(power.k)
}

# Main execution
n <- 20
m <- 1000
correlations <- c(seq(-0.9, -0.1, 0.1), seq(0.1, 0.9, 0.1))

# Calculate powers
power.p <- calculate_pearson_power(n, m, correlations)
power.s <- calculate_spearman_power(n, m, correlations)
power.k <- calculate_kendall_power(n, m, correlations)

# Combine results into a data frame
power_results <- cbind(power.p, power.s, power.k)
rownames(power_results) <- c(seq(-0.9, -0.1, 0.1), seq(0.1, 0.9, 0.1))
print(power_results)

## ----ppp----------------------------------------------------------------------
n <- 20
m <- 1000
set.seed(123)
p.p1 <- replicate(m, expr = {
  x <- rexp(n, 4)
  y <- runif(n, 0, 1)
  cortest <- cor.test(x, y, alternative = "two.sided", method = "pearson")
  cortest$p.value } )
power.p1 <- mean(  p.p1 <= 0.05 )
power.p1

## ----lll----------------------------------------------------------------------
n <- 20
m <- 1000
set.seed(123)
p.s1 <- replicate(m, expr = {
  x <- rexp(n, 4)
  y <- runif(n, 0, 1)
  cortest <- cor.test(x, y, alternative = "two.sided", method = "spearman")
  cortest$p.value } )
power.s1 <- mean(  p.s1 <= 0.05 )
power.s1

