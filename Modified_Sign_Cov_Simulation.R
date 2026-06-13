# ---------------------------------------------
# Residual-based Bergsma--Dassios type test
# Simulation for different r values
# ---------------------------------------------

rm(list = ls())

library(splines)

set.seed(123)

# Sign covariance kernel (Bergsma-type)
a_kernel <- function(x1, x2, x3, x4){
  sign(abs(x1 - x2) + abs(x3 - x4) -
         abs(x1 - x3) - abs(x2 - x4))
}

# r-th order difference
r_diff <- function(e, r){
  n <- length(e)
  if(r == 0) return(e)
  
  for(k in 1:r){
    e <- diff(c(0, e))  # simple padded difference
  }
  return(e)
}

# B-spline regression estimator
bspline_fit <- function(x, y, df = 6){
  fit <- lm(y ~ bs(x, df = df))
  return(fitted(fit))
}

# Test statistic
tau_star_modified <- function(x, eps_hat){
  n <- length(x)
  comb <- combn(n, 4)
  
  stat <- 0
  m <- ncol(comb)
  
  for(i in 1:m){
    idx <- comb[, i]
    
    x1 <- x[idx[1]]; x2 <- x[idx[2]]
    x3 <- x[idx[3]]; x4 <- x[idx[4]]
    
    e1 <- eps_hat[idx[1]]; e2 <- eps_hat[idx[2]]
    e3 <- eps_hat[idx[3]]; e4 <- eps_hat[idx[4]]
    
    stat <- stat +
      a_kernel(x1,x2,x3,x4) *
      a_kernel(e1,e2,e3,e4)
  }
  
  return(stat / choose(n,4))
}

# One simulation run
one_run <- function(n, delta, r = 1, model = "sin"){
  
  x <- runif(n)
  
  # alternative models
  if(model == "linear"){
    g <- delta * x
  } else if(model == "sin"){
    g <- delta * sin(2*pi*x)
  } else if(model == "quad"){
    g <- delta * (x^2 - x)
  }
  
  eps <- g + rnorm(n, 0, 1)
  y <- 3*x + eps   # m(x)=3x (can be arbitrary)
  
  # estimate regression via B-splines
  mhat <- bspline_fit(x, y, df = 6)
  
  eps_hat <- y - mhat
  
  # r-th order transformation
  eps_r <- r_diff(eps_hat, r)
  
  # align lengths after differencing
  x_r <- x[(length(x) - length(eps_r) + 1):length(x)]
  
  stat <- tau_star_modified(x_r, eps_r)
  
  return(stat)
}

# Monte Carlo experiment
mc_power <- function(n, delta_vals, r_vals, B = 1000, model = "sin"){
  
  results <- array(0, dim = c(length(delta_vals), length(r_vals)))
  
  for(i in seq_along(delta_vals)){
    for(j in seq_along(r_vals)){
      
      delta <- delta_vals[i]
      r <- r_vals[j]
      
      stat_vec <- replicate(B, one_run(n, delta, r, model))
      
      # critical value via permutation-free normal approximation (simplified)
      crit <- quantile(stat_vec, 0.95)
      
      results[i, j] <- mean(stat_vec > crit)
    }
  }
  
  return(results)
}

# ---------------------------------------------
# Example run
# ---------------------------------------------

delta_vals <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
r_vals <- c(1, 2, 4, 6, 8, 10)

power_table <- mc_power(
  n = 200,
  delta_vals = delta_vals,
  r_vals = r_vals,
  B = 1000,
  model = "sin"
)

print(power_table)