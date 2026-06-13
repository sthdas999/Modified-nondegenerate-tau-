############################################################
# Modified Bergsma–Dassios tau* with residuals
# Application: mtcars dataset
############################################################

# Load dataset
data(mtcars)

# Response and covariate
Y <- mtcars$mpg
X <- mtcars$disp
n <- length(Y)

############################################################
# Step 1: Penalized B-spline smoothing (via smooth.spline)
############################################################

df <- 6  # degrees of freedom (order parameter)

fit <- smooth.spline(x = X, y = Y, df = df)

m_hat <- predict(fit, X)$y

# Residuals
eps_hat <- Y - m_hat

############################################################
# Step 2: r-th order difference of residuals
############################################################

diff_r <- function(eps, r) {
  n <- length(eps)
  eps_r <- rep(NA, n)
  
  for (t in 1:n) {
    if (t + r - 1 <= n) {
      s <- 0
      for (j in 0:r) {
        if ((t + r - 1 - j) >= 1 && (t + r - 1 - j) <= n) {
          s <- s + (-1)^j * choose(r, j) * eps[t + r - 1 - j]
        }
      }
      eps_r[t] <- s
    }
  }
  
  return(eps_r)
}

############################################################
# Step 3: Kernel function a(p,q,r,s)
############################################################

a_kernel <- function(p, q, r, s) {
  sign(abs(p - q) + abs(r - s) - abs(p - r) - abs(q - s))
}

############################################################
# Step 4: Modified tau* statistic
############################################################

tau_star_modified <- function(X, eps) {
  
  n <- length(X)
  stat <- 0
  count <- 0
  
  for (i1 in 1:(n-3)) {
    for (i2 in (i1+1):(n-2)) {
      for (i3 in (i2+1):(n-1)) {
        for (i4 in (i3+1):n) {
          
          kx <- a_kernel(X[i1], X[i2], X[i3], X[i4])
          ke <- a_kernel(eps[i1], eps[i2], eps[i3], eps[i4])
          
          stat <- stat + kx * ke
          count <- count + 1
        }
      }
    }
  }
  
  return(stat / count)
}

############################################################
# Step 5: Compute statistic for different r values
############################################################

r_values <- c(1, 2, 4, 6, 8, 10)

results <- numeric(length(r_values))

for (k in seq_along(r_values)) {
  
  r <- r_values[k]
  
  eps_r <- diff_r(eps_hat, r)
  
  valid_idx <- which(!is.na(eps_r))
  
  results[k] <- tau_star_modified(X[valid_idx], eps_r[valid_idx])
  
}

############################################################
# Step 6: Output results
############################################################

names(results) <- paste0("r=", r_values)
print(results)