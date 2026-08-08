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
# ============================================================
# POWER AND TYPE-I ERROR STUDY FOR EXISTING MEASURES
# Kendall's tau, Bergsma-Dassios tau-star and dCov
# ============================================================


# ------------------------------------------------------------
# Kendall's tau statistic
# ------------------------------------------------------------

kendall_stat <- function(x, eps_hat){
  
  return(
    abs(
      cor(
        x,
        eps_hat,
        method = "kendall"
      )
    )
  )
}


# ------------------------------------------------------------
# Bergsma-Dassios tau-star statistic
# ------------------------------------------------------------

tau_star_stat <- function(x, eps_hat){
  
  n <- length(x)
  
  comb <- combn(n, 4)
  
  stat <- 0
  
  m <- ncol(comb)
  
  
  for(i in 1:m){
    
    idx <- comb[, i]
    
    x1 <- x[idx[1]]
    x2 <- x[idx[2]]
    x3 <- x[idx[3]]
    x4 <- x[idx[4]]
    
    e1 <- eps_hat[idx[1]]
    e2 <- eps_hat[idx[2]]
    e3 <- eps_hat[idx[3]]
    e4 <- eps_hat[idx[4]]
    
    
    stat <- stat +
      a_kernel(x1, x2, x3, x4) *
      a_kernel(e1, e2, e3, e4)
    
  }
  
  
  return(
    stat / choose(n, 4)
  )
}


# ------------------------------------------------------------
# Distance covariance statistic
# ------------------------------------------------------------

dcov_stat <- function(x, eps_hat){
  
  n <- length(x)
  
  
  # Distance matrices
  
  A <- abs(
    outer(x, x, "-")
  )
  
  B <- abs(
    outer(eps_hat, eps_hat, "-")
  )
  
  
  # Double-centering
  
  A_centered <- matrix(0, n, n)
  
  B_centered <- matrix(0, n, n)
  
  
  for(i in 1:n){
    
    for(j in 1:n){
      
      A_centered[i, j] <-
        A[i, j] -
        mean(A[i, ]) -
        mean(A[, j]) +
        mean(A)
      
      B_centered[i, j] <-
        B[i, j] -
        mean(B[i, ]) -
        mean(B[, j]) +
        mean(B)
      
    }
    
  }
  
  
  dcov_squared <-
    sum(
      A_centered * B_centered
    ) / n^2
  
  
  dcov_squared <-
    max(
      dcov_squared,
      0
    )
  
  
  return(
    sqrt(dcov_squared)
  )
}


# ------------------------------------------------------------
# Generate residuals for one data set
# ------------------------------------------------------------

generate_residuals <- function(
    n,
    delta,
    model = "sin"
){
  
  x <- runif(n)
  
  
  # Alternative models
  
  if(model == "linear"){
    
    g <- delta * x
    
  } else if(model == "sin"){
    
    g <- delta * sin(2 * pi * x)
    
  } else if(model == "quad"){
    
    g <- delta * (x^2 - x)
    
  }
  
  
  # Error term
  
  eps <- g + rnorm(
    n,
    mean = 0,
    sd = 1
  )
  
  
  # Response
  
  y <- 3 * x + eps
  
  
  # B-spline estimate
  
  mhat <- bspline_fit(
    x,
    y,
    df = 6
  )
  
  
  # Estimated residuals
  
  eps_hat <- y - mhat
  
  
  return(
    list(
      x = x,
      eps_hat = eps_hat
    )
  )
}


# ------------------------------------------------------------
# Calculate all three existing statistics
# ------------------------------------------------------------

existing_statistics <- function(
    n,
    delta,
    model = "sin"
){
  
  data <- generate_residuals(
    n = n,
    delta = delta,
    model = model
  )
  
  
  x <- data$x
  
  eps_hat <- data$eps_hat
  
  
  # Kendall's tau
  
  tau_value <- kendall_stat(
    x,
    eps_hat
  )
  
  
  # Bergsma-Dassios tau-star
  
  tau_star_value <- tau_star_stat(
    x,
    eps_hat
  )
  
  
  # Distance covariance
  
  dcov_value <- dcov_stat(
    x,
    eps_hat
  )
  
  
  return(
    c(
      Kendall_tau = tau_value,
      tau_star = tau_star_value,
      dCov = dcov_value
    )
  )
}


# ------------------------------------------------------------
# Generate null distributions
# ------------------------------------------------------------

null_distribution_existing <- function(
    n,
    model = "sin",
    B = 1000
){
  
  null_values <- matrix(
    NA,
    nrow = B,
    ncol = 3
  )
  
  
  colnames(null_values) <- c(
    "Kendall_tau",
    "tau_star",
    "dCov"
  )
  
  
  for(b in 1:B){
    
    null_values[b, ] <-
      existing_statistics(
        n = n,
        delta = 0,
        model = model
      )
    
  }
  
  
  return(null_values)
}


# ------------------------------------------------------------
# Calculate critical values from the null distribution
# ------------------------------------------------------------

existing_critical_values <- function(
    null_values,
    alpha = 0.05
){
  
  critical_values <- numeric(3)
  
  
  critical_values[1] <-
    quantile(
      null_values[, "Kendall_tau"],
      1 - alpha,
      na.rm = TRUE
    )
  
  
  critical_values[2] <-
    quantile(
      null_values[, "tau_star"],
      1 - alpha,
      na.rm = TRUE
    )
  
  
  critical_values[3] <-
    quantile(
      null_values[, "dCov"],
      1 - alpha,
      na.rm = TRUE
    )
  
  
  names(critical_values) <- c(
    "Kendall_tau",
    "tau_star",
    "dCov"
  )
  
  
  return(
    critical_values
  )
}


# ------------------------------------------------------------
# Type-I error study
# ------------------------------------------------------------

existing_type1_error <- function(
    null_values,
    critical_values
){
  
  type1 <- numeric(3)
  
  
  type1[1] <-
    mean(
      null_values[, "Kendall_tau"] >
        critical_values["Kendall_tau"],
      na.rm = TRUE
    )
  
  
  type1[2] <-
    mean(
      null_values[, "tau_star"] >
        critical_values["tau_star"],
      na.rm = TRUE
    )
  
  
  type1[3] <-
    mean(
      null_values[, "dCov"] >
        critical_values["dCov"],
      na.rm = TRUE
    )
  
  
  names(type1) <- c(
    "Kendall_tau",
    "tau_star",
    "dCov"
  )
  
  
  return(
    type1
  )
}


# ------------------------------------------------------------
# Power study for one value of delta
# ------------------------------------------------------------

existing_power_one_delta <- function(
    n,
    delta,
    model,
    critical_values,
    B = 1000
){
  
  alternative_values <- matrix(
    NA,
    nrow = B,
    ncol = 3
  )
  
  
  colnames(alternative_values) <- c(
    "Kendall_tau",
    "tau_star",
    "dCov"
  )
  
  
  for(b in 1:B){
    
    alternative_values[b, ] <-
      existing_statistics(
        n = n,
        delta = delta,
        model = model
      )
    
  }
  
  
  power <- numeric(3)
  
  
  power[1] <-
    mean(
      alternative_values[, "Kendall_tau"] >
        critical_values["Kendall_tau"],
      na.rm = TRUE
    )
  
  
  power[2] <-
    mean(
      alternative_values[, "tau_star"] >
        critical_values["tau_star"],
      na.rm = TRUE
    )
  
  
  power[3] <-
    mean(
      alternative_values[, "dCov"] >
        critical_values["dCov"],
      na.rm = TRUE
    )
  
  
  names(power) <- c(
    "Kendall_tau",
    "tau_star",
    "dCov"
  )
  
  
  return(
    power
  )
}


# ============================================================
# COMPLETE EXISTING-METHODS STUDY
# ============================================================

existing_methods_study <- function(
    n = 200,
    delta_vals = c(
      0,
      0.2,
      0.4,
      0.6,
      0.8,
      1.0
    ),
    B = 1000,
    model = "sin",
    alpha = 0.05
){
  
  cat("\n")
  cat("=============================================\n")
  cat("EXISTING METHODS: TYPE-I ERROR AND POWER\n")
  cat("=============================================\n")
  cat("n     =", n, "\n")
  cat("model =", model, "\n")
  cat("B     =", B, "\n")
  cat("alpha =", alpha, "\n")
  cat("=============================================\n")
  
  
  # ----------------------------------------------------------
  # Step 1: Null simulation
  # ----------------------------------------------------------
  
  cat("\nGenerating null distribution...\n")
  
  
  null_values <-
    null_distribution_existing(
      n = n,
      model = model,
      B = B
    )
  
  
  # ----------------------------------------------------------
  # Step 2: Critical values
  # ----------------------------------------------------------
  
  critical_values <-
    existing_critical_values(
      null_values,
      alpha = alpha
    )
  
  
  cat("\n")
  cat("Critical values:\n")
  print(critical_values)
  
  
  # ----------------------------------------------------------
  # Step 3: Type-I error
  # ----------------------------------------------------------
  
  type1 <-
    existing_type1_error(
      null_values,
      critical_values
    )
  
  
  cat("\n")
  cat("Empirical Type-I error:\n")
  print(type1)
  
  
  # ----------------------------------------------------------
  # Step 4: Power
  # ----------------------------------------------------------
  
  power_table <- matrix(
    NA,
    nrow = length(delta_vals),
    ncol = 4
  )
  
  
  colnames(power_table) <- c(
    "delta",
    "Kendall_tau",
    "tau_star",
    "dCov"
  )
  
  
  power_table[, "delta"] <-
    delta_vals
  
  
  for(i in seq_along(delta_vals)){
    
    delta <- delta_vals[i]
    
    
    cat(
      "\nCalculating power for delta =",
      delta,
      "\n"
    )
    
    
    if(delta == 0){
      
      power_table[i, 2:4] <-
        type1
      
    } else {
      
      power_table[i, 2:4] <-
        existing_power_one_delta(
          n = n,
          delta = delta,
          model = model,
          critical_values = critical_values,
          B = B
        )
      
    }
    
  }
  
  
  return(
    list(
      critical_values = critical_values,
      type1_error = type1,
      power = power_table
    )
  )
}


# ============================================================
# RUN EXISTING-METHODS STUDY
# ============================================================

existing_results <- existing_methods_study(
  n = 200,
  delta_vals = delta_vals,
  B = 1000,
  model = "sin",
  alpha = 0.05
)


# ============================================================
# DISPLAY TYPE-I ERROR
# ============================================================

cat("\n")
cat("=============================================\n")
cat("TYPE-I ERROR: EXISTING METHODS\n")
cat("=============================================\n")

print(
  existing_results$type1_error
)


# ============================================================
# DISPLAY POWER TABLE
# ============================================================

cat("\n")
cat("=============================================\n")
cat("EMPIRICAL POWER: EXISTING METHODS\n")
cat("=============================================\n")

print(
  existing_results$power
)


# ============================================================
# SAVE RESULTS
# ============================================================

write.csv(
  existing_results$power,
  "existing_methods_power_n200_sin.csv",
  row.names = FALSE
)


write.csv(
  as.data.frame(
    t(existing_results$type1_error)
  ),
  "existing_methods_type1_n200_sin.csv",
  row.names = FALSE
)


# ============================================================
# END
# ============================================================