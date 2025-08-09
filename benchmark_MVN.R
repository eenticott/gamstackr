# Benchmark script for MVN_weights vs MVN_weights_fast vs MVN_weights_faster
# Compares speed and output equivalence

library(microbenchmark)
source("R/MVN_weights_fast.R")
source("R/MVN_weights_faster.R")
source("R/check_derivatives_mvn.R")

# Simulate data
set.seed(123)
N <- 500
K <- 5
dim_num <- 3
x <- matrix(rnorm(K * dim_num), nrow = dim_num)
eta <- matrix(rnorm(N * dim_num), nrow = N, ncol = dim_num)
theta <- log(runif(dim_num, 0.5, 2))

# Get functions
MVN_orig <- MVN_weights(x)
MVN_fast <- MVN_weights_fast(x)
MVN_faster <- MVN_weights_faster(x)

# Check output equivalence (f_eval)
orig_out <- MVN_orig(eta, theta, deriv = 2)  # Using deriv=2 to test second derivatives
fast_out <- MVN_fast(eta, theta, deriv = 2)
faster_out <- MVN_faster(eta, theta, deriv = 2)

cat('\nComparing f_eval:\n')
cat('Max abs diff orig vs fast:', max(abs(orig_out$f_eval - fast_out$f_eval)), '\n')
cat('Max abs diff orig vs faster:', max(abs(orig_out$f_eval - faster_out$f_eval)), '\n')

# Validate against numerical derivatives
cat('\nValidating derivatives against numDeriv:\n')
cat('\nOriginal implementation:\n')
orig_check <- check_derivatives_mvn(MVN_orig, eta, theta, dim_num)
print(orig_check$max_differences)
cat('\nAll derivatives within tolerance:', all(unlist(orig_check$checks)), '\n')

cat('\nFast implementation:\n')
fast_check <- check_derivatives_mvn(MVN_fast, eta, theta, dim_num)
print(fast_check$max_differences)
cat('\nAll derivatives within tolerance:', all(unlist(fast_check$checks)), '\n')

cat('\nFaster implementation:\n')
faster_check <- check_derivatives_mvn(MVN_faster, eta, theta, dim_num)
print(faster_check$max_differences)
cat('\nAll derivatives within tolerance:', all(unlist(faster_check$checks)), '\n')

# Benchmark speed
bm <- microbenchmark(
  orig = MVN_orig(eta, theta, deriv = 2),
  fast = MVN_fast(eta, theta, deriv = 2),
  faster = MVN_faster(eta, theta, deriv = 2),
  times = 20
)
print(bm)

# Check derivatives
cat('\nComparing first derivatives (first dimension):\n')
cat('Max abs diff orig vs fast (f_eta_eval):', 
    max(abs(orig_out$f_eta_eval[[1]] - fast_out$f_eta_eval[[1]])), '\n')
cat('Max abs diff orig vs faster (f_eta_eval):', 
    max(abs(orig_out$f_eta_eval[[1]] - faster_out$f_eta_eval[[1]])), '\n')

cat('\nComparing first derivatives (theta):\n')
cat('Max abs diff orig vs fast (f_theta_eval):', 
    max(abs(orig_out$f_theta_eval[[1]] - fast_out$f_theta_eval[[1]])), '\n')
cat('Max abs diff orig vs faster (f_theta_eval):', 
    max(abs(orig_out$f_theta_eval[[1]] - faster_out$f_theta_eval[[1]])), '\n')

# Check second derivatives if available
if (!is.null(faster_out$f_eta2_eval)) {
  cat('\nComparing second derivatives (first dimension):\n')
  cat('Max abs diff orig vs faster (f_eta2_eval):', 
      max(abs(orig_out$f_eta2_eval[[1]][[1]] - faster_out$f_eta2_eval[[1]][[1]])), '\n')
  cat('Max abs diff orig vs faster (f_theta2_eval):', 
      max(abs(orig_out$f_theta2_eval[[1]][[1]] - faster_out$f_theta2_eval[[1]][[1]])), '\n')
}

# Test with larger dataset
cat('\nTesting with larger dataset:\n')
N_large <- 2000
K_large <- 10
x_large <- matrix(rnorm(K_large * dim_num), nrow = dim_num)
eta_large <- matrix(rnorm(N_large * dim_num), nrow = N_large, ncol = dim_num)
theta_large <- log(runif(dim_num, 0.5, 2))

MVN_orig_large <- MVN_weights(x_large)
MVN_fast_large <- MVN_weights_fast(x_large)
MVN_faster_large <- MVN_weights_faster(x_large)

bm_large <- microbenchmark(
  orig = MVN_orig_large(eta_large, theta_large, deriv = 2),
  fast = MVN_fast_large(eta_large, theta_large, deriv = 2),
  faster = MVN_faster_large(eta_large, theta_large, deriv = 2),
  times = 20
)
print(bm_large)

# Visualize timings
if (requireNamespace('ggplot2', quietly = TRUE)) {
  library(ggplot2)
  print(autoplot(bm))
  print(autoplot(bm_large))
}
