# Benchmark script for MVN_weights vs MVN_weights_fast vs MVN_weights_faster
# Compares speed and output equivalence, with numerical derivative validation

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

# Check analytical derivatives against numerical approximations
cat('\nValidating derivatives against numerical approximations:\n')

# Test original implementation
cat('\nOriginal implementation:\n')
orig_check <- check_derivatives_mvn(MVN_orig, eta[1:2,], theta, dim_num, k=1)
cat('First derivatives correct:', orig_check$grad_correct, '\n')
cat('Second derivatives correct:', 
    all(c(orig_check$hess_etaeta_correct, 
          orig_check$hess_etatheta_correct,
          orig_check$hess_thetatheta_correct)), '\n')
cat('Max differences:\n')
print(orig_check[grep("max_diff", names(orig_check))])

# Test fast implementation
cat('\nFast implementation:\n')
tryCatch({
    fast_check <- check_derivatives_mvn(MVN_fast, eta[1:2,], theta, dim_num, k=1)
    cat('First derivatives correct:', fast_check$grad_correct, '\n')
    if (!is.null(fast_check$hess_etaeta_correct)) {
        cat('Second derivatives correct:', 
            all(c(fast_check$hess_etaeta_correct, 
                fast_check$hess_etatheta_correct,
                fast_check$hess_thetatheta_correct)), '\n')
        cat('Max differences:\n')
        print(fast_check[grep("max_diff", names(fast_check))])
    } else {
        cat('Note: Second derivatives not available for fast implementation\n')
    }
}, error = function(e) {
    cat('Error in derivative check:', e$message, '\n')
})

# Test faster implementation
cat('\nFaster implementation:\n')
tryCatch({
    faster_check <- check_derivatives_mvn(MVN_faster, eta[1:2,], theta, dim_num, k=1)
    cat('First derivatives correct:', faster_check$grad_correct, '\n')
    if (!is.null(faster_check$hess_etaeta_correct)) {
        cat('Second derivatives correct:', 
            all(c(faster_check$hess_etaeta_correct, 
                faster_check$hess_etatheta_correct,
                faster_check$hess_thetatheta_correct)), '\n')
        cat('Max differences:\n')
        print(faster_check[grep("max_diff", names(faster_check))])
    } else {
        cat('Note: Second derivatives not available for faster implementation\n')
    }
}, error = function(e) {
    cat('Error in derivative check:', e$message, '\n')
})

# Check output equivalence
cat('\nComparing implementations against original:\n')
orig_out <- MVN_orig(eta, theta, deriv = 2)
fast_out <- MVN_fast(eta, theta, deriv = 2)
faster_out <- MVN_faster(eta, theta, deriv = 2)

cat('Max abs diff orig vs fast (f_eval):', max(abs(orig_out$f_eval - fast_out$f_eval)), '\n')
cat('Max abs diff orig vs faster (f_eval):', max(abs(orig_out$f_eval - faster_out$f_eval)), '\n')

# Benchmark speed with small dataset
cat('\nBenchmarking with small dataset (N=500, K=5):\n')
bm <- microbenchmark(
  orig = MVN_orig(eta, theta, deriv = 2),
  fast = MVN_fast(eta, theta, deriv = 2),
  faster = MVN_faster(eta, theta, deriv = 2),
  times = 20
)
print(bm)

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

cat('\nBenchmarking with large dataset (N=2000, K=10):\n')
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
