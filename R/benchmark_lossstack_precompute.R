# Benchmarking precomputing and reusing intermediate results (point 2)

library(microbenchmark)
library(Rfast)

# Simulate data
set.seed(123)
n <- 50000
K <- 5
neta <- 3

preds <- matrix(rnorm(n * K), n, K)
list_of_X <- lapply(1:neta, function(i) matrix(rnorm(n * 4), n, 4))
W_f_eta_eval <- lapply(1:neta, function(i) matrix(rnorm(n * K), n, K))
loss_eval_l1 <- rnorm(n)
v <- 1.5

# Original: recompute Rfast::rowsums(W_f_eta_eval[[i]] * preds) in each loop
grad_loop <- function() {
  leta <- list()
  for (i in seq_len(neta)) {
    # rowsums is recomputed for each i
    leta[[i]] <- -v * Rfast::colsums(loss_eval_l1 * Rfast::rowsums(W_f_eta_eval[[i]] * preds) * list_of_X[[i]])
  }
  leta
}

# Optimized: precompute rowsums for all i, reuse in loop
grad_precompute <- function() {
  leta <- list()
  # Precompute all rowsums for each i
  rowsums_list <- lapply(W_f_eta_eval, function(mat) Rfast::rowsums(mat * preds))
  for (i in seq_len(neta)) {
    leta[[i]] <- -v * Rfast::colsums(loss_eval_l1 * rowsums_list[[i]] * list_of_X[[i]])
  }
  leta
}

# Benchmark
bm <- microbenchmark(
  loop = grad_loop(),
  precompute = grad_precompute(),
  times = 100L
)

print(bm)
