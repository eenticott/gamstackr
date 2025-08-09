# Test and benchmark for Latent vs Latent_fast
library(microbenchmark)
source("R/Latent_fast.R")

# Assume Latent and Latent_fast are both loaded in the environment
K <- 5
q <- 2
N <- 100
set.seed(1)
eta <- matrix(rnorm(N*q), N, q)
theta <- rnorm((K-1) + (K-1-q)*q)
Latent_orig <- Latent(K, q)
Latent_opt <- Latent_fast(K, q)

# Check output equivalence for value and gradient
out_orig <- Latent_orig(eta, theta, deriv=1)
out_opt  <- Latent_opt(eta, theta, deriv=1)

identical_eval <- all.equal(out_orig$f_eval, out_opt$f_eval, tolerance=1e-8)
identical_grad <- all.equal(
  do.call(cbind, out_orig$f_eta_eval),
  do.call(cbind, out_opt$f_eta_eval),
  tolerance=1e-8
)

cat("f_eval identical:", identical_eval, "\n")
cat("f_eta_eval identical:", identical_grad, "\n")

# Benchmark
bm <- microbenchmark(
  orig = Latent_orig(eta, theta, deriv=1),
  opt  = Latent_opt(eta, theta, deriv=1),
  times = 50
)
print(bm)
