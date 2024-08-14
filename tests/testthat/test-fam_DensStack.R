library(numDeriv)
source("R/fam_DensStack.R")

check_dens_deriv <- function(weight_func, p, N, tol = 1e-8) {
  neta <- attr(weight_func, "neta")
  ntheta <- attr(weight_func, "ntheta")
  K <- attr(weight_func, "num_weights")

  dens <- matrix(rnorm(N*K), nrow = N, ncol = K)

  # Input check for p
  if (length(p) > 1) {
    if (length(p) != neta) {
      stop("Length of p must match number of eta.")
    }
  } else {
    p <- rep(p, neta)
  }
  nbeta <- sum(p)

  list_of_beta <- lapply(p, function(i) rnorm(i))
  list_of_X <- lapply(p, function(i) matrix(rnorm(i*N), nrow = N, ncol = i))
  list_of_eta <- get_list_of_eta(list_of_X, list_of_beta)

  theta <- rnorm(ntheta)

  tmp_fun <- function(pars) {
    list_of_beta <- split(pars[1:nbeta], rep(1:neta, p))
    theta <- pars[(nbeta+1):(nbeta + ntheta)]
    get_ll_dens_derivs(list_of_beta, list_of_X, theta, weight_func, (dens))$l
  }

  out <- get_ll_dens_derivs(list_of_beta, list_of_X, theta, weight_func, (dens), 2)

  calc_grad <- out$lb
  calc_hess <- out$lbb
  num_grad <- grad(tmp_fun, c(unlist(list_of_beta), theta))
  num_hess <- hessian(tmp_fun, c(unlist(list_of_beta), theta))

  # Check gradient
  grad_correct <- mean(abs(num_grad-calc_grad)) < tol

  # Check hess
  hess_correct <- mean(abs(num_hess - calc_hess)) < tol

  return(all(grad_correct, hess_correct))
}

test_that("Dens derivs work", {
  expect_true(check_dens_deriv(ordinal(5), 3, 100))
  expect_true(check_dens_deriv(nested(multinomial(2), list(ordinal(5), ordinal(3))), 3, 100))
  expect_true(check_dens_deriv(nested(ordinal(3), list(ordinal(5), ordinal(3), ordinal(4))), 3, 100))
})




