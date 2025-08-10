library(numDeriv)
#source("R/fam_LossStack.R")
set.seed(1)
# Test derivatives for log_dens stacking

check_loss_deriv <- function(weight_func, loss, p, N, tol = 1e-5) {
  # Get relevant attributes
  neta <- attr(weight_func, "neta")
  ntheta <- attr(weight_func, "ntheta")
  K <- attr(weight_func, "num_weights")

  # Initialise prediction matrix
  preds <- matrix(rnorm(N*K), nrow = N, ncol = K)
  # Input check for p
  # If its a vector it should match eta length
  if (length(p) > 1) {
    if (length(p) != neta) {
      stop("Length of p must match number of eta.")
    }
  } else {
    p <- rep(p, neta)
  }

  # total number of betas is sum of p
  nbeta <- sum(p)

  # Create random test variables to calculate derivs on
  list_of_beta <- lapply(p, function(i) rnorm(i))
  list_of_X <- lapply(p, function(i) matrix(rnorm(i*N), nrow = N, ncol = i))
  list_of_eta <- get_list_of_eta(list_of_X, list_of_beta)

  theta <- rnorm(ntheta)

  y <- rnorm(N)

  # Create function to be used inside num deriv
  tmp_fun <- function(pars) {
    list_of_beta <- split(pars[1:nbeta], rep(1:neta, p))
    theta <- pars[(nbeta+1):(nbeta + ntheta)]
    logv <- pars[(nbeta+ntheta+1):(nbeta+ntheta+1)]
    get_loss_derivs(list_of_beta, list_of_X, theta, logv, loss, weight_func, preds, y, deriv = 0)$l
  }

  # Calculate derivatives
  out <- get_loss_derivs(list_of_beta, list_of_X, theta, logv = 1, loss, weight_func, preds, y, deriv = 2)

  calc_grad <- out$lb
  calc_hess <- out$lbb

  num_grad <- pracma::grad(tmp_fun, c(unlist(list_of_beta), theta, 1))
  num_hess <- pracma::hessian(tmp_fun, c(unlist(list_of_beta), theta, 1))

  # Check gradient
  grad_correct <- mean(abs(num_grad-calc_grad)) < tol

  # Check hess
  hess_correct <- mean(abs(num_hess - calc_hess)) < tol
  print(abs(num_hess - calc_hess) < tol)
  # Check both gradient and hessian match numerical derivatives
  return(all(grad_correct, hess_correct))
}

# Test different combinations of inner functions
test_that("Dens derivs work", {
  expect_true(check_loss_deriv(ordinal(5), square_loss, 3, 100))
  expect_true(check_loss_deriv(nested(multinomial(2), list(ordinal(5), ordinal(3))), square_loss, 3, 100))
  expect_true(check_loss_deriv(nested(ordinal(3), list(ordinal(5), ordinal(3), ordinal(4))), square_loss, 3, 100))
  b <- pinball_loss(0.5)
  expect_true(check_loss_deriv(ordinal(5), pinball_loss(0.5), 3, 100, tol = 1e-5))
  expect_true(check_loss_deriv(nested(multinomial(2), list(ordinal(5), ordinal(3))), pinball_loss(0.5), 3, 100,tol=1e-5))

})
