library(numDeriv)

# Test whether inner function derivatives are correct

check_deriv <- function(weight_func, k, deriv = 1, tol = 1e-8, print_derivs = FALSE) {
  # weight_func - Inner function to test
  # k - which weight index to check
  # deriv - check first or first and second derivs
  # tol - Max error rate between true and numerical
  # print_derivs - Print out deriv values when calculated (for debugging)
  neta <- attr(weight_func, "neta")
  ntheta <- attr(weight_func, "ntheta")

  # Randomise eta and theta to evaluate derivs at
  eta <- matrix(rnorm(neta), nrow = 1)
  theta <- rnorm(ntheta)

  # Calculate weights and derivatives
  out <- weight_func(eta, theta, deriv = deriv)

  tmp_fun <- function(pars) {
    if (neta > 0) {
      eta = matrix(pars[1:neta], nrow = 1)
    }
    if (ntheta > 0) {
      theta = pars[(neta+1):(ntheta+neta)]
    }
    weight_func(eta, theta, deriv = 0)$f_eval[k]
  }

  num_grad <- grad(tmp_fun, c(eta[1,], theta))

  eta_grad <- NULL
  theta_grad <- NULL

  if (neta > 0) {
    eta_grad <- do.call("c", lapply(out$f_eta_eval, "[[", k))
    if (print_derivs) {
      print(paste("Numerical eta grad",paste(num_grad[1:neta], collapse = " ")))
      print(paste("Calculated eta grad",paste(eta_grad, collapse = " ")))
    }
  }
  if (ntheta > 0) {
    theta_grad <- do.call("c", lapply(out$f_theta_eval, "[[", k))
    if (print_derivs) {
      print(paste("Numerical theta grad", paste(num_grad[(neta+1):(ntheta+1+neta)], collapse = " ")))
      print(paste("Calculated theta grad",paste(theta_grad, collapse = " ")))

    }
  }
  calc_grad <- c(eta_grad, theta_grad)

  grad_correct <- mean(abs(calc_grad - num_grad)) < tol

  if (deriv > 1) {
    num_hess <- hessian(tmp_fun, c(eta[1,], theta))
    hess_etaeta_correct <- NULL
    hess_etatheta_correct <- NULL
    hess_thetatheta_correct <- NULL

    if (neta > 0) {
      calc_etaeta <- matrix(unlist(lapply(out$f_eta2_eval, function(i) {lapply(i, "[[", k)})), neta, neta)
      num_etaeta <- num_hess[1:neta, 1:neta]
      hess_etaeta_correct <- mean(abs(calc_etaeta - num_etaeta)) < tol
      if (print_derivs) {
        print(cat("Numerical etaeta block", num_etaeta))
        print(cat("Calculated etaeta block", calc_etaeta))
      }
      if (ntheta > 0) {
        calc_etatheta <- matrix(unlist(lapply(out$f_eta_theta_eval, function(i) {lapply(i, "[[", k)})), neta, ntheta, byrow = TRUE)
        num_etatheta <- num_hess[1:neta, (neta+1):(neta+ntheta)]
        hess_etatheta_correct <- mean(abs(calc_etatheta - num_etatheta)) < tol
        if (print_derivs) {
          print(cat("Numerical etatheta block", num_etatheta))
          print(cat("Calculated etatheta block", calc_etatheta))
        }
      }
    }
    if (ntheta > 0) {
      calc_thetatheta <- matrix(unlist(lapply(out$f_theta2_eval, function(i) {lapply(i, "[[", k)})), ntheta, ntheta)
      num_thetatheta <- num_hess[(neta+1):(neta+ntheta), (neta+1):(neta+ntheta)]
      hess_thetatheta_correct <- mean(abs(calc_thetatheta - num_thetatheta)) < tol
      if (print_derivs) {
        print(cat("Numerical thetatheta block", num_thetatheta))
        print(cat("Calculated thetatheta block", calc_thetatheta))
      }
    }
    return(all(grad_correct,
               hess_etaeta_correct,
               hess_etatheta_correct,
               hess_thetatheta_correct))
  } else {
    return(grad_correct)
  }
}

# Test each of the inner functions
test_that("ordinal derivs", {
  wf <- ordinal(5)
  expect_true(check_deriv(wf, 1, deriv = 2))
  expect_true(check_deriv(wf, 3, deriv = 2))
  expect_true(check_deriv(wf, 5, deriv = 2))

  wf <- ordinal(9)
  expect_true(check_deriv(wf, 1, deriv = 2))
  expect_true(check_deriv(wf, 3, deriv = 2))
  expect_true(check_deriv(wf, 5, deriv = 2))

  wf <- nested(multinomial(3), list(ordinal(3), ordinal(5), ordinal(3)))
  expect_true(check_deriv(wf, 2, deriv = 2))

  wf <- Latent(5, 3)
  expect_true(check_deriv(wf, 2, deriv = 2))
})

test_that("mvn derivs", {
  x <- matrix(rnorm(15), nrow = 3, ncol = 5)
  wf <- MVN_weights_cpp(x)
  expect_true(check_deriv(wf, 1, deriv = 2))

  x <- matrix(rnorm(5), nrow = 1, ncol = 5)
  wf <- MVN_weights_cpp(x)
  expect_true(check_deriv(wf, 2, deriv = 2))

  x <- matrix(rnorm(15), nrow = 3, ncol = 5)
  wf <- MVN_weights_cpp(x)
  expect_true(check_deriv(wf, 5, deriv = 2))
})

test_that("mvn derivs", {
  x <- matrix(rnorm(15), nrow = 3, ncol = 5)
  wf <- MVN_weights(x)
  expect_true(check_deriv(wf, 1, deriv = 2))

  x <- matrix(rnorm(5), nrow = 1, ncol = 5)
  wf <- MVN_weights(x)
  expect_true(check_deriv(wf, 2, deriv = 2))

  x <- matrix(rnorm(15), nrow = 3, ncol = 5)
  wf <- MVN_weights(x)
  expect_true(check_deriv(wf, 5, deriv = 2))
})

