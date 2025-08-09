#' Exponential link function and its derivatives
#'
#' This function computes the exponential link function and its first and second derivatives.
#'
#' @param x Numeric input value
#'
#' @return A list containing the function value (f), first derivative (f1), and second derivative (f2)
#'
log_link <- function(x) {
  out <- list()
  out$f <- exp(x)
  out$f1 <- exp(x)
  out$f2 <- exp(x)
  return(out)
}

#' Identity link function and its derivatives
#'
#' This function computes the identity link function and its first and second derivatives.
#'
#' @param x Numeric input value
#'
#' @return A list containing the function value (f), first derivative (f1), and second derivative (f2)
#'
id_link <- function(x) {
  out <- list()
  out$f <- x
  out$f1 <- 1
  out$f2 <- 0
  return(out)
}

#' Apply link function to transform parameters
#'
#' This function applies a link function to transform parameters.
#'
#' @param f Function to apply
#' @param link_func String specifying the link function to use (default: "log")
#' @param eta Linear predictor values
#' @param tau Parameter values on the link scale
#'
#' @return Transformed parameter values
#'
link_theta <- function(f, link_func = "log", eta, tau) {
  if (link_func == "log") {
    theta_eval <- f(eta, exp(tau))
    tau_eval <- theta_eval * exp(tau)
    return(tau_eval)
  }
}
