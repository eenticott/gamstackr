##' Log link function and derivatives
#'
#' Computes the log link and its first and second derivatives.
#'
#' @param x Numeric vector or matrix.
#'
#' @return List with elements: f (link), f1 (first derivative), f2 (second derivative).
#' @export
log_link <- function(x) {
  out <- list()
  out$f <- exp(x)
  out$f1 <- exp(x)
  out$f2 <- exp(x)
  return(out)
}

##' Identity link function and derivatives
#'
#' Computes the identity link and its first and second derivatives.
#'
#' @param x Numeric vector or matrix.
#'
#' @return List with elements: f (link), f1 (first derivative), f2 (second derivative).
#' @export
id_link <- function(x) {
  out <- list()
  out$f <- x
  out$f1 <- 1
  out$f2 <- 0
  return(out)
}

##' Transform theta using a link function
#'
#' Applies a link function to parameters for model transformation.
#'
#' @param f Function. Link function to apply.
#' @param link_func Character. Name of the link function (default: "log").
#' @param eta Numeric. Linear predictor.
#' @param tau Numeric. Additional parameter(s).
#'
#' @return Transformed parameter value(s).
#' @export
link_theta <- function(f, link_func = "log", eta, tau) {
  if (link_func == "log") {
    theta_eval <- f(eta, exp(tau))
    tau_eval <- theta_eval * exp(tau)
    return(tau_eval)
  }
  return(NULL)
}
