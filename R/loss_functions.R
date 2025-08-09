# square_loss <- function(y, param, deriv = 0, ...) {
#   # For square loss need derivatives of gaussian distribution
#   # derivatives are in terms of t where exp(t) = sig2
#   if (is.list(param) ) param <- do.call("cbind", param)
#   if (is.vector(param)) param <- matrix(param, nrow = 1)
#   if (ncol(param) != 2) stop("Wrong number of parameters provided")
#
#   p <- ncol( param )
#   mu <- param[ , 1, drop = TRUE]
#   sig2 <- exp(param[ , 2, drop = TRUE]) # sigma2 (force positive)
#   tau2 <- 1/sig2
#   tau <- sqrt(tau2)
#
#   n <- length(y)
#
#   if (length(mu) == 1) {
#     mu <- rep(mu, n)
#     tau <- rep(tau, n)
#     tau2 <- rep(tau2, n)
#   }
#
#   ymu <- y - mu
#   ymu2 <- ymu ^ 2
#   d0 <- - .5 * log(2 * pi) + log(tau) - .5 * tau2 * ymu2
#   out <- list()
#   out$d0 <- d0
#
#   if( deriv > 0 )
#   {
#     d1 <- tau2 * ymu
#     d2 <- (1/tau - tau*ymu2) * -1/2 * sig2^(-3/2)
#     out[["d1"]] <- list(d1, d2*sig2)
#
#     if( deriv > 1 ){
#       d11 <- - tau2
#       d12 <- (2 * d1 / tau) * -1/2 * sig2^(-3/2)
#       d22 <- ((- ymu2 - 1/tau2) * 1/4 * sig2^(-3)) + ((1/tau - tau*ymu2) * 3/4 * sig2^(-5/2))
#       out[["d2"]] <- list(d11, d12*sig2, sig2 * d2 + sig2^2 * d22)
#
#       # Haven't changed later ones to sig2
#       if( deriv > 2 ){
#         zeros <- rep(0, n)
#         d111 <- zeros
#         d112 <- - 2 * tau
#         d122 <- 2 * ymu
#         d222 <- 2 / tau^3
#         out[["d3"]] <- list(d111, d112, d122, d222)
#
#         if( deriv > 3){
#           d1111 <- d1112 <- d1222 <- zeros
#           d1122 <- rep(-2, n)
#           d2222 <- -6 / tau2^2
#           out[["d3"]] <- list(d1111, d1112, d1122, d1222, d2222)
#         }
#       }
#     }
#   }
#
#   return( out )
# }
#
# test_f <- function(pars) {
#   square_loss(y = c(1), param = list(pars[1], pars[2]), deriv = 0)$d0
# }
#
# square_loss(y = c(1), param = list(2, -2), deriv = 2)
#
# numDeriv::grad(test_f, c(2,-2))
# numDeriv::hessian(test_f, c(2,-2))

##' Squared loss function
#'
#' Computes squared loss and its derivatives for regression tasks.
#'
#' @param y Numeric vector. True values.
#' @param yhat Numeric vector. Predicted values.
#' @param deriv Integer. Order of derivative to compute (0 = value, 1 = gradient, 2 = Hessian).
#'
#' @return List with elements: loss (numeric vector), l1 (first derivative), l2 (second derivative), C (function for normalization constant and its derivatives).
#' @export
square_loss <- function(y, yhat, deriv = 0) {
  # Calculate loss
  loss <- (y-yhat)^2
  l1 <- NULL
  l2 <- NULL
  if (deriv > 0) {
    # Calculate first derivative w.r.t. yhat
    l1 <- -2*(y-yhat)
    if (deriv > 1) {
      # Calculate second derivative w.r.t. yhat
      l2 <- 2
    }
  }

  # Specify a function to calculate constant C and its derivs w.r.t. v
  C <- function(v, deriv = 0) {
    C <- sqrt(pi/v)
    Cv <- NULL
    Cvv <- NULL
    if (deriv > 0) {
      Cv = -0.5*sqrt(pi)*v^(-3/2)
      if (deriv > 1) {
        Cvv = 0.75 * sqrt(pi)*v^(-5/2)
      }
    }
    return(list(C = C, Cv = Cv, Cvv = Cvv))
  }
  return(list(loss=loss, l1 = l1, l2 = l2, C = C))
}
# Taken from qgam package, comp stable for log(1+exp(x))
##' Numerically stable log(1 + exp(x))
#'
#' Computes log(1 + exp(x)) in a numerically stable way, adapted from the qgam package.
#'
#' @param x Numeric vector.
#'
#' @return Numeric vector of log(1 + exp(x)).
#' @keywords internal
log1px <- function (x)
{
  indx <- .bincode(x, c(-Inf, -37, 18, 33.3, Inf), right = TRUE,
                   include.lowest = TRUE)
  kk <- which(indx == 1)
  if (length(kk)) {
    x[kk] <- exp(x[kk])
  }
  kk <- which(indx == 2)
  if (length(kk)) {
    x[kk] <- log1p(exp(x[kk]))
  }
  kk <- which(indx == 3)
  if (length(kk)) {
    x[kk] <- x[kk] + exp(-x[kk])
  }
  return(x)
}


##' Pinball loss function generator
#'
#' Returns a function that computes the pinball (quantile) loss and its derivatives for a given quantile tau.
#'
#' @param tau Numeric. Quantile level (between 0 and 1).
#'
#' @return Function(y, yhat, deriv = 0) that computes pinball loss and derivatives.
#' @export
pinball_loss <- function(tau) {
  assign(".sigma", NULL, env=environment())
  function(y, yhat, deriv = 0) {
    # Calculate loss
    # Fix scale at sd of residuals
    putSigma <- function(sigma) {assign(".sigma", sigma, env=parent.env(environment(sys.function())))}
    getSigma <- function() {get(".sigma", env=parent.env(environment(sys.function())))}
    if (is.null(getSigma())) {
      if (length(y) == 1) {
        putSigma(1)
      } else {
        putSigma(sd(y-yhat))
      }
    }
    sigma <- getSigma()
    lambda <- 0.01*sigma
    sigma <- 1
    loss <- (tau-1) * (y-yhat)/sigma + lambda*(log1px((y-yhat)/(lambda*sigma)))
    l1 <- NULL
    l2 <- NULL
    if (deriv > 0) {
      # Calculate first derivative w.r.t. yhat
      d = (y-yhat)/(lambda*sigma)
      s = pracma::sigmoid(d)
      l1 <- -(tau-1)/sigma - s/sigma
      if (deriv > 1) {
        # Calculate second derivative w.r.t. yhat
        l2 <- 1/(lambda*sigma^2) * s * (1 - s)
      }
    }

    # Specify a function to calculate constant C and its derivs w.r.t. v
    C <- function(v, deriv = 0) {
      # Define the Beta function parameter values
      x <- v * lambda * (1 - tau)
      y <- v * lambda * tau

      # Calculate I using the Beta function
      C <- sigma * lambda * beta(x, y)
      Cv <- NULL
      Cvv <- NULL

      if (deriv > 0) {
        # Calculate the first derivative
        digamma_x <- Rfast::Digamma(x)
        digamma_y <- Rfast::Digamma(y)
        digamma_sum <- Rfast::Digamma(x + y)
        Cv <- C * lambda * ((1 - tau) * digamma_x + tau * digamma_y - digamma_sum)

        if (deriv > 1) {
          # Calculate the second derivative
          trigamma_x <- Rfast::Trigamma(x)
          trigamma_y <- Rfast::Trigamma(y)
          trigamma_sum <- Rfast::Trigamma(x + y)
          Cvv <- (Cv * lambda * ((1 - tau) * digamma_x + tau * digamma_y - digamma_sum)) +
            C * (lambda^2 * ((1-tau)^2 * trigamma_x + tau^2 * trigamma_y - trigamma_sum))
        }
      }
      # Return the results as a list
      return(list(C = C, Cv = Cv, Cvv = Cvv))
    }
    return(list(loss=loss, l1 = l1, l2 = l2, C = C))
  }
}
