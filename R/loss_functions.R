#
#' Squared error loss function and its derivatives
#'
#' This function calculates the squared error loss and its derivatives with respect to predictions.
#'
#' @param y Vector of observed values
#' @param yhat Vector of predicted values
#' @param deriv Integer indicating derivative order (0=value, 1=gradient, 2=Hessian)
#'
#' @return A list containing the loss, derivatives, and a function C for calculating normalization constants
#' @export
#'
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
#' Numerically stable computation of log(1+exp(x))
#'
#' This function computes log(1+exp(x)) in a numerically stable way by using different
#' approximations depending on the range of x.
#'
#' @param x Numeric vector
#'
#' @return Numeric vector containing log(1+exp(x))
#'
#' @export
#'
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


sigmoid <- function(x, a = 1, b = 0) {
  if (length(x) == 0) {return(c())}

  stopifnot(is.numeric(x), is.numeric(a), is.numeric(b))
  a <- a[1]
  b <- b[1]
  return(1/(1 + exp(-a * (x - b))))
}

#' Create a pinball loss function for quantile regression
#'
#' This function creates a pinball loss function for a specified quantile level.
#'
#' @param tau Quantile level (between 0 and 1)
#'
#' @return A function that calculates the pinball loss and its derivatives
#'
#' @importFrom Rfast Digamma Trigamma
#'
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
      s = sigmoid(d)
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
      C <- sigma * lambda * base::beta(x, y)
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
