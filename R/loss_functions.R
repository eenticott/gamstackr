square_loss <- function(y, param, deriv = 0, ...) {
  # For square loss need derivatives of gaussian distribution
  # derivatives are in terms of t where exp(t) = sig2
  if (is.list(param) ) param <- do.call("cbind", param)
  if (is.vector(param)) param <- matrix(param, nrow = 1)
  if (ncol(param) != 2) stop("Wrong number of parameters provided")

  p <- ncol( param )
  mu <- param[ , 1, drop = TRUE]
  sig2 <- exp(param[ , 2, drop = TRUE]) # sigma2 (force positive)
  tau2 <- 1/sig2
  tau <- sqrt(tau2)

  n <- length(y)

  if (length(mu) == 1) {
    mu <- rep(mu, n)
    tau <- rep(tau, n)
    tau2 <- rep(tau2, n)
  }

  ymu <- y - mu
  ymu2 <- ymu ^ 2
  d0 <- - .5 * log(2 * pi) + log(tau) - .5 * tau2 * ymu2
  out <- list()
  out$d0 <- d0

  if( deriv > 0 )
  {
    d1 <- tau2 * ymu
    d2 <- (1/tau - tau*ymu2) * -1/2 * sig2^(-3/2)
    out[["d1"]] <- list(d1, d2*sig2)

    if( deriv > 1 ){
      d11 <- - tau2
      d12 <- (2 * d1 / tau) * -1/2 * sig2^(-3/2)
      d22 <- ((- ymu2 - 1/tau2) * 1/4 * sig2^(-3)) + ((1/tau - tau*ymu2) * 3/4 * sig2^(-5/2))
      out[["d2"]] <- list(d11, d12*sig2, sig2 * d2 + sig2^2 * d22)

      # Haven't changed later ones to sig2
      if( deriv > 2 ){
        zeros <- rep(0, n)
        d111 <- zeros
        d112 <- - 2 * tau
        d122 <- 2 * ymu
        d222 <- 2 / tau^3
        out[["d3"]] <- list(d111, d112, d122, d222)

        if( deriv > 3){
          d1111 <- d1112 <- d1222 <- zeros
          d1122 <- rep(-2, n)
          d2222 <- -6 / tau2^2
          out[["d3"]] <- list(d1111, d1112, d1122, d1222, d2222)
        }
      }
    }
  }

  return( out )
}

test_f <- function(pars) {
  square_loss(y = c(1), param = list(pars[1], pars[2]), deriv = 0)$d0
}

square_loss(y = c(1), param = list(2, 1), deriv = 2)

numDeriv::grad(test_f, c(2,1))
numDeriv::hessian(test_f, c(2,1))

