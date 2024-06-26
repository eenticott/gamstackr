llk_gaussian <- function(y, param, deriv = 0, ...) {
  # derivatives are sigma2 not tau
  if (is.list(param) ) param <- do.call("cbind", param)
  if (is.vector(param)) param <- matrix(param, nrow = 1)
  if (ncol(param) != 2) stop("Wrong number of parameters provided")

  p <- ncol( param )
  mu <- param[ , 1, drop = TRUE]
  tau <- param[ , 2, drop = TRUE] # 1 / sigma
  tau2 <- tau^2                   # 1 / sigma^2
  n <- length(y)

  if (length(mu) == 1) {
    mu <- rep(mu, n)
    tau <- rep(tau, n)
    tau2 <- rep(tau2, n)
  }

  ymu <- y - mu
  ymu2 <- ymu ^ 2
  theta = 1/tau2
  d0 <- - .5 * log(2 * pi) + log(tau) - .5 * tau2 * ymu2
  out <- list()
  out$d0 <- d0

  if( deriv > 0 )
  {
    d1 <- tau2 * ymu
    d2 <- (1/tau - tau*ymu2) * -1/2 * theta^(-3/2)
    out[["d1"]] <- list(d1, d2)

    if( deriv > 1 ){
      d11 <- - tau2
      d12 <- (2 * d1 / tau) * -1/2 * theta^(-3/2)
      d22 <- ((- ymu2 - 1/tau2) * 1/4 * theta^(-3)) + ((1/tau - tau*ymu2) * 3/4 * theta^(-5/2))
      out[["d2"]] <- list(d11, d12, d22)

      # Haven't changed later ones to theta
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
