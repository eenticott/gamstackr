# Wrap jacobian function provided by family
# o is gamObject, dat is data.frame, jj is index of response of interest
.jacobian_wrap <- function(o, dat, jj, ...){
  # NB stats:: needed otherwise we end up using the predict() defined
  # in the family!!
  eta <- as.matrix(predict(o, newdata = dat, ...))
  # D sigma / D eta
  DmuDeta <- o$family$jacobian(eta = eta, jj = jj)
  # Dalpha/Dbeta = Dalpha/Deta * Deta/Dbeta
  X <- model.matrix(o, newdata = dat)
  lpi <- attr(X, "lpi")
  if( is.null(lpi) ) { lpi <- list(1:ncol(X)) }
  J <- lapply(seq_along(lpi), function(ii) {
    X[ , lpi[[ii]], drop = FALSE] * DmuDeta[, ii]
  })
  J <- do.call("cbind", J)
  return( J )
}
