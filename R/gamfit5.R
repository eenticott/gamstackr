gam.fit5 <- function (x, y, lsp, Sl, weights = NULL, offset = NULL, deriv = 2,
          family, control = gam.control(), Mp = -1, start = NULL,
          gamma = 1)
{
  print("Using this gamfit")
  penalized <- if (length(Sl) > 0)
    TRUE
  else FALSE
  warn <- list()
  nSp <- length(lsp)
  q <- ncol(x)
  nobs <- length(y)
  if (penalized) {
    Eb <- attr(Sl, "E")
    rp <- ldetS(Sl, rho = lsp, fixed = rep(FALSE, length(lsp)),
                np = q, root = TRUE)
    x <- Sl.repara(rp$rp, x)
    Eb <- Sl.repara(rp$rp, Eb)
    St <- crossprod(rp$E)
    E <- rp$E
    attr(E, "use.unscaled") <- TRUE
    if (!is.null(start))
      start <- Sl.repara(rp$rp, start)
  }
  else {
    deriv <- 0
    rp <- list(ldetS = 0, rp = list())
    St <- matrix(0, q, q)
    E <- matrix(0, 0, q)
  }
  start0 <- start
  eval(family$initialize)
  if (!is.null(start0))
    start <- start0
  coef <- as.numeric(start)
  if (is.null(weights))
    weights <- rep.int(1, nobs)
  if (is.null(offset))
    offset <- rep.int(0, nobs)
  llf <- family$ll
  ll <- llf(y, x, coef, weights, family, offset = offset,
            deriv = 1)
  ll0 <- ll$l - drop(t(coef) %*% St %*% coef)/2
  grad <- ll$lb - St %*% coef
  iconv <- max(abs(grad)) < control$epsilon * abs(ll0)
  Hp <- -ll$lbb + St
  rank.checked <- FALSE
  rank <- q
  drop <- NULL
  eigen.fix <- FALSE
  converged <- FALSE
  check.deriv <- FALSE
  eps <- 1e-05
  drop <- NULL
  bdrop <- rep(FALSE, q)
  perturbed <- 0
  for (iter in 1:(2 * control$maxit)) {
    if (check.deriv) {
      fdg <- ll$lb * 0
      fdh <- ll$lbb * 0
      for (k in 1:length(coef)) {
        coef1 <- coef
        coef1[k] <- coef[k] + eps
        ll.fd <- llf(y, x, coef1, weights, family, offset = offset,
                     deriv = 1)
        fdg[k] <- (ll.fd$l - ll$l)/eps
        fdh[, k] <- (ll.fd$lb - ll$lb)/eps
      }
    }
    D <- diag(Hp)
    if (sum(!is.finite(D)) > 0)
      stop("non finite values in Hessian")
    if (min(D) <= 0) {
      Dthresh <- max(D) * sqrt(.Machine$double.eps)
      if (-min(D) < Dthresh) {
        indefinite <- FALSE
        D[D < Dthresh] <- Dthresh
      }
      else indefinite <- TRUE
    }
    else indefinite <- FALSE
    if (indefinite) {
      if (eigen.fix) {
        eh <- eigen(Hp, symmetric = TRUE)
        ev <- abs(eh$values)
        Hp <- eh$vectors %*% (ev * t(eh$vectors))
      }
      else {
        Ib <- diag(rank) * abs(min(D))
        Ip <- diag(rank) * abs(max(D) * .Machine$double.eps^0.5)
        Hp <- Hp + Ip + Ib
      }
      D <- rep(1, ncol(Hp))
      indefinite <- TRUE
    }
    else {
      D <- D^-0.5
      Hp <- D * t(D * Hp)
      Ip <- diag(rank) * .Machine$double.eps^0.5
    }
    L <- suppressWarnings(chol(Hp, pivot = TRUE))
    mult <- 1
    while (attr(L, "rank") < rank) {
      if (eigen.fix) {
        eh <- eigen(Hp, symmetric = TRUE)
        ev <- eh$values
        thresh <- max(min(ev[ev > 0]), max(ev) * 1e-06) *
          mult
        mult <- mult * 10
        ev[ev < thresh] <- thresh
        Hp <- eh$vectors %*% (ev * t(eh$vectors))
        L <- suppressWarnings(chol(Hp, pivot = TRUE))
      }
      else {
        L <- suppressWarnings(chol(Hp + Ip, pivot = TRUE))
        Ip <- Ip * 100
      }
      indefinite <- TRUE
    }
    piv <- attr(L, "pivot")
    ipiv <- piv
    ipiv[piv] <- 1:ncol(L)
    if (converged)
      break
    step <- D * (backsolve(L, forwardsolve(t(L), (D * grad)[piv]))[ipiv])
    c.norm <- sum(coef^2)
    if (c.norm > 0) {
      s.norm <- sqrt(sum(step^2))
      c.norm <- sqrt(c.norm)
      if (s.norm > 0.1 * c.norm)
        step <- step * 0.1 * c.norm/s.norm
    }
    s.norm <- sqrt(sum(step^2))
    if (length(coef) != length(step)) {
      browser()
    }
    coef1 <- coef + step
    ll <- llf(y, x, coef1, weights, family, offset = offset,
              deriv = 1)
    ll1 <- ll$l - drop(t(coef1) %*% St %*% coef1)/2
    khalf <- 0
    fac <- 2
    while ((!is.finite(ll1) || ll1 <= ll0) && khalf < 25) {
      step <- step/fac
      coef1 <- coef + step
      ll <- llf(y, x, coef1, weights, family, offset = offset,
                deriv = 0)
      ll1 <- ll$l - (t(coef1) %*% St %*% coef1)/2
      if (is.finite(ll1) && ll1 >= ll0) {
        ll <- llf(y, x, coef1, weights, family, offset = offset,
                  deriv = 1)
      }
      if (max(abs(coef - coef1)) < max(abs(coef)) * .Machine$double.eps)
        khalf <- 100
      khalf <- khalf + 1
      if (khalf > 5)
        fac <- 5
    }
    if (!is.finite(ll1) || (ll1 <= ll0 && !iconv)) {
      step <- drop(grad) * s.norm/sqrt(sum(grad^2))
      khalf <- 0
    }
    while ((!is.finite(ll1) || (ll1 <= ll0 && !iconv)) &&
           khalf < 25) {
      step <- step/10
      coef1 <- coef + step
      ll <- llf(y, x, coef1, weights, family, offset = offset,
                deriv = 0)
      ll1 <- ll$l - (t(coef1) %*% St %*% coef1)/2
      if (is.finite(ll1) && ll1 >= ll0) {
        ll <- llf(y, x, coef1, weights, family, offset = offset,
                  deriv = 1)
      }
      if (max(abs(coef - coef1)) < max(abs(coef)) * .Machine$double.eps)
        khalf <- 100
      khalf <- khalf + 1
    }
    if ((is.finite(ll1) && ll1 >= ll0 && khalf < 25) ||
        iter == control$maxit) {
      coef <- coef + step
      grad <- ll$lb - St %*% coef
      Hp <- -ll$lbb + St
      ok <- (iter == control$maxit || max(abs(grad)) <
               control$epsilon * abs(ll0))
      if (ok) {
        if (indefinite) {
          if (perturbed == 5)
            stop("indefinite penalized likelihood in gam.fit5 ")
          if (iter < 4 || rank.checked) {
            perturbed <- perturbed + 1
            coef <- coef * (1 + (runif(length(coef)) *
                                   0.02 - 0.01) * perturbed) + (runif(length(coef)) -
                                                                  0.5) * mean(abs(coef)) * 1e-05 * perturbed
            ll <- llf(y, x, coef, weights, family, offset = offset,
                      deriv = 1)
            ll0 <- ll$l - (t(coef) %*% St %*% coef)/2
          }
          else {
            rank.checked <- TRUE
            if (penalized) {
              Sb <- crossprod(Eb)
              Hb <- -ll$lbb/norm(ll$lbb, "F") + Sb/norm(Sb,
                                                        "F")
            }
            else Hb <- -ll$lbb/norm(ll$lbb, "F")
            D <- abs(diag(Hb))
            D[D < 1e-50] <- 1
            D <- D^-0.5
            Hb <- t(D * Hb) * D
            qrh <- qr(Hb, LAPACK = TRUE)
            rank <- Rrank(qr.R(qrh))
            if (rank < q) {
              drop <- sort(qrh$pivot[(rank + 1):q])
              bdrop <- 1:q %in% drop
              lpi <- attr(x, "lpi")
              xat <- attributes(x)
              xat$dim <- xat$dimnames <- NULL
              coef <- coef[-drop]
              St <- St[-drop, -drop]
              x <- x[, -drop]
              if (!is.null(lpi)) {
                ii <- (1:q)[!bdrop]
                ij <- rep(NA, q)
                ij[ii] <- 1:length(ii)
                for (i in 1:length(lpi)) {
                  lpi[[i]] <- ij[lpi[[i]][!(lpi[[i]] %in%
                                              drop)]]
                }
              }
              if (length(xat) > 0)
                for (i in 1:length(xat)) attr(x, names(xat)[i]) <- xat[[i]]
              attr(x, "lpi") <- lpi
              attr(x, "drop") <- drop
              ll <- llf(y, x, coef, weights, family,
                        offset = offset, deriv = 1)
              ll0 <- ll$l - (t(coef) %*% St %*% coef)/2
            }
          }
        }
        else {
          converged <- TRUE
        }
      }
      else ll0 <- ll1
    }
    else {
      if (is.null(drop))
        bdrop <- rep(FALSE, q)
      if (iconv && iter == 1) {
        converged <- TRUE
        coef <- start
      }
      else {
        converged <- FALSE
        warn[[length(warn) + 1]] <- paste("gam.fit5 step failed: max magnitude relative grad =",
                                          max(abs(grad/drop(ll0))))
      }
      break
    }
  }
  if (iter == 2 * control$maxit && converged == FALSE)
    warn[[length(warn) + 1]] <- gettextf("gam.fit5 iteration limit reached: max abs grad = %g",
                                         max(abs(grad)))
  ldetHp <- 2 * sum(log(diag(L))) - 2 * sum(log(D))
  if (!is.null(drop)) {
    fcoef <- rep(0, length(bdrop))
    fcoef[!bdrop] <- coef
  }
  else fcoef <- coef
  dVkk <- d1l <- d2l <- d1bSb <- d2bSb <- d1b <- d2b <- d1ldetH <- d2ldetH <- d1b <- d2b <- NULL
  if (deriv > 0) {
    m <- nSp
    d1b <- matrix(0, rank, m)
    Sib <- Sl.termMult(rp$Sl, fcoef, full = TRUE)
    if (nSp)
      for (i in 1:m) d1b[, i] <- -D * (backsolve(L, forwardsolve(t(L),
                                                                 (D * Sib[[i]][!bdrop])[piv]))[ipiv])
    dVkk <- crossprod(L[, ipiv] %*% (d1b/D))
    if (!is.null(drop)) {
      fd1b <- matrix(0, q, m)
      fd1b[!bdrop, ] <- d1b
    }
    else fd1b <- d1b
    ll <- llf(y, x, coef, weights, family, offset = offset,
              deriv = 3, d1b = d1b)
    if (deriv > 1) {
      d2b <- matrix(0, rank, m * (m + 1)/2)
      k <- 0
      for (i in 1:m) for (j in i:m) {
        k <- k + 1
        v <- -ll$d1H[[i]] %*% d1b[, j] + Sl.mult(rp$Sl,
                                                 fd1b[, j], i)[!bdrop] + Sl.mult(rp$Sl, fd1b[,
                                                                                             i], j)[!bdrop]
        d2b[, k] <- -D * (backsolve(L, forwardsolve(t(L),
                                                    (D * v)[piv]))[ipiv])
        if (i == j)
          d2b[, k] <- d2b[, k] + d1b[, i]
      }
      llr <- llf(y, x, coef, weights, family, offset = offset,
                 deriv = 4, d1b = d1b, d2b = d2b, Hp = Hp, rank = rank,
                 fh = L, D = D)
      d2l <- matrix(0, m, m)
      for (i in 1:m) for (j in i:m) {
        d2l[j, i] <- d2l[i, j] <- t(d1b[, i]) %*% ll$lbb %*%
          d1b[, j]
      }
    }
  }
  if (deriv > 0) {
    d1ldetH <- rep(0, m)
    d1Hp <- list()
    for (i in 1:m) {
      A <- -ll$d1H[[i]] + Sl.mult(rp$Sl, diag(q), i)[!bdrop,
                                                     !bdrop]
      d1Hp[[i]] <- D * (backsolve(L, forwardsolve(t(L),
                                                  (D * A)[piv, ]))[ipiv, , drop = FALSE])
      d1ldetH[i] <- sum(diag(d1Hp[[i]]))
    }
  }
  if (deriv > 1) {
    d2ldetH <- matrix(0, m, m)
    k <- 0
    for (i in 1:m) for (j in i:m) {
      k <- k + 1
      d2ldetH[i, j] <- -sum(d1Hp[[i]] * t(d1Hp[[j]])) -
        llr$trHid2H[k]
      if (i == j) {
        A <- Sl.mult(rp$Sl, diag(q), i, full = TRUE)[!bdrop,
                                                     !bdrop]
        bind <- rowSums(abs(A)) != 0
        A <- A[, bind, drop = FALSE]
        A <- D * (backsolve(L, forwardsolve(t(L), (D *
                                                     A)[piv, ]))[ipiv, , drop = FALSE])
        d2ldetH[i, j] <- d2ldetH[i, j] + sum(diag(A[bind,
                                                    , drop = FALSE]))
      }
      else d2ldetH[j, i] <- d2ldetH[i, j]
    }
  }
  if (deriv > 0) {
    Skb <- Sl.termMult(rp$Sl, fcoef, full = TRUE)
    d1bSb <- rep(0, m)
    for (i in 1:m) {
      Skb[[i]] <- Skb[[i]][!bdrop]
      d1bSb[i] <- sum(coef * Skb[[i]])
    }
  }
  if (deriv > 1) {
    d2bSb <- matrix(0, m, m)
    for (i in 1:m) {
      Sd1b <- St %*% d1b[, i]
      for (j in i:m) {
        d2bSb[j, i] <- d2bSb[i, j] <- 2 * sum(d1b[,
                                                  i] * Skb[[j]] + d1b[, j] * Skb[[i]] + d1b[,
                                                                                            j] * Sd1b)
      }
      d2bSb[i, i] <- d2bSb[i, i] + sum(coef * Skb[[i]])
    }
  }
  REML <- -as.numeric((ll$l - drop(t(coef) %*% St %*% coef)/2)/gamma +
                        rp$ldetS/2 - ldetHp/2 + Mp * (log(2 * pi)/2) - log(gamma)/2)
  REML1 <- if (deriv < 1)
    NULL
  else -as.numeric(-d1bSb/(2 * gamma) + rp$ldet1/2 - d1ldetH/2)
  if (control$trace) {
    cat("\niter =", iter, "  ll =", ll$l, "  REML =", REML,
        "  bSb =", t(coef) %*% St %*% coef/2, "\n")
    cat("log|S| =", rp$ldetS, "  log|H+S| =", ldetHp, "  n.drop =",
        length(drop), "\n")
    if (!is.null(REML1))
      cat("REML1 =", REML1, "\n")
  }
  REML2 <- if (deriv < 2)
    NULL
  else -((d2l - d2bSb/2)/gamma + rp$ldet2/2 - d2ldetH/2)
  lpi <- attr(x, "lpi")
  if (is.null(lpi)) {
    linear.predictors <- if (is.null(offset))
      as.numeric(x %*% coef)
    else as.numeric(x %*% coef + offset)
    fitted.values <- family$linkinv(linear.predictors)
  }
  else {
    fitted.values <- linear.predictors <- matrix(0, nrow(x),
                                                 length(lpi))
    if (!is.null(offset))
      offset[[length(lpi) + 1]] <- 0
    for (j in 1:length(lpi)) {
      linear.predictors[, j] <- as.numeric(x[, lpi[[j]],
                                             drop = FALSE] %*% coef[lpi[[j]]])
      if (!is.null(offset[[j]]))
        linear.predictors[, j] <- linear.predictors[,
                                                    j] + offset[[j]]
      fitted.values[, j] <- family$linfo[[j]]$linkinv(linear.predictors[,
                                                                        j])
    }
  }
  coef <- Sl.repara(rp$rp, fcoef, inverse = TRUE)
  if (!is.null(drop) && !is.null(d1b)) {
    db.drho <- matrix(0, length(bdrop), ncol(d1b))
    db.drho[!bdrop, ] <- d1b
  }
  else db.drho <- d1b
  if (!is.null(d1b))
    db.drho <- Sl.repa(rp$rp, db.drho, l = -1)
  ret <- list(coefficients = coef, family = family, y = y,
              prior.weights = weights, fitted.values = fitted.values,
              linear.predictors = linear.predictors, scale.est = 1,
              REML = REML, REML1 = REML1, REML2 = REML2, rank = rank,
              aic = -2 * ll$l, l = ll$l, lbb = ll$lbb, L = L, bdrop = bdrop,
              D = D, St = St, rp = rp$rp, db.drho = db.drho, S1 = rp$ldet1,
              iter = iter, H = ll$lbb, dH = ll$d1H, dVkk = dVkk, warn = warn)
  ret
}
