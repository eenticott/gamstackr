orth_mat <- function(vec, qq) {
  rbind(diag(qq), matrix(vec, ncol = qq))
}

LatentStack <- function(logP, q, ridgePen = 1e-5) {
  # mgcv latent stacking family
  # q, int, number of latent predictors to estimate

  # Prep
  link <- "identity"
  P <- exp(logP)
  K <- ncol(logP)
  nv <- (K-1-q)*q

  nlp <- q # Number of linear predictors
  # Prepare link functions
  link <- lapply(1:nlp, function(x) "identity")
  stats <- list()
  for (ii in 1:nlp) {
    stats[[ii]] <- stats::make.link(link[[ii]])
    fam <- structure(list(link=link[[ii]],canonical="none",linkfun=stats[[ii]]$linkfun,
                          mu.eta=stats[[ii]]$mu.eta), class="family")
    fam <- mgcv::fix.family.link(fam)
    stats[[ii]]$d2link <- fam$d2link
    stats[[ii]]$d3link <- fam$d3link
    stats[[ii]]$d4link <- fam$d4link
  }

  ### Saving extra parameters in .GlobalEnv environment
  ## Add anything you might want to retrieve/store internally
  # logP
  assign(".logP", logP, envir = environment())
  getLogP <- function() get(".logP")
  putLogP <- function(.x) assign(".logP", .x, envir = environment(sys.function()))
  # P
  assign(".P", P, envir = environment())
  getP <- function() get(".P")
  putP <- function(.x) assign(".P", .x, envir = environment(sys.function()))
  # ridgePen
  assign(".ridgePen", ridgePen, envir = environment())
  getRidgePen <- function() get(".ridgePen")
  putRidgePen <- function(.x) assign(".ridgePen", .x, envir = environment(sys.function()))
  # P
  assign(".p", NULL, envir = environment())
  getp <- function() get(".p")
  putp <- function(.x) assign(".p", .x, envir = environment(sys.function()))
  # q
  assign(".q", q, envir = environment())
  getq <- function() get(".q")
  putq <- function(.x) assign(".q", .x, envir = environment(sys.function()))
  # nv
  assign(".nv", nv, envir = environment())
  getnv <- function() get(".nv")
  putnv <- function(.x) assign(".nv", .x, envir = environment(sys.function()))
  # lpi
  assign(".lpi", NULL, envir = environment())
  getlpi <- function() get(".lpi")
  putlpi <- function(.x) assign(".lpi", .x, envir = environment(sys.function()))
  # ncoef
  assign(".ncoef", NULL, envir = environment())
  getncoef <- function() get(".ncoef")
  putncoef <- function(.x) assign(".ncoef", .x, envir = environment(sys.function()))

  residuals <- function(object, type=c("deviance","pearson","response")) {
    return(as.matrix(object$y)[, 1])
  }

  preinitialize <- function(G) {
    ## G is a gam pre-fit object. Pre-initialize can manipulate some of its
    ## elements, returning a named list of the modified ones.

    # Some elements that may need adjusted
    # G$X
    # G$term.names
    # G$S1
    # attributes(G$X, "lpi")
    nv <- getnv()
    atr <- attributes(G$X)
    lpi <- atr$lpi

    # Have lpi operate over same columns
    Xids <- unique(colnames(G$X))
    Xid <- 1:length(Xids)
    names(Xid) <- Xids
    newlpi <- lapply(lpi, function(ii) as.numeric(Xid[colnames(G$X[,ii])]))
    newX <- unique(G$X, MARGIN = 2)

    G$family$newX <- newX
    G$family$newlpi <- newlpi

    p <- ncol(newX) # number of x parameters
    putp(p) # store

    putlpi(lpi)
    nbeta <- lapply(lpi, function(lpi_ii) length(lpi_ii))
    nbeta <- do.call("sum", nbeta)
    nparams <- nv + K-1 + nbeta
    putncoef(nparams)
    G$X <- cbind(G$X,matrix(0,nrow(G$X), nv + K-1)) # add dummy cols
    attr(G$X, "lpi") <- lpi
    attr(G$X, "dim") <- c(nrow(G$X), ncol(G$X))
    if (!is.null(G$Sl)) {
      ## pad out sqrt of balanced penalty matrix to account for extra params
      attr(G$Sl, "E") <- cbind(attr(G$Sl, "E"),
                               matrix(0, p, nv + K - 1))
      dfc <- expand.grid((q+1):(K-1),(q+1):(K-1))
      G$term.names <- c(G$term.names,
                        paste0("V", paste(dfc[,1], dfc[,2], sep = ",")),
                        paste0("B",1:(K-1),"0"))
    }
    list(X=G$X,term.names=G$term.names,family=G$family, Sl = G$Sl)
  } ## preinitialize

  initialize <- expression({
    my_init_fun <- function(y, nobs, E, x, family, offset){
      # Initialize at 0 for now
      ncoef <- family$getncoef()
      lpi <- family$getlpi()

      p <- lapply(lpi, function(lpi_ii) length(lpi_ii))
      augp <- cumsum(p)
      starting_coefs <- rep(0, ncoef)


      for (ii in 1:length(p)) {
        if (ii == 1) {
          starting_coefs[1:augp[ii]] <- 1:augp[ii]
        } else {
          starting_coefs[(augp[ii-1]+1):augp[ii]] <- (augp[ii-1]+1):augp[ii]
        }
      }

      return(starting_coefs)
    }
    if(is.null(start)){
      start <- my_init_fun(y = y, nobs = nobs, E = E, x = x, family = family, offset = offset)
    }
  })

  ll <- function(y,x,coef,wt,family,offset=NULL,deriv=0,
                 dlb=0,d2b=0,Hp=NULL,rank=0,fh=NULL,D=NULL) {
    ## deriv: 0 - eval
    ##        1 - grad and Hess
    ##        2 - diagonal of first deriv of Hess
    ##        3 - first deriv of Hess
    l <- NULL
    lb <- NULL
    lbb <- NULL

    lpi <- family$newlpi
    x <- family$newX
    p <- lapply(lpi, function(lpi_ii) length(lpi_ii))
    nbeta <- do.call("sum", p)
    max_beta <- do.call("max", lpi)
    beta_matrix <- matrix(0, nrow = max_beta, ncol = length(lpi))

    augp <- cumsum(p)
    for (ii in 1:length(p)) {
      if (ii == 1) {
        beta_matrix[lpi[[ii]],ii] <- coef[1:augp[ii]]
      } else {
        beta_matrix[lpi[[ii]],ii] <- coef[(augp[ii-1]+1):augp[ii]]
      }
    }

    nu0 <- sapply(1:q, function(ii) {
      x[,lpi[[ii]],drop=FALSE] %*% coef[lpi[[ii]]]
    })

    KK <- K
    pp <- family$getp()
    qq <- q

    par <- coef[-(1:nbeta)]
    beta0 <- par[1:(KK-1)]
    delta <- par[-(1:(KK-1))]

    delta_matrix <- orth_mat(vec = delta, qq)
    X <- logP - rowMaxs(logP)
    expX <- exp(X)

    xb <- x[,1:max_beta]
    nu <- xb %*% beta_matrix %*% t(delta_matrix)
    nu <- t(t(nu) + beta0)

    nu1 <- cbind(0, nu)
    nuCen <- nu1 - rowMaxs(nu1)
    a <- exp(nuCen)/rowSums(exp(nuCen))
    l <- sum(log(rowSums(a * P)))

    if (deriv >= 1) {
      keep_idx <- c(as.vector(sapply(1:length(lpi), function(i) 1:max_beta %in% lpi[[i]])), rep(TRUE, nv + K-1))
      # Gradient
      nupX <- nu + X[, -1, drop = FALSE]
      w <- exp(nupX - log(expX[, 1] + rowSums(exp(nupX))))
      ln <- w - a[, -1, drop = FALSE]

      gr_beta0 <- colSums(ln)

      txln <- t(xb) %*% ln
      gr_beta <- txln %*% delta_matrix
      gr_V <- t(txln[, -(1:qq), drop = FALSE]) %*% beta_matrix
      lb <- c(as.numeric(gr_beta), gr_beta0, as.numeric(gr_V))
      lb <- lb[keep_idx]
      # Hessian
      lnn <- list()
      for (jj in 1:(KK - 1)) {
        lnn[[jj]] <- list()
        for (kk in 1:(KK - 1)) {
          lnn[[jj]][[kk]] <- list()
        }
      }
      # coun <- 0
      for (jj in 1:(KK - 1)) {
        for (kk in jj:(KK - 1)) {
          # coun <- coun + 1
          # lnn[[coun]] <- ln[, jj] * (as.numeric(jj == kk) - w[, kk]) - a[, jj + 1] * ln[, kk]
          lnn[[jj]][[kk]] <- lnn[[kk]][[jj]] <-
            ln[, jj] * (as.numeric(jj == kk) - w[, kk]) - a[, jj + 1] * ln[, kk]
        }
        lnn[[jj]] <- do.call(cbind, lnn[[jj]])
        # lnn <- do.call(cbind, lnn)
      }
      lnn <- do.call(abind, args = list(lnn, along = 3))

      lbb0 <- apply(lnn, 2:3, sum)

      D <- array.mult(t(delta_matrix), delta_matrix, lnn)
      lbb <- list()
      for (QQ1 in 1:qq) {
        lbb[[QQ1]] <- list()
        for (QQ2 in 1:qq) {
          lbb[[QQ1]][[QQ2]] <- list()
        }
      }
      for (QQ1 in 1:qq) {
        for (QQ2 in QQ1:qq) {
          lbb[[QQ1]][[QQ2]] <- lbb[[QQ2]][[QQ1]] <- t(xb) %*% (xb * D[, QQ1, QQ2])
        }
        lbb[[QQ1]] <- do.call(cbind, lbb[[QQ1]])
      }
      lbb <- do.call(rbind, lbb)

      l_b0_b <- array.mult(t(xb), delta_matrix, aperm(lnn, c(2, 1, 3)))
      l_b0_b <- matrix(l_b0_b, nrow = KK-1, ncol = pp*qq)

      ldd <- list()
      for (QQ1 in 1:qq) {
        ldd[[QQ1]] <- list()
        for (QQ2 in 1:qq) {
          ldd[[QQ1]][[QQ2]] <- list()
        }
      }
      for (QQ1 in 1:qq) {
        for (QQ2 in QQ1:qq) {
          lnn12 <- lnn[ , (qq+1):(KK-1), (qq+1):(KK-1)]
          ldd[[QQ1]][[QQ2]] <- ldd[[QQ2]][[QQ1]] <-
            apply(as.numeric(xb %*% beta_matrix[, QQ1]) *
                    as.numeric(xb %*% beta_matrix[, QQ2]) *
                    lnn12, 2:3, sum)
        }
        ldd[[QQ1]] <- do.call(cbind, ldd[[QQ1]])
      }
      ldd <- do.call(rbind, ldd)

      l_b0_d <- t(xb %*% beta_matrix) %*% matrix(lnn, nrow = n, ncol = (KK-1) * (KK-1))
      l_b0_d <- array(l_b0_d, dim = c(qq, KK-1, KK-1))
      l_b0_d <- aperm(l_b0_d, c(2, 3, 1))
      l_b0_d <- matrix(l_b0_d[, (qq+1):(KK-1), ], nrow = KK-1)


      l_b_d <- list()
      for (QQ1 in 1:qq) {
        l_b_d[[QQ1]] <- list()
        for (QQ2 in 1:qq) {
          l_b_d[[QQ1]][[QQ2]] <- list()
        }
      }
      QQ1 <- 1
      for (QQ1 in 1:qq) {
        for (QQ2 in 1:qq) {
          l_b_d_qq <- t(xb * as.numeric(xb %*% beta_matrix[, QQ2])) %*%
            matrix(lnn[ , (qq+1):(KK-1), ], nrow = n)
          l_b_d_qq <- array(l_b_d_qq, dim = c(pp, KK-1-qq, KK - 1))
          l_b_d_qq <- matrix(l_b_d_qq, ncol = KK-1) %*% delta_matrix[, QQ1]
          l_b_d_qq <- matrix(l_b_d_qq, nrow = pp)
          if (QQ1 == QQ2) l_b_d_qq <- l_b_d_qq + t(xb) %*% ln[,(qq+1):(KK-1)]
          l_b_d[[QQ1]][[QQ2]] <- l_b_d_qq

        }
        l_b_d[[QQ1]] <- do.call(cbind, l_b_d[[QQ1]])
      }
      l_b_d <- do.call(rbind, l_b_d)

      hes <- rbind(
        cbind(lbb0, l_b0_b, l_b0_d),
        cbind(t(l_b0_b), lbb, l_b_d),
        cbind(t(l_b0_d), t(l_b_d), ldd)
      )

      hes <- rbind(
        cbind(lbb, t(l_b0_b), l_b_d),
        cbind(l_b0_b, lbb0, l_b0_d),
        cbind(t(l_b_d), t(l_b0_d), ldd)
      )

      lbb <- hes[keep_idx, keep_idx]
    }
    return(list(l = l, lb = lb, lbb = lbb))
  }

  jacobian <- function(eta, jj, ...) {
    return(NULL)
  }

  predict <- function(family,se=FALSE,eta=NULL,y=NULL,X=NULL,
                      beta=NULL,off=NULL,Vb=NULL)  {

    return(list("fit" = NULL))
  }

  structure(list(family = "NestedStack",ll = ll,nlp = q,
                 link = "identity",
                 getLogP = getLogP,
                 getP = getP,
                 getRidgePen = getRidgePen,
                 putLogP = putLogP,
                 putP = putP,
                 getlpi = getlpi,
                 putlpi = putlpi,
                 putp = putp,
                 getp = getp,
                 getncoef = getncoef,
                 putncoef = putncoef,
                 preinitialize = preinitialize,
                 initialize = initialize,
                 jacobian = jacobian,
                 mu.eta = stats[[1]]$mu.eta,
                 # postproc=postproc,
                 tri = mgcv::trind.generator(max(1, K - 1)), ## symmetric indices for accessing deriv. arrays
                 residuals=residuals,
                 linfo = stats,
                 #rd=rd,
                 # dev.resids = dev.resids,
                 linkinv = stats$linkinv, # MAYBE IT'S NEEDED IN gam.fit5
                 d2link=1,d3link=1,d4link=1, ## signals to fix.family.link that all done,
                 predict = predict,
                 ls=1, ## signals that ls not needed here
                 available.derivs = 0, ## signal only first derivatives available...
                 discrete.ok = FALSE
  ), class = c("general.family", "extended.family","family"))

}
