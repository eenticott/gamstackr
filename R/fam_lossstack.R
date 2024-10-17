get_loss_derivs <- function(list_of_beta,
                            list_of_X,
                            theta,
                            tau,
                            loss,
                            weight,
                            preds,
                            y) {

  K <- ncol(preds)
  Ws <- weight

  n_eta <- attr(Ws, "neta")
  n_theta <- attr(Ws, "ntheta")

  eta <- get_list_of_eta(list_of_X, list_of_beta)
  W <- Ws(eta, theta, 2)
  mu <- rowSums(W$f_eval * preds)

  llk_deriv <- loss(y, param = list(mu, tau), deriv = 2)

  ll <- sum(llk_deriv$d0)


  # First derivs
  db <- NULL
  dt <- NULL
  if (n_eta > 0) {
    dmudeta <- lapply(1:n_eta, function(k) {rowSums(W$f_eta_eval[[k]] * preds)})
    db <- sapply(1:n_eta,function(k) colSums(dmudeta[[k]] * llk_deriv$d1[[1]] * list_of_X[[k]]))
  }
  if (n_theta > 0) {
    dmudtheta <- lapply(1:n_theta, function(k) {rowSums(W$f_theta_eval[[k]] * preds)})
    dt <- sapply(1:n_theta,function(k) sum(dmudtheta[[k]] * llk_deriv$d1[[1]])) # could make more efficient by combining theta into matrix then multiplying by llk_deriv
  }

  dtau <- sum(llk_deriv$d1[[2]])

  grad <- c(db, dt, dtau)

  # Second derivs
  # beta beta deriv block
  dbb <- NULL
  if (n_eta > 0) {
    dbb <- list()
    for (k in 1:n_eta) {
      dbb[[k]] <- list()
      for (j in 1:n_eta) {
        dbb[[k]][[j]] <- (t(list_of_X[[k]]) %*% (dmudeta[[k]]*dmudeta[[j]] * list_of_X[[j]] * llk_deriv$d2[[1]])) +
          t(list_of_X[[k]]) %*% (rowSums(W$f_eta2_eval[[k]][[j]] * preds) * llk_deriv$d1[[1]] * list_of_X[[j]])
      }
      dbb[[k]] <- do.call("cbind", dbb[[k]])
    }
    dbb <- do.call("rbind", dbb)
  }

  # beta theta deriv block
  dbt <- NULL
  if (n_eta > 0 & n_theta > 0) {
    dbt <- list()
    for (k in 1:n_eta) {
      dbt[[k]] <- list()
      for (j in 1:n_theta) {
        dbt[[k]][[j]] <- colSums(llk_deriv$d2[[1]] * dmudeta[[k]] * dmudtheta[[j]] * list_of_X[[k]]) +
          colSums(llk_deriv$d1[[1]] * rowSums(W$f_eta_theta_eval[[k]][[j]]*preds) * list_of_X[[k]])
      }
      dbt[[k]] <- do.call("cbind", dbt[[k]])
    }
    dbt <- do.call("rbind", dbt)
  }

  # beta tau deriv block
  dbtau <- NULL
  if (n_eta > 0) {
    dbtau <- do.call("c", lapply(1:n_eta, function(k) colSums(llk_deriv$d2[[2]] * dmudeta[[k]] * list_of_X[[k]])))
  }

  # theta theta deriv block
  dtt <- NULL
  if (n_theta > 0) {
    dtt <- list()
    for (k in 1:n_theta) {
      dtt[[k]] <- list()
      for (j in 1:n_theta) {
        dtt[[k]][[j]] <- sum(llk_deriv$d1[[1]]*rowSums(W$f_theta2_eval[[k]][[j]]*preds) +
                               llk_deriv$d2[[1]]*dmudtheta[[k]]*dmudtheta[[j]])
      }
      dtt[[k]] <- do.call("cbind", dtt[[k]])
    }
    dtt <- do.call("rbind", dtt)
  }


  # theta tau deriv block
  dttau <- NULL
  if (n_theta > 0) {
    dttau <- do.call("c", lapply(1:n_theta, function(k) sum(llk_deriv$d2[[2]] * dmudtheta[[k]])))
  }

  # tau tau deriv block
  dtautau <- matrix(sum(llk_deriv$d2[[3]]))

  my_t <- function(x) {
    if (is.null(x)) return(NULL)
    return(t(x))
  }

  hess <- rbind(cbind(dbb, dbt, dbtau),
  cbind(my_t(dbt), dtt, dttau),
  cbind(my_t(dbtau), my_t(dttau), dtautau))

  return(list(l = ll, lb = grad, lbb = hess))
}

LossStack <- function(preds, loss, weights, ridgePen = 1e-5) {
  # mgcv point forecast stacking family

  # Prep
  link <- "identity"
  K <- ncol(preds)
  n_theta <- attr(weights, "ntheta")
  nlp <- attr(weights, "neta") # Number of linear predictors

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
  # preds
  assign(".preds", preds, envir = environment())
  getPreds <- function() get(".preds")
  putPreds <- function(.x) assign(".preds", .x, envir = environment(sys.function()))

  # ridgePen
  assign(".ridgePen", ridgePen, envir = environment())
  getRidgePen <- function() get(".ridgePen")
  putRidgePen <- function(.x) assign(".ridgePen", .x, envir = environment(sys.function()))

  # lpi
  assign(".lpi", NULL, envir = environment())
  getlpi <- function() get(".lpi")
  putlpi <- function(.x) assign(".lpi", .x, envir = environment(sys.function()))

  # weightf
  assign(".weightf",weights, envir = environment())
  getweightf <- function() get(".weightf")

  residuals <- function(object, type=c("deviance","pearson","response")) {
    return(as.matrix(object$y)[, 1])
  }

  preinitialize <- function(G) {
    ## G is a gam pre-fit object. Pre-initialize can manipulate some of its
    ## elements, returning a named list of the modified ones.
    nbeta <- ncol(G$X)
    ntheta <- attr(G$family$weights, "ntheta") + 1
    atr <- attributes(G$X)
    G$X <- cbind(G$X,matrix(0,nrow(G$X),ntheta)) ## add dummy columns to G$X

    attr(G$X, "lpi") <- atr$lpi
    attr(G$X, "dim") <- c(nrow(G$X), ncol(G$X))

    putlpi(atr$lpi)
    if (!is.null(G$Sl))
      attr(G$Sl, "E") <- cbind(attr(G$Sl, "E"),
                               matrix(0, nbeta, ntheta))
    if (ntheta != 0) {
      G$term.names <- c(G$term.names, paste("sig2", 1:ntheta,
                                            sep = "."))
    }
    ## pad out sqrt of balanced penalty matrix to account for extra params

    list(X=G$X,term.names=G$term.names,family=G$family, Sl = G$Sl)
  } ## preinitialize

  initialize <- expression({
    my_init_fun <- function(y, nobs, E, x, family, offset){
      W_f <- family$getweightf()

      # Initialize at 0 for now
      coefs <- rep(0, ncol(x))
      lpi <- attr(x,"lpi")
      p <- lapply(lpi, function(lpi_ii) length(lpi_ii))
      neta <- attr(W_f, "neta")
      nbeta <- do.call("sum", p)
      ntheta <- attr(W_f, "ntheta")
      ntau <- 1

      preds <- family$getPreds()
      scores <- -(preds - y)^2 # Want to score based on loss function, max score should be best model
      init_coefs <- attr(W_f, "init_fun")(scores)

      ridgePen <- family$getRidgePen()
      use.unscaled <- if (!is.null(attr(E,"use.unscaled"))) TRUE else FALSE

      # Penreg func
      start_coef <- function(x1, e1, y1) {
        if (!is.null(ridgePen)) {
          sqrt(ridgePen)
          e1rp <- matrix(0, nrow = nrow(e1),
                         ncol = ncol(e1))
          diag(e1rp) <- sqrt(ridgePen)
          e1 <- e1 + e1rp
        }
        if (use.unscaled) {
          qrx <- qr(rbind(x1,e1))
          x1 <- rbind(x1,e1)
          startji <- qr.coef(qr(x1),c(y1,rep(0,nrow(E))))
          startji[!is.finite(startji)] <- 0
        } else startji <- penReg(x1,e1,y1)

        return(startji)
      }

      list_of_X <- list()
      list_of_beta <- list()

      for (jj in 1:neta) {
        y1 <- init_coefs$init_mu[,jj]
        x1 <- x[, lpi[[jj]], drop = FALSE]
        e1 <- E[, lpi[[jj]], drop = FALSE]

        coefs[lpi[[jj]]] <- start_coef(x1, e1, y1)
        list_of_beta[[jj]] <- coefs[lpi[[jj]]]

        # exp fix for multinomial
        if (attr(family$weights, "name") == "multinomial") {
          y1 <- x1 %*% coefs[lpi[[jj]]]
          coefs[lpi[[jj]]] <- start_coef(x1, e1, log(y1))
          list_of_beta[[jj]] <- coefs[lpi[[jj]]]
          list_of_X[[jj]] <- x1
        }
      }
      list_of_eta <- .internals()[["get_list_of_eta"]](list_of_X, list_of_beta)

      # Initial theta
      if (ntheta > 0) {
        coefs[nbeta + 1:ntheta] <- init_coefs$init_theta
      }

      # logscale of residual variance under initial weights
      coefs[nbeta + ntheta + 1] <- log(var(rowSums(W_f(list_of_eta, theta, deriv = 0)$f_eval * preds)-y))

      return(coefs)
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

    if (deriv == 1) {
      deriv = 2 # My deriv function returns hessian on 2.
    }
    lpi <- family$getlpi()
    p <- lapply(lpi, function(lpi_ii) length(lpi_ii))

    list_of_beta <- list()
    list_of_X <- list()

    for (i in 1:length(lpi)) {
      list_of_beta[[i]] <- coef[lpi[[i]]]
      list_of_X[[i]] <- x[,lpi[[i]], drop=FALSE]
    }
    theta <- NULL
    if (family$n_theta > 0) {
      theta <- coef[do.call("sum", p):(length(coef)-1)]
    }
    tau <- coef[length(coef)]
    preds <- family$getPreds()
    out <- get_loss_derivs(list_of_beta = list_of_beta,
                           list_of_X = list_of_X,
                           theta = theta,
                           tau = tau,
                           loss = loss,
                           weight = weights,
                           preds = getPreds())

    return(out)
  }

  jacobian <- function(eta, jj, ...) {
    return(NULL)
  }

  predict <- function(family,se=FALSE,eta=NULL,y=NULL,X=NULL,
                      beta=NULL,off=NULL,Vb=NULL)  {
    lpi <- family$getlpi()
    p <- lapply(lpi, function(lpi_ii) length(lpi_ii))
    K <- length(lpi)
    N <- nrow(X)
    list_of_beta <- list()
    list_of_X <- list()

    for (i in 1:K) {
      list_of_beta[[i]] <- beta[lpi[[i]]]
      list_of_X[[i]] <- X[,lpi[[i]], drop=FALSE]
    }

    list_of_eta <- .internals()[["get_list_of_eta"]](list_of_X, list_of_beta)
    theta <- NULL
    if (attr(family$getweightf(), "ntheta") > 0) {theta <- beta[do.call("sum", p):(length(beta)-1)]
}

    W <- family$getweightf()(list_of_eta, theta, 0)

    return(list("fit" = W$f_eval))
  }

  structure(list(family = "LossStack",ll = ll,nlp = attr(weights, "neta"),
                 link = "identity",
                 getPreds = getPreds,
                 getRidgePen = getRidgePen,
                 getlpi = getlpi,
                 putlpi = putlpi,
                 getweightf = getweightf,
                 weights = weights,
                 preinitialize = preinitialize,
                 initialize = initialize,
                 jacobian = jacobian,
                 n_theta = n_theta,
                 mu.eta = stats[[1]]$mu.eta,
                 # postproc=postproc,
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

# pre_fam <- LossStack(preds, square_loss, Weighter, ridgePen = 1)
#
# x <- seq(0, 10, len = 100)
# y <- 2*x + rnorm(100)
# df <- cbind(x, y)
# preds <- cbind(1.8*x, 2.2*x, 2.05*x)
# X <- cbind(1, x)
#
#
# Weighter <- MVN_weights4(matrix(c(1,2,3,4,5,6),nrow=2),dim_num=2)
# attr(Weighter, "init_func")(preds)
#
#
# get_loss_derivs(list_of_beta = list(c(1,2), c(1,2)),
#            list_of_X = list(X, X),
#            theta = c(1,1),
#            tau = c(1),
#            loss = square_loss,
#            weight = Weighter)
#
# ll <- function(pars) {
#   list_of_beta <- list(pars[1:2], pars[3:4])
#   theta <- pars[5:6]
#   tau <- pars[7]
#   get_loss_derivs(list_of_beta, list_of_X, theta, tau, square_loss, Weighter)$ll
# }
#
# numDeriv::grad(ll, c(1,2,1,2,1,1,1))
# numDeriv::hessian(ll, c(1,2,1,2,1,1,1))
