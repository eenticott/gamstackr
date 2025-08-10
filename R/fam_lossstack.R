#' Calculate loss function derivatives for stacking
#'
#' This function calculates the derivatives of the loss function for stacking models.
#' It computes the gradient and Hessian matrix needed for model fitting.
#'
#' @param list_of_beta List of beta parameters
#' @param list_of_X List of design matrices
#' @param theta Theta parameters
#' @param logv Log of the variance parameter
#' @param loss Loss function
#' @param weight Weight function
#' @param preds Matrix of predictions from different models
#' @param y Response variable
#' @param deriv Integer indicating derivative order (0=value, 1=gradient, 2=Hessian)
#'
#' @return List containing log-likelihood (l), gradient (lb), and Hessian (lbb)
#'
get_loss_derivs <- function(list_of_beta, list_of_X, theta, logv, loss, weight, preds, y, deriv) {
  neta <- attr(weight, "neta")
  ntheta <- attr(weight, "ntheta")
  eta <- get_list_of_eta(list_of_X, list_of_beta)
  W <- weight(eta, theta, deriv)
  mu <- Rfast::rowsums(W$f_eval*preds)
  v = exp(logv)
  loss_eval <- loss(y, mu, deriv)
  C_fun <- loss_eval$C
  C_eval <- C_fun(v, deriv)

  ll <- sum(-v * loss_eval$loss - log(C_eval$C))

  if (is.nan(ll)) {
    browser()
  }
  if (is.na(ll)) {
    browser()
  }
  if (is.infinite(ll)) {
    browser()
  }
  grad <- NULL
  hess <- NULL
  if (deriv >= 1) {
    leta <- list()
    if (neta > 0) {
      for (i in 1:neta) {
        leta[[i]] <- -v*Rfast::colsums(loss_eval$l1*Rfast::rowsums(W$f_eta_eval[[i]] * preds)*list_of_X[[i]])
      }
    }

    ltheta <- list()
    if (ntheta > 0) {
      for (i in 1:ntheta) {
        ltheta[[i]] <- -v*(sum(loss_eval$l1*Rfast::rowsums(W$f_theta_eval[[i]] * preds)))
      }
    }

    lv <- v*sum(-loss_eval$loss - C_eval$Cv/C_eval$C)
    grad <- c(unlist(leta), unlist(ltheta), lv)
    if (deriv >= 2) {
      lbetabeta <- NULL
      if (neta > 0) {
        for (i in 1:neta) {
          lbetabeta[[i]] <- list()
          for (j in 1:neta) {
            lbetabeta[[i]][[j]] <- t(list_of_X[[i]]) %*% ((-v*loss_eval$l1*Rfast::rowsums(W$f_eta2_eval[[i]][[j]] * preds)-
                                                             v*loss_eval$l2*Rfast::rowsums(W$f_eta_eval[[i]]*preds)*Rfast::rowsums(W$f_eta_eval[[j]]*preds))*list_of_X[[j]])
          }
          lbetabeta[[i]] <- do.call("cbind", lbetabeta[[i]])
        }
        lbetabeta <- do.call("rbind", lbetabeta)
      }

      lbetatheta <- NULL
      if (neta > 0 & ntheta > 0) {
        for (i in 1:neta) {
          lbetatheta[[i]] <- list()
          for (j in 1:ntheta) {
            lbetatheta[[i]][[j]] <- Rfast::colsums((-v*loss_eval$l1*Rfast::rowsums(W$f_eta_theta_eval[[i]][[j]] * preds)-
                                                      v*loss_eval$l2*Rfast::rowsums(W$f_eta_eval[[i]]*preds)*Rfast::rowsums(W$f_theta_eval[[j]]*preds))*list_of_X[[i]])
          }
          lbetatheta[[i]] <- do.call("cbind", lbetatheta[[i]])
        }
        lbetatheta <- do.call("rbind", lbetatheta)
      }


      lbetav <- NULL
      if (neta > 0) {
        for (i in 1:neta) {
          lbetav[[i]] <- Rfast::colsums(-loss_eval$l1*Rfast::rowsums(W$f_eta_eval[[i]] * preds)*list_of_X[[i]])
        }
        lbetav <- do.call("c", lbetav)
      }

      lthetatheta <- NULL
      if (ntheta > 0) {
        for (i in 1:ntheta) {
          lthetatheta[[i]] <- list()
          for (j in 1:ntheta) {
            lthetatheta[[i]][[j]] <- sum(-v*loss_eval$l1*Rfast::rowsums(W$f_theta2_eval[[i]][[j]] * preds)-
                                           v*loss_eval$l2*Rfast::rowsums(W$f_theta_eval[[i]]*preds)*Rfast::rowsums(W$f_theta_eval[[j]]*preds))
          }
          lthetatheta[[i]] <- do.call("cbind", lthetatheta[[i]])
        }
        lthetatheta <- do.call("rbind", lthetatheta)
      }

      lthetav <- NULL
      if (ntheta > 0) {
        for (i in 1:ntheta) {
          lthetav[[i]] <- -sum(loss_eval$l1*Rfast::rowsums(W$f_theta_eval[[i]] * preds))
        }
        lthetav <- do.call("rbind", lthetav)
      }

      lvv <- nrow(preds) * (-C_eval$Cvv/C_eval$C + C_eval$Cv*C_eval$Cv/C_eval$C^2)

      my_t <- function(x) {
        if (is.null(x)) return(NULL)
        return(t(x))
      }
      # Assemble Hessian block-wise depending on availability of beta/theta
      have_beta <- !is.null(lbetabeta) || !is.null(lbetav) || !is.null(lbetatheta)
      have_theta <- !is.null(lthetatheta) || !is.null(lthetav) || !is.null(lbetatheta)

      if (have_beta && have_theta) {
        hess <- rbind(
          cbind(lbetabeta,           lbetatheta,        v * lbetav),
          cbind(my_t(lbetatheta),    lthetatheta,       v * lthetav),
          cbind(v * my_t(lbetav),    my_t(v * lthetav), lvv * v^2 + lv)
        )
      } else if (have_beta && !have_theta) {
        hess <- rbind(
          cbind(lbetabeta,        v * lbetav),
          cbind(v * my_t(lbetav), lvv * v^2 + lv)
        )
      } else if (!have_beta && have_theta) {
        hess <- rbind(
          cbind(lthetatheta,       v * lthetav),
          cbind(my_t(v * lthetav), lvv * v^2 + lv)
        )
      } else {
        hess <- matrix(lvv * v^2 + lv, nrow = 1, ncol = 1)
      }
    }
  }
  return(list(l = ll, lb = grad, lbb = hess))
}

#' Create a loss stacking family for mgcv
#'
#' This function creates an mgcv family object for loss-based stacking models.
#' It allows stacking multiple prediction models using flexible weight functions.
#'
#' @param preds Matrix of predictions from different models
#' @param loss Loss function to evaluate predictions
#' @param weights Weight function to combine predictions
#' @param RidgePen Ridge penalty parameter for regularization
#'
#' @return An mgcv family object for loss-based stacking
#' @export
#'
LossStack <- function(preds, loss, weights, RidgePen = 1e-5) {
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
  assign(".ridgePen", RidgePen, envir = environment())
  getRidgePen <- function() get(".ridgePen")
  putRidgePen <- function(.x) assign(".ridgePen", .x, envir = environment(sys.function()))

  # lpi
  assign(".lpi", NULL, envir = environment())
  getlpi <- function() get(".lpi")
  putlpi <- function(.x) assign(".lpi", .x, envir = environment(sys.function()))

  # weightf
  assign(".weightf",weights, envir = environment())
  getweightf <- function() get(".weightf")

  # loss
  assign(".loss",loss, envir = environment())
  getLoss <- function() get(".loss")

  assign(".coef", NULL, envir = environment())
  getCoef <- function() get(".coef")
  putCoef <- function(.x) assign(".coef", .x, envir = environment(sys.function()))

  residuals <- function(object, type=c("deviance","pearson","response")) {
    return(as.matrix(object$y)[, 1])
  }

  preinitialize <- function(G) {
    ## G is a gam pre-fit object. Pre-initialize can manipulate some of its
    ## elements, returning a named list of the modified ones.
    nbeta <- ncol(G$X)
    ntheta <- attr(G$family$weights, "ntheta")
    nv <- 1
    atr <- attributes(G$X)
    G$X <- cbind(G$X, matrix(0, nrow(G$X), ntheta + nv)) ## add dummy columns to G$X

    # Ensure lpi is present: default to a single block covering beta columns
    lpi0 <- atr$lpi
    if (is.null(lpi0)) lpi0 <- list(seq_len(nbeta))
    attr(G$X, "lpi") <- lpi0
    attr(G$X, "dim") <- c(nrow(G$X), ncol(G$X))

    putlpi(lpi0)
    if (!is.null(G$Sl))
      attr(G$Sl, "E") <- cbind(attr(G$Sl, "E"),
                               matrix(0, nbeta, ntheta+nv))
    if (ntheta != 0) {
      G$term.names <- c(G$term.names, paste("theta", 1:(ntheta),
                                            sep = "."), "logv")
    }
    ## pad out sqrt of balanced penalty matrix to account for extra params

    list(X=G$X,term.names=G$term.names,family=G$family, Sl = G$Sl)
  } ## preinitialize

  initialize <- expression({
    my_init_fun <- function(y, nobs, E, x, family, offset){
      W_f <- family$getweightf()
      loss <- family$getLoss()

      # Initialize at 0 for now
      coefs <- rep(0, ncol(x))
      lpi <- attr(x,"lpi")
      p <- lapply(lpi, function(lpi_ii) length(lpi_ii))
      neta <- attr(W_f, "neta")
      nbeta <- do.call("sum", p)
      ntheta <- attr(W_f, "ntheta")
      nv <- 1
      preds <- family$getPreds()
      scores <- -loss(y, preds)$loss # Want to score based on loss function, max score should be best model
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
        list_of_X[[jj]] <- x1
        # exp fix for multinomial
        if (attr(family$weights, "name") == "multinomial") {
          y1 <- x1 %*% coefs[lpi[[jj]]]
          coefs[lpi[[jj]]] <- start_coef(x1, e1, log(y1))
          list_of_beta[[jj]] <- coefs[lpi[[jj]]]
          list_of_X[[jj]] <- x1
        }
      }
      # Initial theta
      if (ntheta > 0) {
        coefs[(nbeta+1):(nbeta+ntheta)] <- init_coefs$init_theta
      }

      list_of_eta <- .internals()[["get_list_of_eta"]](list_of_X, list_of_beta)
      theta <- init_coefs$init_theta

      # logscale of residual variance under initial weights
      coefs[nbeta + ntheta + 1] <- log(1/var(Rfast::rowsums(W_f(list_of_eta, theta, deriv = 0)$f_eval * preds)-y))

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
    orig_lpi <- family$getlpi()
    given_lpi <- attr(x,"lpi")
    # Fallback to original lpi if current design matrix lacks lpi attribute
    lpi <- if (is.null(given_lpi)) orig_lpi else given_lpi
    # Guard against any zero/invalid indices (mgcv sometimes carries placeholders)
    if (!is.null(lpi)) lpi <- lapply(lpi, function(ix) ix[ix > 0])
    # Final fallback: single block over all beta columns
    if (is.null(lpi) || length(lpi) == 0) lpi <- list(seq_len(ncol(x)))
    if (!is.null(attr(x, "drop"))) {
      drop_idx <- attr(x, "drop")
      orig_nx <-  Reduce("+",lapply(orig_lpi, function(lpi_ii) length(lpi_ii)))
      cur_nx <- Reduce("+",lapply(lpi, function(lpi_ii) length(lpi_ii)))
      theta_pad_idx <- drop_idx[drop_idx > orig_nx] - orig_nx
      lpi_pad_idx <- theta_pad_idx + cur_nx
      for (idx in lpi_pad_idx) {
        coef <- c(coef[1:(idx-1)], 0, coef[idx:length(coef)])
      }
    }


    p <- lapply(lpi, function(lpi_ii) length(lpi_ii))

    list_of_beta <- list()
    list_of_X <- list()

    for (i in seq_along(lpi)) {
      list_of_beta[[i]] <- coef[lpi[[i]]]
      list_of_X[[i]] <- x[,lpi[[i]], drop=FALSE]
    }
    theta <- NULL
    if (family$n_theta > 0) {
      theta <- coef[(do.call("sum", p)+1):(length(coef)-1)]
    }
    logv <- coef[length(coef)]
    preds <- family$getPreds()
    out <- get_loss_derivs(list_of_beta = list_of_beta,
                           list_of_X = list_of_X,
                           theta = theta,
                           logv = logv,
                           loss = loss,
                           weight = weights,
                           preds = getPreds(),
                           y = y,
                           deriv = deriv)

    l <- out$l
    lb <- out$lb
    lbb <- out$lbb

    nbeta <- do.call("sum", lapply(list_of_beta, length))
    ntheta <- length(theta)

    # theta penalty
    t_pen <- attr(weights, "theta_pen")(theta, deriv = deriv)
    theta_pen <- t_pen$p
    if (is.null(theta_pen)) {
      theta_pen <- 0
    }
    theta_pen_d1 <- t_pen$pt
    theta_pen_d2 <- t_pen$ptt

    beta_pen <- sum(unlist(list_of_beta)^2)
    beta_pen_d1 <- 2 * unlist(list_of_beta)
    beta_pen_d2 <- 2 * diag(nbeta)

    l_pen <- RidgePen * (beta_pen + theta_pen)
    l <- l - l_pen

    if (deriv > 0) {
      lb_pen <- RidgePen * c(beta_pen_d1, theta_pen_d1)
      d_penlist <- c(list(beta_pen_d2, theta_pen_d2))
      lbb_pen <- RidgePen * as.matrix(Matrix::bdiag(d_penlist[!sapply(d_penlist, is.null)]))

      lb[1:length(lb_pen)] <- lb[1:length(lb_pen)] - lb_pen
      lbb[1:length(lb_pen),1:length(lb_pen)] <- lbb[1:length(lb_pen),1:length(lb_pen)] - lbb_pen
    }


    if (exists("lpi_pad_idx")) {
      return(list(l=l, lb = lb[-lpi_pad_idx], lbb = lbb[-lpi_pad_idx, -lpi_pad_idx]))
    }

    return(list(l=l,lb=lb,lbb=lbb))
  }

  jacobian <- function(eta, jj, ...) {
    # Values we need
    # Inner functions
    weight <- getweightf()
    # K
    K <- attr(weight, "num_weights")
    # lpi
    lpi <- getlpi()
    coef <- getCoef()
    neta <- attr(weight, "neta")
    ntheta <- attr(weight, "ntheta")
    theta <- tail(coef, ntheta)
    if ( !is.matrix(eta) ){
      eta <- as.matrix(eta)
    }
    #browser()

    # Consider outer and inner weights separately
    # --------------------------------------------------------------------------
    # Outer weight jacobian (jj in (1:K-1))

    # Special case, only one multinomial weight equal to one
    if (jj == 1 && K == 1){
      return(eta * 0)
    }

    store <- weight(eta, theta, deriv = 1)

    eta_deriv <- sapply(store$f_eta_eval, function(mat) mat[,jj])
    theta_deriv <- sapply(store$f_theta_eval, function(mat) mat[,jj])
    eta_idx <- 1:neta
    theta_idx <- (neta+1):(neta+ntheta)

    return(list(DmuDeta = eta_deriv, DmuDtheta = theta_deriv, eta_idx = eta_idx, theta_idx = theta_idx))
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
    if (attr(family$getweightf(), "ntheta") > 0) {theta <- beta[(do.call("sum", p)+1):(length(beta)-1)]
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
                 getLoss = getLoss,
                 weights = weights,
                 preinitialize = preinitialize,
                 initialize = initialize,
                 jacobian = jacobian,
                 n_theta = n_theta,
                 getCoef = getCoef,
                 putCoef = putCoef,
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
