log_rowSums_a_times_b <- function(log_a, log_b) {
  log_c <- log_a + log_b
  log_c_maxs <- matrixStats::rowMaxs(log_c)
  log_c_maxs + log(rowSums(exp(log_c - log_c_maxs)))
}

rs_AB <- function(logA, B) {
  # Compute log(A * B) as logA + log(|B|)
  log_A_times_B <- logA + log(abs(B))

  # Compute the element-wise product in logarithm space
  # This is effectively log(A * B)
  exp_log_A_times_B <- exp(log_A_times_B) * sign(B)

  # Compute row sums
  row_sums <- rowSums(exp_log_A_times_B)

  return(row_sums)
}

AoB <- function(A, logB) {
  # Compute log(A/B) as log(|A|) - log(B)
  log_A_over_B <- log(abs(A)) - logB
  exp_log_A_over_B <- exp(log_A_over_B) * sign(A)
  return(exp_log_A_over_B)
}

rs_AB_o_C <- function(logA, B, logC) {
  # Computes rowSums(A*B)/C stable for logA and logC
  return(AoB(rs_AB(logA, B), logC))
}

get_ll_dens_derivs <- function(list_of_beta,
                            list_of_X,
                            theta,
                            weight,
                            log_dens,
                            deriv = 0) {

  K <- ncol(log_dens)
  Ws <- weight
  dens <- exp(log_dens)
  n_eta <- attr(Ws, "neta")
  n_theta <- attr(Ws, "ntheta")

  eta <- get_list_of_eta(list_of_X, list_of_beta)
  W <- Ws(eta, theta, 2)

  llk_deriv <- log_rowSums_a_times_b(log(W$f_eval), log_dens)

  #llk_deriv <- (rowSums(W$f_eval*dens))

  ll <- sum(llk_deriv)

  grad <- NULL
  hess <- NULL

  if (deriv > 0) {
    db <- NULL
    dt <- NULL
    # First derivs
    if (n_eta > 0) {
      db <- (do.call("c", lapply(1:n_eta, function(k) colSums(rs_AB_o_C(log_dens, W$f_eta_eval[[k]], llk_deriv) *
                                                                list_of_X[[k]]))))
    }
    if (n_theta > 0) {
      dt <- (do.call("c", lapply(1:n_theta, function(k) sum(rs_AB_o_C(log_dens, W$f_theta_eval[[k]], llk_deriv)))))
    }

    grad <- c(db, dt)

    if (deriv > 1) {
      dbb <- NULL
      dbt <- NULL
      dtt <- NULL
      if (n_eta > 0) {
        dbb <- list()
        for (i in 1:n_eta) {
          dbb[[i]] <- list()
          for (j in 1:n_eta) {
            dbb[[i]][[j]] <- t(list_of_X[[i]]) %*%
              ((rs_AB_o_C(log_dens, W$f_eta2_eval[[i]][[j]], llk_deriv) -
                  (rs_AB_o_C(log_dens, W$f_eta_eval[[i]], llk_deriv) *
                     rs_AB_o_C(log_dens, W$f_eta_eval[[j]], llk_deriv))) *
                 list_of_X[[j]])
          }
          dbb[[i]] <- do.call("cbind", dbb[[i]])
        }
        dbb <- do.call("rbind", dbb)
        if (n_theta > 0) {
          dbt <- list()
          for (i in 1:n_eta) {
            dbt[[i]] <- list()
            for (j in 1:n_theta) {
              dbt[[i]][[j]] <- colSums((rs_AB_o_C(log_dens, W$f_eta_theta_eval[[i]][[j]], llk_deriv) -
                                          (rs_AB_o_C(log_dens, W$f_eta_eval[[i]], llk_deriv)*
                                             rs_AB_o_C(log_dens, W$f_theta_eval[[j]], llk_deriv))) *
                                         list_of_X[[i]])
            }
            dbt[[i]] <- do.call("cbind", dbt[[i]])
          }
          dbt <- do.call("rbind", dbt)

        }
      }
      if (n_theta > 0) {
        dtt <- list()
        for (i in 1:n_theta) {
          dtt[[i]] <- list()
          for (j in 1:n_theta) {
            dtt[[i]][[j]] <- sum((rs_AB_o_C(log_dens, W$f_theta2_eval[[i]][[j]], llk_deriv) -
                                    (rs_AB_o_C(log_dens, W$f_theta_eval[[i]], llk_deriv) *
                                       rs_AB_o_C(log_dens, W$f_theta_eval[[j]], llk_deriv))))
          }
          dtt[[i]] <- do.call("cbind", dtt[[i]])
        }
        dtt <- do.call("rbind", dtt)
      }
      my_t <- function(x) {
        if (is.null(x)) return(NULL)
        return(t(x))
      }

      hess <- rbind(cbind(dbb, dbt),
                    cbind(my_t(dbt), dtt))
    }
  }
  return(list(l = ll, lb = grad, lbb = hess))
}

DensStack <- function(logP, weight, RidgePen = 1e-5) {
  ### mgcv family for nested stacking
  ### inner_funcs is list of lists of inner weight functions

  ## Prep
  link <- "identity"

  neta <- attr(weight, "neta")
  ntheta <- attr(weight, "ntheta")
  K <- attr(weight, "num_weights")

  P <- lapply(logP, exp)

  # Prepare link functions
  link <- lapply(1:(K - 1 + neta), function(x) "identity")
  stats <- list()
  for (ii in 1:(K - 1 + neta)) {
    stats[[ii]] <- stats::make.link(link[[ii]])
    fam <- structure(list(link=link[[ii]],canonical="none",linkfun=stats[[ii]]$linkfun,
                          mu.eta=stats[[ii]]$mu.eta), class="family")
    fam <- mgcv::fix.family.link(fam)
    stats[[ii]]$d2link <- fam$d2link
    stats[[ii]]$d3link <- fam$d3link
    stats[[ii]]$d4link <- fam$d4link
  }

  ## Save parameters to global env
  # logP
  assign(".logP", logP, envir = environment())
  getLogP <- function() get(".logP")
  putLogP <- function(.x) assign(".logP", .x, envir = environment(sys.function()))

  # P
  assign(".P", P, envir = environment())
  getP <- function() get(".P")
  putP <- function(.x) assign(".P", .x, envir = environment(sys.function()))

  # lpi
  assign(".lpi", NULL, envir = environment())
  getlpi <- function() get(".lpi")
  putlpi <- function(.x) assign(".lpi", .x, envir = environment(sys.function()))

  # weight
  assign(".weight", weight, envir = environment())
  getWeight <- function() get(".weight")
  putWeight <- function(.x) assign(".weight", .x, envir = environment(sys.function()))

  # Ridgepen
  assign(".RidgePen", RidgePen, envir = environment())
  getRidgePen <- function() get(".RidgePen")

  # n_params
  assign(".nparams", ntheta + neta, envir = environment())
  getNparams <- function() get(".nparams")

  # num_weights
  assign(".nweights", K, envir = environment())
  getNweights <- function() get(".nweights")

  # coef
  assign(".coef", NULL, envir = environment())
  getCoef <- function() get(".coef")
  putCoef <- function(.x) assign(".coef", .x, envir = environment(sys.function()))

  # copied residuals from fam_stackProb
  residuals <- function(object, type=c("deviance","pearson","response")) {
    return(as.matrix(object$y)[, 1])
  }

  #######
  preinitialize <- function(G) {
    ## G is a gam pre-fit object. Pre-initialize can manipulate some of its
    ## elements, returning a named list of the modified ones.
    ## A) extends model matrix with dummy
    nbeta <- ncol(G$X)
    atr <- attributes(G$X)
    G$X <- cbind(G$X,matrix(0,nrow(G$X),ntheta)) ## add dummy columns to G$X
    attr(G$X, "lpi") <- atr$lpi
    attr(G$X, "dim") <- c(nrow(G$X), ncol(G$X))
    putlpi(atr$lpi)
    if (!is.null(G$Sl))
      attr(G$Sl, "E") <- cbind(attr(G$Sl, "E"),
                               matrix(0, nbeta, ntheta))
    if (ntheta != 0) {
      G$term.names <- c(G$term.names, paste("my_thetas", 1:ntheta,
                                            sep = "."))
    }
    ## pad out sqrt of balanced penalty matrix to account for extra params

    list(X=G$X,term.names=G$term.names,family=G$family, Sl = G$Sl)
  } ## preinitialize

  # Initialize starting parameters
  initialize <- expression({
    my_init_fun <- function(y, nobs, E, x, family, offset){
      # Get required parameters
      logP <- family$getLogP()
      P <- family$getP()
      weight <- family$getWeight()
      lpi <- attr(x,"lpi")
      p <- lapply(lpi, function(lpi_ii) length(lpi_ii))
      nbeta <- do.call("sum", p)
      neta <- attr(weight, "neta")
      ntheta <- attr(weight, "ntheta")

      coefs <- rep(0, ncol(x))
      use.unscaled <- if (!is.null(attr(E,"use.unscaled"))) TRUE else FALSE
      ridgePen <- family$getRidgePen()

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

      # Get intial eta and theta from weight function
      init_pars <- attr(weight, "init_func")(logP)

      # Initial betas
      list_of_X <- list()
      list_of_beta <- list()
      for (jj in 1:neta) { # Multi eta case
        y <- init_pars$init_mu[,jj]
        x1 <- x[, lpi[[jj]], drop = FALSE]
        e1 <- E[, lpi[[jj]], drop = FALSE]
        coefs[lpi[[jj]]] <- start_coef(x1, e1, y)
        list_of_beta[[jj]] <- coefs[lpi[[jj]]]
        list_of_X[[jj]] <- x1
      }

      # Initial thetas
      coefs[(nbeta+1):(nbeta+ntheta)] <- init_pars$init_theta
      return(coefs)
    }
    if(is.null(start)){
      start <- my_init_fun(y = y, nobs = nobs, E = E, x = x, family = family, offset = offset)
    }
  }) ## initialize

  ll <- function(y,x,coef,wt,family,offset=NULL,deriv=0,
                 dlb=0,d2b=0,Hp=NULL,rank=0,fh=NULL,D=NULL) {
    ## deriv: 0 - eval
    ##        1 - grad and Hess
    ##        2 - diagonal of first deriv of Hess
    given_lpi <- attr(x,"lpi")
    orig_lpi <- getlpi()
    lpi <- given_lpi
    if (deriv == 1) {deriv <- 2}
    # Handle dropped coefficients
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

    if (!identical(given_lpi, lpi)) {
      stop("Mismatch in design matrices.")
    }

    if (is.null(lpi)) {
      stop("Missing design matrix.")
    }

    # Retrieve parameters from family
    weight <- family$getWeight()
    log_dens <- family$getLogP()

    p <- lapply(lpi, function(lpi_ii) length(lpi_ii))

    neta <- attr(weight, "neta")
    ntheta <- attr(weight, "ntheta")
    nbeta = do.call("sum",p)

    list_of_beta <- list()
    list_of_X <- list()

    for (ii in 1:neta) {
      list_of_beta[[ii]] <- coef[lpi[[ii]]]
      list_of_X[[ii]] <- x[, lpi[[ii]], drop = FALSE]
    }

    theta <- coef[(nbeta+1):(nbeta+ntheta)]

    derivs <- get_ll_dens_derivs(list_of_beta,
                              list_of_X,
                              theta,
                              weight,
                              exp(log_dens),
                              deriv = deriv)

    #derivs$llbb <- as.numeric(derivs$llbb) # not sure why it isn't already numeric
    # happened after adding NULLs in for ID inner function

    if (any(is.nan(derivs$lbb))) {
      print(derivs$lbb)}
    # browser()
    #   print(list_of_betaT)
    #   print(list_of_theta)
    # }

    # theta penalty
    t_pen <- attr(weight, "theta_pen")(theta, deriv = deriv)
    theta_pen <- t_pen$p
    theta_pen_d1 <- t_pen$pt
    theta_pen_d2 <- t_pen$ptt


    # beta penalty - Penalise size of beta or eta?
    N <- nrow(log_dens)

    # # beta penalties (penalise squared average eta)
    # x_eta_hat <- list()
    # beta_pens <- list()
    # beta_pens_d1 <- list()
    # beta_pens_d2 <- list()
    # for (ii in 1:neta) {
    #   x_eta_hat[[ii]] <- colSums(list_of_X_eta[[ii]])
    #   beta_pens[[ii]] <- ((x_eta_hat[[ii]] %*% list_of_beta[[ii]]) / N)^2
    #   if (deriv > 0) {
    #     beta_pens_d1[[ii]] <- 2/(N^2) * x_eta_hat[[ii]] %*% t(x_eta_hat[[ii]]) %*% list_of_beta[[ii]]
    #     beta_pens_d2[[ii]] <- 2/(N^2) * x_eta_hat[[ii]] %*% t(x_eta_hat[[ii]])
    #   }
    # }
    # beta_pen <- sum(unlist(beta_pens))
    # if (deriv > 0) {
    #   beta_pen_d1 <- unlist(beta_pens_d1)
    #   beta_pen_d2 <- as.matrix(Matrix::bdiag(c(beta_pens_d2)))
    # }

    # beta penalties - penalise square of beta
    # Consider dividing by standard deviation of X cols
    beta_pen <- sum(unlist(list_of_beta)^2)
    beta_pen_d1 <- 2 * unlist(list_of_beta)
    beta_pen_d2 <- 2 * diag(nbeta)

    l <- derivs$l
    l_pen <- RidgePen * (beta_pen + theta_pen)
    l <- l - l_pen

    lb <- NULL
    lbb <- NULL
    if (deriv > 0) {
      lb <- derivs$lb
      lb_pen <- RidgePen * c(beta_pen_d1, theta_pen_d1)
      lbb <- derivs$lbb
      d_penlist <- c(list(beta_pen_d2, theta_pen_d2))
      lbb_pen <- RidgePen * as.matrix(Matrix::bdiag(d_penlist[!sapply(d_penlist, is.null)]))

      lb <- lb - lb_pen
      lbb <- lbb - lbb_pen
    }

    if (exists("lpi_pad_idx")) {
      return(list(l = l, lb = lb[-lpi_pad_idx], lbb = lbb[-lpi_pad_idx, -lpi_pad_idx]))
    }
    return(list(l = l, lb = lb, lbb = lbb))
  }

  jacobian <- function(eta, jj, ...){
    # Values we need
    # Inner functions
    weight <- getWeight()
    # K
    K <- attr(weight, "num_weights")
    # lpi
    lpi <- getlpi()

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
    if (jj < K) {
      alpha <- cbind(1, exp(eta[,1:(K-1)])) / rowSums(cbind(1, exp(eta[,1:(K-1)])))
      # D alpha / D eta
      DaDe <- sapply(1:(K - 1), function(.kk) {
        alpha[, jj] * (as.numeric(jj == .kk + 1) - alpha[, .kk + 1])
      })
      if(nrow(alpha) == 1) { DaDe <- matrix(DaDe, nrow = 1) }
      return(list(DmuDeta = DaDe, eta_idx = 1:(K-1)))
    }
    # Inner weight jacobian (jj >= K)
    k <- jj - K # get index excluding outer etas
    coefs <- getCoef()
    n_each_weight <- getNweights()
    # get which branch of inner weights we need
    w <- min(which(k <= cumsum(n_each_weight)))
    w_idx <- k - c(0, cumsum(n_each_weight))[w]
    # use for indexing correctly
    each_eta <- c(0,cumsum(n_each_eta))
    each_theta <- c(0,cumsum(n_each_theta))
    ntheta <- do.call("sum", n_each_theta)
    neta <- do.call("sum", n_each_eta)
    ncoef <- length(coefs)

    # Need list_of_etaT[[w]] and list_of_theta[[w]]

    # Get indexes of inner etas (along lpi) and thetas (along coef vectors)
    eta_idx <- (K-1 + each_eta[w] + 1):(K-1 + each_eta[w+1])
    theta_idx <- (ncoef - ntheta + each_theta[w] + 1):(ncoef - ntheta + each_theta[w+1])

    etaT <- eta[,eta_idx]
    theta <- coefs[theta_idx]

    derivs <- .internals()[["eval_deriv"]](inners[[w]], etaT, theta, deriv = 1)

    if (is.list(derivs$f_eta_eval)) {
      eta_deriv <- matrix(nrow = nrow(eta),ncol = length(derivs$f_eta_eval))
      for (i in 1:length(derivs$f_eta_eval)) {
        eta_deriv[,i] <- derivs$f_eta_eval[[i]][,w_idx]
      }
    } else {
      eta_deriv <- as.matrix(derivs$f_eta_eval[,w_idx])
    }

    if (is.list(derivs$f_theta_eval)) {
      theta_deriv <- matrix(nrow = nrow(eta),ncol = length(derivs$f_theta_eval))
      for (i in 1:length(derivs$f_theta_eval)) {
        theta_deriv[,i] <- derivs$f_theta_eval[[i]][,w_idx]
      }
    } else {
      theta_deriv <- as.matrix(derivs$f_theta_eval[,w_idx])
    }

    return(list(DmuDeta = eta_deriv, DmuDtheta = theta_deriv, eta_idx = eta_idx, theta_idx = theta_idx))
  }

  predict <- function(family,se=FALSE,eta=NULL,y=NULL,X=NULL,
                      beta=NULL,off=NULL,Vb=NULL) {
    ## optional function to give predicted values - idea is that
    ## predict.gam(...,type="response") will use this, and that
    ## either eta will be provided, or {X, beta, off, Vb}. family$data
    ## contains any family specific extra information.
    ## if se = FALSE returns one item list containing matrix otherwise
    ## list of two matrices "fit" and "se.fit"...

    weight <- family$getWeight()
    lpi <- family$getlpi()
    K <- attr(weight, "num_weights")
    N <- nrow(X)
    p <- lapply(lpi, function(lpi_ii) length(lpi_ii))
    nbeta <- do.call("sum", p)
    ntheta <- attr(weight, "ntheta")

    # get eta parameters
    list_of_beta <- list()
    list_of_X <- list()

    for (ii in 1:neta) {
      list_of_beta[[ii]] <- beta[lpi[[ii]]]
      list_of_X[[ii]] <- X[, lpi[[ii]]]
    }

    list_of_eta <-.internals()[["get_list_of_eta"]](list_of_X, list_of_beta)

    theta <- beta[(nbeta+1):(nbeta+ntheta)]

    W <- weight(list_of_eta, theta)

    return(list("fit" = W$f_eval))
  }

  structure(list(family = "DensStack",ll = ll,nlp = neta,
                 link = "identity",
                 getLogP = getLogP,
                 getP = getP,
                 getRidgePen = getRidgePen,
                 putLogP = putLogP,
                 putP = putP,
                 getlpi = getlpi,
                 putlpi = putlpi,
                 getNparams = getNparams,
                 getCoef = getCoef,
                 putCoef = putCoef,
                 getWeight = getWeight,
                 n.theta = ntheta,
                 preinitialize = preinitialize,
                 initialize = initialize,
                 jacobian = jacobian,
                 mu.eta = stats[[1]]$mu.eta,
                 # postproc=postproc,
                 tri = mgcv::trind.generator(max(1, neta)), ## symmetric indices for accessing deriv. arrays
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
