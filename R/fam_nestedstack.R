# TODO: Make more efficient by saving X_eta, X_etaT in preinit
# TODO: Add function that takes outer_then_inner and gives multiplied_weights

get_derivatives <- function(list_of_beta,
                            list_of_betaT,
                            list_of_theta,
                            list_of_X_eta,
                            list_of_X_etaT,
                            list_of_inner_functions,
                            list_of_densities,
                            derivs = 1, elementwise = FALSE) {
  list_of_log_densities <- list_of_densities
  list_of_densities <- lapply(list_of_log_densities, exp)

  N <- nrow(list_of_densities[[1]])
  neta <- unlist(lapply(list_of_inner_functions, function(x) attr(x, "neta")))
  list_of_etaT <- get_list_of_eta(list_of_X_etaT, list_of_betaT)
  list_of_eta <- get_list_of_eta(list_of_X_eta, list_of_beta)
  if (length(list_of_eta) == 0) {
    list_of_eta <- list(matrix(rep(1, N), nrow = N))
    alpha_matrix <- list_of_eta[[1]]
  } else {
    alpha_matrix <- eta_to_alpha(list_of_eta)
  }

  if (derivs == 1) {
    derivs = 2
  } else if (derivs == 0) {
    derivs = 0
  }

  K <-length(list_of_theta)
  # Store evaluated inputs
  eval_store <- list()
  for (k in 1:length(list_of_inner_functions)) {
    eval_store[[k]] <- eval_deriv(list_of_inner_functions[[k]], list_of_etaT[[k]], list_of_theta[[k]],deriv = derivs)
    if (attr(list_of_inner_functions[[k]], "name") == "id") {
      eval_store[[k]]$f_eval <- matrix(1, nrow = nrow(alpha_matrix))
    }
  }

  ll_eval <- ll(get_eval("f_eval", eval_store), alpha_matrix, list_of_log_densities)

  l_calc <- function(pars) {
    list_of_betaT <- list(list(pars[1:(length(pars)-attr(list_of_inner_functions[[k]], "ntheta"))]))
    list_of_theta <- list(pars[(length(pars)-attr(list_of_inner_functions[[k]], "ntheta")+1):length(pars)])
    N <- nrow(list_of_densities[[1]])
    neta <- unlist(lapply(list_of_inner_functions, function(x) attr(x, "neta")))
    list_of_etaT <- get_list_of_eta(list_of_X_etaT, list_of_betaT)
    list_of_eta <- get_list_of_eta(list_of_X_eta, list_of_beta)
    if (length(list_of_eta) == 0) {
      list_of_eta <- list(matrix(rep(1, N), nrow = N))
      alpha_matrix <- list_of_eta[[1]]
    } else {
      alpha_matrix <- eta_to_alpha(list_of_eta)
    }

    if (derivs == 1) {
      derivs = 2
    } else if (derivs == 0) {
      derivs = 0
    }

    K <-length(list_of_theta)
    # Store evaluated inputs
    eval_store <- list()
    for (k in 1:length(list_of_inner_functions)) {
      eval_store[[k]] <- eval_deriv(list_of_inner_functions[[k]], list_of_etaT[[k]], list_of_theta[[k]],deriv = derivs)
      if (attr(list_of_inner_functions[[k]], "name") == "id") {
        eval_store[[k]]$f_eval <- matrix(1, nrow = nrow(alpha_matrix))
      }
    }
    ll_eval <- ll(get_eval("f_eval", eval_store), alpha_matrix, list_of_log_densities)
    return(sum(ll_eval))
  }

  grad <- NULL
  hessian <- NULL

  if (derivs >= 1) {
    # First derivatives
    ll_alpha <- alpha_matrix/exp(ll_eval)
    log_ll_alpha <- log(alpha_matrix) - ll_eval

    ll_etaT_deriv <- ll_etaT(log_ll_alpha, get_eval("f_eta_eval", eval_store), list_of_densities)

    ll_theta_deriv <- unlist(ll_theta(log_ll_alpha, get_eval("f_theta_eval", eval_store), list_of_densities))

    if (K == 1) {
      ll_eta_deriv <- NULL
    } else {
      ll_eta_deriv <- ll_eta(get_eval("f_eval", eval_store), list_of_eta, list_of_densities, alpha_matrix)
    }

    ll_betaT_deriv <- eta_to_beta(ll_etaT_deriv, list_of_X_etaT, elementwise)

    ll_beta_deriv <- eta_to_beta(ll_eta_deriv, list_of_X_eta, elementwise)
    if (elementwise) {
      grad <- cbind(ll_beta_deriv, ll_betaT_deriv)
    } else {
      grad <- c(unlist(ll_beta_deriv), unlist(ll_betaT_deriv), unlist(ll_theta_deriv))
    }

    # Second derivatives
    f_eta2T_eval <- get_eval("f_eta2_eval", eval_store)
    f_theta2_eval <- get_eval("f_theta2_eval", eval_store)
    f_eta_theta_eval <- get_eval("f_eta_theta_eval", eval_store)
    f_theta_eval <- get_eval("f_theta_eval", eval_store)
    f_eta_eval <- get_eval("f_eta_eval", eval_store)

    if (K == 1) {
      ll_eta_theta_eval <- NULL
      ll_eta_etaT_eval <- NULL
      eta2 <- NULL
    } else {
      ll_eta_theta_eval <- ll_eta_theta(get_eval("f_theta_eval", eval_store), get_eval("f_eval", eval_store), list_of_densities, log_ll_alpha, list_of_X_eta)
      ll_eta_etaT_eval <- ll_eta_etaT(ll_etaT_deriv, get_eval("f_eval", eval_store),list_of_densities, log_ll_alpha, list_of_X_eta, list_of_X_etaT)
      eta2 <- ll_eta2(alpha_matrix, list_of_eta, list_of_densities, get_eval("f_eval", eval_store), ll_eta_deriv, list_of_X_eta)
    }


    betaT2_derivs <- list()
    theta2_derivs <- list()
    betaT_theta_derivs <- list()
    beta_theta_derivs <- list()
    beta_betaT_derivs <- list()

    for (k1 in (1:K)[neta != 0]) {
      betaT2_derivs[[k1]] <- list()
      theta2_derivs[[k1]] <- list()
      betaT_theta_derivs[[k1]] <- list()
      beta_theta_derivs[[k1]] <- list()
      beta_betaT_derivs[[k1]] <- list()

      for (k2 in (1:K)[neta != 0]) {
        betaT2_derivs[[k1]][[k2]] <- ll_betaT2(log_ll_alpha,
                                               f_eta2T_eval,
                                               ll_etaT_deriv,
                                               k1, k2,
                                               list_of_densities,
                                               list_of_X_etaT)
        theta2_derivs[[k1]][[k2]] <- ll_theta2(log_ll_alpha, f_theta2_eval, list_of_densities, f_theta_eval, k1, k2)
        betaT_theta_derivs[[k1]][[k2]] <- ll_etaT_theta(log_ll_alpha, f_eta_theta_eval, f_eta_eval, f_theta_eval, list_of_densities, k1, k2, list_of_X_etaT)
      }

      betaT2_derivs[[k1]] <- do.call("cbind", betaT2_derivs[[k1]])
      theta2_derivs[[k1]] <- do.call("cbind", theta2_derivs[[k1]])
      betaT_theta_derivs[[k1]] <- do.call("cbind", betaT_theta_derivs[[k1]])
    }

    for (k1 in (1:(K-1))[(1:(K-1))>0]) {
      if (K == 1) {
        beta_theta_derivs[[k1]] <- NULL
        beta_betaT_derivs[[k1]] <- NULL
      } else {
        a <- do.call("cbind", ll_eta_theta_eval[[k1]])
        if (!is.null(a)){a <- -a}
        b <- do.call("cbind", ll_eta_etaT_eval[[k1]])
        if (!is.null(b)){b <- -b}
        beta_theta_derivs[[k1]] <- a
        beta_betaT_derivs[[k1]] <- b
      }
    }

    betaT2_derivs <- do.call("rbind", betaT2_derivs)
    theta2_derivs <- do.call("rbind", theta2_derivs)
    betaT_theta_derivs <- do.call("rbind", betaT_theta_derivs)
    beta_theta_derivs <- do.call("rbind", beta_theta_derivs)
    beta_betaT_derivs <- do.call("rbind", beta_betaT_derivs)

    if (K==1) {
      beta_betaT_derivs <- NULL
      beta_theta_derivs <- NULL
    }

    t_alt <- function(mat) {
      if (is.null(mat)) {
        return(NULL)
      } else {
        return(t(mat))
      }
    }
    hessian <- rbind(
      cbind(eta2, beta_betaT_derivs, beta_theta_derivs),
      cbind(t_alt(beta_betaT_derivs), betaT2_derivs, betaT_theta_derivs),
      cbind(t_alt(beta_theta_derivs), t_alt(betaT_theta_derivs), theta2_derivs)
    )
    hessian <- matrix(as.numeric(hessian), nrow = nrow(hessian))
  }
  return(list(ll = sum(ll_eval), llb = grad, llbb = hessian))
}

#' mgcv family for nested stacking
#'
#' @param P List of Matrices of model densities
#' @param inner_funcs List of inner function objects
#' @param RidgePen Positive numeric the scalar applied to the ridge penalty in regression
#'
#' @return An initialised mgcv family object
#' @export
#'
#' @examples
#' inners = list(ordinal(3), id())
#' P = list(matrix(rnorm(300), nrow = 100), matrix(rnorm(100), nrow = 100))
#' NestedStack(P, inners)
NestedStack <- function(logP, inner_funcs, RidgePen = 1e-5) {
  ### mgcv family for nested stacking
  ### inner_funcs is list of lists of inner weight functions

  ## Prep
  inners <- inner_funcs
  link <- "identity"
  neta <- do.call("sum",(lapply(inners, function(x) attr(x, "neta"))))
  ntheta <- do.call("sum",(lapply(inners, function(x) attr(x, "ntheta"))))
  n.theta = ntheta
  K <- length(inner_funcs)
  n_each_eta <- lapply(inner_funcs, function(x) attr(x, "neta"))
  n_each_theta <- lapply(inner_funcs, function(x) attr(x, "ntheta"))
  n_each_weight <- lapply(inner_funcs, function(x) attr(x, "num_weights"))
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

  # inner
  assign(".inner", inner_funcs, envir = environment())
  getInner <- function() get(".inner")
  putInner <- function(.x) assign(".inner", .x, envir = environment(sys.function()))

  # Ridgepen
  assign(".RidgePen", RidgePen, envir = environment())
  getRidgePen <- function() get(".RidgePen")

  # n_params
  assign(".nparams", ntheta + neta + K - 1, envir = environment())
  getNparams <- function() get(".nparams")

  # num_weights
  assign(".nweights", n_each_weight, envir = environment())
  getNweights <- function() get(".nweights")

  # coef
  assign(".coef", NULL, envir = environment())
  getCoef <- function() get(".coef")
  putCoef <- function(.x) assign(".coef", .x, envir = environment(sys.function()))

  # MWF
  assign(".MWF", TRUE, envir = environment())
  getMWF <- function() get(".MWF")
  putMWF <- function(.x) assign(".MWF", .x, envir = environment(sys.function()))

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
    inner_f <- G$family$getInner()
    ntheta <- Reduce("+", lapply(inner_f, function(i) attr(i, "ntheta")))
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

  initialize <- expression({
    my_init_fun <- function(y, nobs, E, x, family, offset){
      # Get required parameters
      logP <- family$getLogP()
      P <- family$getP()
      lpi <- attr(x,"lpi")
      p <- lapply(lpi, function(lpi_ii) length(lpi_ii))
      ntheta <- family$n.theta
      inners <- family$getInner()
      neta <- lapply(inners, function(x) attr(x, "neta"))
      ntheta <- lapply(inners, function(x) attr(x, "ntheta"))
      K <- length(inners)
      coefs <- rep(0, ncol(x))
      use.unscaled <- if (!is.null(attr(E,"use.unscaled"))) TRUE else FALSE
      ridgePen <- family$getRidgePen()
      # Find initial eta and theta
      init_pars <- list()
      for (i in 1:K) {
        init_pars[[i]] <- attr(inners[[i]], "init_func")(logP[[i]])
      }

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



      # INITIAL BETATS
      list_of_X_etaT <- list()
      list_of_betaT <- list()
      for (ii in 1:K) {
        ij <- K + c(0, cumsum(neta))[ii]
        n_eta <- neta[[ii]]
        if (n_eta == 0) { # NULL CASE
          list_of_betaT[ii] <- NULL
          list_of_X_etaT[ii] <- NULL
        } else
          if (n_eta == -1) { # Single eta case
            y <- init_pars[[ii]]$init_mu
            x1 <- x[, lpi[[ij]], drop = FALSE]
            e1 <- E[, lpi[[ij]], drop = FALSE]
            coefs[lpi[[ij]]] <- start_coef(x1, e1, y)
            list_of_X_etaT[[ii]] <- x1
            list_of_betaT[[ii]] <- coefs[lpi[[ij]]]
          } else {
            list_of_X_etaT[[ii]] <- list()
            list_of_betaT[[ii]] <- list()
            for (jj in 1:n_eta) { # Multi eta case
              y <- init_pars[[ii]]$init_mu[,jj]
              x1 <- x[, lpi[[ij + jj - 1]], drop = FALSE]
              e1 <- E[, lpi[[ij + jj - 1]], drop = FALSE]
              coefs[lpi[[ij + jj - 1]]] <- start_coef(x1, e1, y)
              list_of_betaT[[ii]][[jj]] <- coefs[lpi[[ij + jj - 1]]]
              list_of_X_etaT[[ii]][[jj]] <- x1
            }
          }
      }

      list_of_etaT <- .internals()[["get_list_of_eta"]](list_of_X_etaT, list_of_betaT)

      # INITIAL THETAS
      list_of_theta <- list()
      ij = do.call("sum",p)
      ti <- c(0, cumsum(ntheta))
      for (ii in 1:K) {
        if (is.null(init_pars[[ii]]$init_theta)) {
          list_of_theta[ii] <- list(NULL)
        } else {
          coefs[(ij + ti[ii] + 1):(ij + ti[ii+1])] <- init_pars[[ii]]$init_theta
          list_of_theta[[ii]] <- init_pars[[ii]]$init_theta
        }
      }

      list_of_densities <- family$getP()
      N <- nrow(list_of_densities[[1]])
      init_dens <- matrix(nrow = N, ncol = K)
      # INITIAL BETAS
      eval_store <- list()
      for (k in 1:length(inners)) {
        f_eval <- .internals()[["eval_deriv"]](inners[[k]], list_of_etaT[[k]], list_of_theta[[k]],deriv = 0)$f_eval
        p <- list_of_densities[[k]]
        if(nrow(f_eval) != nrow(p)) {
          f_eval = matrix(1, nrow = nrow(p))
        }
        init_dens[,k] <- rowSums(f_eval * p)
      }

      y <- max.col(init_dens)

      id_mat <- matrix(0, nrow = N, ncol = K)

      id_mat[matrix(c(1:N, y), ncol = 2)] <- 1

      if (K > 1) {
        for (k in 1:(K-1)) {
          y1 <- id_mat[,k+1]
          x1 <- x[,lpi[[k]], drop = FALSE]
          e1 <- E[,lpi[[k]], drop = FALSE]
          coefs[lpi[[k]]] <- start_coef(x1, e1, y1)
        }
      }


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
    ##        3 - first deriv of Hess

    ## coef: 1 : K-1 coefs for outer models
    ##       K+c(0,cumsum(neta))[i] : K+c(0,cumsum(neta))[i+1] coefs for inner model i
    #print(str(x))
    given_lpi <- attr(x,"lpi")
    orig_lpi <- getlpi()
    lpi <- given_lpi
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
    inners <- family$getInner()
    list_of_densities <- family$getLogP()

    neta <- lapply(inners, function(x) attr(x, "neta"))
    ntheta <- lapply(inners, function(x) attr(x, "ntheta"))

    p <- lapply(lpi, function(lpi_ii) length(lpi_ii))

    list_of_beta <- list()
    list_of_X_eta <- list()

    if (K > 1) {
      for (ii in 1:(K-1)) {
        list_of_beta[[ii]] <- coef[lpi[[ii]]]
        list_of_X_eta[[ii]] <- x[, lpi[[ii]], drop = FALSE]
      }
    }

    list_of_betaT <- list()
    list_of_X_etaT <- list()

    # Need to create nested list structures
    for (ii in 1:K) {
      ij <- K + c(0, cumsum(neta))[ii]
      n_eta <- neta[[ii]]
      if (n_eta == 0) {
        list_of_betaT[ii] <- list(NULL)
        list_of_X_etaT[ii] <- list(NULL)
      } else
        if (n_eta == -1) {
          list_of_betaT[[ii]] <- coef[lpi[[ij]]]
          list_of_X_etaT[[ii]] <- x[, lpi[[ij]], drop = FALSE]
        } else {
          list_of_betaT[[ii]] <- list()
          list_of_X_etaT[[ii]] <- list()
          for (jj in 1:n_eta) {
            list_of_betaT[[ii]][[jj]] <- coef[lpi[[ij + jj - 1]]]
            list_of_X_etaT[[ii]][[jj]] <- x[, lpi[[ij + jj - 1]], drop = FALSE]
          }
        }
    }

    # Investigate how lpi deals with theta
    list_of_theta <- list()
    ij = do.call("sum",p)
    ti <- c(0, cumsum(ntheta))
    for (ii in 1:K) {
      if (ntheta[[ii]] == 0) {
        list_of_theta[ii] <- list(NULL)
      } else {
        list_of_theta[[ii]] <- coef[(ij + ti[ii] + 1):(ij + ti[ii + 1])]
      }
    }

    derivs <- get_derivatives(list_of_beta,
                              list_of_betaT,
                              list_of_theta,
                              list_of_X_eta,
                              list_of_X_etaT,
                              inners,
                              list_of_densities,
                              derivs = deriv)

    derivs$llbb <- as.numeric(derivs$llbb) # not sure why it isn't already numeric
    # happened after adding NULLs in for ID inner function

    if (any(is.nan(derivs$llbb))) {
       print(derivs$llbb)}
      # browser()
    #   print(list_of_betaT)
    #   print(list_of_theta)
    # }

    no_beta <- length(unlist(list_of_beta))
    no_betaT <- length(unlist(list_of_betaT))
    no_theta <- length(unlist(list_of_theta))

    t_pen <- list()
    for (i in 1:K) {
      t_pen[[i]] <- attr(inners[[i]], "theta_pen")(list_of_theta[[i]], deriv = deriv)
    }

    N <- nrow(list_of_densities[[1]])

    # beta penalties (penalise squared average eta)
    if (K > 1) {
      x_eta_hat <- list()
      beta_pens <- list()
      beta_pens_d1 <- list()
      beta_pens_d2 <- list()
      for (ii in 1:(K-1)) {
        x_eta_hat[[ii]] <- colSums(list_of_X_eta[[ii]])
        beta_pens[[ii]] <- ((x_eta_hat[[ii]] %*% list_of_beta[[ii]]) / N)^2

        if (deriv > 0) {
          beta_pens_d1[[ii]] <- 2/(N^2) * x_eta_hat[[ii]] %*% t(x_eta_hat[[ii]]) %*% list_of_beta[[ii]]
          beta_pens_d2[[ii]] <- 2/(N^2) * x_eta_hat[[ii]] %*% t(x_eta_hat[[ii]])
        }
      }

      beta_pen <- sum(unlist(beta_pens))
      if (deriv > 0) {
        beta_pen_d1 <- unlist(beta_pens_d1)
        beta_pen_d2 <- as.matrix(Matrix::bdiag(c(beta_pens_d2)))
      }

    } else {
      beta_pen <- 0
      beta_pen_d1 <- NULL
      beta_pen_d2 <- NULL
    }

    betaT_pen <- sum(unlist(list_of_betaT)^2)
    theta_pen <- sum(unlist(sapply(t_pen, "[[", "p")))

    l <- derivs$ll
    l_pen <- RidgePen * (beta_pen + betaT_pen + theta_pen)
    l <- l - l_pen

    lb <- NULL
    lbb <- NULL
    if (deriv > 0) {
      lb <- derivs$llb
      lb_pen <- RidgePen * c(beta_pen_d1, 2 * unlist(list_of_betaT), unlist(lapply(t_pen, "[[", "pt")))
      lbb <- derivs$llbb
      d_penlist <- c(list(beta_pen_d2, 2 * diag(no_betaT)), lapply(t_pen, "[[", "ptt"))
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
    inners <- getInner()
    # K
    K <- length(inners)
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

    list_of_inner_functions <- family$getInner()

    lpi <- family$getlpi()

    K <- length(list_of_inner_functions)
    N <- nrow(X)
    # get eta parameters
    list_of_beta <- list()
    list_of_X_eta <- list()

    if (K > 1) {
      for (ii in (1:(K-1))) {
        list_of_beta[[ii]] <- beta[lpi[[ii]]]
        list_of_X_eta[[ii]] <- X[, lpi[[ii]]]
      }
    }

    list_of_eta <-.internals()[["get_list_of_eta"]](list_of_X_eta, list_of_beta)
    list_of_betaT <- list()
    list_of_X_etaT <- list()
    ####
    neta <- family$n_each_eta
    # get etaT parameters
    for (ii in 1:K) {
      ij <- K + c(0, cumsum(neta))[ii]
      n_eta <- neta[[ii]]
      if (n_eta == 0) {
        list_of_betaT[ii] <- list(NULL)
        list_of_X_etaT[ii] <- list(NULL)
      } else
        if (n_eta == -1) {
          list_of_betaT[[ii]] <- beta[lpi[[ij]]]
          list_of_X_etaT[[ii]] <- X[, lpi[[ij]], drop = FALSE]
        } else {
          list_of_betaT[[ii]] <- list()
          list_of_X_etaT[[ii]] <- list()
          for (jj in 1:n_eta) {
            list_of_betaT[[ii]][[jj]] <- beta[lpi[[ij + jj - 1]]]
            list_of_X_etaT[[ii]][[jj]] <- X[, lpi[[ij + jj - 1]], drop = FALSE]
          }
        }
    }

    list_of_etaT <- .internals()[["get_list_of_eta"]](list_of_X_etaT, list_of_betaT)
    #####
    # get theta parameters
    p <- lapply(lpi, function(lpi_ii) length(lpi_ii))
    ntheta = family$n_each_theta
    list_of_theta <- list()
    ij = do.call("sum",p)
    ti <- c(0, cumsum(ntheta))
    for (ii in 1:K) {
      if (ntheta[[ii]] == 0) {
        list_of_theta[ii] <- list(NULL)
      } else {

        list_of_theta[[ii]] <- beta[(ij + ti[ii] + 1):(ij + ti[ii + 1])]
      }
    }

    inner_weights <- list()
    for (k in 1:length(list_of_inner_functions)) {
      attr(list_of_inner_functions[[k]], "init_func")(list_of_densities[[k]])
      inner_weights[[k]] <- .internals()[["eval_deriv"]](list_of_inner_functions[[k]], list_of_etaT[[k]], list_of_theta[[k]], deriv = 0)
    }

    outer_weights <- eta_to_alpha(list_of_eta)
    if (K == 1) {
      outer_weights <- matrix(1, nrow = N)
    }

    for (i in 1:length(inner_weights)) {
      inner_weights[[i]] <- inner_weights[[i]]$f_eval
      if (nrow(inner_weights[[i]]) != N) {
        inner_weights[[i]] <- matrix(1, nrow = N)
      }
    }

    # Do we return outer weights then inner weights, or multiplied weights?
    all_inner_weights <- do.call("cbind", inner_weights)

    outer_then_inner <- cbind(outer_weights, all_inner_weights)

    multiplied_weights <- list_times_list(matrix_to_lov(outer_weights), inner_weights)
    multiplied_weights <- do.call("cbind", multiplied_weights)

    MWF <- family$getMWF()
    if (MWF) {
      return(list("fit" = multiplied_weights))
    }

    return(list("fit" = outer_then_inner))
  }

  structure(list(family = "NestedStack",ll = ll,nlp = K - 1 + neta,
                 link = "identity",
                 getLogP = getLogP,
                 getP = getP,
                 getRidgePen = getRidgePen,
                 putLogP = putLogP,
                 putP = putP,
                 getlpi = getlpi,
                 putlpi = putlpi,
                 getInner = getInner,
                 putInner = putInner,
                 getNparams = getNparams,
                 getCoef = getCoef,
                 putCoef = putCoef,
                 getMWF = getMWF,
                 putMWF = putMWF,
                 n.theta = n.theta,
                 n_each_theta = n_each_theta,
                 n_each_eta = n_each_eta,
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
