# get_loss_derivs <- function(list_of_beta,
#                             list_of_X,
#                             theta,
#                             tau,
#                             loss,
#                             weight,
#                             preds,
#                             y,
#                             deriv) {
#
#   K <- ncol(preds)
#   Ws <- weight
#
#   n_eta <- attr(Ws, "neta")
#   n_theta <- attr(Ws, "ntheta")
#
#   eta <- get_list_of_eta(list_of_X, list_of_beta)
#   W <- Ws(eta, theta, deriv)
#   mu <- Rfast::rowsums(W$f_eval * preds)
#
#   llk_deriv <- loss(y, param = list(mu, tau), deriv = deriv)
#
#   ll <- sum(llk_deriv$d0)
#   grad <- NULL
#   hess <- NULL
#
#   # First derivs
#   if (deriv > 0) {
#     db <- NULL
#     dt <- NULL
#     if (n_eta > 0) {
#       dmudeta <- lapply(1:n_eta, function(k) {Rfast::rowsums(W$f_eta_eval[[k]] * preds)})
#       db <- sapply(1:n_eta,function(k) colSums(dmudeta[[k]] * llk_deriv$d1[[1]] * list_of_X[[k]]))
#     }
#     if (n_theta > 0) {
#       dmudtheta <- lapply(1:n_theta, function(k) {Rfast::rowsums(W$f_theta_eval[[k]] * preds)})
#       dt <- sapply(1:n_theta,function(k) sum(dmudtheta[[k]] * llk_deriv$d1[[1]])) # could make more efficient by combining theta into matrix then multiplying by llk_deriv
#     }
#
#     dtau <- sum(llk_deriv$d1[[2]])
#
#     grad <- c(db, dt, dtau)
#
#     # Second derivs
#     if (deriv > 1) {
#       # beta beta deriv block
#       dbb <- NULL
#       if (n_eta > 0) {
#         dbb <- list()
#         for (k in 1:n_eta) {
#           dbb[[k]] <- list()
#           for (j in 1:n_eta) {
#             dbb[[k]][[j]] <- (t(list_of_X[[k]]) %*% (dmudeta[[k]]*dmudeta[[j]] * list_of_X[[j]] * llk_deriv$d2[[1]])) +
#               t(list_of_X[[k]]) %*% (Rfast::rowsums(W$f_eta2_eval[[k]][[j]] * preds) * llk_deriv$d1[[1]] * list_of_X[[j]])
#           }
#           dbb[[k]] <- do.call("cbind", dbb[[k]])
#         }
#         dbb <- do.call("rbind", dbb)
#       }
#
#       # beta theta deriv block
#       dbt <- NULL
#       if (n_eta > 0 & n_theta > 0) {
#         dbt <- list()
#         for (k in 1:n_eta) {
#           dbt[[k]] <- list()
#           for (j in 1:n_theta) {
#             dbt[[k]][[j]] <- colSums(llk_deriv$d2[[1]] * dmudeta[[k]] * dmudtheta[[j]] * list_of_X[[k]]) +
#               colSums(llk_deriv$d1[[1]] * Rfast::rowsums(W$f_eta_theta_eval[[k]][[j]]*preds) * list_of_X[[k]])
#           }
#           dbt[[k]] <- do.call("cbind", dbt[[k]])
#         }
#         dbt <- do.call("rbind", dbt)
#       }
#
#       # beta tau deriv block
#       dbtau <- NULL
#       if (n_eta > 0) {
#         dbtau <- do.call("c", lapply(1:n_eta, function(k) colSums(llk_deriv$d2[[2]] * dmudeta[[k]] * list_of_X[[k]])))
#       }
#
#       # theta theta deriv block
#       dtt <- NULL
#       if (n_theta > 0) {
#         dtt <- list()
#         for (k in 1:n_theta) {
#           dtt[[k]] <- list()
#           for (j in 1:n_theta) {
#             dtt[[k]][[j]] <- sum(llk_deriv$d1[[1]]*Rfast::rowsums(W$f_theta2_eval[[k]][[j]]*preds) +
#                                    llk_deriv$d2[[1]]*dmudtheta[[k]]*dmudtheta[[j]])
#           }
#           dtt[[k]] <- do.call("cbind", dtt[[k]])
#         }
#         dtt <- do.call("rbind", dtt)
#       }
#
#
#       # theta tau deriv block
#       dttau <- NULL
#       if (n_theta > 0) {
#         dttau <- do.call("c", lapply(1:n_theta, function(k) sum(llk_deriv$d2[[2]] * dmudtheta[[k]])))
#       }
#
#       # tau tau deriv block
#       dtautau <- matrix(sum(llk_deriv$d2[[3]]))
#
#       my_t <- function(x) {
#         if (is.null(x)) return(NULL)
#         return(t(x))
#       }
#
#       hess <- rbind(cbind(dbb, dbt, dbtau),
#                     cbind(my_t(dbt), dtt, dttau),
#                     cbind(my_t(dbtau), my_t(dttau), dtautau))
#     }
#   }
#
#   return(list(l = ll, lb = grad, lbb = hess))
# }
#
# LossStack <- function(preds, loss, weights, RidgePen = 1e-5, fix_scale = FALSE) {
#   # mgcv point forecast stacking family
#
#   # Prep
#   link <- "identity"
#   K <- ncol(preds)
#   n_theta <- attr(weights, "ntheta")
#   nlp <- attr(weights, "neta") # Number of linear predictors
#
#   # Prepare link functions
#   link <- lapply(1:nlp, function(x) "identity")
#   stats <- list()
#   for (ii in 1:nlp) {
#     stats[[ii]] <- stats::make.link(link[[ii]])
#     fam <- structure(list(link=link[[ii]],canonical="none",linkfun=stats[[ii]]$linkfun,
#                           mu.eta=stats[[ii]]$mu.eta), class="family")
#     fam <- mgcv::fix.family.link(fam)
#     stats[[ii]]$d2link <- fam$d2link
#     stats[[ii]]$d3link <- fam$d3link
#     stats[[ii]]$d4link <- fam$d4link
#   }
#
#
#
#   ### Saving extra parameters in .GlobalEnv environment
#   ## Add anything you might want to retrieve/store internally
#   # preds
#   assign(".preds", preds, envir = environment())
#   getPreds <- function() get(".preds")
#   putPreds <- function(.x) assign(".preds", .x, envir = environment(sys.function()))
#
#   # ridgePen
#   assign(".ridgePen", ridgePen, envir = environment())
#   getRidgePen <- function() get(".ridgePen")
#   putRidgePen <- function(.x) assign(".ridgePen", .x, envir = environment(sys.function()))
#
#   # scale
#   assign(".scale", NULL, envir = environment())
#   getScale <- function() get(".scale")
#   putScale <- function(.x) assign(".scale", .x, envir = environment(sys.function()))
#
#   # lpi
#   assign(".lpi", NULL, envir = environment())
#   getlpi <- function() get(".lpi")
#   putlpi <- function(.x) assign(".lpi", .x, envir = environment(sys.function()))
#
#   # weightf
#   assign(".weightf",weights, envir = environment())
#   getweightf <- function() get(".weightf")
#
#   residuals <- function(object, type=c("deviance","pearson","response")) {
#     return(as.matrix(object$y)[, 1])
#   }
#
#   preinitialize <- function(G) {
#     ## G is a gam pre-fit object. Pre-initialize can manipulate some of its
#     ## elements, returning a named list of the modified ones.
#     nbeta <- ncol(G$X)
#     ntheta <- attr(G$family$weights, "ntheta")
#     ntau <- as.numeric(!G$family$fix_scale)
#     atr <- attributes(G$X)
#     G$X <- cbind(G$X,matrix(0,nrow(G$X),ntheta+ntau)) ## add dummy columns to G$X
#
#     attr(G$X, "lpi") <- atr$lpi
#     attr(G$X, "dim") <- c(nrow(G$X), ncol(G$X))
#
#     putlpi(atr$lpi)
#     if (!is.null(G$Sl))
#       attr(G$Sl, "E") <- cbind(attr(G$Sl, "E"),
#                                matrix(0, nbeta, ntheta+ntau))
#     if (ntheta != 0) {
#       G$term.names <- c(G$term.names, paste("theta", 1:(ntheta+ntau),
#                                             sep = "."))
#     }
#     ## pad out sqrt of balanced penalty matrix to account for extra params
#
#     list(X=G$X,term.names=G$term.names,family=G$family, Sl = G$Sl)
#   } ## preinitialize
#
#   initialize <- expression({
#     my_init_fun <- function(y, nobs, E, x, family, offset){
#       W_f <- family$getweightf()
#
#       # Initialize at 0 for now
#       coefs <- rep(0, ncol(x))
#       lpi <- attr(x,"lpi")
#       p <- lapply(lpi, function(lpi_ii) length(lpi_ii))
#       neta <- attr(W_f, "neta")
#       nbeta <- do.call("sum", p)
#       ntheta <- attr(W_f, "ntheta")
#       ntau <- 1
#
#       preds <- family$getPreds()
#       scores <- loss()$d0 # Want to score based on loss function, max score should be best model
#       init_coefs <- attr(W_f, "init_fun")(scores)
#
#       ridgePen <- family$getRidgePen()
#       use.unscaled <- if (!is.null(attr(E,"use.unscaled"))) TRUE else FALSE
#
#       # Penreg func
#       start_coef <- function(x1, e1, y1) {
#         if (!is.null(ridgePen)) {
#           sqrt(ridgePen)
#           e1rp <- matrix(0, nrow = nrow(e1),
#                          ncol = ncol(e1))
#           diag(e1rp) <- sqrt(ridgePen)
#           e1 <- e1 + e1rp
#         }
#         if (use.unscaled) {
#           qrx <- qr(rbind(x1,e1))
#           x1 <- rbind(x1,e1)
#           startji <- qr.coef(qr(x1),c(y1,rep(0,nrow(E))))
#           startji[!is.finite(startji)] <- 0
#         } else startji <- penReg(x1,e1,y1)
#
#         return(startji)
#       }
#
#       list_of_X <- list()
#       list_of_beta <- list()
#
#       for (jj in 1:neta) {
#         y1 <- init_coefs$init_mu[,jj]
#         x1 <- x[, lpi[[jj]], drop = FALSE]
#         e1 <- E[, lpi[[jj]], drop = FALSE]
#
#         coefs[lpi[[jj]]] <- start_coef(x1, e1, y1)
#         list_of_beta[[jj]] <- coefs[lpi[[jj]]]
#         list_of_X[[jj]] <- x1
#         # exp fix for multinomial
#         if (attr(family$weights, "name") == "multinomial") {
#           y1 <- x1 %*% coefs[lpi[[jj]]]
#           coefs[lpi[[jj]]] <- start_coef(x1, e1, log(y1))
#           list_of_beta[[jj]] <- coefs[lpi[[jj]]]
#           list_of_X[[jj]] <- x1
#         }
#       }
#       # Initial theta
#       if (ntheta > 0) {
#         coefs[(nbeta+1):(nbeta+ntheta)] <- init_coefs$init_theta
#       }
#
#       list_of_eta <- .internals()[["get_list_of_eta"]](list_of_X, list_of_beta)
#       theta <- init_coefs$init_theta
#
#       # logscale of residual variance under initial weights
#       coefs[nbeta + ntheta + 1] <- log(var(Rfast::rowsums(W_f(list_of_eta, theta, deriv = 0)$f_eval * preds)-y))
#
#       return(coefs)
#     }
#     if(is.null(start)){
#       start <- my_init_fun(y = y, nobs = nobs, E = E, x = x, family = family, offset = offset)
#     }
#   })
#
#   ll <- function(y,x,coef,wt,family,offset=NULL,deriv=0,
#                  dlb=0,d2b=0,Hp=NULL,rank=0,fh=NULL,D=NULL) {
#     ## deriv: 0 - eval
#     ##        1 - grad and Hess
#     ##        2 - diagonal of first deriv of Hess
#     ##        3 - first deriv of Hess
#
#     if (deriv == 1) {
#       deriv = 2 # My deriv function returns hessian on 2.
#     }
#     lpi <- family$getlpi()
#     p <- lapply(lpi, function(lpi_ii) length(lpi_ii))
#
#     list_of_beta <- list()
#     list_of_X <- list()
#
#     for (i in 1:length(lpi)) {
#       list_of_beta[[i]] <- coef[lpi[[i]]]
#       list_of_X[[i]] <- x[,lpi[[i]], drop=FALSE]
#     }
#     theta <- NULL
#     if (family$n_theta > 0) {
#       theta <- coef[(do.call("sum", p)+1):(length(coef)-1)]
#     }
#     tau <- coef[length(coef)]
#     preds <- family$getPreds()
#     out <- get_loss_derivs(list_of_beta = list_of_beta,
#                            list_of_X = list_of_X,
#                            theta = theta,
#                            tau = tau,
#                            loss = loss,
#                            weight = weights,
#                            preds = getPreds(),
#                            y = y,
#                            deriv = deriv)
#
#     return(out)
#   }
#
#   jacobian <- function(eta, jj, ...) {
#     return(NULL)
#   }
#
#   predict <- function(family,se=FALSE,eta=NULL,y=NULL,X=NULL,
#                       beta=NULL,off=NULL,Vb=NULL)  {
#     lpi <- family$getlpi()
#     p <- lapply(lpi, function(lpi_ii) length(lpi_ii))
#     K <- length(lpi)
#     N <- nrow(X)
#     list_of_beta <- list()
#     list_of_X <- list()
#
#     for (i in 1:K) {
#       list_of_beta[[i]] <- beta[lpi[[i]]]
#       list_of_X[[i]] <- X[,lpi[[i]], drop=FALSE]
#     }
#
#     list_of_eta <- .internals()[["get_list_of_eta"]](list_of_X, list_of_beta)
#     theta <- NULL
#     if (attr(family$getweightf(), "ntheta") > 0) {theta <- beta[(do.call("sum", p)+1):(length(beta)-1)]
# }
#
#     W <- family$getweightf()(list_of_eta, theta, 0)
#
#     return(list("fit" = W$f_eval))
#   }
#
#   structure(list(family = "LossStack",ll = ll,nlp = attr(weights, "neta"),
#                  link = "identity",
#                  getPreds = getPreds,
#                  getRidgePen = getRidgePen,
#                  getlpi = getlpi,
#                  putlpi = putlpi,
#                  getweightf = getweightf,
#                  weights = weights,
#                  preinitialize = preinitialize,
#                  initialize = initialize,
#                  fix_scale = fix_scale,
#                  jacobian = jacobian,
#                  n_theta = n_theta,
#                  mu.eta = stats[[1]]$mu.eta,
#                  # postproc=postproc,
#                  residuals=residuals,
#                  linfo = stats,
#                  #rd=rd,
#                  # dev.resids = dev.resids,
#                  linkinv = stats$linkinv, # MAYBE IT'S NEEDED IN gam.fit5
#                  d2link=1,d3link=1,d4link=1, ## signals to fix.family.link that all done,
#                  predict = predict,
#                  ls=1, ## signals that ls not needed here
#                  available.derivs = 0, ## signal only first derivatives available...
#                  discrete.ok = FALSE
#   ), class = c("general.family", "extended.family","family"))
#
# }
#
# # pre_fam <- LossStack(preds, square_loss, Weighter, ridgePen = 1)
# #
# # x <- seq(0, 10, len = 100)
# # y <- 2*x + rnorm(100)
# # df <- cbind(x, y)
# # preds <- cbind(1.8*x, 2.2*x, 2.05*x)
# # X <- cbind(1, x)
# #
# #
# # Weighter <- MVN_weights4(matrix(c(1,2,3,4,5,6),nrow=2),dim_num=2)
# # attr(Weighter, "init_func")(preds)
# #
# #
# # get_loss_derivs(list_of_beta = list(c(1,2), c(1,2)),
# #            list_of_X = list(X, X),
# #            theta = c(1,1),
# #            tau = c(1),
# #            loss = square_loss,
# #            weight = Weighter)
# #
# # ll <- function(pars) {
# #   list_of_beta <- list(pars[1:2], pars[3:4])
# #   theta <- pars[5:6]
# #   tau <- pars[7]
# #   get_loss_derivs(list_of_beta, list_of_X, theta, tau, square_loss, Weighter)$ll
# # }
# #
# # numDeriv::grad(ll, c(1,2,1,2,1,1,1))
# # numDeriv::hessian(ll, c(1,2,1,2,1,1,1))


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
      lbetabeta <- list()
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

      lbetatheta <- list()
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


      lbetav <- list()
      if (neta > 0) {
        for (i in 1:neta) {
          lbetav[[i]] <- Rfast::colsums(-loss_eval$l1*Rfast::rowsums(W$f_eta_eval[[i]] * preds)*list_of_X[[i]])
        }
        lbetav <- do.call("c", lbetav)
      }

      lthetatheta <- list()
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

      lthetav <- list()
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
      hess <- rbind(cbind(lbetabeta, lbetatheta, v*lbetav),
                    cbind(my_t(lbetatheta), lthetatheta, v*lthetav),
                    cbind(v*my_t(lbetav), my_t(v*lthetav), lvv*v^2 + lv))
    }
  }
  return(list(l = ll, lb = grad, lbb = hess))
}

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
    G$X <- cbind(G$X,matrix(0,nrow(G$X),ntheta+nv)) ## add dummy columns to G$X

    attr(G$X, "lpi") <- atr$lpi
    attr(G$X, "dim") <- c(nrow(G$X), ncol(G$X))

    putlpi(atr$lpi)
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


    p <- lapply(lpi, function(lpi_ii) length(lpi_ii))

    list_of_beta <- list()
    list_of_X <- list()

    for (i in 1:length(lpi)) {
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
