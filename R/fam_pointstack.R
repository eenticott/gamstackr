get_point_derivatives <- function(list_of_beta,
                            list_of_X,
                            sigma2,
                            y,
                            predictions,
                            deriv = 1) {
  K <- ncol(predictions)
  mu <- predictions
  sigma2 <- exp(sigma2)
  nbeta <- sapply(list_of_beta, length)
  bigX <- do.call("cbind", list_of_X)

  Xb <- list()
  for (k in 1:K) {
    Xb[[k]] <- list_of_X[[k]]%*%list_of_beta[[k]]
  }

  y_est <- rowSums(sapply(1:K, function(k) exp(Xb[[k]])*mu[,k]))

  ll <- sum(dnorm(x=y, mean=y_est, sd = sqrt(sigma2), log=TRUE))

  ll1 <- NULL
  ll2 <- NULL

  if (deriv >= 1) {
    llk_deriv <- llk_gaussian(y, list(y_est, 1/sqrt(sigma2)), deriv = deriv)
    d1 <- (sapply(1:K, function(k) exp(Xb[[k]])*mu[,k])[,rep(1:K, each = nbeta[k])]*bigX)
    lls <- (llk_deriv$d1[[2]]*sigma2)
    ll1 <- c(colSums(d1*(llk_deriv$d1[[1]])), sum(lls))

    if (deriv >= 2) {
      d2 <- Matrix::bdiag(lapply(1:K, function(k)t(list_of_X[[k]]) %*% ((exp(Xb[[k]])*mu[,k])[,rep(1, nbeta[k])] * -llk_deriv$d1[[1]] * list_of_X[[k]])))
      llbb <- (t(d1) %*% (-llk_deriv$d2[[1]] * (d1))) + d2
      llbs <- -colSums(llk_deriv$d2[[2]] * d1)*sigma2
      llss <- sum(-(llk_deriv$d2[[3]]*sigma2^2) - lls)
      ll2 <- -as.matrix(rbind(cbind(llbb, llbs), cbind(t(llbs), llss)))
    }
  }
  return(list(l = ll, lb = ll1, lbb = ll2))
}

LossStack <- function(preds, ridgePen = 1e-5) {
  # mgcv point forecast stacking family

  # Prep
  link <- "identity"
  K <- ncol(preds)

  nlp <- K # Number of linear predictors
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

  residuals <- function(object, type=c("deviance","pearson","response")) {
    return(as.matrix(object$y)[, 1])
  }

  preinitialize <- function(G) {
    ## G is a gam pre-fit object. Pre-initialize can manipulate some of its
    ## elements, returning a named list of the modified ones.
    nbeta <- ncol(G$X)
    ntheta <- 1
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
      # Initialize at 0 for now
      nbeta <- ncol(x)-1
      ntheta <- 1

      preds <- family$getPreds()
      theta_init <- log(var(rowMeans(preds)-y)) # logscale of error variance with even weights

      starting_coefs <- c(rep(0, nbeta), theta_init)
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
    sig2 <- coef[length(coef)]
    preds <- family$getPreds()
    out <- get_point_derivatives(list_of_beta = list_of_beta,
                          list_of_X = list_of_X,
                          sigma2 = sig2,
                          y = y,
                          predictions = preds,
                          deriv = deriv)

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
    print(eta)

    W <- matrix(nrow = N, ncol = K)
    for (k in 1:K) {
      W[,k] <- exp(list_of_X[[k]]%*%list_of_beta[[k]])
    }
    return(list("fit" = W))
  }

  structure(list(family = "LossStack",ll = ll,nlp = K,
                 link = "identity",
                 getPreds = getPreds,
                 getRidgePen = getRidgePen,
                 getlpi = getlpi,
                 putlpi = putlpi,
                 preinitialize = preinitialize,
                 initialize = initialize,
                 jacobian = jacobian,
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

