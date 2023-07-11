TemplateStack <- function(logP, Ridgepen = 1e-5) {
  #

  # Prep
  link <- "identity"
  P <- lapply(logP, exp)
  K <- ncol(logP)

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


    list(X=G$X,term.names=G$term.names,family=G$family, Sl = G$Sl)
  } ## preinitialize

  initialize <- expression({
    my_init_fun <- function(y, nobs, E, x, family, offset){
      startings_coefs <- rep(0, ncol(x))
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

    return(list(l = l, lb = lb, lbb = lbb))
  }

  jacobian <- function(eta, jj, ...) {
    return(NULL)
  }

  predict <- function(family,se=FALSE,eta=NULL,y=NULL,X=NULL,
                      beta=NULL,off=NULL,Vb=NULL)  {
    return(list("fit" = NULL))
  }
}
