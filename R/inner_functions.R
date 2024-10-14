# Inner functions
# id
# ordinal
# MVN_weights
# MVN_weights_cpp
# latent
# multinomial
# multinomial


# Functions that can be used for internal model weights and derivatives
rowsums <- Cpp_rowSums

#' Matrix of max(i,j) where i,j are matrix coords
#'
#' @param n Integer: Size of n x n matrix
#'
#' @return Matrix of size n x n
max_index_mat <- function(n) { # used in ordinal penalty
  out = matrix(1:n, ncol = n, nrow = n, byrow = T)
  out[lower.tri(out)] <- t(out)[lower.tri(out)]
  return((n+1) - out)
}


# Inner functions ==============================================================

# Ordinal inner functions ------------------------------------------------------
#' Produces an object containing derivatives for ordered stacking
#'
#' @param K Integer: How many weights in the ordering
#'
#' @return List of functions for finding ordinal weights
#' @export
#'
#' @examples
#' ordinal(3)
ordinal <- function(K) { # nolint
  if (K <= 2) {
    stop("Ordinal weights require at least 3 experts.")
  }

  # Initialise storage
  store <- list()
  force(store)
  assign(".store", NULL, envir = environment())
  getstore <- function() get(".store")
  putstore <- function(.x) assign(".store", .x, envir = environment(sys.function()))

  # Cumulative distribution function
  F <- function(x) {
    1 / (1 + exp(- (x)))
  }

  # weight function at j
  f_j <- function(eta, theta, j) {
    if (j == 1) {
      return(F(theta[1] - eta))
    } else if (j == length(theta) + 1) {
      return(1 - F(theta[j - 1] - eta))
    } else {
      return(F(theta[j] - eta) - F(theta[j - 1] - eta))
    }
  }

  f <- function(eta, theta) {
    return(do.call("cbind",
              lapply(1:(length(theta) + 1), function(x) {f_j(eta, theta, x)}))
      )
  }

  # Computationally stable calc of diff in F
  F_diff <- function(x) {
    abs_x <- abs(x)
    exp_x <- exp(-abs_x)
    exp_2x <- exp_x^2
    exp_x_p1_2 <- (1 + exp_x)^2
    out <- exp_x / exp_x_p1_2
    return(out)
  }

  # First derivative of f_j w.r.t. eta
  f_eta_j <- function(eta, theta, j) {
    if (j == 1) {
      return(-F_diff(theta[1] - eta))
    } else if (j == length(theta) + 1) {
      return(F_diff(theta[j - 1] - eta))
    } else {
      return(F_diff(theta[j - 1] - eta) - F_diff(theta[j] - eta))
    }
  }

  # First derivative of f w.r.t. eta (as vector)
  f_eta <- function(eta, theta) {
    return(do.call("cbind",
            lapply(1:(length(theta) + 1), function(x) {f_eta_j(eta, theta, x)})))

  }

  # First derivative of f w.r.t. theta_n
  f_theta_n <- function(eta, theta, n) {
    store <- matrix(0, length(eta), length(theta) + 1)
    store[, n] <- F_diff(theta[n] - eta)
    store[, n + 1] <- -F_diff(theta[n] - eta)
    return(store)
  }

  # First derivative of f w.r.t. theta
  f_theta <- function(eta, theta) {
    out <- list()

    for (i in seq_len(length(theta))) {
      out[[i]] <- f_theta_n(eta, theta, i)
    }
    return(out)
  }

  # Second derivatives follow similiarly
  F_diff2 <- function(x) {
    abs_x <- abs(x)
    exp_x <- exp(-abs_x)
    exp_2x <- exp_x^2
    exp_3x <- exp_x * exp_2x
    exp_x_p1_3 <- (1 + exp_3x + 3 * exp_2x + 3 * exp_x)
    out <- sign(x) * (exp_2x - exp_x) / exp_x_p1_3
  }

  # for f_theta2 f_thetaj_theta_notj will always be 0
  # still store this as representative structure
  f_theta2_ij <- function(eta, theta, i, j) {
    store <- matrix(0, length(eta), length(theta) + 1)
    if (i != j) {
      return(store)
    }
    store[, i] <- F_diff2(theta[i] - eta)
    store[, i + 1] <- -F_diff2(theta[i] - eta)
    return(store)
  }

  f_theta2 <- function(eta, theta) {
    out <- list()
    for (i in seq_len(length(theta))) {
      out[[i]] <- list()
      for (j in seq_len(length(theta))) {
        out[[i]][[j]] <- f_theta2_ij(eta, theta, i, j)
      }
    }
    return(out)
  }

  f_eta2_j <- function(eta, theta, j) {
    if (j == 1) {
      return(F_diff2(theta[1] - eta))
    } else if (j == length(theta) + 1) {
      return(-F_diff2(theta[j - 1] - eta))
    } else {
      return(F_diff2(theta[j] - eta) - F_diff2(theta[j - 1] - eta))
    }
  }

  f_eta2 <- function(eta, theta) {
    return(do.call("cbind",
            lapply(1:(length(theta) + 1), function(x) {f_eta2_j(eta, theta, x)})))
  }

  f_eta_theta_n <- function(eta, theta, n) {
    store <- matrix(0, length(eta), length(theta) + 1)
    store[, n] <- -F_diff2(theta[n] - eta)
    store[, n + 1] <- F_diff2(theta[n] - eta)
    return(store)
  }

  f_eta_theta <- function(eta, theta) {
    out <- list()
    for (i in seq_len(length(theta))) {
      out[[i]] <- f_eta_theta_n(eta, theta, i)
    }
    return(out)
  }

  get_derivs <- function(eta, tau, deriv) {
    # eta - list or vector of eta
    # tau - vector of extra parameters (on untransformed scale)
    # deriv - integer 0, 1, 2, what level of derivatives needed
    if (is.list(eta)) {
      eta <- do.call("cbind", eta)
    }
    store <- list()
    # transform tau
    tau <- c(stats::qlogis(1/K), tau)
    theta <- c(tau[1], tau[1] + cumsum(exp(tau[-1])))
    store$f_eval <- f(eta, theta)

    if (deriv >= 1) {
      # first derivatives
      store$f_eta_eval <- list(f_eta(eta, theta))

      # Take cumulative sum of reverse order, then reverse so sum at front and just last term at end
      store$f_theta_eval <- list_by_vector(rev(mat_cumsum(rev(f_theta(eta, theta)))), c(tau[1], exp(tau[-1])))[-1]
      if (K == 2) {
        store$f_theta_eval <- NULL

      }
      if (deriv >= 2) {
        # second derivatives
        store$f_eta2_eval <- list(list(f_eta2(eta, theta)))
        store$f_eta_theta_eval <- list(list_by_vector(rev(mat_cumsum(rev(f_eta_theta(eta, theta)))),
                                                 c(tau[1], exp(tau[-1])))[-1])
        if (K > 2) {
          f_theta2_orig <- f_theta2(eta, theta)
          f_theta2_d <- get_diag(f_theta2_orig)
          sums <- rev(mat_cumsum(rev(f_theta2_d[-1])))
          f_theta2_eval <- list()
          exp_taus <- exp(tau[-1])
          d <- length(theta)
          for (i in 2:d) {
            f_theta2_eval[[i - 1]] <- list()
            for (j in 2:d) {
              f_theta2_eval[[i - 1]][[j - 1]] <- exp_taus[i - 1] * exp_taus[j - 1] *
                sums[[intersect((i:d) - 1, (j:d) - 1)[1]]]
              if (i == j) {
                f_theta2_eval[[i - 1]][[j - 1]] <- f_theta2_eval[[i - 1]][[j - 1]] + store$f_theta_eval[[i - 1]]
              }

            }
          }
        } else {
          f_theta2_eval <- list(NULL)
        }

        store$f_theta2_eval <- f_theta2_eval
      }
    }

    return(store)
  }

  init_func <- function(densities) {
    # initialise thetas at equal spacing
    # eta intialised to centre between thresholds of best expert
    y1 <- max.col(densities)
    init_tau <- log(diff(stats::qlogis(seq(1, K - 1) / K)))
    fixed_tau <- stats::qlogis(1 / K)
    tau <- c(fixed_tau, init_tau, 1)
    theta <- c(fixed_tau - 1, tau[1], tau[1] + cumsum(exp(tau[-1])))
    mustart <- (theta[y1 + 1] + theta[y1]) / 2
    if (K == 2) {
      init_tau <- NULL
    }
    return(list(init_mu = matrix(mustart), init_theta = init_tau))
  }

  pen <- function(tau, deriv = 0) {
    # penalise the distance between thresholds
    n <- length(tau) + 2
    theta <- c(stats::qlogis(1 / n), stats::qlogis(1 / n) + cumsum(exp(tau)))

    pen <- sum((theta - stats::qlogis(seq(1, n - 1) / n)) ^ 2)
    pen_grad <- NULL
    pen_hess <- NULL

    if (deriv > 0) {
      pen_grad <- exp(tau) * 2 * rev(cumsum(rev(theta[-1] - stats::qlogis(seq(2, n - 1) / n))))

      pen_hess <- outer(exp(tau), exp(tau)) * 2 * max_index_mat(n - 2) + diag(pen_grad,)
    }
    return(list(p = pen, pt = pen_grad, ptt = pen_hess))
  }
  name = "ordinal"
  return(structure(get_derivs,
                   init_func = init_func,
                   theta_pen = pen,
                   neta = 1,
                   ntheta = K - 2,
                   num_weights = K,
                   name = "ordinal"))
}

# Multivariate normal inner functions ------------------------------------------
#' Produces an object containing derivatives for multivariate normal stacking
#'
#' @param x Array of (dim_num, n_k) specifying locations for each of the models
#' @param dim_num Integer specifying number of dimensions required
#'
#' @return Object containing derivatives for multivariate normal stacking
#' @export
#'
#' @examples
#' x <- matrix(c(1:6), nrow = 3)
#' MVN_weights(x)
MVN_weights <- function(ex_coords) {
  # Input
  # -----
  # ex_coords - Array of (dim_num, K) specifying coordinates of each expert
  x <- ex_coords
  force(x)
  n_k <- ncol(x)
  dim_num <- nrow(x)

  x <- t(x)

  get_derivs <- function(eta, theta, deriv = 1) {
    # eta - list or vector of eta
    # tau - vector of extra parameters (on real line)
    # deriv - integer 0, 1, 2, what level of derivatives needed
    if (is.list(eta)) {eta <- do.call("cbind", eta)}
    tau <- exp(theta)
    dens_matrix <- matrix(nrow = nrow(eta) * n_k, ncol = dim_num)

    # find density of each coordinate in each dimension
    for (n in 1:dim_num) {
      dens_matrix[, n] <-  llk_gaussian(rep(x[, n], each = nrow(eta)),
                                        param = list(rep(eta[, n], n_k), 1 / sqrt(tau[n])),
                                        deriv = 0)$d0
    }

    # computationally stable calc of exp(dens)/rowsums(dens)
    log_dens <- matrix(rowsums(dens_matrix), ncol = n_k)
    shift_mat <- log_dens - matrixStats::rowMaxs(log_dens)

    store <- list()

    expshift <- exp(shift_mat)
    rs_expshift <- rowsums(expshift)

    store$f_eval <- expshift / rs_expshift

    # derivatives follow from formulas
    if (deriv >= 1) {
      dens_list <- list()
      for (n in 1:dim_num) {
        dens_list[[n]] <-  llk_gaussian(rep(x[, n], each = nrow(eta)),
                                        param = list(rep(eta[, n], n_k), 1 / sqrt(tau[n])),
                                        deriv = 2)
      }

      detas <- sapply(dens_list, "[[", "d1")[1, ]

      dtaus <- sapply(dens_list, "[[", "d1")[2, ]

      f_eta_out <- list()
      f_tau_out <- list()

      for (i in 1:dim_num) {
        l1 <- matrix(detas[[i]], ncol = n_k)
        l2 <- matrix(dtaus[[i]], ncol = n_k)

        f_eta_out[[i]] <- store$f_eval * (l1 - (rowsums(l1 * expshift)) / rs_expshift)
        f_tau_out[[i]] <- store$f_eval * (l2 - (rowsums(l2 * expshift)) / rs_expshift)
      }
      # remember to delete eval_dens from func defs
      store$f_eta_eval <- f_eta_out
      store$f_tau_eval <- f_tau_out
      store$f_theta_eval <- list_by_vector(f_tau_out, tau)

      if (deriv >= 2) {
        out <- list()

        for (alpha in 1:dim_num) {
          out[[alpha]] <- list()
          dalpha2 <- matrix(dens_list[[alpha]]$d2[[1]], ncol = n_k)
          dalpha1 <- matrix(dens_list[[alpha]]$d1[[1]], ncol = n_k)

          for (beta in 1:dim_num) {
            dbeta1 <- matrix(dens_list[[beta]]$d1[[1]], ncol = n_k)
            if (alpha == beta) {
              d2 <- (dalpha2 + dalpha1**2)
            } else {
              d2 <- dalpha1 * dbeta1
            }
            p1 <- d2 * store$f_eval
            p2 <- - rowsums(dbeta1 * expshift) * store$f_eval * dalpha1 / rs_expshift
            p3 <- - (rowsums(d2 * expshift) * store$f_eval) / rs_expshift
            p35 <- - rowsums(dalpha1 * expshift) * store$f_eval * dbeta1 / rs_expshift
            p4 <- rowsums(dbeta1 * expshift) * rowsums(dalpha1 * expshift) * 2 * store$f_eval / rs_expshift^2
            out[[alpha]][[beta]] <- p1 + p2 + p3 + p4 + p35
          }
        }
        store$f_eta2_eval <- out

        out <- list()
        for (alpha in 1:dim_num) {
          out[[alpha]] <- list()
          dalpha2 <- matrix(dens_list[[alpha]]$d2[[2]], ncol = n_k)
          dalpha1 <- matrix(dens_list[[alpha]]$d1[[1]], ncol = n_k)

          for (beta in 1:dim_num) {
            dbeta1 <- matrix(dens_list[[beta]]$d1[[2]], ncol = n_k)
            if (alpha == beta) {
              d2 <- (dalpha2 + (dalpha1 * dbeta1))
            } else {
              d2 <- dalpha1 * dbeta1
            }

            p1 <- d2 * store$f_eval
            p2 <- - rowsums(dbeta1 * expshift) * store$f_eval * dalpha1 / rs_expshift
            p3 <- - (rowsums(d2 * expshift) * store$f_eval) / rs_expshift
            p35 <- - rowsums(dalpha1 * expshift) * store$f_eval * dbeta1 / rs_expshift
            p4 <- rowsums(dbeta1 * expshift) * rowsums(dalpha1 * expshift) * 2 * store$f_eval / rs_expshift^2
            out[[alpha]][[beta]] <- (p1 + p2 + p3 + p4 + p35) * tau[beta]
          }
        }

        store$f_eta_theta_eval <- out

        out <- list()
        for (alpha in 1:dim_num) {
          out[[alpha]] <- list()
          dalpha2 <- matrix(dens_list[[alpha]]$d2[[3]], ncol = n_k)
          dalpha1 <- matrix(dens_list[[alpha]]$d1[[2]], ncol = n_k)

          for (beta in 1:dim_num) {
            dbeta1 <- matrix(dens_list[[beta]]$d1[[2]], ncol = n_k)
            if (alpha == beta) {
              d2 <- (dalpha2 + dalpha1**2)
            } else {
              d2 <-  dalpha1 * dbeta1
            }
            p1 <- d2 * store$f_eval
            p2 <- - rowsums(dbeta1 * expshift) * store$f_eval * dalpha1 / rs_expshift
            p3 <- - (rowsums(d2 * expshift) * store$f_eval) / rs_expshift
            p35 <- - rowsums(dalpha1 * expshift) * store$f_eval * dbeta1 / rs_expshift
            p4 <- rowsums(dbeta1 * expshift) * rowsums(dalpha1 * expshift) * 2 * store$f_eval / rs_expshift^2

            if (alpha == beta) {
              out[[alpha]][[beta]] = (p1 + p2 + p3 + p4 + p35) * tau[alpha]^2 + tau[alpha] * store$f_tau_eval[[alpha]]
            } else {
              out[[alpha]][[beta]] = (p1 + p2 + p3 + p4 + p35) * tau[alpha] * tau[beta]
            }
          }
        }

        store$f_theta2_eval <- out

      }
    }
    return(store)
  }

  init_func <- function(densities) {
    # initialise theta at variance in each dimension
    # initialise eta at best expert each point
    theta <- (apply(t(x), 1, (stats::var)))
    y1 <- max.col(densities)
    mustart <- (x[y1, ,drop=FALSE])
    return(list(init_theta = log(theta), init_mu = mustart))
  }

  pen <- function(tau, deriv = 0) {
    # penalise size of theta to avoid flat densities
    theta <- exp(tau)
    v <- (apply(t(x), 1, (stats::var)))
    pen <- sum((theta - v)^2)
    pen_grad <- NULL
    pen_hess <- NULL

    if (deriv > 0) {
      pen_grad <- theta * 2 * (theta - v)
      if (length(pen_grad) == 1) {
        pen_hess <- matrix((4 * theta - 2 * v) * theta)
      } else {
        pen_hess <- diag((4 * theta - 2 * v) * theta)
      }
    }

    return(list(p = pen, pt = pen_grad, ptt = pen_hess))
  }

  name = "MVN"
  return(structure(
    get_derivs,
    init_func = init_func,
    arg_list = list("x" = x, dim_num = dim_num),
    neta = dim_num,
    ntheta = dim_num,
    theta_pen = pen,
    num_weights = n_k,
    name = name
  ))
}

# Multivariate normal inner functions cpp --------------------------------------
#' Produces an object containing derivatives for multivariate normal stacking
#' Uses cpp for faster computation
#'
#' @param ex_coords Array of (dim_num, K) specifying locations for each of the experts, dim_num is dimension of expert space
#'
#' @return Object containing derivatives for multivariate normal stacking
#' @export
#'
#' @examples
#' x <- matrix(c(1:6), nrow = 3)
#' MVN_weights(x, 3)
MVN_weights_cpp <- function(ex_coords) {
  # Input
  # -----
  # ex_coords - Array of (dim_num, K) specifying locations for each of the models
  # Suppose we have N data points
  # n_k inner functions
  # inner functions f_j depend on:
  #   etaT_1,..., etaT_p
  #   etaT denotes tilde(eta)
  #   theta_1,..., theta_q
  x <- ex_coords
  force(x)
  n_k <- ncol(x)
  dim_num <- nrow(x)
  x <- t(x)
  store <- list()
  force(store)
  assign(".store", NULL, envir = environment())
  getstore <- function() get(".store")
  putstore <- function(.x) assign(".store", .x, envir = environment(sys.function()))

  assign(".dens_mat", NULL, envir = environment())
  getdens <- function() get(".dens_mat")
  putdens <- function(.x) assign(".dens_mat", .x, envir = environment(sys.function()))

  assign(".n_k", n_k, envir = environment())
  getn_k <- function() get(".n_k")

  assign(".x", x, envir = environment())
  getx <- function() get(".x")

  assign(".dim_num", dim_num, envir = environment())
  getdim_num <- function() get(".dim_num")

  # For memory efficiency it is helpful to save list of matrices of correct size
  # This initialises with correct size, also gets called in predict
  init_store <- function(N, K) {
    dim_num <- getdim_num()
    store <- list()
    store$f_eval <- matrix(nrow = N, ncol = K)
    store$f_eta_eval <- list()
    store$f_theta_eval <- list()
    store$f_tau_eval <- list()
    store$f_eta2_eval <- list()
    store$f_theta2_eval <- list()
    for (i in 1:dim_num) {
      store$f_eta_eval[[i]] <- matrix(nrow = N, ncol = K)
      store$f_theta_eval[[i]] <- matrix(nrow = N, ncol = K)
      store$f_tau_eval[[i]] <- matrix(nrow = N, ncol = K)
      store$f_eta2_eval[[i]] <- lapply(1:dim_num, function(x) matrix(nrow = N, ncol = K))
      store$f_theta2_eval[[i]] <- lapply(1:dim_num, function(x) matrix(nrow = N, ncol = K))
      store$f_eta_theta_eval[[i]] <- lapply(1:dim_num, function(x) matrix(nrow = N, ncol = K))
    }
    putstore(store)
     # Init dens storage
    putdens(matrix(nrow=N*K, ncol=dim_num))
  }

  get_derivs <- function(eta, theta, deriv = 1) {
    # eta - list or vector of eta
    # tau - vector of extra parameters (on real line)
    # deriv - integer 0, 1, 2, what level of derivatives needed
    if (is.list(eta)) {eta <- do.call("cbind", eta)}
    store <- getstore()
    K <- getn_k()
    dim_num <- getdim_num()
    N <- nrow(eta)
    x <- getx()
    # If store is null or not correct structure need to re-initialise
    # init_store calls putstore with correct object
    if (is.null(store)) {
      init_store(N, K)
      store <- getstore()
    } else if (nrow(eta) != nrow(store$f_eval)) {
      init_store(N, K)
      store <- getstore()
    }
    get_derivs_cpp(eta, theta, deriv, K, dim_num, x, getdens(), store)
    return(store)
  }

  init_func <- function(densities) {
    x <- getx()
    theta <- (apply(t(x), 1, (stats::var)))
    y1 <- max.col(densities)
    mustart <- (x[y1, ,drop=FALSE])
    N <- nrow(densities)
    dim_num <- getdim_num()
    n_k <- getn_k()
    return(list(init_theta = log(theta), init_mu = mustart))
  }

  pen <- function(tau, deriv = 0) {
    theta <- exp(tau)
    v <- (apply(t(x), 1, (stats::var)))
    pen <- sum((theta - v)^2)
    pen_grad <- NULL
    pen_hess <- NULL

    if (deriv > 0) {
      pen_grad <- theta * 2 * (theta - v)
      if (length(pen_grad) == 1) {
        pen_hess <- matrix((4 * theta - 2 * v) * theta)
      } else {
        pen_hess <- diag((4 * theta - 2 * v) * theta)
      }
    }

    return(list(p = pen, pt = pen_grad, ptt = pen_hess))
  }

  name = "MVN"
  return(structure(
    get_derivs,
    init_func = init_func,
    arg_list = list("x" = x, dim_num = dim_num),
    neta = dim_num,
    ntheta = dim_num,
    theta_pen = pen,
    num_weights = n_k,
    name = name,
    putdens = putdens,
    putstore = putstore
  ))
}


# ========================================
#' Produces equal weights without any dependence on parameters
#'
#' @return Object of derivatives for equally weighted models
#' @export
#'
#' @examples
#' id()
id <- function() {
  # Suppose we have N data points
  # n_k inner functions
  # inner functions f_j depend on:
  #   etaT_1,..., etaT_p
  #   etaT denotes tilde(eta)
  #   theta_1,..., theta_q
  get_derivs <- function(eta, theta, deriv) {
    store <- list()
    if (deriv >= 0) {
      store$f_eval <- matrix(1)
      if (deriv >= 1) {
        zero_mat <- NULL

        store$f_eta_eval <- NULL

        store$f_theta_eval <- NULL

        if (deriv >= 2) {
          store$f_eta2_eval <- NULL

          store$f_eta_theta_eval <- NULL

          store$f_theta2_eval <- NULL
        }
      }
    }
    return(store)
  }

  pen <- function(tau, deriv = 0) {
    pen <- NULL
    pen_grad <- NULL
    pen_hess <- NULL

    return(list(p = pen, pt = pen_grad, ptt = pen_hess))
  }

  init_func <- function(densities) {
    return(list(init_theta = NULL, init_mu = NULL))
  }

  name <- "id"

  return(structure(get_derivs, neta = 0, ntheta = 0, theta_pen = pen, init_func = init_func, num_weights = 1, name = name))
}

# ========================================
#' Multinomial weights
#'
#' @param K - Integer specifying how many weights are used
#' @param normalised - Boolean, should weights be normalised
#'
#' @return - Weight object
#' @export
#'
#' @examples
multinomial <- function(K, normalised = TRUE) {
  # If normalised we don't need K eta
  neta <- K - (1*normalised)

  get_derivs <- function(eta, theta, deriv) {
    # eta - list or vector of eta
    # tau - vector of extra parameters (on real line)
    # deriv - integer 0, 1, 2, what level of derivatives needed
    store <- list()
    if (is.list(eta)) {
      eta <- do.call("cbind", eta)
    }
    if (normalised) {
      eta <- cbind(0, eta)
      rs_eta <- rowSums(exp(eta))
      exp_eta <- exp(eta)
      if (deriv >= 0) {
        store$f_eval <- exp_eta/rs_eta
        if (deriv >= 1) {
          store$f_eta_eval <- list()
          for (i in 2:K) {
            f_eta_eval <- matrix(-exp_eta*exp_eta[,i]/rs_eta^2, nrow = nrow(eta), ncol = K)
            f_eta_eval[,i] <-  f_eta_eval[,i] + exp_eta[,i]/rs_eta
            store$f_eta_eval[[i-1]] <- f_eta_eval
          }
          store$f_theta_eval <- NULL

          if (deriv >= 2) {
            store$f_eta2_eval <- list()
            for (i in 2:K) {
              dki <- matrix(1:K == i, nrow = nrow(eta), ncol = K, byrow = TRUE)
              store$f_eta2_eval[[i-1]] <- list()
              for (j in 2:K) {
                store$f_eta2_eval[[i-1]][[j-1]] <- store$f_eta_eval[[j-1]]*dki -
                  store$f_eta_eval[[j-1]]*store$f_eval[,i] -
                  store$f_eta_eval[[j-1]][,i] * store$f_eval
              }
            }
            store$f_eta_theta_eval <- NULL
            store$f_theta2_eval <- NULL
          }
        }
      }
      return(store)
    } else {
      if (deriv >= 0) {
        store$f_eval <- exp(eta)
        if (deriv >= 1) {
          store$f_eta_eval <- list()
          for (i in 1:neta) {
            store$f_eta_eval[[i]] <- matrix(0, nrow = nrow(eta), ncol = neta)
            store$f_eta_eval[[i]][,i] <- exp(eta[,i])
          }
          store$f_theta_eval <- NULL

          if (deriv >= 2) {
            store$f_eta2_eval <- list()
            for (i in 1:neta) {
              store$f_eta2_eval[[i]] <- list()
              for (j in 1:neta) {
                store$f_eta2_eval[[i]][[j]] <- matrix(0, nrow = nrow(eta), ncol = neta)
                if (i == j) {
                  store$f_eta2_eval[[i]][[j]][,j] <- exp(eta[,j])
                }
              }
            }

            store$f_eta_theta_eval <- NULL
            store$f_theta2_eval <- NULL
          }
        }
      }
      return(store)
    }
  }

  pen <- function(tau, deriv = 0) {
    pen <- NULL
    pen_grad <- NULL
    pen_hess <- NULL

    return(list(p = pen, pt = pen_grad, ptt = pen_hess))
  }

  init_func <- function(scores) {
    # Initialise with binary score indicating best expert
    y1 <- max.col(scores)
    x <- diag(ncol(scores))
    mustart <- x[y1,,drop=FALSE]
    #mustart <- (mustart + 0.0000000001)/rowSums(mustart + 0.0000000001) # add pertubation for stability when taking log
    #mustart <- log(mustart)
    return(list(init_theta = NULL, init_mu = mustart))
  }

  name <- "multinomial"

  return(structure(get_derivs,
                   neta = neta,
                   ntheta = 0,
                   theta_pen = pen,
                   init_func = init_func,
                   num_weights = K, name = name))
}


# =================================================
# Need to deal with fringe cases, 0 eta, 0 theta in inner and outer and both
#' Nested weights, weights a set of inner weighting functions by function given
#' in outer weight
#'
#' @param outer_weight - Weighting function to multiply inner weights
#' @param inner_weight - List of inner weighting functions
#'
#' @return A weighting object that multiplies set of outer weights by inner weights
#' @export
#'
#' @examples
nested <- function(outer_weight, inner_weight) {
  if (attr(outer_weight, "num_weights") != length(inner_weight)) {
    stop("Number of outer weights should match number of inner functions.")
  }
  K <- length(inner_weight)
  n_k <- sapply(inner_weight, function(w_f) attr(w_f, "num_weights"))
  n_weights <- sum(n_k)

  n_outer_eta <- attr(outer_weight, "neta")
  n_outer_theta <- attr(outer_weight, "ntheta")

  n_inner_eta <- sapply(inner_weight, function(w_f) attr(w_f, "neta"))
  n_inner_theta <- sapply(inner_weight, function(w_f) attr(w_f, "ntheta"))

  total_eta <- sum(n_inner_eta) + n_outer_eta
  total_theta <- sum(n_inner_theta) + n_outer_theta

  cumsum_eta <- c(0, cumsum(n_inner_eta))
  cumsum_theta <- c(0, cumsum(n_inner_theta))

  # Create flag variables for outer, inner parameters
  have_outer_eta <- n_outer_eta > 0
  have_outer_theta <- n_outer_theta > 0
  have_inner_eta <- sum(n_inner_eta) > 0
  have_inner_theta <- sum(n_inner_theta) > 0

  outer_set <- rep(1:K, times = n_k)

  eta_idx <- rep(1:(K+1), times = c(n_outer_eta, n_inner_eta))
  theta_idx <- rep(1:(K+1), times = c(n_outer_theta, n_inner_theta))

  get_derivs <- function(eta, theta, deriv = 0) {
    # eta - list or vector of eta
    # tau - vector of extra parameters (on real line)
    # deriv - integer 0, 1, 2, what level of derivatives needed
    if (is.list(eta)) {
      eta <- do.call("cbind", eta)
    }
    store <- list()
    outer_eta <- eta[,eta_idx == 1,drop=FALSE]
    outer_theta <- theta[theta_idx == 1]
    outer_w <- outer_weight(outer_eta, outer_theta, deriv)

    inner_eta <- list()
    inner_theta <- list()
    inner_w <- list()

    f_eval <- list()


    for (i in 1:K) {
      inner_eta[[i]] <- eta[,eta_idx == (i+1),drop=FALSE]
      inner_theta[[i]] <- theta[theta_idx == (i+1)]
      inner_w[[i]] <- inner_weight[[i]](inner_eta[[i]], inner_theta[[i]], deriv = deriv)
      f_eval[[i]] <- outer_w$f_eval[,i] * inner_w[[i]]$f_eval
    }

    store$f_eval <- do.call("cbind", f_eval)

    if (deriv >= 1) {
      ow_idx <- rep(1:attr(outer_weight, "num_weights"), times = n_k)
      eta_ele_idx <- do.call("c", lapply(n_inner_eta, function(n) 1:n))
      theta_ele_idx <- do.call("c", lapply(n_inner_theta, function(n) 1:n))
      eta_set_idx <- rep(1:length(inner_weight), times = n_inner_eta)
      theta_set_idx <- rep(1:length(inner_weight), times = n_inner_theta)

      inner_eval <- do.call("cbind", lapply(inner_w, "[[", "f_eval"))
      outer_eval <- outer_w$f_eval
      store$f_eta_eval <- list()
      store$f_theta_eval <- list()

      if (have_outer_eta) {
        for (i in 1:n_outer_eta) {
          store$f_eta_eval[[i]] <- inner_eval * outer_w$f_eta_eval[[i]][,ow_idx]
        }
      }

      if (have_inner_eta) {
        for (i in (1:sum(n_inner_eta))) {
          eta_idx <- i + n_outer_eta
          store$f_eta_eval[[eta_idx]] <- matrix(0, nrow = nrow(store$f_eval), ncol = n_weights)
          store$f_eta_eval[[eta_idx]][, ow_idx == eta_set_idx[i]] <- outer_w$f_eval[,eta_set_idx[i]] *
            inner_w[[eta_set_idx[i]]]$f_eta_eval[[eta_ele_idx[i]]]
        }
      }

      if (have_outer_theta) {
        for (i in 1:n_outer_theta) {
          store$f_theta_eval[[i]] <- inner_eval * outer_w$f_theta_eval[[i]][,ow_idx]
        }
      }

      if (have_inner_theta) {
        for (i in 1:sum(n_inner_theta)) {
          theta_idx <- i + n_outer_theta
          store$f_theta_eval[[theta_idx]] <- matrix(0, nrow = nrow(store$f_eval), ncol = n_weights)
          store$f_theta_eval[[theta_idx]][, ow_idx == theta_set_idx[i]] <- outer_w$f_eval[,theta_set_idx[i]] *
            inner_w[[theta_set_idx[i]]]$f_theta_eval[[theta_ele_idx[i]]]
        }
      }
      if (deriv > 1) {
        store$f_eta2_eval <- list()
        store$f_eta_theta_eval <- list()
        store$f_theta2_eval <- list()
        if (have_outer_eta) {
          for (i in 1:n_outer_eta) {
            store$f_eta2_eval[[i]] <- list()
            store$f_eta_theta_eval[[i]] <- list()

            for (j in 1:n_outer_eta) {
              store$f_eta2_eval[[i]][[j]] <- inner_eval * outer_w$f_eta2_eval[[i]][[j]][,ow_idx]
            }

            if (have_inner_eta) {
              for (j in 1:sum(n_inner_eta)) {
                jx <- j + n_outer_eta
                store$f_eta2_eval[[i]][[jx]] <- matrix(0, nrow = nrow(store$f_eval), ncol = n_weights)
                store$f_eta2_eval[[i]][[jx]][,ow_idx == eta_set_idx[j]] <- outer_w$f_eta_eval[[i]][,eta_set_idx[j]] *
                  inner_w[[eta_set_idx[j]]]$f_eta_eval[[eta_ele_idx[j]]]
              }
            }

            if (have_outer_theta) {
              for (j in 1:n_outer_theta) {
                store$f_eta_theta_eval[[i]][[j]] <- inner_eval * outer_w$f_eta_theta_eval[[i]][[j]][,ow_idx]
              }
            }

            if (have_inner_theta) {
              for (j in 1:sum(n_inner_theta)) {
                jx <- j + n_outer_theta
                store$f_eta_theta_eval[[i]][[jx]] <- matrix(0, nrow = nrow(store$f_eval), ncol = n_weights)
                store$f_eta_theta_eval[[i]][[jx]][,ow_idx == theta_set_idx[j]] <- outer_w$f_eta_eval[[i]][,theta_set_idx[j]] *
                  inner_w[[theta_set_idx[j]]]$f_theta_eval[[theta_ele_idx[j]]]
              }
            }

          }
        }

        if (have_inner_eta) {
          for (i in 1:sum(n_inner_eta)) {
            ix <- i + n_outer_eta
            store$f_eta2_eval[[ix]] <- list()
            store$f_eta_theta_eval[[ix]] <- list()
            if (have_outer_eta) {
              for (j in 1:n_outer_eta) {
                store$f_eta2_eval[[ix]][[j]] <- store$f_eta2_eval[[j]][[ix]]
              }
            }

            for (j in 1:sum(n_inner_eta)) {
              jx <- j + n_outer_eta
              store$f_eta2_eval[[ix]][[jx]] <-  matrix(0, nrow = nrow(store$f_eval), ncol = n_weights)
              if (eta_set_idx[i] == eta_set_idx[j]) {
                store$f_eta2_eval[[ix]][[jx]][,ow_idx == eta_set_idx[i]] <- outer_w$f_eval[,eta_set_idx[i]] *
                  inner_w[[eta_set_idx[i]]]$f_eta2_eval[[eta_ele_idx[i]]][[eta_ele_idx[j]]]
              }
            }

            if (have_outer_theta) {
              for (j in 1:n_outer_theta) {
                store$f_eta_theta_eval[[ix]][[j]] <- matrix(0, nrow = nrow(store$f_eval), ncol = n_weights)
                store$f_eta_theta_eval[[ix]][[j]][,ow_idx==eta_set_idx[i]] <-
                  outer_w$f_theta_eval[[j]][,eta_set_idx[i]] * inner_w[[eta_set_idx[i]]]$f_eta_eval[[eta_ele_idx[i]]]
              }
            }

            if (have_inner_theta) {
              for (j in 1:sum(n_inner_theta)) {
                jx <- j + n_outer_theta
                store$f_eta_theta_eval[[ix]][[jx]] <- matrix(0, nrow = nrow(store$f_eval), ncol = n_weights)
                if (eta_set_idx[i] == theta_set_idx[j]) {
                  store$f_eta_theta_eval[[ix]][[jx]][,ow_idx==eta_set_idx[i]] <- outer_w$f_eval[,eta_set_idx[i]] *
                    inner_w[[eta_set_idx[i]]]$f_eta_theta_eval[[eta_ele_idx[i]]][[theta_ele_idx[j]]]
                }
              }
            }
          }
        }

        if (have_outer_theta) {
          for (i in 1:n_outer_theta) {
            store$f_theta2_eval[[i]] <- list()
            for (j in 1:n_outer_theta) {
              store$f_theta2_eval[[i]][[j]] <- inner_eval * outer_w$f_theta2_eval[[i]][[j]][,ow_idx]
            }
            if (have_inner_theta) {
              for (j in 1:sum(n_inner_theta)) {
                jx <- j + n_outer_theta
                store$f_theta2_eval[[i]][[jx]] <- matrix(0, nrow = nrow(store$f_eval), ncol = n_weights)
                store$f_theta2_eval[[i]][[jx]][,ow_idx == theta_set_idx[j]] <- outer_w$f_theta_eval[[i]][,theta_set_idx[j]] *
                  inner_w[[theta_set_idx[j]]]$f_theta_eval[[theta_ele_idx[j]]]
              }
            }
          }
        }

        if (have_inner_theta) {
          for (i in 1:sum(n_inner_theta)) {
            ix <- i + n_outer_theta
            store$f_theta2_eval[[ix]] <- list()
            if (have_outer_theta) {
              for (j in 1:n_outer_theta) {
                store$f_theta2_eval[[ix]][[j]] <- store$f_theta2_eval[[j]][[ix]]
              }
            }
            for (j in 1:sum(n_inner_theta)) {
              jx <- j + n_outer_theta
              store$f_theta2_eval[[ix]][[jx]] <- matrix(0, nrow = nrow(store$f_eval), ncol = n_weights)
              if (theta_set_idx[i] == theta_set_idx[j]) {
                store$f_theta2_eval[[ix]][[jx]][,ow_idx == theta_set_idx[i]] <- outer_w$f_eval[,theta_set_idx[i]] *
                  inner_w[[theta_set_idx[i]]]$f_theta2_eval[[theta_ele_idx[i]]][[theta_ele_idx[j]]]
              }
            }
          }
        }

      }
    }
    return(store)
  }

  pen <- function(tau, deriv = 0) {
    all_weights <- c(outer_weight, inner_weight)
    theta_pens <- lapply(1:(K+1), function(i) attr(all_weights[[i]], "theta_pen")(tau[i==theta_idx], deriv = deriv))
    pen <- NULL
    pen_grad <- NULL
    pen_hess <- NULL

    pen <- sum(unlist(sapply(theta_pens, "[[", "p")))

    if (deriv > 0) {
      pen_grad <-  unlist(lapply(theta_pens, "[[", "pt"))
      if (deriv > 0) {
        d_penlist <- lapply(theta_pens, "[[", "ptt")
        pen_hess <- as.matrix(Matrix::bdiag(d_penlist[!sapply(d_penlist, is.null)]))
      }
    }
    return(list(p = pen, pt = pen_grad, ptt = pen_hess))
  }

  init_func <- function(scores) {
    # Initialise with binary score indicating best expert
    k_score <- matrix(nrow = nrow(scores), ncol = K)
    theta_init <- list()
    mu_init <- list()
    score_idx <- rep(1:K, times = n_k)
    for (k in 1:K) {
      inner_score <- scores[,k==score_idx]
      inner_init <- attr(inner_weight[[k]], "init_func")(inner_score)
      theta_init[[k]] <- inner_init$init_theta
      mu_init[[k]] <- inner_init$init_mu
      init_weights <- inner_weight[[k]](mu_init[[k]], theta_init[[k]],deriv=0)$f_eval
      k_score[,k] <- rowSums(init_weights * inner_score)
    }

    outer_init <- attr(outer_weight, "init_func")(k_score)
    outer_mu <- outer_init$init_mu
    outer_theta <- outer_init$init_theta

    init_theta <- c(outer_theta, do.call("c", theta_init))
    init_mu <- cbind(outer_mu, do.call("cbind", mu_init))

    return(list(init_theta = init_theta, init_mu = init_mu))
  }

  return(structure(get_derivs,
                   neta = total_eta,
                   ntheta = total_theta,
                   theta_pen = pen,
                   init_func = init_func,
                   num_weights = n_weights, name = "nested"))
}


#' Latent weights aim to recreate multinomial with less linear predictors by
#' using 'latent' linear predictors
#'
#' @param K integer - Number of experts to weight
#' @param q integer - Number of latent factors to use
#'
#' @return
#' @export
#'
#' @examples
Latent <- function(K, q) {
  nv <- (K-1-q)*q

  neta <- q
  ntheta <- nv + K-1

  get_derivs <- function(eta, theta, deriv = 0) {
    # eta - list or vector of eta
    # tau - vector of extra parameters (on real line)
    # deriv - integer 0, 1, 2, what level of derivatives needed
    if (is.list(eta)) {
      eta <- do.call("cbind", eta)
    }
    beta0 <- theta[1:(K-1)]
    N <- nrow(eta)
    v_mat <- rbind(diag(q), matrix(theta[K:length(theta)], ncol = q))

    eta_aug <- (eta %*% t(v_mat)) + matrix(beta0, ncol = (K-1), nrow = N, byrow = TRUE)

    aug_w <- multinomial(K, normalised = TRUE)
    aug_deriv <- aug_w(eta_aug, 0, deriv = deriv)

    store <- list()
    store$f_eval <- aug_deriv$f_eval

    if (deriv >= 1) {
      store$f_eta_eval <- lapply(1:neta, function(i) matrix(nrow = N, ncol = K))
      for (k in 1:K) {
        M <- do.call("cbind", lapply(aug_deriv$f_eta_eval, function(i) i[,k])) %*%  v_mat
        for (j in 1:neta) {
          store$f_eta_eval[[j]][,k] <- M[,j]
        }
      }
      store$f_theta_eval <- lapply(1:ntheta, function(i) matrix(nrow = N, ncol = K))

      for (i in 1:(K-1)) {
        store$f_theta_eval[[i]] <- aug_deriv$f_eta_eval[[i]]
      }

      vj <- rep(1:q, times = (K-1-q))
      for (i in (K:ntheta)) {
        vx <- i - (K-1)
        ix <- floor(q+1+(i-K)/q)
        jx <- vj[vx]
        store$f_theta_eval[[i]] <- aug_deriv$f_eta_eval[[ix]]*eta[,jx]
      }

      if (deriv >= 2) {
        store$f_eta2_eval <- list()
        for (i in 1:neta) {
          store$f_eta2_eval[[i]] <- list()
          for (j in 1:neta) {
            out = 0
            for (h in 1:(K-1)) {
              for (l in 1:(K-1)) {
                out = out + aug_deriv$f_eta2_eval[[h]][[l]] * v_mat[h,i] * v_mat[l,j]
              }
            }
            store$f_eta2_eval[[i]][[j]] <- out
          }
        }
        store$f_eta_theta_eval <- list()
        for (i in 1:neta) {
          store$f_eta_theta_eval[[i]] <- list()
          for (j in 1:(K-1)) {#eta bj0 derivs
            out = 0
            for (h in 1:(K-1)) {
              out = out + aug_deriv$f_eta2_eval[[h]][[j]] * v_mat[h,i]
            }
            store$f_eta_theta_eval[[i]][[j]] <- out
          }
          for (j in K:ntheta) {
            vx <- j - (K-1) # which element of v
            ix <- floor(q+1+(j-K)/q) # which row of v
            jx <- vj[vx] # which col of v
            if (jx == i) {
              out = aug_deriv$f_eta_eval[[ix]]
            } else {
              out = 0
            }
            for (h in 1:(K-1)) {
              out = out + aug_deriv$f_eta2_eval[[h]][[ix]]*v_mat[h,i]*eta[,jx]
            }
            store$f_eta_theta_eval[[i]][[j]] <- out
          }
        }
        for (i in 1:(K-1)) {
          store$f_theta2_eval[[i]] <- list()
          for (j in 1:(K-1)) {
            store$f_theta2_eval[[i]][[j]] <- aug_deriv$f_eta2_eval[[i]][[j]]
          }
          for (j in K:ntheta) {
            vx <- j - (K-1) # which element of v
            c <- floor(q+1+(j-K)/q) # which row of v
            d <- vj[vx] # which col of v

            store$f_theta2_eval[[i]][[j]] <- aug_deriv$f_eta2_eval[[i]][[c]]*eta[,d]
          }
        }
        for (i in K:ntheta) {
          store$f_theta2_eval[[i]] <- list()
          for (j in 1:(K-1)) {
            store$f_theta2_eval[[i]][[j]] <- store$f_theta2_eval[[j]][[i]]
          }
          for (j in K:ntheta) {
            ix <- i - (K-1)
            jx <- j - (K-1)
            a <- floor(q+1+(i-K)/q)
            b <- vj[ix]
            c <- floor(q+1+(j-K)/q)
            d <- vj[jx]

            store$f_theta2_eval[[i]][[j]] <- aug_deriv$f_eta2_eval[[a]][[c]]*eta[,b]*eta[,d]
          }
        }
      }
    }
    return(store)
  }

  init_func <- function(scores) {
    aug_weights <- multinomial(K)
    init_aug_eta <- attr(aug_weights, "init_func")(scores)
    init_eta <- init_aug_eta[,1:neta]
    b0s <- rep(0, K-1)
    vs <- list()
    for (j in (neta+1):(K-1)) {
      y1 <- init_aug_eta[,j]
      X1 <- init_eta
      coefs <- coef(lm(y1~X1))
      b0s[j] <- coefs[1]
      vs[[j-neta]] <- coefs[-1]
    }
    vs <- do.call("c", vs)
    init_theta <- c(b0s, vs)
    return(list(init_theta=init_theta,init_mu=init_eta))
  }

  pen <- function(tau, deriv = 0) {
    p <- 0
    pt <- rep(0, ntheta)
    ptt <- matrix(0, nrow = ntheta, ncol = ntheta)
    return(list(p = p, pt = pt, ptt = ptt))
  }

  return(structure(get_derivs,
                   neta=neta,
                   ntheta=ntheta,
                   theta_pen=pen,
                   init_func=init_func,
                   num_weights=K, name = "Latent"))
}

# Only used for theoretical results
ordinal_alt <- function(K) {
  # If normalised we don't need K eta
  neta <- 1

  get_derivs <- function(eta, theta, deriv) {
    store <- list()
    if (is.list(eta)) {
      eta <- do.call("cbind", eta)
    }
    theta <- c(0, theta, 1)

    N <- nrow(eta)
    K <- length(theta)
    phi_mat <- matrix(theta, nrow = N, ncol = K, byrow = TRUE)
    etam <- eta%*%theta
    W <- exp(etam)/rowSums(exp(etam))
    store$f_eval <- W
    if (deriv >= 1) {
      store$f_eta_eval <- list()
      store$f_eta_eval[[1]] <- W * (phi_mat - rowSums(W*phi_mat))
      store$f_theta_eval <- list()
      for (i in 1:(K-2)) {
        tmp_mat <- matrix(0, nrow = N, ncol = K)
        tmp_mat[,i+1] <- 1
        store$f_theta_eval[[i]] <- (as.vector(eta) * W) * (tmp_mat - W[,i+1])
      }
      if (deriv >= 2) {
        store$f_eta2_eval <- list()
        store$f_eta2_eval[[1]] <- list()
        store$f_eta2_eval[[1]][[1]] <- W * (phi_mat - rowSums(W*phi_mat))^2 - W*rowSums(W*phi_mat^2) + W*(rowSums(phi_mat*W)^2)
          #W * (phi_mat - rowSums(W*phi_mat))^2 - W*sapply(1:K, FUN = function(k) {rowSums(W*(phi_mat-phi_mat[,k])^2)})
        store$f_theta2_eval <- list()
        for (i in 1:(K-2)) {
          store$f_theta2_eval[[i]] <- list()
          for (j in 1:(K-2)) {
            tmp_mat1 <- matrix(0, nrow = N, ncol = K)
            tmp_mat1[,i+1] <- 1
            tmp_mat2 <- matrix(0, nrow = N, ncol = K)
            tmp_mat2[,j+1] <- 1
            store$f_theta2_eval[[i]][[j]] <- (as.vector(eta^2) * W)* ((tmp_mat1 - W[,i+1]) * (tmp_mat2 - W[,j+1]) - W[,i+1]*((i==j)-W[,j+1]))
          }
        }
        store$f_eta_theta_eval <- list()
        store$f_eta_theta_eval[[1]] <- list()
        for (j in 1:(K-2)) {
          tmp_mat <- matrix(0, nrow = N, ncol = K)
          tmp_mat[,j+1] <- 1
          store$f_eta_theta_eval[[1]][[j]] <- as.vector(eta)*W*(tmp_mat - W[,j+1])*(phi_mat-rowSums(phi_mat*W)) + W*(tmp_mat + as.vector(rowSums(W*phi_mat)*eta)*W[,j+1] - W[,j+1] - as.vector(eta*phi_mat[,j+1]*W[,j+1]))
            #W * (as.vector(eta)*(tmp_mat-W[,j+1])*(phi_mat - rowSums(W*phi_mat)) + (tmp_mat-W[,j+1]))
        }
      }
    }
    return(store)
  }

  name <- "ordinal_alt"

  return(structure(get_derivs,
                   neta = neta,
                   ntheta = K-2,
                   #theta_pen = pen,
                   #init_func = init_func,
                   num_weights = K, name = name))
}








