# Functions that can be used for internal model weights and derivatives
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

# Ordinal inner functions ------------------------------------------------------
#' Produces an object containing derivatives for ordered stacking
#'
#' @param num_weights Integer: How many weights in the ordering
#'
#' @return List of functions for finding ordinal weights
#' @export
#'
#' @examples
#' ordinal(3)
ordinal <- function(num_weights) { # nolint
  F <- function(x) {
    1 / (1 + exp(- (x)))
  }

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
    return(sapply(1:(length(theta) + 1), function(x) {f_j(eta, theta, x)}))
  }

  F_diff <- function(x) {
    abs_x <- abs(x)
    exp_x <- exp(-abs_x)
    exp_2x <- exp_x^2
    exp_x_p1_2 <- (1 + exp_x)^2
    out <- exp_x / exp_x_p1_2
    return(out)
  }

  f_eta_j <- function(eta, theta, j) {
    if (j == 1) {
      return(-F_diff(theta[1] - eta))
    } else if (j == length(theta) + 1) {
      return(F_diff(theta[j - 1] - eta))
    } else {
      return(F_diff(theta[j - 1] - eta) - F_diff(theta[j] - eta))
    }
  }

  f_eta <- function(eta, theta) {
    return(sapply(1:(length(theta) + 1), function(x) {f_eta_j(eta, theta, x)}))
  }

  f_theta_n <- function(eta, theta, n) {
    store <- matrix(0, length(eta), length(theta) + 1)
    store[, n] <- F_diff(theta[n] - eta)
    store[, n + 1] <- -F_diff(theta[n] - eta)
    return(store)
  }

  f_theta <- function(eta, theta) {
    out <- list()

    for (i in seq_len(length(theta))) {
      out[[i]] <- f_theta_n(eta, theta, i)
    }
    return(out)
  }

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
    return(sapply(1:(length(theta) + 1), function(x) {f_eta2_j(eta, theta, x)}))
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
    store <- list()
    tau <- c(stats::qlogis(1/num_weights), tau)
    theta <- c(tau[1], tau[1] + cumsum(exp(tau[-1])))

    store$f_eval <- f(eta, theta)

    if (deriv >= 1) {
      store$f_eta_eval <- f_eta(eta, theta)

      # Take cumulative sum of reverse order, then reverse so sum at front and just last term at end
      store$f_theta_eval <- list_by_vector(rev(mat_cumsum(rev(f_theta(eta, theta)))), c(tau[1], exp(tau[-1])))[-1]

      if (deriv >= 2) {
        store$f_eta2_eval <- f_eta2(eta, theta)

        store$f_eta_theta_eval <- list_by_vector(rev(mat_cumsum(rev(f_eta_theta(eta, theta)))),
                                                 c(tau[1], exp(tau[-1])))[-1]

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
        store$f_theta2_eval <- f_theta2_eval
      }
    }

    return(store)
  }

  init_func <- function(densities) {
    y1 <- max.col(densities)
    init_tau <- log(diff(stats::qlogis(seq(1, num_weights - 1) / num_weights)))
    fixed_tau <- stats::qlogis(1 / num_weights)
    tau <- c(fixed_tau, init_tau, 1)
    theta <- c(fixed_tau - 1, tau[1], tau[1] + cumsum(exp(tau[-1])))
    mustart <- (theta[y1 + 1] + theta[y1]) / 2
    return(list(init_mu = mustart, init_theta = init_tau))
  }

  pen <- function(tau, deriv = 0) {
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

  return(structure(get_derivs, init_func = init_func, theta_pen = pen, neta = 1, ntheta = num_weights - 2))
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
#' MVN_weights(x, 3)
MVN_weights <- function(x, dim_num) {
  # Input
  # -----
  # x - Array of (dim_num, n_k) specifying locations for each of the models
  # Suppose we have N data points
  # n_k inner functions
  # inner functions f_j depend on:
  #   etaT_1,..., etaT_p
  #   etaT denotes tilde(eta)
  #   theta_1,..., theta_q
  force(x)
  n_k <- ncol(x)
  dim_num <- nrow(x)

  x <- t(x)

  get_derivs <- function(eta, theta, deriv = 1) {
    eta <- do.call("cbind", eta)
    tau <- exp(theta)
    dens_matrix <- matrix(nrow = nrow(eta) * n_k, ncol = dim_num)

    for (n in 1:dim_num) {
      dens_matrix[, n] <-  llk_gaussian(rep(x[, n], each = nrow(eta)),
                                        param = list(rep(eta[, n], n_k), 1 / sqrt(tau[n])),
                                        deriv = 0)$d0
    }

    log_dens <- matrix(rowSums(dens_matrix), ncol = n_k)
    shift_mat <- log_dens - matrixStats::rowMaxs(log_dens)

    store <- list()
    store$f_eval <- exp(shift_mat) / rowSums(exp(shift_mat))

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

        f_eta_out[[i]] <- store$f_eval * (l1 - (rowSums(l1 * exp(shift_mat))) / rowSums(exp(shift_mat)))
        f_tau_out[[i]] <- store$f_eval * (l2 - (rowSums(l2 * exp(shift_mat))) / rowSums(exp(shift_mat)))
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
            p2 <- - rowSums(dbeta1 * exp(shift_mat)) * store$f_eval * dalpha1 / rowSums(exp(shift_mat))
            p3 <- - (rowSums(d2 * exp(shift_mat)) * store$f_eval) / rowSums(exp(shift_mat))
            p35 <- - rowSums(dalpha1 * exp(shift_mat)) * store$f_eval * dbeta1 / rowSums(exp(shift_mat))
            p4 <- rowSums(dbeta1 * exp(shift_mat)) * rowSums(dalpha1 * exp(shift_mat)) * 2 * store$f_eval / rowSums(exp(shift_mat))^2
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
            p2 <- - rowSums(dbeta1 * exp(shift_mat)) * store$f_eval * dalpha1 / rowSums(exp(shift_mat))
            p3 <- - (rowSums(d2 * exp(shift_mat)) * store$f_eval) / rowSums(exp(shift_mat))
            p35 <- - rowSums(dalpha1 * exp(shift_mat)) * store$f_eval * dbeta1 / rowSums(exp(shift_mat))
            p4 <- rowSums(dbeta1 * exp(shift_mat)) * rowSums(dalpha1 * exp(shift_mat)) * 2 * store$f_eval / rowSums(exp(shift_mat))^2
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
            p2 <- - rowSums(dbeta1 * exp(shift_mat)) * store$f_eval * dalpha1 / rowSums(exp(shift_mat))
            p3 <- - (rowSums(d2 * exp(shift_mat)) * store$f_eval) / rowSums(exp(shift_mat))
            p35 <- - rowSums(dalpha1 * exp(shift_mat)) * store$f_eval * dbeta1 / rowSums(exp(shift_mat))
            p4 <- rowSums(dbeta1 * exp(shift_mat)) * rowSums(dalpha1 * exp(shift_mat)) * 2 * store$f_eval / rowSums(exp(shift_mat))^2

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
    theta <- (apply(t(x), 1, (stats::var)))
    y1 <- max.col(densities)
    mustart <- (x[y1, ])
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

      pen_hess <- diag((4 * theta - 2 * v) * theta)
    }

    return(list(p = pen, pt = pen_grad, ptt = pen_hess))
  }


  return(structure(
    get_derivs,
    init_func = init_func,
    arg_list = list("x" = x, dim_num = dim_num),
    neta = dim_num,
    ntheta = dim_num,
    theta_pen = pen
  ))
}

# ========================================
#' Produces equal weights without any dependence on parameters
#'
#' @param num_weights Integer: How many equally weighted parts
#' @param N  Integer: Number of data points
#'
#' @return Object of derivatives for equally weighted models
#' @export
#'
#' @examples
#' id(3, 1000)
id <- function(num_weights, N) {
  # Suppose we have N data points
  # n_k inner functions
  # inner functions f_j depend on:
  #   etaT_1,..., etaT_p
  #   etaT denotes tilde(eta)
  #   theta_1,..., theta_q
  get_derivs <- function(eta, theta, deriv) {
    store <- list()
    if (deriv >= 0) {
      store$f_eval <- matrix(1/num_weights, nrow = N, ncol = num_weights)
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

  return(structure(get_derivs, neta = 0, ntheta = 0, theta_pen = pen, init_func = init_func))
}
