#' Fast Multivariate Normal weights (vectorized)
#'
#' Produces an object containing derivatives for multivariate normal stacking (vectorized for speed)
#'
#' @param ex_coords Array of (dim_num, K) specifying locations for each of the experts
#' @return Object containing derivatives for multivariate normal stacking
#' @export
#'
#' @examples
#' x <- matrix(c(1:6), nrow = 3)
#' MVN_weights_fast(x)

MVN_weights_fast <- function(ex_coords) {
  x <- ex_coords
  n_k <- ncol(x)
  dim_num <- nrow(x)
  x <- t(x)

  get_derivs <- function(eta, theta, deriv = 1) {
    if (is.list(eta)) eta <- do.call("cbind", eta)
    tau <- exp(theta)
    N <- nrow(eta)
    # Use llk_gaussian2 for all density and derivative calculations, as in the original
    dens_list <- vector("list", dim_num)
    if (deriv == 0) {
      for (n in 1:dim_num) {
        dens_list[[n]] <- llk_gaussian2(matrix(x[,n], ncol = n_k, nrow = N),
                                        param = list(eta[,n], 1/sqrt(tau[n])),
                                        deriv = 0)$d0
      }
      log_dens <- Reduce("+", dens_list)
      shift_mat <- log_dens - matrixStats::rowMaxs(log_dens)
      expshift <- exp(shift_mat)
      rs_expshift <- rowsums(expshift)
      store <- list()
      store$f_eval <- expshift / rs_expshift
      return(store)
    }
    # For deriv >= 1, compute all at once to avoid redundant calls
    dens_list <- vector("list", dim_num)
    for (n in 1:dim_num) {
      dens_list[[n]] <- llk_gaussian2(matrix(x[,n], ncol = n_k, nrow = N),
                                      param = list(eta[,n], 1/sqrt(tau[n])),
                                      deriv = deriv)
    }
    # d0: log densities, d1: list(eta, tau), d2: list(eta, tau, eta2, ...)
    log_dens <- Reduce("+", lapply(dens_list, "[[", "d0"))
    shift_mat <- log_dens - matrixStats::rowMaxs(log_dens)
    expshift <- exp(shift_mat)
    rs_expshift <- rowsums(expshift)
    store <- list()
    store$f_eval <- expshift / rs_expshift

    # Derivatives
    detas <- sapply(dens_list, "[[", "d1")[1, ]
    dtaus <- sapply(dens_list, "[[", "d1")[2, ]
    f_eta_out <- vector("list", dim_num)
    f_tau_out <- vector("list", dim_num)
    for (i in 1:dim_num) {
      l1 <- detas[[i]]
      l2 <- dtaus[[i]]
      f_eta_out[[i]] <- store$f_eval * (l1 - (rowsums(l1 * expshift)) / rs_expshift)
      f_tau_out[[i]] <- store$f_eval * (l2 - (rowsums(l2 * expshift)) / rs_expshift)
    }
    store$f_eta_eval <- f_eta_out
    store$f_tau_eval <- f_tau_out
    store$f_theta_eval <- list_by_vector(f_tau_out, tau)

    # Second derivatives (optional, as in original)
    if (deriv >= 2) {
      # Copy from original if needed
      store$f_eta2_eval <- NULL
      store$f_eta_theta_eval <- NULL
      store$f_theta2_eval <- NULL
    }
    return(store)
  }

  init_func <- function(densities) {
    theta <- apply(x, 2, stats::var)
    y1 <- max.col(densities)
    mustart <- x[y1, , drop = FALSE]
    return(list(init_theta = log(theta), init_mu = mustart))
  }

  pen <- function(tau, deriv = 0) {
    theta <- exp(tau)
    v <- apply(x, 2, stats::var)
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

  name = "MVN_fast"
  structure(
    get_derivs,
    init_func = init_func,
    arg_list = list("x" = x, dim_num = dim_num),
    neta = dim_num,
    ntheta = dim_num,
    theta_pen = pen,
    num_weights = n_k,
    name = name
  )
}
