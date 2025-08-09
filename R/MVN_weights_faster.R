#' Faster Multivariate Normal weights (fully vectorized)
#'
#' Produces an object containing derivatives for multivariate normal stacking
#' with additional optimizations and full vectorization
#'
#' @param ex_coords Array of (dim_num, K) specifying locations for each of the experts
#' @return Object containing derivatives for multivariate normal stacking
#' @export
#'
#' @examples
#' x <- matrix(c(1:6), nrow = 3)
#' MVN_weights_faster(x)
MVN_weights_faster <- function(ex_coords) {
  x <- ex_coords
  n_k <- ncol(x)
  dim_num <- nrow(x)
  x <- t(x)
  
  # Pre-compute variance for initialization and penalty
  x_var <- apply(x, 2, stats::var)
  
  get_derivs <- function(eta, theta, deriv = 1) {
    if (is.list(eta)) eta <- do.call("cbind", eta)
    tau <- exp(theta)
    N <- nrow(eta)
    
    # Pre-allocate storage for efficiency
    store <- list()
    dens_array <- array(0, dim = c(N, n_k, dim_num))
    deta_array <- if(deriv >= 1) array(0, dim = c(N, n_k, dim_num)) else NULL
    dtau_array <- if(deriv >= 1) array(0, dim = c(N, n_k, dim_num)) else NULL
    
    # Vectorized computation for all dimensions
    for (n in 1:dim_num) {
      # Reuse matrix for efficiency
      x_mat <- matrix(x[,n], ncol = n_k, nrow = N)
      gauss_result <- llk_gaussian2(x_mat,
                                  param = list(eta[,n], 1/sqrt(tau[n])),
                                  deriv = deriv)
      dens_array[,,n] <- gauss_result$d0
      if (deriv >= 1) {
        deta_array[,,n] <- gauss_result$d1[[1]]
        dtau_array[,,n] <- gauss_result$d1[[2]]
      }
    }
    
    # Sum log densities across dimensions
    log_dens <- rowSums(dens_array, dims = 2)
    shift_mat <- log_dens - matrixStats::rowMaxs(log_dens)
    expshift <- exp(shift_mat)
    rs_expshift <- rowSums(expshift)
    
    # Compute f_eval with matrix operations
    store$f_eval <- expshift / rs_expshift
    
    if (deriv >= 1) {
      # Pre-compute common terms for derivatives
      expshift_rs <- expshift / rs_expshift  # normalized weights
      weighted_sum_cache <- list()
      
      # Vectorized first derivatives
      store$f_eta_eval <- vector("list", dim_num)
      store$f_tau_eval <- vector("list", dim_num)
      
      for (i in 1:dim_num) {
        # Cache weighted sums for reuse
        weighted_sum_deta <- rowSums(deta_array[,,i] * expshift)
        weighted_sum_dtau <- rowSums(dtau_array[,,i] * expshift)
        weighted_sum_cache[[i]] <- list(deta = weighted_sum_deta, dtau = weighted_sum_dtau)
        
        # Compute derivatives using cached values
        store$f_eta_eval[[i]] <- store$f_eval * (deta_array[,,i] - weighted_sum_deta / rs_expshift)
        store$f_tau_eval[[i]] <- store$f_eval * (dtau_array[,,i] - weighted_sum_dtau / rs_expshift)
      }
      
      # Compute theta derivatives
      store$f_theta_eval <- list_by_vector(store$f_tau_eval, tau)
      
      if (deriv >= 2) {
        # Initialize second derivative storage
        store$f_eta2_eval <- vector("list", dim_num)
        store$f_eta_theta_eval <- vector("list", dim_num)
        store$f_theta2_eval <- vector("list", dim_num)
        
        for (alpha in 1:dim_num) {
          store$f_eta2_eval[[alpha]] <- vector("list", dim_num)
          store$f_eta_theta_eval[[alpha]] <- vector("list", dim_num)
          store$f_theta2_eval[[alpha]] <- vector("list", dim_num)
          
          dalpha1 <- deta_array[,,alpha]
          dalpha2 <- dtau_array[,,alpha]
          ws_alpha_deta <- weighted_sum_cache[[alpha]]$deta
          ws_alpha_dtau <- weighted_sum_cache[[alpha]]$dtau
          
          for (beta in 1:dim_num) {
            dbeta1 <- deta_array[,,beta]
            dbeta2 <- dtau_array[,,beta]
            ws_beta_deta <- weighted_sum_cache[[beta]]$deta
            ws_beta_dtau <- weighted_sum_cache[[beta]]$dtau
            
            if (alpha == beta) {
              d2_eta <- (dalpha1^2)
              d2_tau <- (dalpha2^2)
            } else {
              d2_eta <- dalpha1 * dbeta1
              d2_tau <- dalpha2 * dbeta2
            }
            
            # Compute second derivatives using cached values
            store$f_eta2_eval[[alpha]][[beta]] <- store$f_eval * (
              d2_eta - 
              ws_beta_deta * dalpha1 / rs_expshift -
              ws_alpha_deta * dbeta1 / rs_expshift +
              2 * ws_alpha_deta * ws_beta_deta * store$f_eval / rs_expshift
            )
            
            store$f_eta_theta_eval[[alpha]][[beta]] <- store$f_eval * (
              d2_tau - 
              ws_beta_dtau * dalpha1 / rs_expshift -
              ws_alpha_deta * dbeta2 / rs_expshift +
              2 * ws_alpha_deta * ws_beta_dtau * store$f_eval / rs_expshift
            ) * tau[beta]
            
            if (alpha == beta) {
              store$f_theta2_eval[[alpha]][[beta]] <- (
                store$f_eta_theta_eval[[alpha]][[beta]] * tau[alpha] +
                tau[alpha] * store$f_tau_eval[[alpha]]
              )
            } else {
              store$f_theta2_eval[[alpha]][[beta]] <- 
                store$f_eta_theta_eval[[alpha]][[beta]] * tau[alpha]
            }
          }
        }
      }
    }
    return(store)
  }
  
  init_func <- function(densities) {
    y1 <- max.col(densities)
    mustart <- x[y1, , drop = FALSE]
    return(list(init_theta = log(x_var), init_mu = mustart))
  }
  
  pen <- function(tau, deriv = 0) {
    theta <- exp(tau)
    pen <- sum((theta - x_var)^2)
    pen_grad <- NULL
    pen_hess <- NULL
    if (deriv > 0) {
      pen_grad <- theta * 2 * (theta - x_var)
      if (length(pen_grad) == 1) {
        pen_hess <- matrix((4 * theta - 2 * x_var) * theta)
      } else {
        pen_hess <- diag((4 * theta - 2 * x_var) * theta)
      }
    }
    return(list(p = pen, pt = pen_grad, ptt = pen_hess))
  }
  
  name = "MVN_faster"
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
