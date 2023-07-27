# Implements derivatives given the inner function derivatives.
#rowSums <- matrixStats::rowSums2
rowsums <- Cpp_rowSums
# h is a commonly seen function in the derivatives
h <- function(dens, f_deriv, ll_alpha) {
  if (is.list(f_deriv)) {
    out <- list()
    for (i in 1:length(f_deriv)) {
      #out[[i]] <- exp(log_rowSums_a_times_b(dens * log(d_deriv[[i]]))) * ll_alpha
      out[[i]] <- rowsums(dens * f_deriv[[i]]) * ll_alpha
    }
  } else {
    #out <- exp(log_rowSums_a_times_b(dens, log(f_deriv))) * ll_alpha
    out <-  rowsums(dens * f_deriv) * ll_alpha
  }
  return(out)
}

# takes cross product around diagonal matrix D (but D inputted as vector)
XtDX <- function(X1, D, X2) {
  as.matrix(crossprod(X1, D*X2))
}


log_rowSums_a_times_b <- function(log_a, log_b) {
  log_c <- log_a + log_b
  log_c_maxs <- matrixStats::rowMaxs(log_c)
  log_c_maxs + log(rowSums(exp(log_c - log_c_maxs)))
}


# Likelihood -------------------------------------------------------------------
ll <- function(list_of_f_eval, alpha_matrix, list_of_densities) {
  N <- nrow(list_of_densities[[1]])
  K <- length(list_of_f_eval)
  out <- matrix(nrow = N, ncol = K)
  log_f_d <- sapply(1:K, function(i) {
    log_rowSums_a_times_b(log(list_of_f_eval[[i]]), list_of_densities[[i]])})
  out <- log_rowSums_a_times_b(log(alpha_matrix), log_f_d)
  return(out)
}

# First derivatives
# First derivatives ------------------------------------------------------------
gd <- function(indexes, list_of_eta, list_of_densities, list_of_f_eval) {
  # SHOULD RENAME d, VERY CONFUSING GIVEN MY NOTATION IN WRITTEN DERIVATIVES
  N <- nrow(list_of_densities[[1]])
  K <- length(indexes)
  out <- matrix(nrow = N, ncol = K)
  list_of_eta2 <- c(list(matrix(0, nrow = N)), list_of_eta)
  for (i in indexes) {
    #out[,i] = exp(list_of_eta2[[i]] + log_rowSums_a_times_b(list_of_densities[[i]] * log(list_of_f_eval[[i]])))
    out[,i] = exp(list_of_eta2[[i]] + log(rowsums(list_of_densities[[i]] * list_of_f_eval[[i]])))
  }
  return(out)
}

ll_eta <- function(list_of_f_eval, list_of_eta, list_of_densities, alpha_mat) {
  K = length(list_of_densities)
  ds <- gd(1:(K), list_of_eta, list_of_densities, list_of_f_eval)
  out <- (ds/rowsums(ds) - alpha_mat/rowsums(alpha_mat))[,-1, drop = FALSE]
  return(lapply(seq_len(ncol(out)), function(i) out[,i]))
}


ll_etaTK <- function(ll_alpha, f_eta_eval, p) {
  # likelihood: evaluted likelihood
  # f_eta: list of evaluated f_etas
  # alpha: Evaluated vector of alphas
  # eta: current eta parameters
  # theta: current theta parameters
  # p: densities
  if (is.list(f_eta_eval)) {
    out = list()
    for (i in 1:length(f_eta_eval)) {
      if (is.null(f_eta_eval)) {
        out[i] <- list(NULL)
      } else {
        out[[i]] = h(p, f_eta_eval[[i]], ll_alpha)
      }
    }
    #N = length(f_eta_eval)
    #f_eta_eval = do.call(rbind, f_eta_eval)
  } else {
    out = h(p, f_eta_eval, ll_alpha)
  }

  #out = h(repmat(p, N, 1), f_eta_eval, rep(ll_alpha, N))
  return(out)
}

ll_etaT <- function(ll_alpha, list_of_f_eta_eval, list_of_densities) {
  out <- list()
  for (i in 1:length(list_of_f_eta_eval)) {
    if (is.null(list_of_f_eta_eval[[i]])) {
      out[i] <- list(NULL)
    } else {
      out[[i]] <- ll_etaTK(ll_alpha[,i], list_of_f_eta_eval[[i]], list_of_densities[[i]])
    }
  }
  return(out)
}



ll_thetaK <- function(ll_alpha, f_theta_eval, p) {
  # likelihood: evaluted likelihood
  # f_theta: list of evaluated f_thetas
  # alpha: Evaluated vector of alphas
  # eta: current eta parameters
  # theta: current theta parameters
  # p: densities
  f_theta = f_theta_eval
  out = vector(length = length(f_theta))
  for (n in 1:length(f_theta)) {
    out[n] = sum(h(p, f_theta[[n]], ll_alpha))
  }
  return(out)
}

ll_theta <- function(ll_alpha, list_of_f_theta_eval, list_of_densities) {
  out <- list()
  for (i in 1:length(list_of_f_theta_eval)) {
    if (is.null(list_of_f_theta_eval[[i]])) {
      out[i] <- list(NULL)
    } else {
      out[[i]] <- ll_thetaK(ll_alpha[,i], list_of_f_theta_eval[[i]], list_of_densities[[i]])
    }
  }
  return(out)
}

# Second derivatives -----------------------------------------------------------
## Second derivatives for eta
ll_eta2 <- function(alpha, list_of_eta, list_of_densities, list_of_f_eval, list_of_l_eta, list_of_X_eta) {
  K = ncol(alpha)
  N = nrow(list_of_densities[[1]])
  # SHOULD RENAME d, VERY CONFUSING GIVEN MY NOTATION IN WRITTEN DERIVATIVES
  ds <- gd(1:K, list_of_eta, list_of_densities, list_of_f_eval)
  ds <- ds/rowsums(ds)
  out <- list()
  for (i in 2:K) {
    out[[i]] <- list()
    for (j in 2:K) {
      if (i == j) {
        eta2_deriv <- list_of_l_eta[[i-1]] * (1 - ds[,i]) - (alpha[,i] * list_of_l_eta[[i-1]])
        out[[i-1]][[j-1]] <- XtDX(list_of_X_eta[[i-1]], eta2_deriv, list_of_X_eta[[j-1]])
      } else {
        eta2_deriv <- (list_of_l_eta[[i-1]] * (-ds[,j])) - (alpha[,i] * list_of_l_eta[[j-1]])
        out[[i-1]][[j-1]] <- XtDX(list_of_X_eta[[i-1]], eta2_deriv, list_of_X_eta[[j-1]])
      }
    }
    out[[i-1]] <- do.call("cbind", out[[i-1]])
  }

  out <- do.call("rbind", out)
  return(out)
}

## Second derivatives for eta
ll_betaT2 <- function(ll_alpha, f_eta2T_eval, list_of_ll_etaT, k1, k2, list_of_densities, list_of_X_etaT) {
  no_alpha <- if (is.list(f_eta2T_eval[[k1]])) {
    length(f_eta2T_eval[[k1]])
  } else if (is.null(f_eta2T_eval[[k1]])) {
    return(NULL)
  } else {
    f_eta2T_eval[[k1]] <- list(list(f_eta2T_eval[[k1]]))
    list_of_X_etaT[[k1]] <- list(list_of_X_etaT[[k1]])
    list_of_ll_etaT[[k1]] <- list(list_of_ll_etaT[[k1]])
    1
  }
  no_beta <- if (is.list(f_eta2T_eval[[k2]])) {
    length(f_eta2T_eval[[k2]])
  } else if (is.null(f_eta2T_eval[[k2]])) {
    return(NULL)
  } else {
    f_eta2T_eval[[k1]] <- list(list(f_eta2T_eval[[k1]]))# added list(list())
    list_of_X_etaT[[k2]] <- list(list_of_X_etaT[[k2]])
    list_of_ll_etaT[[k2]] <- list(list_of_ll_etaT[[k2]])
    1
  }

  out <- list()
  for (alpha in 1:no_alpha) {
    out[[alpha]] <- list()
    for (beta in 1:no_beta) {
      if (k1 == k2) {
        part1 = h(list_of_densities[[k1]], f_eta2T_eval[[k1]][[alpha]][[beta]], ll_alpha[,k1]) # don't want beta argument when only single list
      } else {
        part1 = 0
      }
      part2 = list_of_ll_etaT[[k1]][[alpha]] * list_of_ll_etaT[[k2]][[beta]]

      output <- XtDX(list_of_X_etaT[[k1]][[alpha]], part1 - part2, list_of_X_etaT[[k2]][[beta]])

      out[[alpha]][[beta]] <- output
    }
    out[[alpha]] <- do.call("cbind", out[[alpha]])
  }
  out <- do.call("rbind", out)
  out
}


## Second derivatives for theta
ll_theta2_nq <- function(ll_alpha, f_theta2_eval, list_of_densities, f_theta_eval1, f_theta_eval2, k1, k2, n, q) {
  p1 = list_of_densities[[k1]]
  p2 = list_of_densities[[k2]]

  s1 = h(p1, f_theta_eval1[[n]], ll_alpha[,k1])
  s2 = h(p2, f_theta_eval2[[q]], ll_alpha[,k2])

  if (is.null(f_theta_eval1) | is.null(f_theta_eval2)) {
    return(NULL)
  }

  if (k1 == k2) {
    part1 = h(p1, f_theta2_eval[[k1]][[n]][[q]], ll_alpha[,k1])
  } else {
    part1 = 0
  }
  return(sum(part1 - s1*s2))
}

ll_theta2 <- function(ll_alpha, f_theta2_eval, list_of_densities, f_theta_eval, k1, k2) {
  f_theta_eval1 <- f_theta_eval[[k1]]
  f_theta_eval2 <- f_theta_eval[[k2]]

  no_n = length(f_theta_eval1)
  no_q = length(f_theta_eval2)

  out = matrix(ncol = no_q, nrow = no_n)

  for (i in 1:no_n) {
    for (j in 1:no_q) {
      out[i, j] = ll_theta2_nq(ll_alpha, f_theta2_eval, list_of_densities, f_theta_eval1, f_theta_eval2, k1, k2, i, j)
    }
  }
  return(out)
}

## Mixed derivatives
### eta tilde and theta
ll_etaT_theta <- function(ll_alpha, f_eta_theta_eval, f_eta_eval, f_theta_eval, list_of_densities, k1, k2, list_of_X_etaT) {
  N <- nrow(ll_alpha)
  if (is.list(f_eta_eval[[k1]])) {
    no_etas = length(f_eta_eval[[k1]])
  } else {
    no_etas = 1
  }
  p1 <- (list_of_densities[[k1]])
  p2 <- (list_of_densities[[k2]])

  part1 <- list()
  if (k1 == k2) {
    for (q in 1:length(f_eta_theta_eval[[k1]])) {
      if (no_etas != 1) {
        part1[[q]] <- h(p1, lapply(1:no_etas, function(x) {f_eta_theta_eval[[k1]][[x]][[q]]}), ll_alpha[,k1])
      } else {
        part1[[q]] <- h(p1, f_eta_theta_eval[[k1]][[q]], ll_alpha[,k1])
      }

    }
  } else {
    for (q in 1:length(f_theta_eval[[k2]])) {
      if (is.list(f_eta_eval[[k1]])) {
        part1[[q]] <- as.list(rep(0, length(f_eta_eval[[k1]])))
      } else {
        part1[[q]] <- 0
      }
    }
  }
  #ll_eta_k1 = h(repmat(p1, no_etas, 1), unlist(f_eta_eval[[k1]]), ll_alpha[,k1])
  ll_eta_k1 = h(p1, f_eta_eval[[k1]], ll_alpha[,k1])
  out = list()
  for (q in 1:length(f_theta_eval[[k2]])) {
    if (!is.list(ll_eta_k1)) {
      out[[q]] <-  t(list_of_X_etaT[[k1]]) %*% (part1[[q]] - (ll_eta_k1 * h(p2, f_theta_eval[[k2]][[q]], ll_alpha[,k2])))

    } else {
      #out[[q]] <-  do.call("rbind", list_take_list(part1[[q]], lapply(ll_eta_k1, "*", h(p2, f_theta_eval[[k2]][[q]], ll_alpha[,k2]))))
      out[[q]] <- do.call("rbind",list_by_list(list_of_X_etaT[[k1]],
                                               list_take_list(part1[[q]],
                                                              lapply(ll_eta_k1, "*", h(p2, f_theta_eval[[k2]][[q]], ll_alpha[,k2])))))
    }
  }
  return(do.call("cbind", out))
}

### Eta eta tilde derivatives
ll_eta_etaT <- function(list_of_etaT_derivs, list_of_f_eval, list_of_densities, ll_alpha, list_of_X_eta, list_of_X_etaT) {
  K <- length(list_of_etaT_derivs)

  out = list()
  for (k1 in 1:(K-1)) {
    out[[k1]] <- list()
    if (is.null(list_of_etaT_derivs[[k1]])) {
      next
    }
    for (k2 in 1:K) {
      if (is.null(list_of_etaT_derivs[[k2]])) {
        next
      }
      if (k1 == k2) {
        part1 = list_of_etaT_derivs[[k1]]
      } else {
        part1 = 0
      }

      if (is.list(list_of_etaT_derivs[[k2]])) {
        part2 = lapply(list_of_etaT_derivs[[k2]], "*", h(list_of_densities[[k1]], list_of_f_eval[[k1]], ll_alpha[,k1]))

        if (is.list(part1)) {
          Xk1k2 <- list_take_list(part1, part2)
        } else {
          Xk1k2 <- lapply(part2, "-")
        }
        out_inner <- list()
        for (j in 1:length(part2)) {
          out_inner[[j]] <- XtDX(list_of_X_eta[[(k1)]], Xk1k2[[j]],list_of_X_etaT[[k2]][[j]])
        }
        out[[k1]][[k2]] <- do.call("cbind", out_inner)

      } else {
        part2 = h(list_of_densities[[k1]], list_of_f_eval[[k1]], ll_alpha[,k1]) * list_of_etaT_derivs[[k2]]

        out[[k1]][[k2]] <- XtDX(list_of_X_eta[[k1]], part1 - part2, list_of_X_etaT[[k2]])
      }
    }
  }
  return(out)
}

ll_eta_theta <- function(list_of_f_theta_eval, list_of_f_eval, list_of_densities, ll_alpha, list_of_X_eta) {
  # TODO: Make more efficient by removing eta1 calc
  K <- length(list_of_densities)
  N <- nrow(list_of_densities[[1]])
  out <- list()

  for (k1 in 1:(K-1)) {
    out[[k1]] <- list()
    if (is.null(list_of_f_theta_eval[[k1]])) {
      out[[k1]] <- list(NULL)
      next
    }
    for (k2 in 1:K) {
      out[[k1]][[k2]] <- list()
      if (is.null(list_of_f_theta_eval[[k2]])) {
        out[[k1]][k2] <- list(NULL)
        next
      }
      for (q in 1:length(list_of_f_theta_eval[[k2]])) {
        if (k1 == k2) {
          part1 = h(list_of_densities[[k1]], list_of_f_theta_eval[[k1]][[q]], ll_alpha[,k1])
        } else {
          part1 = 0
        }
        part2 =  h(list_of_densities[[k1]], list_of_f_eval[[k1]], ll_alpha[,k1]) * h(list_of_densities[[k2]], list_of_f_theta_eval[[k2]][[q]], ll_alpha[,k2])

        Xk1k2q = part1 - part2

        out[[k1]][[k2]][[q]] = t(list_of_X_eta[[(k1)]]) %*% Xk1k2q
      }
      out[[k1]][[k2]] <- do.call("cbind", out[[k1]][[k2]])
    }
    #out[[k1]] <- do.call("cbind", out[[k1]])
  }

  #out <- do.call("rbind", out)
  return(out)
}


re_order_eta <- function(ll_eta_theta, ll_eta_etaT) {
  K = length(ll_eta_theta)
  out = list()
  for (k in 1:K) {
    out[[k]] <- cbind(t(sapply(ll_eta_etaT, "[[", k)), t(sapply(ll_eta_theta, "[[", k)))
  }

  out <- do.call("cbind", out)
  return(out)
}

# ## Model derivatives, calcualtes all derivatives for parameters in models k1 and k2.
# ll_model_ij <- function(ll_alpha, list_of_densities, eval_store, k1, k2, list_of_X_etaT, list_of_ll_etaT) {
#   f_eta2T_eval <- get_eval("f_eta2_eval", eval_store)
#   f_theta2_eval <- get_eval("f_theta2_eval", eval_store)
#   f_eta_theta_eval <- get_eval("f_eta_theta_eval", eval_store)
#   f_theta_eval <- get_eval("f_theta_eval", eval_store)
#   f_eta_eval <- get_eval("f_eta_eval", eval_store)
#
#   ll_betaT2(ll_alpha, f_eta2T_eval, list_of_ll_etaT, 1, 1, list_of_densities, list_of_X_etaT)
#   etaTetaT <- diag(ll_etaT2(ll_alpha[,k1], f_eta2T_eval[[k1]], list_of_ll_etaT, k1, k2, list_of_densities[[k1]])[,1])
#   betaTbetaT <- etaeta_to_betabeta(etaTetaT, list_of_X_etaT, k1, k2)
#
#   etaTtheta <- ll_etaT_theta(ll_alpha, f_eta_theta_eval, f_eta_eval, f_theta_eval, list_of_densities, k1, k2, list_of_X_etaT)
#   betaTtheta <- eta_to_beta(etaTtheta, list_of_X_etaT[[k1]])[[1]]
#
#   thetaetaT <- ll_etaT_theta(ll_alpha, f_eta_theta_eval, f_eta_eval, f_theta_eval, list_of_densities, k2, k1)
#   thetabetaT <- t(eta_to_beta(list(thetaetaT), list(list_of_X_etaT[[k2]]))[[1]])
#
#   thetatheta <- ll_theta2(ll_alpha, f_theta2_eval, list_of_densities, f_theta_eval, k1, k2)
#
#   #print(ncol(cbind(betaTbetaT, betaTtheta)))
#   #print(ncol(cbind(thetabetaT, thetatheta)))
#   return(rbind(cbind(betaTbetaT, betaTtheta), cbind(thetabetaT, thetatheta)))
# }
