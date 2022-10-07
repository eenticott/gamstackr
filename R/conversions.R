eta_to_alpha <- function(list_of_eta) {
  if (length(list_of_eta) == 0) {
    return(1)
  }
  # exp_eta_matrix <- cbind(1, exp(matrix(unlist(list_of_eta), ncol = length(list_of_eta))))
  # out <- exp_eta_matrix / rowSums(exp_eta_matrix)
  eta <- matrix(unlist(list_of_eta), ncol = length(list_of_eta)) #turn eta into matrix
  eta1 <- cbind(0, eta) # add 0 column for fixed first eta
  etaCen <- eta1 - apply(eta1, 1, max) # Centre the matrix
  exp_etaCen <- exp(etaCen)
  out <- exp_etaCen / rowSums(exp_etaCen)
  return(out)
}

eta_to_beta <- function(list_of_ll_eta, list_of_X_eta) {
  K <- length(list_of_ll_eta)
  if (K == 0) {
    return(NULL)
  }
  out <- list()
  for (i in 1:K) {
    if (is.null(list_of_ll_eta[[i]])) {
      out[i] <- list(NULL)
    } else
    if (is.list(list_of_ll_eta[[i]])) {
      out[[i]] <- list_by_list(list_of_X_eta[[i]], list_of_ll_eta[[i]])
    } else {
      out[[i]] <- t(list_of_X_eta[[i]]) %*% list_of_ll_eta[[i]]
    }
  }
  return(out)
}

etaeta_to_betabeta <- function(ll_etaeta, list_of_X_eta, k1, k2) {
  if (is.null(list_of_X_eta[[k1]]) | is.null(list_of_X_eta[[k2]])) {
    return(NULL)
  }
  return(as.matrix(t(list_of_X_eta[[k1]])) %*% ll_etaeta %*% as.matrix(list_of_X_eta[[k2]]))
}
