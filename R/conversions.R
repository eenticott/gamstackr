#' Convert eta parameters to alpha (weight) parameters
#'
#' This function converts a list of eta parameters to alpha parameters using the softmax function.
#' The transformation ensures that the resulting weights sum to 1.
#'
#' @param list_of_eta A list of eta parameters
#'
#' @return A matrix of alpha parameters (weights) that sum to 1 across rows
#'
eta_to_alpha <- function(list_of_eta) {
  if (length(list_of_eta) == 0) {
    return(matrix(1))
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

#' Convert eta derivatives to beta derivatives
#'
#' This function converts derivatives with respect to eta to derivatives with respect to beta
#' by applying the chain rule.
#'
#' @param list_of_ll_eta A list of derivatives with respect to eta
#' @param list_of_X_eta A list of design matrices
#' @param elementwise Logical; if TRUE, performs elementwise multiplication instead of matrix multiplication
#'
#' @return A list of derivatives with respect to beta
#'
eta_to_beta <- function(list_of_ll_eta, list_of_X_eta, elementwise = FALSE) {
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
      out[[i]] <- list_by_list(list_of_X_eta[[i]], list_of_ll_eta[[i]], elementwise)
    } else {
      if (elementwise) {
        out[[i]] <- list_of_X_eta[[i]] * list_of_ll_eta[[i]]
      } else {
        out[[i]] <- t(list_of_X_eta[[i]]) %*% list_of_ll_eta[[i]]
      }
    }
  }
  return(out)
}

#' Convert second derivatives with respect to eta to second derivatives with respect to beta
#'
#' This function converts second derivatives (Hessian) with respect to eta to second derivatives
#' with respect to beta by applying the chain rule.
#'
#' @param ll_etaeta A matrix of second derivatives with respect to eta
#' @param list_of_X_eta A list of design matrices
#' @param k1 Index of the first eta parameter
#' @param k2 Index of the second eta parameter
#'
#' @return A matrix of second derivatives with respect to beta
#'
etaeta_to_betabeta <- function(ll_etaeta, list_of_X_eta, k1, k2) {
  if (is.null(list_of_X_eta[[k1]]) | is.null(list_of_X_eta[[k2]])) {
    return(NULL)
  }
  return(as.matrix(t(list_of_X_eta[[k1]])) %*% ll_etaeta %*% as.matrix(list_of_X_eta[[k2]]))
}
