
#' Convert a list of eta values to alpha weights (softmax transformation)
#'
#' @param list_of_eta List of eta vectors (linear predictors)
#'
#' @return Matrix of alpha weights (rows sum to 1)
#' @export
#'
#' @examples
#' eta_to_alpha(list(matrix(rnorm(10), ncol=1)))
eta_to_alpha <- function(list_of_eta) {
  if (length(list_of_eta) == 0) {
    return(matrix(1))
  }
  # Convert list of eta to matrix
  eta <- matrix(unlist(list_of_eta), ncol = length(list_of_eta))
  # Add 0 column for fixed first eta (identifiability)
  eta1 <- cbind(0, eta)
  # Centre the matrix for numerical stability
  etaCen <- eta1 - apply(eta1, 1, max)
  exp_etaCen <- exp(etaCen)
  out <- exp_etaCen / rowSums(exp_etaCen)
  return(out)
}


#' Convert eta derivatives to beta derivatives
#'
#' @param list_of_ll_eta List of derivatives with respect to eta
#' @param list_of_X_eta List of design matrices for eta
#' @param elementwise Logical, if TRUE multiplies elementwise, else matrix product
#'
#' @return List of beta derivatives
#' @export
eta_to_beta <- function(list_of_ll_eta, list_of_X_eta, elementwise = FALSE) {
  K <- length(list_of_ll_eta)
  if (K == 0) {
    return(NULL)
  }
  out <- list()
  for (i in 1:K) {
    if (is.null(list_of_ll_eta[[i]])) {
      out[i] <- list(NULL)
    } else if (is.list(list_of_ll_eta[[i]])) {
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


#' Convert eta-eta second derivatives to beta-beta second derivatives
#'
#' @param ll_etaeta Second derivative matrix with respect to eta
#' @param list_of_X_eta List of design matrices for eta
#' @param k1 Index for first eta
#' @param k2 Index for second eta
#'
#' @return Matrix of beta-beta second derivatives
#' @export
etaeta_to_betabeta <- function(ll_etaeta, list_of_X_eta, k1, k2) {
  if (is.null(list_of_X_eta[[k1]]) || is.null(list_of_X_eta[[k2]])) {
    return(NULL)
  }
  return(as.matrix(t(list_of_X_eta[[k1]])) %*% ll_etaeta %*% as.matrix(list_of_X_eta[[k2]]))
}
