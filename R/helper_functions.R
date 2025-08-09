# TODO Move these to utils.R
##' Evaluate a set of inner derivatives on given parameters
#'
#' Evaluates a set of inner derivatives for a given set of parameters.
#'
#' @param inner_deriv_object Function. Output of any inner function.
#' @param eta Numeric vector or list. Linear predictor values.
#' @param theta Numeric vector. Extra parameters.
#' @param deriv Integer. Level of derivatives desired (0 = value, 1 = gradient, etc.).
#'
#' @return Named list of derivatives.
#' @export
eval_deriv <- function(inner_deriv_object, eta, theta, deriv = 0) {
  inner_deriv_object(eta, theta, deriv)
}

##' Calculate eta for given design matrices and coefficients
#'
#' Calculates eta (linear predictors) for a list of design matrices and coefficients.
#'
#' @param list_of_X_eta List of design matrices for eta.
#' @param list_of_beta List of coefficients for eta.
#'
#' @return List of eta vectors.
#' @export
get_list_of_eta <- function(list_of_X_eta, list_of_beta) {
  eta_list <- list()
  if (length(list_of_X_eta) == 0) {
    return(eta_list)
  }
  for (k in seq_along(list_of_X_eta)) {
    if (is.list(list_of_X_eta[[k]])) {
      eta_list[[k]] <- list()
  for (k_2 in seq_along(list_of_X_eta[[k]])) {
        # TODO: Check that beta has same list_of_list structure as eta
        eta_list[[k]][[k_2]] <- as.matrix(list_of_X_eta[[k]][[k_2]]) %*% list_of_beta[[k]][[k_2]]
      }
    } else if (is.null(list_of_X_eta[[k]])) {
      eta_list[k] <- list(NULL)
    } else {
      eta_list[[k]] <- as.matrix(list_of_X_eta[[k]]) %*% list_of_beta[[k]]
    }
  }
  return(eta_list)
}

##' Retrieve evaluated derivatives from storage object
#'
#' Fetches all evaluated derivatives of the specified type from an eval_store.
#'
#' @param eval_string Character. Name of derivative to retrieve (e.g., 'd0', 'd1').
#' @param eval_store List. Storage object created by eval_deriv.
#'
#' @return List of evaluated derivatives.
#' @export
get_eval <- function(eval_string, eval_store) {
  # Fetches all evaluated derivatives of the specified string from an eval_store
  lapply(eval_store, `[[`, eval_string)
}
