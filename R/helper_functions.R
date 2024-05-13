# TODO Move these to utils.R
#' Evaluates a set of inner derivatives on given parameters
#'
#' @param inner_deriv_object Output of any inner function
#' @param eta Vector of linear predictor values
#' @param theta Vector of extra parameters
#' @param deriv integer dictating what level of derivatives are desired
#'
#' @return Named list of derivatives
eval_deriv <- function(inner_deriv_object, eta, theta, deriv = 0) {
  store <- inner_deriv_object(eta, theta, deriv)
}

#' Calculates eta for given design matrices and coefficients
#'
#' @param list_of_X_eta List of design matrices for eta
#' @param list_of_beta List of coefficients for eta
#'
#' @return A list of eta vectors
get_list_of_eta <- function(list_of_X_eta, list_of_beta) {
  eta_list <- list()
  if (length(list_of_X_eta) == 0) {
    return(eta_list)
  }



  for (k in 1:length(list_of_X_eta)) {
    if (is.list(list_of_X_eta[[k]])) {
      eta_list[[k]] <- list()
      for (k_2 in 1:length(list_of_X_eta[[k]])) {
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

#' Retrieve evaluated derivatives from storage object
#'
#' @param eval_string String for name of derivative you want
#' @param eval_store A storage object created by eval_deriv
#'
#' @return A list of evaluated derivatives
get_eval <- function(eval_string, eval_store) {
  # Fetches all evaluated derivatives of the specified string from an eval_store
  lapply(eval_store, `[[`, eval_string)
}
