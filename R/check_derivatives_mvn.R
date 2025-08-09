# Helper functions for numerical derivative validation
library(numDeriv)

#' Compare analytical and numerical derivatives for MVN weights
#' @param mvn_func MVN weight function to test
#' @param eta Matrix of eta values
#' @param theta Vector of theta values
#' @param dim_num Number of dimensions
#' @param k Weight index to check
#' @param tol Tolerance for comparison
#' @param print_derivs Print derivatives for debugging
#' @return List indicating if derivatives match numerical approximations
check_derivatives_mvn <- function(mvn_func, eta, theta, dim_num, k = 1, tol = 1e-8, print_derivs = FALSE) {
  N <- nrow(eta)
  neta <- attr(mvn_func, "neta")
  ntheta <- attr(mvn_func, "ntheta")
  
  # Calculate weights and derivatives
  out <- mvn_func(eta, theta, deriv = 2)
  
  # Function that returns kth weight for numerical derivatives
  tmp_fun <- function(pars) {
    if (neta > 0) {
      eta_tmp = matrix(pars[1:neta], nrow = 1)
    }
    if (ntheta > 0) {
      theta_tmp = pars[(neta+1):(ntheta+neta)]
    }
    mvn_func(eta_tmp, theta_tmp, deriv = 0)$f_eval[k]
  }
  
  # First derivatives
  num_grad <- grad(tmp_fun, c(eta[1,], theta))
  
  eta_grad <- NULL
  theta_grad <- NULL
  
  if (neta > 0 && !is.null(out$f_eta_eval)) {
    eta_grad <- tryCatch({
      do.call("c", lapply(out$f_eta_eval, "[[", k))
    }, error = function(e) NULL)
    
    if (print_derivs && !is.null(eta_grad)) {
      print(paste("Numerical eta grad", paste(num_grad[1:neta], collapse = " ")))
      print(paste("Calculated eta grad", paste(eta_grad, collapse = " ")))
    }
  }
  
  if (ntheta > 0 && !is.null(out$f_theta_eval)) {
    theta_grad <- tryCatch({
      do.call("c", lapply(out$f_theta_eval, "[[", k))
    }, error = function(e) NULL)
    
    if (print_derivs && !is.null(theta_grad)) {
      print(paste("Numerical theta grad", paste(num_grad[(neta+1):(ntheta+neta)], collapse = " ")))
      print(paste("Calculated theta grad", paste(theta_grad, collapse = " ")))
    }
  }
  
  calc_grad <- c(eta_grad, theta_grad)
  grad_correct <- mean(abs(calc_grad - num_grad)) < tol
  
  # Second derivatives
  num_hess <- hessian(tmp_fun, c(eta[1,], theta))
  
  hess_etaeta_correct <- NULL
  hess_etatheta_correct <- NULL
  hess_thetatheta_correct <- NULL
  
  if (neta > 0 && !is.null(out$f_eta2_eval)) {
    tryCatch({
      calc_etaeta <- matrix(unlist(lapply(out$f_eta2_eval, function(i) {
        if(is.null(i)) return(NULL)
        lapply(i, "[[", k)
      })), neta, neta)
      if(!is.null(calc_etaeta)) {
        num_etaeta <- num_hess[1:neta, 1:neta]
        hess_etaeta_correct <- mean(abs(calc_etaeta - num_etaeta)) < tol
        
        if (print_derivs) {
          print("Eta-Eta Hessian:")
          print("Numerical:")
          print(num_etaeta)
          print("Analytical:")
          print(calc_etaeta)
        }
      }
    }, error = function(e) {
      hess_etaeta_correct <- NULL
    })
    
    if (ntheta > 0 && !is.null(out$f_eta_theta_eval)) {
      tryCatch({
        calc_etatheta <- matrix(unlist(lapply(out$f_eta_theta_eval, function(i) {
          if(is.null(i)) return(NULL)
          lapply(i, "[[", k)
        })), neta, ntheta, byrow = TRUE)
        if(!is.null(calc_etatheta)) {
          num_etatheta <- num_hess[1:neta, (neta+1):(neta+ntheta)]
          hess_etatheta_correct <- mean(abs(calc_etatheta - num_etatheta)) < tol
          
          if (print_derivs) {
            print("Eta-Theta Hessian:")
            print("Numerical:")
            print(num_etatheta)
            print("Analytical:")
            print(calc_etatheta)
          }
        }
      }, error = function(e) {
        hess_etatheta_correct <- NULL
      })
    }
  }
  
  if (ntheta > 0 && !is.null(out$f_theta2_eval)) {
    tryCatch({
      calc_thetatheta <- matrix(unlist(lapply(out$f_theta2_eval, function(i) {
        if(is.null(i)) return(NULL)
        lapply(i, "[[", k)
      })), ntheta, ntheta)
      if(!is.null(calc_thetatheta)) {
        num_thetatheta <- num_hess[(neta+1):(neta+ntheta), (neta+1):(neta+ntheta)]
        hess_thetatheta_correct <- mean(abs(calc_thetatheta - num_thetatheta)) < tol
        
        if (print_derivs) {
          print("Theta-Theta Hessian:")
          print("Numerical:")
          print(num_thetatheta)
          print("Analytical:")
          print(calc_thetatheta)
        }
      }
    }, error = function(e) {
      hess_thetatheta_correct <- NULL
    })
  }
  
  # Return results
  results <- list(
    grad_correct = grad_correct,
    hess_etaeta_correct = hess_etaeta_correct,
    hess_etatheta_correct = hess_etatheta_correct,
    hess_thetatheta_correct = hess_thetatheta_correct,
    max_diff_grad = mean(abs(calc_grad - num_grad))
  )
  
  if (neta > 0) {
    results$max_diff_etaeta <- if(!is.null(hess_etaeta_correct)) mean(abs(calc_etaeta - num_etaeta)) else NULL
    results$max_diff_etatheta <- if(!is.null(hess_etatheta_correct)) mean(abs(calc_etatheta - num_etatheta)) else NULL
  }
  
  if (ntheta > 0) {
    results$max_diff_thetatheta <- if(!is.null(hess_thetatheta_correct)) mean(abs(calc_thetatheta - num_thetatheta)) else NULL
  }
  
  return(results)
}
