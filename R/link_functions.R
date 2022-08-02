log_link <- function(x) {
  out <- list()
  out$f <- exp(x)
  out$f1 <- exp(x)
  out$f2 <- exp(x)
  return(out)
}

id_link <- function(x) {
  out <- list()
  out$f <- x
  out$f1 <- 1
  out$f2 <- 0
  return(out)
}

link_theta <- function(f, link_func = "log", eta, tau) {
  if (link_func == "log") {
    theta_eval <- f(eta, exp(tau))
    tau_eval <- theta_eval * exp(tau)
    return(tau_eval)
  }
}
