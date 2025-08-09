Latent_fast <- function(K, q) {
  nv <- (K-1-q)*q
  neta <- q
  ntheta <- nv + K-1

  get_derivs <- function(eta, theta, deriv = 0) {
    if (is.list(eta)) eta <- do.call("cbind", eta)
    beta0 <- theta[1:(K-1)]
    N <- nrow(eta)
    v_mat <- rbind(diag(q), matrix(theta[K:length(theta)], ncol = q))
    eta_aug <- (eta %*% t(v_mat)) + matrix(beta0, ncol = (K-1), nrow = N, byrow = TRUE)

    aug_w <- multinomial(K, normalised = TRUE)
    aug_deriv <- aug_w(eta_aug, 0, deriv = deriv)

    store <- list()
    store$f_eval <- aug_deriv$f_eval

    if (deriv >= 1) {
      # Stack all aug_deriv$f_eta_eval into a 3D array: [N, K, K-1]
      aug_f_eta_eval <- array(
        unlist(aug_deriv$f_eta_eval, recursive = FALSE),
        dim = c(N, K, K-1)
      )
      # For each j in 1:neta, f_eta_eval[[j]] = sum_h aug_f_eta_eval[,,h] * v_mat[h,j]
      store$f_eta_eval <- lapply(seq_len(neta), function(j) {
        # Matrix multiplication over the 3rd dimension
        # Result: N x K matrix
        mat <- matrix(0, nrow = N, ncol = K)
        for (h in 1:(K-1)) {
          mat <- mat + aug_f_eta_eval[,,h] * v_mat[h,j]
        }
        mat
      })
      # f_theta_eval and higher derivatives as before
      store$f_theta_eval <- lapply(1:ntheta, function(i) matrix(NA, nrow = N, ncol = K))
      for (i in 1:(K-1)) {
        store$f_theta_eval[[i]] <- aug_deriv$f_eta_eval[[i]]
      }
      vj <- rep(1:q, times = (K-1-q))
      for (i in (K:ntheta)) {
        vx <- i - (K-1)
        ix <- floor(q+1+(i-K)/q)
        jx <- vj[vx]
        store$f_theta_eval[[i]] <- aug_deriv$f_eta_eval[[ix]] * eta[,jx]
      }
    }
    store
  }

  init_func <- function(scores) {
    aug_weights <- multinomial(K)
    init_aug_eta <- attr(aug_weights, "init_func")(scores)$init_mu
    init_eta <- init_aug_eta[,1:neta]
    b0s <- rep(0, K-1)
    vs <- list()
    for (j in (neta+1):(K-1)) {
      y1 <- init_aug_eta[,j]
      X1 <- init_eta
      coefs <- coef(lm(y1~X1))
      b0s[j] <- coefs[1]
      vs[[j-neta]] <- coefs[-1]
    }
    vs <- do.call("c", vs)
    init_theta <- c(b0s, vs)
    list(init_theta=init_theta,init_mu=init_eta)
  }

  pen <- function(tau, deriv = 0) {
    p <- 0
    pt <- rep(0, ntheta)
    ptt <- matrix(0, nrow = ntheta, ncol = ntheta)
    list(p = p, pt = pt, ptt = ptt)
  }

  structure(get_derivs,
            neta=neta,
            ntheta=ntheta,
            theta_pen=pen,
            init_func=init_func,
            num_weights=K, name = "Latent_faster")
}