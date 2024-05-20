#include <RcppArmadillo.h>
using namespace Rcpp;

// Returns rowsums of a matrix as a vector
NumericVector my_rowSums(const NumericMatrix& x) {
  int nr = x.nrow(), nc = x.ncol();
  NumericVector ans(nr);
  for (int j = 0; j < nc; j++) {
    for (int i = 0; i < nr; i++) {
      ans[i] += x(i, j);
    }
  }
  return ans;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix vec2mat(Rcpp::NumericVector vec, int n, int k) {
  // turns a vector back into a matrix
  Rcpp::NumericMatrix mat = Rcpp::no_init(n, k);
  for (auto i = 0; i < n * k; i++) mat[i] = vec[i];
  return mat;
}

// [[Rcpp::export]]
NumericVector rowMaxs(NumericMatrix x) {
  // Returns vector of max of each row
  int nrow = x.nrow();
  NumericVector maxs(nrow);

  // Convert NumericMatrix to arma::mat
  arma::mat x_arma = as<arma::mat>(x);

  // Compute row-wise maximums
  maxs = arma::max(x_arma, 1);
  NumericVector maxr = maxs;
  return maxr;
}

// [[Rcpp::export]]
NumericMatrix matbyvec(NumericMatrix x, NumericVector v) {
  // Multiplies each column of a matrix elementwise by a vector
  // Get dimensions required for output
  int nrow = x.nrow();
  int ncol = x.ncol();
  NumericMatrix out(nrow, ncol);

  for (int j = 0; j < ncol; j++) {
    out(_,j) = x(_,j) * v;
  }

  return out;
}

// [[Rcpp::export]]
NumericMatrix mattakevec(NumericMatrix x, NumericVector v) {
  // Takes a vector away from each column of a matrix elementwise
  int nrow = x.nrow();
  int ncol = x.ncol();
  NumericMatrix out(nrow, ncol);

  for (int j = 0; j < ncol; j++) {
    out(_,j) = x(_,j) - v;
  }

  return out;
}

// [[Rcpp::export]]
List get_derivs_cpp(NumericMatrix eta, NumericVector theta, int deriv, int n_k, int dim_num, NumericMatrix x) {
  // Finds up to second derivatives of mvn weight method for given coefficients

  // Retrieve tau from theta vector
  NumericVector tau = exp(theta);
  int N = eta.nrow();
  // Need to evaluate n_k * N * dim_num gaussian densities
  NumericMatrix dens_matrix(n_k * N, dim_num);
  Function llk_gaussian("llk_gaussian");
  NumericMatrix log_dens(N, n_k);

  for (int n = 0; n < dim_num; n++) {
    auto x_n = x(_, n);
    NumericVector res(x_n.size() * N);
    for(int i = 0; i < x_n.size(); ++i) {
      for(int j = 0; j < N; ++j) {
        res[i * N + j] = x_n[i];
      }
    }

    List llk = llk_gaussian(res, List::create(rep(eta(_, n), n_k), 1 / sqrt(tau[n])), 0);
    NumericVector d0 = llk["d0"];
    dens_matrix(_, n) = d0;
  }

  // Turn vector of densities back into matrix
  log_dens = vec2mat(my_rowSums(dens_matrix), N, n_k);

  // Store shift matrix
  NumericMatrix expshift(N, n_k);
  NumericVector rmax = rowMaxs(log_dens);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < n_k; j++) {
      expshift(i, j) = exp(log_dens(i,j) - rmax(i));
    }
  }

  NumericVector rs_expshift = rowSums(expshift);

  // Create output matrix for likelihood evaluation
  NumericMatrix f_out(N, n_k);

  // Evaluate each likelihood column
  for (int j = 0; j < n_k; j++) {
    f_out(_,j) = expshift(_,j) / rs_expshift;
  }

  // Store all outputs in a named list
  List store;

  store["f_eval"] = f_out;

  // If deriv>=1 we evaluate first derivatives
  if (deriv >= 1) {
    // Create all of the lists that we need to store gaussian derivatives
    List detas(dim_num);
    List dtaus(dim_num);
    List detaetas(dim_num);
    List dtautaus(dim_num);
    List detataus(dim_num);

    for (int n = 0; n < dim_num; n++) {
      // rep works different in cpp, this code replicates res = rep(x, each = length(x_n))
      auto x_n = x(_, n);
      NumericVector res(x_n.size() * N);
      for(int i = 0; i < x_n.size(); ++i) {
        for(int j = 0; j < N; ++j) {
          res[i * N + j] = x_n[i];
        }
      }

      // Evaluate likelihood densities using given R function
      List llk = llk_gaussian(res, List::create(rep(eta(_, n), n_k), 1 / sqrt(tau[n])), 2);
      // Store all outputs in a list
      List d1 = llk["d1"];
      List d2 = llk["d2"];
      detas[n] = d1[0];
      dtaus[n] = d1[1];
      detaetas[n] = d2[0];
      detataus[n] = d2[1];
      dtautaus[n] = d2[2];
    }

    List f_eta_out(dim_num);
    List f_tau_out(dim_num);

    for (int n = 0; n < dim_num; n++) {
      NumericMatrix l1 = (vec2mat(detas[n], N, n_k));
      NumericMatrix l2 = (vec2mat(dtaus[n], N, n_k));

      NumericMatrix expl1 = vec2mat(expshift*l1, N, n_k);
      NumericVector rs1ors = my_rowSums(expl1)/rs_expshift;
      NumericMatrix l1diff = mattakevec(l1,rs1ors);

      NumericMatrix expl2 = vec2mat(expshift*l2, N, n_k);
      NumericVector rs2ors = my_rowSums(expl2)/rs_expshift;
      NumericMatrix l2diff = mattakevec(l2,rs2ors);

      NumericMatrix out1(N, n_k);
      NumericMatrix out2(N, n_k);

      for (int j = 0; j < n_k; j++) {
        out1(_,j) = f_out(_,j) * l1diff(_,j);
        out2(_,j) = f_out(_,j) * l2diff(_,j);
      }

      f_eta_out[n] = out1;
      f_tau_out[n] = out2;
    }

    store["f_eta_eval"] = f_eta_out;
    store["f_tau_eval"] = f_tau_out;
    List f_theta_out(dim_num);
    for (int i = 0; i < dim_num; i++) {
      NumericMatrix fti = f_tau_out[i];
      f_theta_out[i] = vec2mat(fti * rep(tau[i], fti.length()), N, n_k);
    }
    store["f_theta_eval"] = f_theta_out;

    if (deriv >= 2) {
      List f_eta2_out(dim_num);
      //
      for (int alpha = 0; alpha < dim_num; alpha++) {
        List f_eta2_out_alpha(dim_num);
        NumericMatrix dalpha1 = vec2mat(detas[alpha], N, n_k);
        NumericMatrix dalpha2 = vec2mat(detaetas[alpha], N, n_k);
        for (int beta = 0; beta < dim_num; beta++) {
          NumericMatrix dbeta1 = vec2mat(detas[beta], N, n_k);
          NumericMatrix d2(N, n_k);
          if (alpha == beta) {
            d2 = vec2mat(dalpha2 + (dalpha1*dalpha1), N, n_k);
          } else {
            d2 = vec2mat(dalpha1 * dbeta1, N, n_k);
          }

          // Split calculation into 5 parts
          NumericMatrix p1 = vec2mat(d2 * f_out, N, n_k);
          NumericMatrix p2 = vec2mat(matbyvec(f_out, my_rowSums(vec2mat(dbeta1 * expshift, N, n_k))) *
            matbyvec(dalpha1, 1/ rs_expshift), N, n_k);
          NumericMatrix p3 = matbyvec(vec2mat(matbyvec(f_out,my_rowSums(vec2mat(d2*expshift, N, n_k))), N, n_k), 1/rs_expshift);
          NumericMatrix p35 = vec2mat(matbyvec(f_out, my_rowSums(vec2mat(dalpha1*expshift, N, n_k))) * matbyvec(dbeta1, 1/rs_expshift), N, n_k);
          NumericMatrix p4 = matbyvec(f_out, 2 * (my_rowSums(vec2mat(dbeta1*expshift, N, n_k)) *
            my_rowSums(vec2mat(dalpha1*expshift, N, n_k)) * 1/(rs_expshift*rs_expshift)));

          f_eta2_out_alpha[beta] = vec2mat(p1 - p2 - p3 + p4 - p35, N, n_k);
        }
        f_eta2_out[alpha] = f_eta2_out_alpha;
      }

      store["f_eta2_eval"] = f_eta2_out;

      List f_eta_theta_out(dim_num);
      for (int alpha = 0; alpha < dim_num; alpha++) {
        List f_eta_theta_out_alpha(dim_num);
        NumericMatrix dalpha1 = vec2mat(detas[alpha], N, n_k);
        NumericMatrix dalpha2 = vec2mat(detataus[alpha], N, n_k);
        for (int beta = 0; beta < dim_num; beta++) {
          NumericMatrix dbeta1 = vec2mat(dtaus[beta], N, n_k);
          NumericMatrix d2(N, n_k);
          if (alpha == beta) {
            d2 = vec2mat(dalpha2 + (dalpha1*dbeta1), N, n_k);
          } else {
            d2 = vec2mat(dalpha1 * dbeta1, N, n_k);
          }
          NumericMatrix p1 = vec2mat(d2 * f_out, N, n_k);
          // Rcout << "At p1" << p1 << std::endl;
          NumericMatrix p2 = vec2mat(matbyvec(f_out, my_rowSums(vec2mat(dbeta1 * expshift, N, n_k))) *
            matbyvec(dalpha1, 1/ rs_expshift), N, n_k);
          // Rcout << "At p2" << p2 << std::endl;
          NumericMatrix p3 = matbyvec(vec2mat(matbyvec(f_out,my_rowSums(vec2mat(d2*expshift, N, n_k))), N, n_k), 1/rs_expshift);
          //Rcout << "At p3" << p3 << std::endl;
          NumericMatrix p35 = vec2mat(matbyvec(f_out, my_rowSums(vec2mat(dalpha1*expshift, N, n_k))) * matbyvec(dbeta1, 1/rs_expshift), N, n_k);
          //vec2mat(vec2mat(my_rowSums(vec2mat(dalpha1*expshift, N, n_k)) * f_out, N, n_k) *
          //matbyvec(dbeta1, 1/rs_expshift), N, n_k);

          //Rcout << "At p35" << p35 << std::endl;
          NumericMatrix p4 = matbyvec(f_out, 2 * (my_rowSums(vec2mat(dbeta1*expshift, N, n_k)) *
            my_rowSums(vec2mat(dalpha1*expshift, N, n_k)) * 1/(rs_expshift*rs_expshift)));

          //vec2mat(vec2mat(my_rowSums(vec2mat(dbeta1*expshift, N, n_k)) *
          //  my_rowSums(vec2mat(dalpha1*expshift, N, n_k)), N, n_k) * 2 * matbyvec(f_out, 1/(rs_expshift*rs_expshift)), N, n_k);
          //Rcout << "At p4" << p4 << std::endl;

          f_eta_theta_out_alpha[beta] = vec2mat(tau[beta] * (p1 - p2 - p3 + p4 - p35), N, n_k);
        }
        f_eta_theta_out[alpha] = f_eta_theta_out_alpha;
        // Rcout << "Stored first" << std::endl;
      }

      store["f_eta_theta_eval"] = f_eta_theta_out;

      List f_theta2_out(dim_num);
      for (int alpha = 0; alpha < dim_num; alpha++) {
        List f_theta2_out_alpha(dim_num);
        NumericMatrix dalpha1 = vec2mat(dtaus[alpha], N, n_k);
        NumericMatrix dalpha2 = vec2mat(dtautaus[alpha], N, n_k);
        for (int beta = 0; beta < dim_num; beta++) {
          NumericMatrix dbeta1 = vec2mat(dtaus[beta], N, n_k);
          NumericMatrix d2(N, n_k);
          if (alpha == beta) {
            d2 = vec2mat(dalpha2 + (dalpha1*dbeta1), N, n_k);
          } else {
            d2 = vec2mat(dalpha1 * dbeta1, N, n_k);
          }
          NumericMatrix p1 = vec2mat(d2 * f_out, N, n_k);
          // Rcout << "At p1" << p1 << std::endl;
          NumericMatrix p2 = vec2mat(matbyvec(f_out, my_rowSums(vec2mat(dbeta1 * expshift, N, n_k))) *
            matbyvec(dalpha1, 1/ rs_expshift), N, n_k);
          // Rcout << "At p2" << p2 << std::endl;
          NumericMatrix p3 = matbyvec(vec2mat(matbyvec(f_out,my_rowSums(vec2mat(d2*expshift, N, n_k))), N, n_k), 1/rs_expshift);
          //Rcout << "At p3" << p3 << std::endl;
          NumericMatrix p35 = vec2mat(matbyvec(f_out, my_rowSums(vec2mat(dalpha1*expshift, N, n_k))) * matbyvec(dbeta1, 1/rs_expshift), N, n_k);
          //vec2mat(vec2mat(my_rowSums(vec2mat(dalpha1*expshift, N, n_k)) * f_out, N, n_k) *
          //matbyvec(dbeta1, 1/rs_expshift), N, n_k);

          //Rcout << "At p35" << p35 << std::endl;
          NumericMatrix p4 = matbyvec(f_out, 2 * (my_rowSums(vec2mat(dbeta1*expshift, N, n_k)) *
            my_rowSums(vec2mat(dalpha1*expshift, N, n_k)) * 1/(rs_expshift*rs_expshift)));

          //vec2mat(vec2mat(my_rowSums(vec2mat(dbeta1*expshift, N, n_k)) *
          //  my_rowSums(vec2mat(dalpha1*expshift, N, n_k)), N, n_k) * 2 * matbyvec(f_out, 1/(rs_expshift*rs_expshift)), N, n_k);
          //Rcout << "At p4" << p4 << std::endl;
          if (alpha == beta) {
            NumericMatrix f1 = f_tau_out[alpha];
            NumericVector f2 = rep(tau[alpha], f1.length());
            NumericMatrix part1 = vec2mat(f1 * f2, N, n_k);
            NumericMatrix part2 = vec2mat((tau[alpha] * tau[alpha] * (p1 - p2 - p3 + p4 - p35)), N, n_k);

            f_theta2_out_alpha[beta] = vec2mat(part1 + part2, N, n_k);
          } else {
            f_theta2_out_alpha[beta] = vec2mat(tau[beta] * tau[alpha] * (p1 - p2 - p3 + p4 - p35), N, n_k);
          }
        }
        f_theta2_out[alpha] = f_theta2_out_alpha;
        // Rcout << "Stored first" << std::endl;
      }

      store["f_theta2_eval"] = f_theta2_out;
    }
  }
  return store;
}
