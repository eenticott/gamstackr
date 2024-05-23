#include <RcppArmadillo.h>
using namespace Rcpp;

List llk_gaussian(arma::vec y, arma::mat param, int deriv = 0) {
  if (param.n_cols != 2) {
    stop("Wrong number of parameters provided");
  }

  arma::vec mu = param.col(0);
  arma::vec tau = param.col(1);
  arma::vec tau2 = tau % tau;
  int n = y.n_elem;

  if (mu.n_elem == 1) {
    mu = NumericVector(n, mu(0));
    tau = NumericVector(n, tau(0));
    tau2 = NumericVector(n, tau2(0));
  }
  arma::vec ymu = y - mu;
  arma::vec ymu2 = ymu % ymu;
  arma::vec theta = 1 / tau2;
  arma::vec d0 = -0.5 * log(2 * M_PI) + log(tau) - 0.5 * tau2 % ymu2;

  List out;
  out["d0"] = d0;

  if (deriv > 0) {
    arma::vec d1 = tau2 % ymu;
    arma::vec d2 = (1 / tau - tau % ymu2) % (-0.5 * pow(theta, -1.5));
    out["d1"] = List::create(d1, d2);

    if (deriv > 1) {
      arma::vec d11 = -tau2;
      arma::vec d12 = (2 * d1 % tau) % (-0.5 * pow(theta, -1.5));
      arma::vec d22 = ((-ymu2 - 1 / tau2) % (0.25 * pow(theta, -3))) + ((1 / tau - tau % ymu2) % (0.75 * pow(theta, -2.5)));
      out["d2"] = List::create(d11, d12, d22);

      if (deriv > 2) {
        arma::vec zeros(n, 0.0);
        arma::vec d111 = zeros;
        arma::vec d112 = -2 * tau;
        arma::vec d122 = 2 * ymu;
        arma::vec d222 = 2 / pow(tau, 3);
        out["d3"] = List::create(d111, d112, d122, d222);

        if (deriv > 3) {
          arma::vec d1111 = zeros;
          arma::vec d1112 = zeros;
          arma::vec d1122 = NumericVector(n, -2.0);
          arma::vec d1222 = zeros;
          arma::vec d2222 = -6 / (tau2 * tau2);
          out["d3"] = List::create(d1111, d1112, d1122, d1222, d2222);
        }
      }
    }
  }
  return out;
}

// [[Rcpp::export]]
List get_derivs_cpp(arma::mat & eta, arma::vec & theta, int deriv, int n_k, int dim_num, arma::mat& x, arma::mat& dens_matrix, List store) {
  // Finds up to second derivatives of mvn weight method for given coefficients

  // Retrieve tau from theta vector
  arma::vec tau = exp(theta);
  int N = eta.n_rows;
  // Need to evaluate n_k * N * dim_num gaussian densities
  // NumericMatrix dens_matrix(n_k * N, dim_num);
   // Function llk_gaussian("llk_gaussian");
  arma::mat log_dens(N, n_k);

  for (int n = 0; n < dim_num; n++) {
    arma::mat xn = x.col(n);
    arma::mat etan = eta.col(n);
    arma::vec taun(1);
    taun(0) = tau(n);
    arma::vec res = arma::repelem(xn, N, 1);
    arma::mat param(res.n_elem, 2);
    param.col(0) = repmat(etan, n_k, 1);
    param.col(1) = arma::repelem(1/sqrt(taun), n_k * N, 1);
    List llk = llk_gaussian(res, param, 0);
    arma::vec d0 = llk["d0"];
    dens_matrix.col(n) = d0;
  }

  // Turn vector of densities back into matrix
  //log_dens = vec2mat(my_rowSums(dens_matrix), N, n_k);
  log_dens = arma::sum(dens_matrix, 1);
  log_dens.reshape(N, n_k);
  // Store shift matrix
  arma::vec rmax = max(log_dens, 1);
  arma::mat expshift = exp(log_dens.each_col() - rmax);

  arma::vec rs_expshift = arma::sum(expshift, 1);

  // Create output matrix for likelihood evaluation
  // NumericMatrix f_out(N, n_k);
  arma::mat f_out = expshift.each_col() / rs_expshift;
  // Evaluate each likelihood column
  // for (int j = 0; j < n_k; j++) {
  //   f_out(_,j) = expshift(_,j) / rs_expshift;
  // }

  // Store all outputs in a named list
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
      // auto x_n = x(_, n);
      // NumericVector res(x_n.size() * N);
      // for(int i = 0; i < x_n.size(); ++i) {
      //   for(int j = 0; j < N; ++j) {
      //     res[i * N + j] = x_n[i];
      //   }
      // }
      arma::mat xn = x.col(n);
      arma::mat etan = eta.col(n);
      arma::vec taun(1);
      taun(0) = tau(n);
      arma::vec res = arma::repelem(xn, N, 1);
      arma::mat param(res.n_elem, 2);
      param.col(0) = repmat(etan, n_k, 1);
      param.col(1) = arma::repelem(1/sqrt(taun), n_k * N, 1);
      List llk = llk_gaussian(res, param, 2);
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
    List f_theta_out(dim_num);

    for (int n = 0; n < dim_num; n++) {
      arma:: vec detasn = detas[n];
      arma:: vec dtausn = dtaus[n];

      arma::mat l1 = reshape(detasn, N, n_k);
      arma::mat l2 = reshape(dtausn, N, n_k);
      arma::mat l1diff = l1.each_col()-(sum(expshift % l1, 1)/rs_expshift);
      arma::mat l2diff = l2.each_col()-(sum(expshift % l2, 1)/rs_expshift);

      f_eta_out(n) = f_out % l1diff;
      f_tau_out(n) = f_out % l2diff;
      f_theta_out(n) = f_out % l2diff * tau(n);
    }

    store["f_eta_eval"] = f_eta_out;
    store["f_tau_eval"] = f_tau_out;
    store["f_theta_eval"] = f_theta_out;

    if (deriv >= 2) {
      List f_eta2_out(dim_num);
      //
      for (int alpha = 0; alpha < dim_num; alpha++) {
        List f_eta2_out_alpha(dim_num);
        arma::mat dalpha1 = detas(alpha);
        arma::mat dalpha2 = detaetas(alpha);
        dalpha1.reshape(N, n_k);
        dalpha2.reshape(N, n_k);

        for (int beta = 0; beta < dim_num; beta++) {
          arma::mat dbeta1 = detas(beta);
          dbeta1.reshape(N, n_k);
          arma::mat d2(N, n_k);
          if (alpha == beta) {
            d2 = dalpha2 + (dalpha1%dalpha1);
          } else {
            d2 = dalpha1 % dbeta1;
          }

          // Split calculation into 5 parts
          arma::mat p1 = d2 % f_out;
          arma::mat p2 = (f_out.each_col() % sum(dbeta1%expshift, 1)) % (dalpha1.each_col() / rs_expshift);
          arma::mat p3 = (f_out.each_col() % sum(d2%expshift, 1)).each_col() / rs_expshift;
          arma::mat p35 = (f_out.each_col() % sum(dalpha1%expshift, 1)) % (dbeta1.each_col() / rs_expshift);
          arma::mat p4 = (f_out.each_col() % (2*sum(dbeta1%expshift, 1) % sum(dalpha1%expshift, 1))).each_col()/(rs_expshift%rs_expshift);

          f_eta2_out_alpha(beta) = p1 - p2 - p3 + p4 - p35;
        }
        f_eta2_out[alpha] = f_eta2_out_alpha;
      }

      store["f_eta2_eval"] = f_eta2_out;

      List f_eta_theta_out(dim_num);

      for (int alpha = 0; alpha < dim_num; alpha++) {
        List f_eta_theta_out_alpha(dim_num);
        arma::mat dalpha1 = detas(alpha);
        arma::mat dalpha2 = detataus(alpha);
        dalpha1.reshape(N, n_k);
        dalpha2.reshape(N, n_k);

        for (int beta = 0; beta < dim_num; beta++) {
          arma::mat dbeta1 = dtaus(beta);
          dbeta1.reshape(N, n_k);
          arma::mat d2(N, n_k);
          if (alpha == beta) {
            d2 = dalpha2 + (dalpha1%dbeta1);
          } else {
            d2 = dalpha1 % dbeta1;
          }

          // Split calculation into 5 parts
          arma::mat p1 = d2 % f_out;
          arma::mat p2 = (f_out.each_col() % sum(dbeta1%expshift, 1)) % (dalpha1.each_col() / rs_expshift);
          arma::mat p3 = (f_out.each_col() % sum(d2%expshift, 1)).each_col() / rs_expshift;
          arma::mat p35 = (f_out.each_col() % sum(dalpha1%expshift, 1)) % (dbeta1.each_col() / rs_expshift);
          arma::mat p4 = (f_out.each_col() % (2*sum(dbeta1%expshift, 1) % sum(dalpha1%expshift, 1))).each_col()/(rs_expshift%rs_expshift);

          f_eta_theta_out_alpha(beta) = (p1 - p2 - p3 + p4 - p35)*tau(beta);
        }
        f_eta_theta_out[alpha] = f_eta_theta_out_alpha;
      }
      store["f_eta_theta_eval"] = f_eta_theta_out;

      List f_theta2_out(dim_num);
      for (int alpha = 0; alpha < dim_num; alpha++) {
        List f_theta2_out_alpha(dim_num);
        arma::mat dalpha1 = dtaus(alpha);
        arma::mat dalpha2 = dtautaus(alpha);
        dalpha1.reshape(N, n_k);
        dalpha2.reshape(N, n_k);
        for (int beta = 0; beta < dim_num; beta++) {
          arma::mat dbeta1 = dtaus(beta);
          dbeta1.reshape(N, n_k);
          arma::mat d2(N, n_k);
          if (alpha == beta) {
            d2 = dalpha2 + (dalpha1%dalpha1);
          } else {
            d2 = dalpha1 % dbeta1;
          }
          // Split calculation into 5 parts
          arma::mat p1 = d2 % f_out;
          arma::mat p2 = (f_out.each_col() % sum(dbeta1%expshift, 1)) % (dalpha1.each_col() / rs_expshift);
          arma::mat p3 = (f_out.each_col() % sum(d2%expshift, 1)).each_col() / rs_expshift;
          arma::mat p35 = (f_out.each_col() % sum(dalpha1%expshift, 1)) % (dbeta1.each_col() / rs_expshift);
          arma::mat p4 = (f_out.each_col() % (2*sum(dbeta1%expshift, 1) % sum(dalpha1%expshift, 1))).each_col()/(rs_expshift%rs_expshift);

          if (alpha == beta) {
            arma::mat fto = f_tau_out(alpha);
            f_theta2_out_alpha[beta] = (fto*tau(alpha))+(tau[alpha] * tau[alpha] * (p1 - p2 - p3 + p4 - p35));
          } else {
            f_theta2_out_alpha[beta] = (tau[beta] * tau[alpha] * (p1 - p2 - p3 + p4 - p35));
          }
        }
        f_theta2_out[alpha] = f_theta2_out_alpha;
      }

      store["f_theta2_eval"] = f_theta2_out;
    }
  }
  return store;
}
