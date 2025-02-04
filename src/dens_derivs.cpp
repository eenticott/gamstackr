#include <Rcpp.h>
using namespace Rcpp;

Environment pkg_env = Environment::namespace_env("gamstackr");

// Reference the function from the package's namespace
Function bind_mats = pkg_env["bind_list"];

// [[Rcpp::export]]
NumericVector log_rowSums_a_times_b_cpp(NumericMatrix log_a, NumericMatrix log_b) {
  int n = log_a.nrow();
  int m = log_a.ncol();
  NumericVector result(n);

  for (int i = 0; i < n; i++) {
    double max_val = -std::numeric_limits<double>::infinity();
    NumericVector temp_row(m);
    for (int j = 0; j < m; j++) {
      temp_row[j] = log_a(i, j) + log_b(i, j);
      if (temp_row[j] > max_val) {
        max_val = temp_row[j];
      }
    }
    double sum_exp = 0.0;
    for (int j = 0; j < m; j++) {
      sum_exp += std::exp(temp_row[j] - max_val);
    }
    result[i] = max_val + std::log(sum_exp);
  }
  return result;
}

// [[Rcpp::export]]
NumericVector rs_AB_cpp(NumericMatrix logA, NumericMatrix B) {
  int n = logA.nrow();
  int m = logA.ncol();
  NumericVector result(n);

  for (int i = 0; i < n; i++) {
    double row_sum = 0.0;
    for (int j = 0; j < m; j++) {
      double log_AB = logA(i, j) + std::log(std::abs(B(i, j)));
      row_sum += std::exp(log_AB) * ((B(i, j) > 0) - (B(i, j) < 0));
    }
    result[i] = row_sum;
  }
  return result;
}

// [[Rcpp::export]]
NumericVector ABoC(NumericVector AB, NumericVector C) {
  int n = AB.length();
  NumericVector result(n);

  for (int i = 0; i < n; i++) {
    result[i] = AB[i]/C[i];
  }
  return result;
}

// [[Rcpp::export]]
NumericVector AoB_cpp(NumericVector A, NumericVector logB) {
  int n = A.length();
  NumericVector result(n);

  for (int i = 0; i < n; i++) {
    double log_val = std::log(std::abs(A(i))) - logB(i);
    result(i) = std::exp(log_val) * ((A(i) > 0) - (A(i) < 0));
  }
  return result;
}

// [[Rcpp::export]]
NumericMatrix outer_prod(NumericVector x, NumericVector y) {
  int m = x.size();
  int n = y.size();
  NumericMatrix result(m, n);

  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      result(i, j) = x[i] * y[j];
    }
  }

  return result;
}

// [[Rcpp::export]]
List get_ll_dens_derivs_cpp(List list_of_beta,
                            List list_of_X,
                            NumericVector theta,
                            Function weight,
                            NumericMatrix log_dens,
                            Function beta_to_eta,
                            int deriv = 0) {
  int K = log_dens.ncol();
  int n = log_dens.nrow();
  int n_eta = as<int>(weight.attr("neta"));
  int n_theta = as<int>(weight.attr("ntheta"));
  NumericMatrix dens(n, K);

  // Pre-compute densities
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < K; j++) {
      dens(i, j) = std::exp(log_dens(i, j));
    }
  }

  // Calculate eta and W
  List eta = beta_to_eta(list_of_X, list_of_beta);
  List W = weight(eta, theta, 2);

  // Extract `f_eval` as NumericMatrix
  NumericMatrix W_f_eval = as<NumericMatrix>(W["f_eval"]);

  NumericMatrix logW(n, K);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < K; j++) {
      logW(i, j) = std::log(W_f_eval(i, j));
    }
  }

  NumericVector llk_deriv = log_rowSums_a_times_b_cpp(logW, log_dens);

  double ll = std::accumulate(llk_deriv.begin(), llk_deriv.end(), 0.0);

  // Initialize gradient and Hessian
  NumericVector grad;
  NumericMatrix hess;
  List dbb(n_eta);
  List dbt(n_eta);
  List dtt(n_theta);
  NumericMatrix dbb_out;
  NumericMatrix dbt_out;
  NumericMatrix dtt_out(n_theta,n_theta);
  if (deriv > 0) {
    NumericVector dt;
    List db(n_eta);

    if (n_eta > 0) {
      List f_eta_eval = as<List>(W["f_eta_eval"]);  // Extract f_eta_eval as a List
      for (int k = 0; k < n_eta; k++) {
        NumericMatrix xk = list_of_X[k];
        NumericVector dbk(xk.ncol());
        NumericMatrix f_eta_eval_k = as<NumericMatrix>(f_eta_eval[k]);  // Access each element as NumericMatrix
        NumericVector rs_AB_o_C_vals = AoB_cpp(rs_AB_cpp(log_dens, f_eta_eval_k), llk_deriv);  // Expecting a NumericVector here
        for (int i = 0; i < n; i++) {
          dbk += rs_AB_o_C_vals[i] * as<NumericMatrix>(list_of_X[k])(i, _);
        }
        db[k] = dbk;
      }
    }

    if (n_theta > 0) {
      dt = NumericVector(n_theta);
      List f_theta_eval = as<List>(W["f_theta_eval"]);  // Extract f_theta_eval as a List
      for (int k = 0; k < n_theta; k++) {
        NumericMatrix f_theta_eval_k = as<NumericMatrix>(f_theta_eval[k]);  // Access each element as NumericMatrix
        NumericVector rs_AB_o_C_vals = AoB_cpp(rs_AB_cpp(log_dens, f_theta_eval_k), llk_deriv);
        for (int i = 0; i < n; i++) {
          dt[k] += rs_AB_o_C_vals[i];
        }
      }
    }

    int n_beta = 0;

    for (int i = 0; i < db.size(); ++i) {
      NumericVector vec = db[i];
      n_beta += vec.size();
    }

    grad = NumericVector(n_beta + n_theta);
    int current_index = 0;
    for (int i = 0; i < db.size(); ++i) {
      NumericVector vec = db[i];
      std::copy(vec.begin(), vec.end(), grad.begin() + current_index);
      current_index += vec.size();  // Update the current index for the next insertion
    }
    for (int i = 0; i < n_theta; ++i) {
      grad[n_beta + i] = dt[i];
    }
    if (deriv > 1) {
      List f_eta_eval = as<List>(W["f_eta_eval"]);  // Extract f_eta_eval as a List
      List f_theta_eval = as<List>(W["f_theta_eval"]);  // Extract f_theta_eval as a List
      List f_eta2_eval = as<List>(W["f_eta2_eval"]);  // Extract f_eta2_eval as a List
      List f_theta2_eval = as<List>(W["f_theta2_eval"]);  // Extract f_eta2_eval as a List
      List f_eta_theta_eval = as<List>(W["f_eta_theta_eval"]);  // Extract f_eta2_eval as a List

      // dbb
      for (int i = 0; i < n_eta; i++) {
        List dbbi(n_eta);
        List f_eta2_eval_i = as<List>(f_eta2_eval[i]);  // Extract each sub-list
        NumericMatrix f_eta_eval_i = as<NumericMatrix>(f_eta_eval[i]);
        for (int j = 0; j < n_eta; j++) {
          NumericMatrix f_eta2_eval_ij = as<NumericMatrix>(f_eta2_eval_i[j]);
          NumericMatrix f_eta_eval_j = as<NumericMatrix>(f_eta_eval[j]);

          // Access as NumericMatrix
          NumericVector rs_AB_o_C_eta2ij = AoB_cpp(rs_AB_cpp(log_dens, f_eta2_eval_ij), llk_deriv);
          NumericVector rs_AB_o_C_etai = AoB_cpp(rs_AB_cpp(log_dens, f_eta_eval_i), llk_deriv);
          NumericVector rs_AB_o_C_etaj = AoB_cpp(rs_AB_cpp(log_dens, f_eta_eval_j), llk_deriv);
          NumericVector multiplier = rs_AB_o_C_eta2ij - (rs_AB_o_C_etai*rs_AB_o_C_etaj);
          NumericMatrix xmx = outer_prod(as<NumericMatrix>(list_of_X[i])(0, _), as<NumericMatrix>(list_of_X[j])(0, _));
          NumericMatrix mat_sum(xmx.nrow(), xmx.ncol());
          for (int k = 0; k < n; k++) {
            mat_sum += multiplier[k] *
              outer_prod(as<NumericMatrix>(list_of_X[i])(k, _), as<NumericMatrix>(list_of_X[j])(k, _));
          }
          dbbi[j] = mat_sum;
        }
        dbb[i] = dbbi;
      }
      dbb_out = bind_mats(dbb);

      // dbt
      for (int i = 0; i < n_eta; i++) {
        List dbti(n_theta);
        List f_eta_theta_eval_i = as<List>(f_eta_theta_eval[i]);  // Extract each sub-list
        NumericMatrix f_eta_eval_i = as<NumericMatrix>(f_eta_eval[i]);

        for (int j = 0; j < n_theta; j++) {

          NumericMatrix f_eta_theta_eval_ij = as<NumericMatrix>(f_eta_theta_eval_i[j]);

          NumericMatrix f_theta_eval_j = as<NumericMatrix>(f_theta_eval[j]);

          // Access as NumericMatrix
          NumericVector rs_AB_o_C_eta_thetaij = AoB_cpp(rs_AB_cpp(log_dens, f_eta_theta_eval_ij), llk_deriv);
          NumericVector rs_AB_o_C_etai = AoB_cpp(rs_AB_cpp(log_dens, f_eta_eval_i), llk_deriv);
          NumericVector rs_AB_o_C_thetaj = AoB_cpp(rs_AB_cpp(log_dens, f_theta_eval_j), llk_deriv);
          NumericVector multiplier = rs_AB_o_C_eta_thetaij - (rs_AB_o_C_etai*rs_AB_o_C_thetaj);

          NumericMatrix xi = list_of_X[i];
          NumericVector mat_sum(xi.ncol());

          for (int k = 0; k < n; k++) {
            mat_sum += multiplier[k] * xi(k, _);
          }
          dbti[j] = mat_sum;
        }
        dbt[i] = dbti;
      }
      dbt_out = bind_mats(dbt);

      // dtt
      for (int i = 0; i < n_theta; i++) {
        List dtti(n_theta);
        List f_theta2_eval_i = as<List>(f_theta2_eval[i]);  // Extract each sub-list
        NumericMatrix f_theta_eval_i = as<NumericMatrix>(f_theta_eval[i]);
        for (int j = 0; j < n_theta; j++) {
          NumericMatrix f_theta2_eval_ij = as<NumericMatrix>(f_theta2_eval_i[j]);
          NumericMatrix f_theta_eval_j = as<NumericMatrix>(f_theta_eval[j]);

          // Access as NumericMatrix
          NumericVector rs_AB_o_C_theta2ij = AoB_cpp(rs_AB_cpp(log_dens, f_theta2_eval_ij), llk_deriv);
          NumericVector rs_AB_o_C_thetai = AoB_cpp(rs_AB_cpp(log_dens, f_theta_eval_i), llk_deriv);
          NumericVector rs_AB_o_C_thetaj = AoB_cpp(rs_AB_cpp(log_dens, f_theta_eval_j), llk_deriv);
          NumericVector multiplier = rs_AB_o_C_theta2ij - (rs_AB_o_C_thetai*rs_AB_o_C_thetaj);
          NumericVector sums(1);
          for (int k = 0; k < n; k++) {
            sums = sums + multiplier[k];
          }
          dtti[j] = sums;
        }
        dtt[i] = dtti;
      }
      dtt_out = bind_mats(dtt);
    }
  }
  return List::create(Named("l") = ll,
                      Named("lb") = grad,
                      Named("dbb") = dbb_out,
                      Named("dbt") = dbt_out,
                      Named("dtt") = dtt_out);
}
