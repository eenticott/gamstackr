#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec rs_AB_cpp2(const arma::vec& exp_log_a_max, const arma::mat& exp_log_a_diff, const arma::mat& b) {
  // Perform row-wise multiplication and summation
  arma::vec row_sums = arma::sum(exp_log_a_diff % b, 1);

  // Element-wise multiplication with exp_log_a_max
  arma::vec result = exp_log_a_max % row_sums;

  return result;
}


// [[Rcpp::export]]
NumericVector cpp_AoB(NumericVector A, NumericVector logB) {
  int n = A.size();
  NumericVector result(n);

  for (int i = 0; i < n; i++) {
    result[i] = exp(log(fabs(A[i])) - logB[i]) * (A[i] < 0 ? -1 : 1);
  }

  return result;
}

// [[Rcpp::export]]
NumericVector stable_rowsum_div_cpp(NumericMatrix logA, NumericMatrix B, NumericVector logC) {
  int n = logA.nrow();  // Number of rows
  int m = logA.ncol();  // Number of columns
  NumericVector result(n);

  for (int i = 0; i < n; i++) {
    double max_log_sum = -std::numeric_limits<double>::infinity();
    NumericVector log_sum_exp_vals(m);

    // Compute log(A * B) and store in a temporary vector
    for (int j = 0; j < m; j++) {
      log_sum_exp_vals[j] = logA(i, j) + log(B(i, j));
      if (log_sum_exp_vals[j] > max_log_sum) {
        max_log_sum = log_sum_exp_vals[j]; // Track the maximum value for stability
      }
    }

    // Calculate stable row-wise log-sum-exp
    double sum_exp = 0.0;
    for (int j = 0; j < m; j++) {
      sum_exp += std::exp(log_sum_exp_vals[j] - max_log_sum);
    }
    double row_log_sum_exp = max_log_sum + std::log(sum_exp);

    // Final result: log(rowsums(A * B) / C)
    result[i] = std::exp(row_log_sum_exp - logC[i]);
  }

  return result;
}
