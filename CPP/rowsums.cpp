#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rowSumsRcpp(NumericMatrix x) {
  int nrow = x.nrow(), ncol = x.ncol();
  NumericVector out(nrow);

  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      out[j] += x(i, j);
    }
  }

  return out;
}