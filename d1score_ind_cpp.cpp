#include <Rcpp.h>
#include <unordered_map>
#include <algorithm>  // for sort

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix d1score_ind_cpp(NumericVector beta, NumericMatrix X, NumericVector y, IntegerVector period) {
  int n = X.nrow();
  int p = X.ncol();
  
  // Calculate ex = exp(X %*% beta)
  NumericVector ex(n);
  for (int i = 0; i < n; ++i) {
    double dot_product = 0.0;
    for (int j = 0; j < p; ++j) {
      dot_product += X(i, j) * beta[j];
    }
    ex[i] = std::exp(dot_product);
  }
  
  // Obtain unique and sorted values from `period` 
  IntegerVector unique_periods = sort_unique(period);  
  int num_periods = unique_periods.size();
  
  
  // Initialize map for cumulative sums (as in `Rfast::group`)
  std::unordered_map<int, double> sum_ex_map;
  std::unordered_map<int, NumericVector> sum_xex_map;
  
  // Initialize `sum_xex_map` with a vector of zeros
  for (int i = 0; i < num_periods; ++i) {
    sum_xex_map[unique_periods[i]] = NumericVector(p, 0.0);
  }
  
  // Calculate sum_ex and sum_xex for each period
  for (int i = 0; i < n; ++i) {
    int current_period = period[i];
    sum_ex_map[current_period] += ex[i];
    for (int j = 0; j < p; ++j) {
      sum_xex_map[current_period][j] += X(i, j) * ex[i];
    }
  }
  
  // Calculate d1i
  NumericMatrix d1i(num_periods, p);
  for (int i = 0; i < num_periods; ++i) {
    int current_period = unique_periods[i];
    double factor = y[i] / sum_ex_map[current_period];
    for (int j = 0; j < p; ++j) {
      d1i(i, j) = sum_xex_map[current_period][j] * factor - sum_xex_map[current_period][j];
    }
  }
  
  return d1i;
}

