#include <Rcpp.h>
#include <algorithm>  // for sort

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix d2func_aggr_cpp(NumericVector beta, NumericMatrix X, NumericVector y, IntegerVector period) {
  int n = X.nrow();
  int p = X.ncol();
  
  // Calculate `ex = exp(X %*% beta)`
  NumericVector ex(n);
  for (int i = 0; i < n; ++i) {
    double dot_product = 0.0;
    for (int j = 0; j < p; ++j) {
      dot_product += X(i, j) * beta[j];
    }
    ex[i] = std::exp(dot_product);
  }
  
  // Obtein unique and sorted values from  `period`
  IntegerVector periods = sort_unique(period);
  int num_periods = periods.size();
  
  // Initialize sum_ex and sum_xex
  NumericVector sum_ex(num_periods);
  NumericMatrix sum_xex(num_periods, p);
  NumericMatrix aa(n, p);
  
  // Calculate sum_ex and sum_xex
  for (int i = 0; i < n; ++i) {
    int current_period = period[i];
    int period_index = std::distance(periods.begin(), std::find(periods.begin(), periods.end(), current_period));
    sum_ex[period_index] += ex[i];
    for (int j = 0; j < p; ++j) {
      aa(i, j) = X(i, j) * ex[i];
      sum_xex(period_index, j) += aa(i, j);
    }
  }
  
  // Initialize hessi
  NumericMatrix hessi(p, p);
  
  // Calculate hessi
  for (int i = 0; i < num_periods; ++i) {
    NumericMatrix sum_xxexi(p, p);
    int count = 0;
    for (int j = 0; j < n; ++j) {
      if (period[j] == periods[i]) {
        NumericMatrix xxex(p, p);
        for (int k = 0; k < p; ++k) {
          for (int l = 0; l < p; ++l) {
            xxex(k, l) = X(j, k) * X(j, l) * ex[j];
          }
        }
        for (int k = 0; k < p; ++k) {
          for (int l = 0; l < p; ++l) {
            sum_xxexi(k, l) += xxex(k, l);
          }
        }
        count++;
      }
    }
    
    double factor = y[i] / (sum_ex[i] * sum_ex[i]);
    for (int k = 0; k < p; ++k) {
      for (int l = 0; l < p; ++l) {
        hessi(k, l) += (sum_xxexi(k, l) * sum_ex[i] - sum_xex(i, k) * sum_xex(i, l)) * factor - sum_xxexi(k, l);
      }
    }
  }
  
  return hessi;
}
