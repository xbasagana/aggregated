#include <Rcpp.h>
#include <unordered_map>
#include <cmath>  // For log and exp

using namespace Rcpp;

// [[Rcpp::export]]
double agr_logLik_cpp(NumericVector beta, NumericMatrix X, NumericVector y, IntegerVector period) {
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
  
  // Initialize a map for cumulative sums (as in `Rfast::group`)
  std::unordered_map<int, double> sum_ex_map;
  
  // Calculate sum_ex for each period
  for (int i = 0; i < n; ++i) {
    int current_period = period[i];
    sum_ex_map[current_period] += ex[i];
  }
  
  // Calculate log-likelihood
  double logLik = 0.0;
  for (int i = 0; i < num_periods; ++i) {
    int current_period = unique_periods[i];
    logLik += y[i] * std::log(sum_ex_map[current_period]) - sum_ex_map[current_period] - lgamma(y[i] + 1);
  }
  
  return logLik;
}
