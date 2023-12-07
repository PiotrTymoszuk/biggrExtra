/*** Cpp functions that sample Random multinomial distibutions
*/

#include <Rcpp.h>
#include "samplers.h"

using namespace Rcpp;

//[[Rcpp::export]]

IntegerVector concatenate(IntegerVector x, IntegerVector y) {

  int n_cmm = x.size() + y.size();
  IntegerVector z(n_cmm);
  CharacterVector names_x = x.names();
  CharacterVector names_y = y.names();
  CharacterVector names_cmm(n_cmm);

  // copying the body

  std::copy(x.begin(), x.end(), z.begin());
  std::copy(y.begin(), y.end(), z.begin() + x.size());

  // copying the names

  std::copy(names_x.begin(), names_x.end(), names_cmm.begin());
  std::copy(names_y.begin(), names_y.end(), names_cmm.begin() + names_x.size());

  z.attr("names") = names_cmm;

  return z;
}

// [[Rcpp::export]]

NumericMatrix drawNorm(NumericVector mu, NumericVector sd, double n) {

  int n_vars = mu.size();
  NumericMatrix out(n, n_vars);

  for(int i = 0; i < n_vars; ++i) {

    NumericMatrix::Column out_col = out(_, i);

    out_col = rnorm(n, mu[i], sd[i]);

  }

  return out;

}

// [[Rcpp::export]]

NumericMatrix vecDraws(CharacterVector x, int size, int n_iter) {

  // returns a numeric matrix with counts of elements in random draws from
  // a character vector with a specified size

  CharacterVector unique_elmts = unique(x);

  CharacterVector res(size);

  IntegerVector freq_vec;
  int freq;

  NumericMatrix out(unique_elmts.length(), n_iter);

  out.fill(0);

  for(int i = 0; i < n_iter; ++i) {

    res = sample(x, size);

    for(int j = 0; j < unique_elmts.length(); ++j) {

      freq = std::count(res.begin(), res.end(), unique_elmts[j]);

      out(j, i) = freq;

    }

  }

  rownames(out) = unique_elmts;

  return out;

}

// [[Rcpp::export]]

NumericMatrix compareCounts(NumericMatrix draws, NumericVector n_observed) {

  // compares the observed number of occurrences of an element
  // with the expected number obtained e.g. by random draws

  int n_features = n_observed.size();
  int n_draws = draws.ncol();

  NumericVector row;
  NumericVector fails;

  int res;

  NumericVector out(n_features);
  NumericVector n_expected(n_features);
  NumericMatrix final_mtx(n_features, 5);

  for(int i = 0; i < n_features; ++i) {

    row = draws(i, _);
    fails = ifelse(row >= n_observed[i], 1, 0);

    res = sum(fails);

    if(res == 0) res = 1;

    out[i] = res;
    n_expected[i] = mean(row);

  }

  final_mtx(_, 0) = n_observed; //observed frequency
  final_mtx(_, 1) = n_expected; //mean expected frequency by chance
  final_mtx(_, 2) = n_observed/n_expected; //Enrichment odds ratio
  final_mtx(_, 3) = out; //numbers of draws with higher expected than observed freq.
  final_mtx(_, 4) = out/n_draws; // raw p value

  rownames(final_mtx) = rownames(draws);

  return final_mtx;

}
