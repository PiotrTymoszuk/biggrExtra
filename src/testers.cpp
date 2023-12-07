/*** Calculation of p values from simulation or bootstrap results*/

#include <Rcpp.h>
#include "ci_functions.h"

using namespace Rcpp;

// [[Rcpp::export]]

NumericVector simPval(NumericMatrix sim,
                      double mu = 1,
                      bool median_estimate = false) {

  // sim is a matrix with features in columns and iterations in rows
  // cells contain simulated values
  // the function counts are occurrences

  int n_iters = sim.nrow();
  int n_features = sim.ncol();

  double est;

  NumericVector fails;
  double n_fails;

  NumericVector p_vals(n_features);

  for(int i = 0; i < n_features; ++i) {

    if(median_estimate) {

      est = Median(sim(_, i));

    } else {

      est = mean(sim(_, i));

    }

    if(est == mu) {

      n_fails = n_iters;

    } else {

      if(est > mu) {

        fails = ifelse(sim(_, i) <= mu, 1.0, 0.0);

      } else {

        fails = ifelse(sim(_, i) >= mu, 1.0, 0.0);

      }

      n_fails = sum(fails);

    }

    if(n_fails == 0) n_fails = 1;

    n_fails = n_fails/n_iters;

    p_vals[i] = n_fails;

  }

  return p_vals;

}


