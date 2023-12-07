/*** Calculation of CI and expected values of MC simulations */

#include <Rcpp.h>
#include <Rmath.h>
#include "ci_functions.h"

using namespace Rcpp;

// [[Rcpp::export]]

double Median(NumericVector x) {

  int x_size = x.size();
  double med;

  x.sort();

  if(x_size % 2 != 0) {

    med = x[x_size/2];

  } else {

    med = (x[x_size/2 - 1] + x[x_size/2])/2;

  }

  return med;

}

// [[Rcpp::export]]

NumericVector Quantile(NumericVector x, NumericVector probs) {

  // calculation of quantiles
  // many thanks to https://github.com/RcppCore/Rcpp/issues/967

  const size_t n = x.size(), np = probs.size();

  if (n==0) return x;
  if (np==0) return probs;

  NumericVector index = (n - 1.) * probs, y = x.sort(), x_hi(np), qs(np);
  NumericVector lo = floor(index), hi = ceiling(index);

  for (size_t i = 0; i < np; ++i) {

    qs[i] = y[lo[i]];
    x_hi[i] = y[hi[i]];

    if ((index[i] > lo[i]) && (x_hi[i] != qs[i])) {

      double h;
      h = index[i] - lo[i];
      qs[i] = (1.- h) * qs[i] + h * x_hi[i];

    }

  }

  return qs;

}

// [[Rcpp::export]]

NumericMatrix estMC(NumericMatrix mc,
                    Function ci_fun,
                    double conf_level = 0.95,
                    double mu = 1) {

  // computes MC estimates: mean and median
  // along with SD and 95% CI with the function defined by the `ci_fun`
  // argument.
  // The function performs also significance testing by counting the number
  // of simulations returning results above or below the value mu

  int n_reacts = mc.ncol();

  // out will store stats: reactions as rows; columns are median, mean, SD,
  // lower and upper CI

  NumericMatrix out(n_reacts, 5);

  NumericVector ci_result;

  for(int i = 0; i < n_reacts; ++i) {

    out(i, 0) = Median(mc(_, i));
    out(i, 1) = mean(mc(_, i));
    out(i, 2) = sd(mc(_, i));

    ci_result = ci_fun(mc(_,i), conf_level);

    out(i, 3) = ci_result[0];
    out(i, 4) = ci_result[1];

  }

  return out;

}

// [[Rcpp::export]]

NumericVector perci(NumericVector theta, double conf_level = 0.95) {

  // computes percentile confidence intervals

  NumericVector ci_probs{(1 - conf_level)/2, (1 + conf_level)/2};

  return Quantile(theta, ci_probs);

}

// [[Rcpp::export]]

NumericVector bca(NumericVector theta, double conf_level = 0.95) {

  // computes BCA confidence intervals based on the R code
  // provided by the coxed package
  // https://rdrr.io/cran/coxed/src/R/bca.R

  double low;
  double high;

  low = (1 - conf_level)/2;
  high = 1 - low;

  int sims = theta.size();

  NumericVector low_theta;

  low_theta = ifelse(theta < mean(theta), 1.0, 0.0);

  double z_inv;

  z_inv = sum(low_theta)/sims;

  double z;

  z = R::qnorm(z_inv, 0, 1, 1, 0);

  NumericVector U;

  U = (sims - 1) * (mean(theta) - theta);

  double top;
  double under;
  double a;

  top = sum(pow(U, 3));

  under = sum(pow(U, 2));

  under = 6 * std::pow(under, 1.5);

  a = top/under;

  double lower_inv;
  double upper_inv;

  double q_low;
  double q_high;

  q_low = R::qnorm(low, 0, 1, 1, 0);

  lower_inv = z + (z + q_low)/(1 - a * (z + q_low));

  lower_inv = R::pnorm(lower_inv, 0, 1, 1, 0);

  q_high = R::qnorm(high, 0, 1, 1, 0);

  upper_inv = z + (z + q_high)/(1 - a * (z + q_low));

  upper_inv = R::pnorm(upper_inv, 0, 1, 1, 0);

  NumericVector probs = {lower_inv, upper_inv};

  return Quantile(theta, probs);

}
