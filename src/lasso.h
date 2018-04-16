#ifndef LASSO_H
#define LASSO_H

#include <RcppArmadillo.h>
#include <Rcpp.h>

Rcpp::List lassoc(const arma::mat &X, const arma::mat &Y, const arma::mat &ind, const double lam = 0.1, std::string crit = "loss", const double tol = 1e-4, const double maxit = 1e4);

#endif
