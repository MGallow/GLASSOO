#ifndef GLASSO_H
#define GLASSO_H

#include <RcppArmadillo.h>
#include <Rcpp.h>

Rcpp::List GLASSOc(const arma::mat &S, const arma::mat &initSigma, const arma::mat &initOmega, const double lam, std::string crit_out = "avg", std::string crit_in = "loss", const double tol_out = 1e-4, const double tol_in = 1e-4, const int maxit_out = 1e4, const int maxit_in = 1e4);

#endif
