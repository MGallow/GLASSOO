#ifndef MISC_H
#define MISC_H

#include <RcppArmadillo.h>
#include <Rcpp.h>


void reducec(const arma::mat &X, arma::mat &Y, const int &p);

void extractc(const arma::mat &X, arma::mat &Y, const int &p);

void updatec(arma::mat &X, const arma::mat &Y, const int &p);

int numzeros(arma::mat &X);

#endif
