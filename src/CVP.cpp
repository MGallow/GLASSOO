// Matt Galloway

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <progress.hpp>
#include "glasso.h"
#include "misc.h"

using namespace Rcpp;




//' @title CV (no folds) penalized precision matrix estimation (c++)
//' @description Cross validation (no folds) function for GLASSO. This function is to be used with CVP_GLASSO.
//'
//' @param n sample size for X_valid (used to calculate crit_cv)
//' @param S_train pxp sample covariance matrix for training data (denominator n).
//' @param S_valid pxp sample covariance matrix for validation data (denominator n).
//' @param lam positive tuning parameters for elastic net penalty. If a vector of parameters is provided, they should be in increasing order.
//' @param diagonal option to penalize the diagonal elements of the estimated precision matrix (\eqn{\Omega}). Defaults to \code{FALSE}.
//' @param crit_out criterion for convergence in outer (blockwise) loop. Criterion \code{avg} will loop until the average absolute parameter change is less than \code{tol_out} times tolerance multiple. Criterion \code{max} will loop until the maximum change in the estimated Sigma after an iteration over the parameter set is less than \code{tol_out}. Defaults to \code{avg}.
//' @param crit_in criterion for convergence in inner (lasso) loop. Criterion for convergence. Criterion \code{loss} will loop until the relative change in the objective for each response after an iteration is less than \code{tol_in}. Criterion \code{avg} will loop until the average absolute change for each response is less than \code{tol_in} times tolerance multiple. Similary, criterion \code{max} will loop until the maximum absolute change is less than \code{tol_in} times tolerance multiple. Defaults to \code{loss}.
//' @param tol_out convergence tolerance for outer (blockwise) loop. Defaults to 1e-4.
//' @param tol_in convergence tolerance for inner (lasso) loop. Defaults to 1e-4.
//' @param maxit_out maximum number of iterations for outer (blockwise) loop. Defaults to 1e4.
//' @param maxit_in maximum number of iterations for inner (lasso) loop. Defaults to 1e4.
//' @param adjmaxit_out adjusted maximum number of iterations. During cross validation this option allows the user to adjust the maximum number of iterations after the first \code{lam} tuning parameter has converged (for each \code{alpha}). This option is intended to be paired with \code{warm} starts and allows for "one-step" estimators. Defaults to 1e4.
//' @param crit_cv cross validation criterion (\code{loglik}, \code{AIC}, or \code{BIC}). Defaults to \code{loglik}.
//' @param start specify \code{warm} or \code{cold} start for cross validation. Default is \code{warm}.
//' @param trace option to display progress of CV. Choose one of \code{progress} to print a progress bar, \code{print} to print completed tuning parameters, or \code{none}.
//' 
//' @return cross validation errors (crit_cv)
//' 
//' @keywords internal
//'
// [[Rcpp::export]]
arma::mat CVP_GLASSOc(const int n, const arma::mat &S_train, const arma::mat &S_valid, const arma::colvec &lam, bool diagonal = false, std::string crit_out = "avg", std::string crit_in = "loss", const double tol_out = 1e-4, const double tol_in = 1e-4, int maxit_out = 1e4, int maxit_in = 1e4, int adjmaxit_out = 1e4, std::string crit_cv = "loglik", std::string start = "warm", std::string trace = "progress") {
  
  // initialization
  int p = S_train.n_rows, l = lam.n_rows;
  double sgn = 0, logdet = 0, lam_, alpha;
  arma::mat Omega, initOmega(p, p, arma::fill::eye), Sigma, initSigma(S_train), CV_error(l, 1, arma::fill::zeros);
  Progress progress(l, trace == "progress");
  
  
  // initial sigma
  if (!diagonal){
    
    // provide estimate that is pd and dual feasible
    initSigma -=arma::diagmat(initSigma); initSigma = arma::abs(initSigma);
    alpha = lam[0]/initSigma.max();
    if (alpha >= 1){
      alpha = 1;
    }
    initSigma = (1 - alpha)*initSigma;
    initSigma.diag() = S_train.diag();
    
  }
  
  // loop over all tuning parameters
  for (int i = 0; i < l; i++){
    
    // set temporary tuning parameters
    lam_ = lam[i];
    
    // update diagonal elements of initSigma, if necessary
    if (diagonal){
      initSigma.diag() = S_train.diag() + lam_;
    }
    
    // compute the penalized likelihood precision matrix estimator at the ith value in lam:
    List GLASSO = GLASSOc(S_train, initSigma, initOmega, lam_, crit_out, crit_in, tol_out, tol_in, maxit_out, maxit_in);
    Omega = as<arma::mat>(GLASSO["Omega"]);
    
    if (start == "warm"){
      
      // option to save initial values for warm starts
      initSigma = as<arma::mat>(GLASSO["Sigma"]);
      initOmega = as<arma::mat>(GLASSO["Omega"]);
      maxit_out = adjmaxit_out;
      
    }
    
    // compute the observed negative validation loglikelihood (close enough)
    arma::log_det(logdet, sgn, Omega);
    CV_error[i] = (n/2)*(arma::accu(Omega % S_valid) - logdet);
    
    // update for crit_cv, if necessary
    if (crit_cv == "AIC"){
      CV_error[i] += numzeros(Omega);
    }
    if (crit_cv == "BIC"){
      CV_error[i] += numzeros(Omega)*log(n)/2;
    }
    
    // update progress bar
    if (trace == "progress"){
      progress.increment();
      
      // if not quiet, then print progress lambda
    } else if (trace == "print"){
      Rcout << "Finished lam = " << lam[i] << "\n";
    }
  }
  
  // return CV errors
  return(CV_error);
}





