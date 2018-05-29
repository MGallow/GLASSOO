// Matt Galloway

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <progress.hpp>
#include "glasso.h"
#include "misc.h"

using namespace Rcpp;




//' @title K fold (c++)
//' @description creates vector of shuffled indices.
//' @param n number of elements.
//' @param K number of folds.
//' @keywords internal
//'
arma::vec kfold(const int &n, const int &K){
  
  // create sequence 1:n
  arma::vec indices = arma::linspace<arma::vec>(1, n, n);
  
  // assign number fold
  for (int i = 0; i < n; i ++){
    indices[i] = i % K;
  }
  
  // shuffle indices
  indices = arma::shuffle(indices);
  
  return indices;
  
}



//--------------------------------------------------------------------------------------------




//' @title CV penalized precision matrix estimation (c++)
//' @description Cross validation function for GLASSO.
//'
//' @param X option to provide a nxp matrix. Each row corresponds to a single observation and each column contains n observations of a single feature/variable.
//' @param S option to provide a pxp sample covariance matrix (denominator n). If argument is \code{NULL} and \code{X} is provided instead then \code{S} will be computed automatically.
//' @param lam positive tuning parameters for elastic net penalty. If a vector of parameters is provided, they should be in increasing order.
//' @param diagonal option to penalize the diagonal elements of the estimated precision matrix (\eqn{\Omega}). Defaults to \code{FALSE}.
//' @param path option to return the regularization path. This option should be used with extreme care if the dimension is large. If set to TRUE, cores will be set to 1 and errors and optimal tuning parameters will based on the full sample. Defaults to FALSE.
//' @param crit_out criterion for convergence in outer (blockwise) loop. Criterion \code{avg} will loop until the average absolute parameter change is less than \code{tol_out} times tolerance multiple. Criterion \code{max} will loop until the maximum change in the estimated Sigma after an iteration over the parameter set is less than \code{tol_out}. Defaults to \code{avg}.
//' @param crit_in criterion for convergence in inner (lasso) loop. Criterion for convergence. Criterion \code{loss} will loop until the relative change in the objective for each response after an iteration is less than \code{tol_in}. Criterion \code{avg} will loop until the average absolute change for each response is less than \code{tol_in} times tolerance multiple. Similary, criterion \code{max} will loop until the maximum absolute change is less than \code{tol_in} times tolerance multiple. Defaults to \code{loss}.
//' @param tol_out convergence tolerance for outer (blockwise) loop. Defaults to 1e-4.
//' @param tol_in convergence tolerance for inner (lasso) loop. Defaults to 1e-4.
//' @param maxit_out maximum number of iterations for outer (blockwise) loop. Defaults to 1e4.
//' @param maxit_in maximum number of iterations for inner (lasso) loop. Defaults to 1e4.
//' @param adjmaxit_out adjusted maximum number of iterations. During cross validation this option allows the user to adjust the maximum number of iterations after the first \code{lam} tuning parameter has converged (for each \code{alpha}). This option is intended to be paired with \code{warm} starts and allows for "one-step" estimators. Defaults to 1e4.
//' @param K specify the number of folds for cross validation.
//' @param crit_cv cross validation criterion (\code{loglik}, \code{AIC}, or \code{BIC}). Defaults to \code{loglik}.
//' @param start specify \code{warm} or \code{cold} start for cross validation. Default is \code{warm}.
//' @param trace option to display progress of CV. Choose one of \code{progress} to print a progress bar, \code{print} to print completed tuning parameters, or \code{none}.
//' 
//' @return list of returns includes:
//' \item{lam}{optimal tuning parameter.}
//' \item{path}{array containing the solution path. Solutions will be ordered by ascending lambda values.}
//' \item{min.error}{minimum average cross validation error (cv.crit) for optimal parameters.}
//' \item{avg.error}{average cross validation error (cv.crit) across all folds.}
//' \item{cv.error}{cross validation errors (cv.crit).}
//' @keywords internal
//'
// [[Rcpp::export]]
List CV_GLASSOc(const arma::mat &X, const arma::mat &S, const arma::colvec &lam, bool diagonal = false, bool path = false, std::string crit_out = "avg", std::string crit_in = "loss", const double tol_out = 1e-4, const double tol_in = 1e-4, int maxit_out = 1e4, int maxit_in = 1e4, int adjmaxit_out = 1e4, int K = 5, std::string crit_cv = "loglik", std::string start = "warm", std::string trace = "progress") {
  
  // initialization
  int n, p = S.n_cols, l = lam.n_rows, initmaxit = maxit_out;
  double sgn = 0, logdet = 0, lam_, alpha;
  arma::mat X_train, X_test, S_train(S), S_test(S), initSigma(S), initOmega;
  arma::mat Omega, identity, CV_errors(l, K, arma::fill::zeros);
  initOmega = identity = arma::eye<arma::mat>(p, p);
  arma::colvec CV_error, zerosl(l, arma::fill::zeros); arma::rowvec X_bar;
  arma::uvec index, index_; arma::vec folds; arma::cube Path;
  Progress progress(l*K, trace == "progress");
  
  
  // no need to create folds if K = 1
  if (K == 1){
    
    // set sample size
    n = S.n_rows;
    
    // initialize Path, if necessary
    if (path){
      Path = arma::zeros<arma::cube>(p, p, l);
    }
    
  } else {
    
    // designate folds and shuffle -- ensures randomized folds
    n = X.n_rows;
    folds = kfold(n, K);
    
  }
  
  // parse data into folds and perform CV
  for (int k = 0; k < K; k++){
    if (K > 1) {
      
      // separate into training and testing data
      index = arma::find(folds != k);
      index_ = arma::find(folds == k);
      
      // training set
      X_train = X.rows(index);
      X_bar = arma::mean(X_train, 0);
      X_train -= arma::ones<arma::colvec>(X_train.n_rows)*X_bar;
      
      // validation set
      X_test = X.rows(index_);
      X_test -= arma::ones<arma::colvec>(X_test.n_rows)*X_bar;
      n = X_test.n_rows;
      
      // sample covariances
      S_train = arma::cov(X_train, 1);
      S_test = arma::cov(X_test, 1);
      
    }
    
    // re-initialize values for each fold
    CV_error = zerosl; maxit_out = initmaxit;
    initSigma = S_train; initOmega = identity;
    
    // initial sigma, if not diagonal
    if (!diagonal){
      
      // provide estimate that is pd and dual feasible
      initSigma -=arma::diagmat(initSigma); initSigma = arma::abs(initSigma);
      alpha = lam[0]/initSigma.max();
      if (alpha >= 1){
        alpha = 1;
      }
      initSigma = (1 - alpha)*S_train;
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
      
      // compute the ridge-penalized likelihood precision matrix estimator at the ith value in lam:
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
      CV_error[i] = (n/2)*(arma::accu(Omega % S_test) - logdet);
      
      // update for crit_cv, if necessary
      if (crit_cv == "AIC"){
        CV_error[i] += numzeros(Omega);
      }
      if (crit_cv == "BIC"){
        CV_error[i] += numzeros(Omega)*log(n)/2;
      }
      
      // save estimate if path = TRUE
      if (path){
        Path.slice(i) = Omega;
      }
      
      // update progress bar
      if (trace == "progress"){
        progress.increment();
        
        // if not quiet, then print progress lambda
      } else if (trace == "print"){
        Rcout << "Finished lam = " << lam[i] << " in fold " << k << "\n";
      }
    }
    
    // if not quiet, then print progress fold
    if (trace == "print"){
      Rcout << "Finished fold" << k << "\n";
    }
    
    // append CV errors
    CV_errors.col(k) = CV_error;
    
  }
  
  
  // determine optimal tuning parameters
  arma::colvec AVG_error = arma::mean(CV_errors, 1);
  double error = AVG_error.min();
  arma::uword ind = AVG_error.index_min();
  int lam_ind = ind % AVG_error.n_rows;
  double best_lam = lam[lam_ind];
  
  
  // return list of coefficients
  return List::create(Named("lam") = best_lam,
                      Named("path") = Path,
                      Named("min.error") = error,
                      Named("avg.error") = AVG_error,
                      Named("cv.error") = CV_errors);
}



