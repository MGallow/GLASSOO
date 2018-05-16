// Matt Galloway

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "lasso.h"

using namespace Rcpp;




//' @title Penalized precision matrix estimation (c++)
//' 
//' @description Penalized precision matrix estimation using the graphical lasso (glasso) algorithm
//' 
//' @details For details on the implementation of 'GLASSOO', see the vignette
//' \url{https://mgallow.github.io/GLASSOO/}.
//'
//' @param S pxp sample covariance matrix (denominator n).
//' @param initSigma initialization matrix for estimated covariance matrix Sigma
//' @param initOmega initialization matrix for estimated precision matrix Omega
//' @param lam tuning parameter for lasso penalty.
//' @param diagonal option to penalize the diagonal elements of the estimated precision matrix (\eqn{\Omega}). Defaults to \code{FALSE}.
//' @param crit1 criterion for convergence in outer (blockwise) loop. Criterion \code{loglik} will loop until the change in the penalized likelihood after an iteration over the parameter set is less than \code{tol1}. Criterion \code{max} will loop until the maximum change in the estimated Sigma after an iteration over the parameter set is less than \code{tol1}. Defaults to \code{max}.
//' @param crit2 criterion for convergence in inner (lasso) loop. Criterion \code{loss} will loop until the change in the lasso objective after an iteration over the parameter set is less than \code{tol2}. Criterion \code{max} will loop until the maximum change in the lasso estimate after an iteration over the parameter set is less than \code{tol2}. Defaults to \code{loss}.
//' @param tol1 convergence tolerance for outer (blockwise) loop. Defaults to 1e-4.
//' @param tol2 convergence tolerance for inner (lasso) loop. Defaults to 1e-4.
//' @param maxit1 maximum number of iterations for outer (blockwise) loop. Defaults to 1e4.
//' @param maxit2 maximum number of iterations for inner (lasso) loop. Defaults to 1e4.
//' 
//' @return returns list of returns which includes:
//' \item{Iterations}{number of iterations.}
//' \item{lam}{optimal tuning parameters.}
//' \item{Omega}{estimated penalized precision matrix.}
//' \item{Sigma}{estimated covariance matrix.}
//' 
//' @references
//' \itemize{
//' \item 
//' For more information on the graphical lasso algorithm, see: \cr
//' Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. "Sparse inverse covariance estimation with the graphical lasso." \emph{Biostatistics} 9.3 (2008): 432-441.\cr
//' \url{http://statweb.stanford.edu/~tibs/ftp/glasso-bio.pdf}
//' }
//' 
//' @author Matt Galloway \email{gall0441@@umn.edu}
//' 
//' @keywords internal
//'
// [[Rcpp::export]]
List GLASSOc(const arma::mat &S, const arma::mat &initSigma, const arma::mat &initOmega, const double lam, bool diagonal = false, std::string crit1 = "max", std::string crit2 = "loss", const double tol1 = 1e-4, const double tol2 = 1e-4, const int maxit1 = 1e4, const int maxit2 = 1e4){

  // allocate memory
  bool criterion = true;
  int P = S.n_cols, iter;
  double lik, lik2, sgn, logdet;
  iter = lik = lik2 = sgn = logdet = 0;
  arma::mat Omega, Sigma, Sigma2, C, B;
  C = arma::ones<arma::mat>(P, P);
  B = arma::zeros<arma::mat>(P - 1, 1);
  Sigma2 = initSigma;
  
  // option to penalize diagonal elements
  if (!diagonal){
    C -= arma::eye<arma::mat>(P, P);
  }

  // loop until convergence
  while (criterion && (iter < maxit1)){
    
    // update values
    iter++;
    Sigma = Sigma2;

    // blockwise coordinate descent
    for (int p = 0; p < P; p++){
      
      
      
    }
    

    // stopping criterion
    if (crit1 == "loglik"){

      // compute penalized likelihood improvement (close enough)
      arma::log_det(logdet, sgn, Omega);
      lik2 = (-P/2)*(arma::accu(Omega % S) - logdet + lam*arma::accu(C % arma::abs(Omega)));
      criterion = (std::abs(lik2 - lik) >= tol1);
      lik = lik2;

    } else {

      // compute estimate change
      criterion = (arma::abs(Sigma2 - Sigma).max() > tol1);

    }

    // R_CheckUserInterrupt
    if (iter % 1000 == 0){
      R_CheckUserInterrupt();
    }
  }

  return List::create(Named("Iterations") = iter,
                      Named("lam") = lam,
                      Named("Omega") = Omega,
                      Named("Sigma") = Sigma2);

}

