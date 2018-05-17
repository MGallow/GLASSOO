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
//' @param lam tuning parameter for lasso penalty.
//' @param crit_out criterion for convergence in outer (blockwise) loop. Criterion \code{avg} will loop until the average absolute parameter change is less than \code{tol_out} times tolerance multiple. Criterion \code{max} will loop until the maximum change in the estimated Sigma after an iteration over the parameter set is less than \code{tol_out}. Defaults to \code{avg}.
//' @param crit_in criterion for convergence in inner (lasso) loop. Criterion for convergence. Criterion \code{loss} will loop until the change in the objective for each response after an iteration is less than \code{tol_in}. Criterion \code{avg} will loop until the average absolute change for each response is less than \code{tol_in} times tolerance multiple. Similary, criterion \code{max} will loop until the maximum absolute change is less than \code{tol_in} times tolerance multiple. Defaults to \code{loss}.
//' @param tol_out convergence tolerance for outer (blockwise) loop. Defaults to 1e-4.
//' @param tol_in convergence tolerance for inner (lasso) loop. Defaults to 1e-4.
//' @param maxit_out maximum number of iterations for outer (blockwise) loop. Defaults to 1e4.
//' @param maxit_in maximum number of iterations for inner (lasso) loop. Defaults to 1e4.
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
List GLASSOc(const arma::mat &S, const arma::mat &initSigma, const double lam, std::string crit_out = "avg", std::string crit_in = "loss", const double tol_out = 1e-4, const double tol_in = 1e-4, const int maxit_out = 1e4, const int maxit_in = 1e4){
  
  // allocate memory
  bool criterion = true;
  int P = S.n_cols, iter = 0;
  double mult;
  arma::mat Sigmatemp(P - 1, P - 1, arma::fill::zeros), Stemp, Omegatemp(P, 1, arma::fill::zeros), ind(P - 1, 1, arma::fill::ones), Beta, Betas, H, Omega(P, P, arma::fill::eye), Sigma(initSigma), Sigma2(initSigma), Sminus(S);
  Beta = Stemp = arma::zeros<arma::mat>(P - 1, 1);
  Betas = H = arma::zeros<arma::mat>(P - 1, P);
  Sminus -= arma::diagmat(Sminus);
  mult = arma::accu(arma::abs(Sminus));
  
  
  // loop until convergence
  while (criterion && (iter < maxit_out)){
    
    // update values
    iter++;
    Sigma = Sigma2;
    
    // blockwise coordinate descent
    for (int p = 0; p < P; p++){
      
      // update Beta with pth column of Betas
      Beta = Betas.col(p);
      
      // set Sigmatemp = Sigma[-p, -p]
      for (int i = 0; i < (P - 1); i++){
        for (int j = 0; j < (P - 1); j++){
          if ((i < p) && (j < p)){
            Sigmatemp(i, j) = Sigma2(i, j);
          } else if ((i < p) && (j >= p)){
            Sigmatemp(i, j) = Sigma2(i, j + 1);
          } else if ((i >= p) && (j < p)){
            Sigmatemp(i, j) = Sigma2(i + 1, j);
          } else {
            Sigmatemp(i, j) = Sigma2(i + 1, j + 1);
          }
        }
      }
      
      // initialize H matrix used in lassoc
      H = Sigmatemp*Beta - Sigmatemp.diag() % Beta;
      
      // set Stemp = S[-p, p]
      for (int i = 0; i < (P - 1); i++){
        if (i < p){
          Stemp(i, 0) = S(i, p);
        } else {
          Stemp(i, 0) = S(i + 1, p);
        }
      }
      
      // execute LASSO
      List LASSO = lassoc(Sigmatemp, Stemp, Beta, H, ind, lam, crit_in, tol_in, maxit_in);
      Betas.col(p) = as<arma::mat>(LASSO["Coefficients"]);
      
      // update Stemp = Sigma12
      Stemp = Sigmatemp*Betas.col(p);
      
      // update Sigma[-p, p] = Sigma[p, -p] = Stemp
      for (int i = 0; i < (P - 1); i++){
        if (i < p){
          Sigma2(i, p) = Stemp(i, 0);
          Sigma2(p, i) = Stemp(i, 0);
        } else {
          Sigma2(i + 1, p) = Stemp(i, 0);
          Sigma2(p, i + 1) = Stemp(i, 0);
        }
      }
      
    }
    
    
    // stopping criterion
    if (crit_out == "avg") {
      
      // compute estimate avg change
      criterion = (arma::accu(arma::abs(Sigma2 - Sigma)) > tol_out*mult);
        
    } else {
      
      // compute estimate max change
      criterion = (arma::abs(Sigma2 - Sigma).max() > tol_out*mult/std::pow(P, 2));
      
    }
    
    // R_CheckUserInterrupt
    if (iter % 1000 == 0){
      R_CheckUserInterrupt();
    }
  }
  
  
  // compute Omega from Sigma
  for (int p = 0; p < P; p++){
    
    // update Stemp = Sigma[-p, p]
    for (int i = 0; i < (P - 1); i++){
      if (i < p){
        Stemp(i, 0) = Sigma2(i, p);
      } else {
        Stemp(i, 0) = Sigma2(i + 1, p);
      }
    }
    
    // update Omega[p, p] and Omegatemp = Omega12
    Omega(p, p) = 1/(Sigma2(p, p) - arma::accu(Stemp % Betas.col(p)));
    Omegatemp = -Omega(p, p)*Betas.col(p);
    
    // set Omega[-p, p] = Omega[p, -p] = Omegatemp
    for (int i = 0; i < (P - 1); i++){
      if (i < p){
        Omega(i, p) = Omegatemp(i, 0);
        Omega(p, i) = Omegatemp(i, 0);
      } else {
        Omega(i + 1, p) = Omegatemp(i, 0);
        Omega(p, i + 1) = Omegatemp(i, 0);
      }
    }
    
  }
  
  return List::create(Named("Iterations") = iter,
                      Named("lam") = lam,
                      Named("Omega") = Omega,
                      Named("Sigma") = Sigma2);
  
}

