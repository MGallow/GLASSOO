// Matt Galloway

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "soft.h"

using namespace Rcpp;





//' @title Lasso Regression (c++)
//' 
//' @description Computes the coefficient estimates for lasso-penalized linear regression.
//' 
//' @details For details on the implementation of 'GLASSO', see the vignette
//' \url{https://mgallow.github.io/GLASSO/}.
//'
//' @param X matrix or data frame
//' @param y matrix or data frame of response values
//' @param ind optional matrix specifying which coefficients will be penalized.
//' @param lam tuning parameter for lasso regularization term. Defaults to 'lam = 0.1'
//' @param crit criterion for convergence. Criterion \code{obj} will loop until the change in the objective after an iteration over the parameter set is less than \code{tol}. Criterion \code{max} will loop until the maximum change in the estimate after an iteration over the parameter set is less than \code{tol}. Defaults to \code{obj}.
//' @param tol tolerance for algorithm convergence. Defaults to 1e-4
//' @param maxit maximum iterations. Defaults to 1e4
//' 
//' @return returns list of returns which includes:
//' \item{Iterations}{number of iterations.}
//' \item{Objective}{value of the objective function.}
//' \item{Coefficients}{estimated regression coefficients.}
//' \item{H}{update H matrix.}
//' 
//' @references
//' \itemize{
//' \item 
//' For more information on the ADMM algorithm, see: \cr
//' Boyd, Stephen, Neal Parikh, Eric Chu, Borja Peleato, Jonathan Eckstein, and others. 2011. 'Distributed Optimization and Statistical Learning via the Alternating Direction Method of Multipliers.' \emph{Foundations and Trends in Machine Learning} 3 (1). Now Publishers, Inc.: 1-122.\cr
//' \url{https://web.stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf}
//' }
//' 
//' @author Matt Galloway \email{gall0441@@umn.edu}
//' 
//' @keywords internal
//'
// [[Rcpp::export]]
List lassoc(const arma::mat &X, const arma::mat &Y, const arma::mat &ind, const double lam = 0.1, std::string crit = "obj", const double tol = 1e-4, const double maxit = 1e4){

  // allocate memory
  int P = X.n_cols, R = Y.n_cols, iter;
  double obj, obj2, temp;
  bool criterion = true;
  arma::mat B, B2, H, H2;
  B = B2 = H = H2 = arma::zeros<arma::mat>(P, R);
  iter = obj = obj2 = temp = 0;
  
  // save values
  arma::mat XX = arma::trans(X)*X;
  arma::mat XY = arma::trans(X)*Y;
  
  
  // loop until convergence
  while (criterion && (iter <= maxit)){

    // keep old values
    B = B2;
    obj = obj2;
    
    // loop over all entries of beta
    for (int r = 0; r < R; r++){
      for (int p = 0; p < P; p++){
        
        // update betas
        B2(p, r) = softc(XY(p, r) - H(p, r), lam*ind(p, r))/XX(p, p);
        temp = 0;
        
        // update H, if necessary
        if (B2(p, r) != B(p, r)){
          
          // update each element in column r, except p
          for (int p_ = 0; p_ < P; p_++){
            if (p_ != p){
              
              // update all rows except row p
              H2(p_, r) = H(p_, r) - XX(p_, p)*(B(p, r) - B2(p, r));
              
              // create temporary sum used in objective, if necessary
              if (crit == "obj"){
                temp += B(p_, r)*(H2(p_, r) - H(p_, r));
              }
            }
          }
          
          // update objective function, if necessary
          if (crit == "obj"){
            obj2 = obj + (B(p, r) - B2(p, r))*(XY(p, r) - H(p, r)/2) + (std::pow(B2(p, r), 2) - std::pow(B(p, r), 2))*XX(p, p)/2 + temp/2 + lam*(std::abs(B2(p, r)) - std::abs(B(p, r)));
          }
        }
      }
    }
    
    // stopping criterion
    iter++;
    if (crit == "obj"){
      
      // compute objective improvement
      criterion = (std::abs(obj2 - obj) > tol);
      
    } else {
      
      // compute estimate change
      criterion = (arma::abs(B2 - B).max() > tol);
      
    }

    // R_CheckUserInterrupt
    if (iter % 1000 == 0){
      R_CheckUserInterrupt();
    }
  }
  
  // compute final objective function value, if necessary
  if (crit != "obj"){
    obj2 = std::pow(arma::norm(Y - X*B2, "fro"), 2)/2 + lam*arma::accu(arma::abs(ind % B2));
  }

  return List::create(Named("Iterations") = iter,
                      Named("Objective") = obj2,
                      Named("Coefficients") = B2,
                      Named("H") = H2);

}

