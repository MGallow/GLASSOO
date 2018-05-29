// Matt Galloway

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "soft.h"

using namespace Rcpp;





//' @title Lasso Regression (c++)
//' 
//' @description Computes the coefficient estimates for lasso-penalized linear regression.
//' 
//' @details For details on the implementation of 'GLASSOO', see the vignette
//' \url{https://mgallow.github.io/GLASSOO/}.
//'
//' @param XX crossproduct of nxp data matrix.
//' @param XY crossproduct of nxp data matrix and nx1 matrix of response values.
//' @param initB initialization for beta regression coefficients.
//' @param lam tuning parameter for lasso regularization term. Defaults to \code{lam = 0.1}.
//' @param crit criterion for convergence. Criterion \code{loss} will loop until the relative change in the objective for each response after an iteration is less than \code{tol}. Criterion \code{avg} will loop until the average absolute change for each response is less than \code{tol} times tolerance multiple. Similary, criterion \code{max} will loop until the maximum absolute change is less than \code{tol} times tolerance multiple. Defaults to \code{loss}.
//' @param tol tolerance for algorithm convergence. Defaults to 1e-4.
//' @param maxit maximum iterations. Defaults to 1e4.
//' 
//' @return returns list of returns which includes:
//' \item{Iterations}{number of iterations.}
//' \item{Coefficients}{estimated regression coefficients.}
//' \item{H}{H matrix}
//' 
//' @references
//' \itemize{
//' \item Friedman, Jerome, et al. "Pathwise coordinate optimization." \emph{The Annals of Applied Statistics} 1.2 (2007): 302-332. \url{https://arxiv.org/pdf/0708.1485.pdf}
//' \item Tibshirani, Robert. 1996. "Regression Shrinkage and Selection via the Lasso." \emph{Journal of the Royal Statistical Society. Series B (Methodological)}. JSTOR: 267-288.
//' \item Tibshirani, Robert, Bien, Jacob, Friedman, Jerome, Hastie, Trevor, Simon, Noah, Jonathan, Taylor, and Tibshirani, Ryan J. "Strong Rules for Discarding Predictors in Lasso-Type Problems." \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}. Wiley Online Library 74 (2): 245-266.
//' \item Ghaoui, Laurent El, Viallon, Vivian, and Rabbani, Tarek. 2010. "Safe Feature Elimination for the Lasso and Sparse Supervised Learning Problems." \emph{arXiv preprint arXiv: 1009.4219}.
//' \item Osborne, Michael R, Presnell, Brett, and Turlach, Berwin A. "On the Lasso and its Dual." \emph{Journal of Computational and Graphical Statistics}. Taylor and Francis 9 (2): 319-337.

//' }
//' 
//' @author Matt Galloway \email{gall0441@@umn.edu}
//' 
//' @keywords internal
//' 
//' @export
//'
// [[Rcpp::export]]
List lassoc(const arma::mat &XX, const arma::mat &XY, const arma::colvec &initB, const arma::colvec &initH, const double lam = 0.1, std::string crit = "loss", const double tol = 1e-4, const double maxit = 1e4){
  
  // allocate memory
  bool criterion = true, check = true;
  int P = XX.n_cols, iter = 0;
  double max, mult, temp = 0, loss2 = 0, loss = 0;
  arma::mat B(initB), B2(initB), H(initH), H2(initH), maxes;
  maxes = arma::abs(XY);
  max = 2*lam - maxes.max();
  mult = arma::accu(arma::abs(XY/XX.diag()));

  // loop until convergence
  while (criterion && (iter < maxit)){
    
    // update values
    iter++;
    B = B2;
    loss = loss2;
    
    // loop over all rows in column r of beta
    for (int p = 0; p < P; p++){
      if ((max > maxes(p)) && check){
        
        // set entry equal to zero if fails STRONG rules for lasso
        B2[p] = 0;
        
      } else {
        
        // otherwise update betas by soft thresholding
        B2[p] = softc(XY[p] - H2[p], lam)/XX(p, p);
        
      }
      
      // if updated beta is different, update H matrix
      if (B2[p] != B[p]){
        
        // update each element in column r of H, except p
        temp = 0;
        H = H2;
        for (int p_ = 0; p_ < P; p_++){
          if (p_ != p){
            
            // update all rows except row p
            H2[p_] -= XX(p_, p)*(B[p] - B2[p]);
            
            // create temporary sum used in loss function, if necessary
            if (crit == "loss"){
              temp += B2[p_]*XX(p_, p);
              
            }
          }
        }
        
        // update loss, if necessary
        if (crit == "loss"){
          loss2 += (B[p] - B2[p])*(XY[p] - temp) + (std::pow(B2[p], 2) - std::pow(B[p], 2))*XX(p, p)/2 + lam*(std::abs(B2[p]) - std::abs(B[p]));
          
        }
      }
    }
    
    // stopping criterion
    if (crit == "loss"){
      
      // compute loss improvement
      criterion = (std::abs((loss2 - loss)/loss) > tol);
      
    } else if (crit == "avg") {
      
      // compute estimate avg change
      criterion = (arma::accu(arma::abs(B2 - B)) > tol*mult);
      
    } else {
      
      // compute estimate max change
      criterion = (arma::abs(B2 - B).max() > tol*mult/P);
      
    }
    
    // R_CheckUserInterrupt
    if (iter % 1000 == 0){
      R_CheckUserInterrupt();
    }
    
    // check KKT conditions for STRONG rule
    if (!criterion){
      if (arma::any(arma::vectorise(arma::abs(-XY + XX*B2 + lam*arma::sign(B2))) > 1)){
        
        // if true, then recalculate with check = false
        criterion = true;
        check = false;
        
      }
    }
    
  }
  
  return List::create(Named("Iterations") = iter,
                      Named("Coefficients") = B2,
                      Named("H") = H2);
  
}

