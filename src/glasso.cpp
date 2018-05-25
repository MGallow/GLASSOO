// Matt Galloway

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "lasso.h"
#include "misc.h"

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
//' @param initOmega initialization matrix for Omega used to initialize the Betas
//' @param lam tuning parameter for lasso penalty.
//' @param crit_out criterion for convergence in outer (blockwise) loop. Criterion \code{avg} will loop until the average absolute parameter change is less than \code{tol_out} times tolerance multiple. Criterion \code{max} will loop until the maximum change in the estimated Sigma after an iteration over the parameter set is less than \code{tol_out}. Defaults to \code{avg}.
//' @param crit_in criterion for convergence in inner (lasso) loop. Criterion for convergence. Criterion \code{loss} will loop until the relative change in the objective for each response after an iteration is less than \code{tol_in}. Criterion \code{avg} will loop until the average absolute change for each response is less than \code{tol_in} times tolerance multiple. Similary, criterion \code{max} will loop until the maximum absolute change is less than \code{tol_in} times tolerance multiple. Defaults to \code{loss}.
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
//' \item Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. 'Sparse inverse covariance estimation with the graphical lasso.' \emph{Biostatistics} 9.3 (2008): 432-441.
//' \item Banerjee, Onureen, Ghauoui, Laurent El, and d'Aspremont, Alexandre. 2008. "Model Selection through Sparse Maximum Likelihood Estimation for Multivariate Gaussian or Binary Data." \emph{Journal of Machine Learning Research} 9: 485-516.
//' \item Tibshirani, Robert. 1996. "Regression Shrinkage and Selection via the Lasso." \emph{Journal of the Royal Statistical Society. Series B (Methodological)}. JSTOR: 267-288.
//' \item Meinshausen, Nicolai and Buhlmann, Peter. 2006. "High-Dimensional Graphs and Variable Selection with the Lasso." \emph{The Annals of Statistics}. JSTOR: 1436-1462.
//' \item Witten, Daniela M, Friedman, Jerome H, and Simon, Noah. 2011. "New Insights and Faster computations for the Graphical Lasso." \emph{Journal of Computation and Graphical Statistics}. Taylor and Francis: 892-900.
//' \item Tibshirani, Robert, Bien, Jacob, Friedman, Jerome, Hastie, Trevor, Simon, Noah, Jonathan, Taylor, and Tibshirani, Ryan J. "Strong Rules for Discarding Predictors in Lasso-Type Problems." \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}. Wiley Online Library 74 (2): 245-266.
//' \item Ghaoui, Laurent El, Viallon, Vivian, and Rabbani, Tarek. 2010. "Safe Feature Elimination for the Lasso and Sparse Supervised Learning Problems." \emph{arXiv preprint arXiv: 1009.4219}.
//' \item Osborne, Michael R, Presnell, Brett, and Turlach, Berwin A. "On the Lasso and its Dual." \emph{Journal of Computational and Graphical Statistics}. Taylor and Francis 9 (2): 319-337.
//' \item Rothman, Adam. 2017. "STAT 8931 notes on an algorithm to compute the Lasso-penalized Gausssian likelihood precision matrix estimator." 
//' }
//' 
//' @author Matt Galloway \email{gall0441@@umn.edu}
//' 
//' @export
//'
//' @keywords internal
//'
// [[Rcpp::export]]
List GLASSOc(const arma::mat &S, const arma::mat &initSigma, const arma::mat &initOmega, const double lam, std::string crit_out = "avg", std::string crit_in = "loss", const double tol_out = 1e-4, const double tol_in = 1e-4, const int maxit_out = 1e4, const int maxit_in = 1e4){
  
  // allocate memory
  bool criterion = true;
  int P = S.n_cols, iter = 0;
  double mult;
  arma::mat Sigmatemp(P - 1, P - 1, arma::fill::zeros), Stemp, Omegatemp(P, 1, arma::fill::zeros), ind(P - 1, 1, arma::fill::ones);
  arma::mat maxes, Beta, Betas, Stemps, zeros, H, Omega(P, P, arma::fill::eye), Sigma(initSigma), Sigma2(initSigma), Sminus(S);
  Beta = zeros = Stemp = arma::zeros<arma::mat>(P - 1, 1);
  Betas = Stemps = H = arma::zeros<arma::mat>(P - 1, P);
  Sminus -= arma::diagmat(Sminus);
  mult = arma::accu(arma::abs(Sminus));
  
  // initialize Betas and Stemps
  extractdividec(initOmega, Betas);
  extractc(S, Stemps);
  maxes = arma::abs(arma::max(Stemps, 0));
  
  // loop until convergence
  while (criterion && (iter < maxit_out)){
    
    // update values
    iter++;
    Sigma = Sigma2;
    
    // blockwise coordinate descent
    for (int p = 0; p < P; p++){
      if (lam > maxes(p)){
        
        // if true, then set Stemp = Betas[, p] = 0 if lam > maxes(p)
        Stemp = zeros;
        Betas.col(p) = zeros;
        
      } else {
        
        // update Beta with pth column of Betas
        Beta = Betas.col(p);
        
        // set Sigmatemp = Sigma[-p, -p]
        reducec(Sigma2, Sigmatemp, p);
        
        // initialize Stemp and H matrix used in lassoc (Sigmatemp.minus*Beta)
        Stemp = Stemps.col(p);
        H = Sigmatemp*Beta - Sigmatemp.diag() % Beta;
        
        // execute LASSO
        List LASSO = lassoc(Sigmatemp, Stemp, Beta, H, ind, lam, crit_in, tol_in, maxit_in);
        Betas.col(p) = as<arma::mat>(LASSO["Coefficients"]);
        
        // update Stemp = Sigma12
        Stemp = Sigmatemp*Betas.col(p);
        
      }
      
      // update Sigma[-p, p] = Sigma[p, -p] = Stemp
      updatec(Sigma2, Stemp, p);
      
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
  extractc(Sigma, Stemps);
  for (int p = 0; p < P; p++){
    
    // update Omega[p, p] and Omegatemp = Omega12
    Omega(p, p) = 1/(Sigma2(p, p) - arma::accu(Stemps.col(p) % Betas.col(p)));
    Omegatemp = -Omega(p, p)*Betas.col(p);
    
    // set Omega[-p, p] = Omega[p, -p] = Omegatemp
    updatec(Omega, Omegatemp, p);
    
  }
  
  return List::create(Named("Iterations") = iter,
                      Named("lam") = lam,
                      Named("Omega") = Omega,
                      Named("Sigma") = Sigma2);
  
}

