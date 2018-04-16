// Matt Galloway

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "lasso.h"

using namespace Rcpp;




//' @title Penalized precision matrix estimation (c++)
//' 
//' @description Penalized precision matrix estimation using the ADMM algorithm
//' 
//' @details For details on the implementation of 'ADMMsigma', see the vignette
//' \url{https://mgallow.github.io/ADMMsigma/}.
//'
//' @param S pxp sample covariance matrix (denominator n).
//' @param lam tuning parameter for elastic net penalty. Defaults to grid of values ???.
//' @param diagonal option to penalize the diagonal elements of the estimated precision matrix (\eqn{\Omega}). Defaults to \code{FALSE}.
//' @param crit criterion for convergence (\code{ADMM}, \code{grad}, or \code{loglik}). If \code{crit != ADMM} then \code{tol1} will be used as the convergence tolerance. Default is \code{ADMM}.
//' @param tol1 absolute convergence tolerance. Defaults to 1e-4.
//' @param tol2 relative convergence tolerance. Defaults to 1e-4.
//' @param maxit maximum number of iterations. Defaults to 1e4.
//' 
//' @return returns list of returns which includes:
//' \item{Iterations}{number of iterations.}
//' \item{lam}{optimal tuning parameters.}
//' \item{Omega}{estimated penalized precision matrix.}
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
List GLASSOc(const arma::mat &S, const double lam, bool diagonal = false, std::string crit = "ADMM", const double tol1 = 1e-4, const double tol2 = 1e-4, const int maxit = 1e4){

  // allocate memory
  bool criterion = true;
  int p = S.n_cols;
  int iter = 0;
  double lik, lik2, sgn, logdet;
  lik = lik2 = sgn = logdet = 0;
  arma::mat Omega, Sigma, Sigma2, C;
  C = arma::ones<arma::mat>(p, p);
  
  // option to penalize diagonal elements
  if (diagonal){
    C -= arma::eye<arma::mat>(p, p);
  }

  // loop until convergence
  while (criterion && (iter <= maxit)){

    // blockwise coordinate descent
    

    // stopping criterion
    iter++;
    if (crit == "loglik"){

      // compute likelihood improvement (close enough)
      arma::log_det(logdet, sgn, Omega);
      lik2 = (-p/2)*(arma::accu(Omega % S) - logdet + lam*arma::accu(C % arma::abs(Omega)));
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

