#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
List reduce(arma::mat X, int p){
  
  int P = X.n_cols;
  arma::mat Y(P - 1, P - 1, arma::fill::zeros);
    
  for (int i = 0; i < (P - 1); i++){
    for (int j = 0; j < (P - 1); j++){
      if ((i < p) && (j < p)){
        Y(i, j) = X(i, j);
      } else if ((i < p) && (j >= p)){
        Y(i, j) = X(i, j + 1);
      } else if ((i >= p) && (j < p)){
        Y(i, j) = X(i + 1, j);
      } else {
        Y(i, j) = X(i + 1, j + 1);
      }
    }
  }
  
  return List::create(Named("X") = X,
                      Named("Y") = Y);
  
}


// [[Rcpp::export]]
List reduce2(arma::mat X, int p){
  
  int P = X.n_cols;
  arma::mat Y(P - 1, 1, arma::fill::zeros);
  
  for (int i = 0; i < (P - 1); i++){
    if (i < p){
      Y(i, 0) = X(i, p);
    } else {
      Y(i, 0) = X(i + 1, p);
    }
  }
  
  return List::create(Named("X") = X,
                      Named("Y") = Y);
  
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
X = tridiag(p = 5)$Omega %>% round(3)
reduce2(X, 0)
reduce2(X, 1)
reduce2(X, 2)
reduce2(X, 3)
reduce2(X, 4)
reduce2(X, 5)

*/
