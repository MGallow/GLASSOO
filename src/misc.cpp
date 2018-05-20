// Matt Galloway

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;


//' @title Reduce (c++)
//' @description This function produces an identical matrix with the pth column and pth row removed (ie: Y = X[-p, -p]).
//'
//' @param X matrix
//' @param Y matrix
//' @param p column number
//' @keywords internal
//'

void reducec(const arma::mat &X, arma::mat &Y, const int &p) {
  
  // loop over all elements of Y and update accordingly
  for (int i = 0; i < Y.n_rows; i++){
    for (int j = 0; j < Y.n_cols; j++){
      
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
  
}



////-----------------------------------------------------



//' @title Extract (c++)
//' @description This function extracts the pth column with the pth row removed from a matrix (ie: Y = X[-p, p]).
//'
//' @param X matrix
//' @param Y matrix
//' @param p column number
//' @keywords internal
//'

void extractc(const arma::mat &X, arma::mat &Y, const int &p) {
  
  // loop over all elements of Y and update accordingly
  for (int i = 0; i < Y.n_rows; i++){
    
    if (i < p){
      Y(i, 0) = X(i, p);
      
    } else {
      Y(i, 0) = X(i + 1, p);
    }
  }
  
}




////-----------------------------------------------------



//' @title Update (c++)
//' @description This function updates an existing matrix by setting Y equal to the pth column with the pth row of X removed and the pth row with the pth column of X removed (ie: X[-p, p] = X[p, -p] = Y).
//'
//' @param X matrix
//' @param Y matrix
//' @param p column number
//' @keywords internal
//'

void updatec(arma::mat &X, const arma::mat &Y, const int &p) {
  
  // loop over all elements of Y and update accordingly
  for (int i = 0; i < Y.n_rows; i++){
    
    if (i < p){
      X(i, p) = Y(i, 0);
      X(p, i) = Y(i, 0);
      
    } else {
      X(i + 1, p) = Y(i, 0);
      X(p, i + 1) = Y(i, 0);
    }
  }
  
}
