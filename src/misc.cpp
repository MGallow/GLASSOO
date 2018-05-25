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
//' @description This function updates an existing matrix by setting the pth columm of Y equal to the pth column with the pth row of X (ie: Y[, p] = X[-p, p] for all p).
//'
//' @param X matrix
//' @param Y matrix
//' @keywords internal
//'

void extractc(const arma::mat &X, arma::mat &Y) {
  
  // loop over all elements of Y and update accordingly
  for (int i = 0; i < Y.n_rows; i++){
    for (int j = 0; j < Y.n_cols; j++){
      
      if (i < j){
        Y(i, j) = X(i, j);
        
      } else {
        Y(i, j) = X(i + 1, j);
      }
    }
  }
  
}




////-----------------------------------------------------



//' @title Extract and Divide (c++)
//' @description This function updates an existing matrix by setting the pth columm of Y equal to the pth column with the pth row of X removed divided the pth diagonal element (ie: Y[, p] = X[-p, p]/X[p, p] for all p).
//'
//' @param X matrix
//' @param Y matrix
//' @keywords internal
//'

void extractdividec(const arma::mat &X, arma::mat &Y) {
  
  // loop over all elements of Y and update accordingly
  for (int i = 0; i < Y.n_rows; i++){
    for (int j = 0; j < Y.n_cols; j++){
      
      if (i < j){
        Y(i, j) = X(i, j)/X(j, j);
        
      } else {
        Y(i, j) = X(i + 1, j)/X(j, j);
      }
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





////-----------------------------------------------------



//' @title Nonzeros (c++)
//' @description This function counts the number of nonzero elements in a matrix
//' @param X matrix
//' @keywords internal
//'

int numzeros(arma::mat &X) {
  
  // loop over all elements of X and count nonzeros
  int num = 0;
  arma::mat::iterator it = X.begin();
  arma::mat::iterator it_end = X.end();
  
  for (; it != it_end; ++it){
    if (*it != 0){
      num++;
    }
  }
  
  return(num);
}
