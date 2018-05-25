## Matt Galloway



#' @title Lasso Regression (c++)
#' 
#' @description Computes the coefficient estimates for lasso-penalized linear regression.
#' 
#' @details For details on the implementation of 'GLASSO', see the vignette
#' \url{https:#mgallow.github.io/GLASSO/}.
#'
#' @param X nxp data matrix.
#' @param Y nxr matrix of response values
#' @param lam tuning parameter for lasso regularization term. Defaults to \code{lam = 0.1}.
#' @param crit criterion for convergence. Criterion \code{loss} will loop until the relative change in the objective for each response after an iteration is less than \code{tol}. Criterion \code{avg} will loop until the average absolute change for each response is less than \code{tol} times tolerance multiple. Similary, criterion \code{max} will loop until the maximum absolute change is less than \code{tol} times tolerance multiple. Defaults to \code{loss}.
#' @param tol tolerance for algorithm convergence. Defaults to 1e-4..
#' @param maxit maximum iterations. Defaults to 1e4
#' @param ind optional matrix specifying which coefficients will be penalized.
#' 
#' @return returns list of returns which includes:
#' \item{Call}{function call.}
#' \item{Iterations}{number of iterations.}
#' \item{Loss}{value of the objective function.}
#' \item{Coefficients}{estimated regression coefficients.}
#' 
#' @references
#' \itemize{
#' \item 
#' For more information on the coordinate descent algorithm, see: \cr
#' Friedman, Jerome, et al. 'Pathwise coordinate optimization.' \emph{The Annals of
#' Applied Statistics} 1.2 (2007): 302-332.\cr
#' \url{https://arxiv.org/pdf/0708.1485.pdf}
#' }
#' 
#' @author Matt Galloway \email{gall0441@@umn.edu}
#' 
#' @keywords internal
#' 
#' @export

# we define the lasso function
LASSO = function(X, Y, lam = 0.1, crit = c("loss", "avg", 
    "max"), tol = 1e-04, maxit = 10000, ind = matrix(1, ncol(X), 
    ncol(Y))) {
    
    # checks
    if (is.null(X) || is.null(Y)) {
        stop("Must provide entry for X and Y!")
    }
    if (lam <= 0) {
        stop("lam must be positive!")
    }
    if (tol <= 0) {
        stop("Entry must be positive!")
    }
    if (maxit%%1 != 0) {
        stop("Entry must be an integer!")
    }
    
    # match values
    X = as.matrix(X)
    Y = as.matrix(Y)
    crit = match.arg(crit)
    lam = sort(lam)
    call = match.call()
    
    # save values
    XX = crossprod(X)
    XY = crossprod(X, Y)
    
    # execute lassoc
    init = matrix(0, nrow = ncol(X), ncol = ncol(Y))
    LASSO = lassoc(XX = XX, XY = XY, initB = init, initH = init, 
        ind = ind, lam = lam, crit = crit, tol = tol, maxit = maxit)
    
    
    # compute loss
    loss = sum((Y - X %*% LASSO$Coefficients)^2)/2 + lam * 
        sum(abs(ind * LASSO$Coefficients))
    
    returns = list(Call = call, Iterations = LASSO$Iterations, 
        Loss = loss, Coefficients = LASSO$Coefficients)
    return(returns)
    
}




##-----------------------------------------------------------------------------------
