## Matt Galloway



#' @title Lasso Regression (c++)
#' 
#' @description Computes the coefficient estimates for lasso-penalized linear regression.
#' 
#' @details For details on the implementation of 'GLASSO', see the vignette
#' \url{https:#mgallow.github.io/GLASSO/}.
#'
#' @param X matrix or data frame
#' @param Y matrix or data frame of response values
#' @param lam tuning parameter for lasso regularization term. Defaults to 'lam = 0.1'
#' @param crit criterion for convergence. Criterion \code{loss} will loop until the change in the objective after an iteration over the parameter set is less than \code{tol}. Criterion \code{sum} will loop until the sum change in the estimate after an interation over the parameter set is less than \code{tol} times tolerance multiple. Similary, criterion \code{max} will loop until the maximum change is less than \code{tol} times tolerance multiple. Defaults to \code{loss}.
#' @param tol tolerance for algorithm convergence. Defaults to 1e-4
#' @param maxit maximum iterations. Defaults to 1e4
#' @param ind optional matrix specifying which coefficients will be penalized.
#' 
#' @return returns list of returns which includes:
#' \item{Iterations}{number of iterations.}
#' \item{Loss}{value of the objective function.}
#' \item{Coefficients}{estimated regression coefficients.}
#' 
#' @references
#' \itemize{
#' \item 
#' For more information on the ADMM algorithm, see: \cr
#' Boyd, Stephen, Neal Parikh, Eric Chu, Borja Peleato, Jonathan Eckstein, and others. 2011. 'Distributed Optimization and Statistical Learning via the Alternating Direction Method of Multipliers.' \emph{Foundations and Trends in Machine Learning} 3 (1). Now Publishers, Inc.: 1-122.\cr
#' \url{https:#web.stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf}
#' }
#' 
#' @author Matt Galloway \email{gall0441@@umn.edu}
#' 
#' @export

# we define the lasso function
LASSO = function(X, Y, lam = 0.1, crit = c("loss", 
    "avg", "max"), tol = 1e-04, maxit = 10000, 
    ind = matrix(1, ncol(X), ncol(Y))) {
    
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
    LASSO = lassoc(XX = XX, XY = XY, initB = init, 
        ind = ind, lam = lam, crit = crit, tol = tol, 
        maxit = maxit)
    
    
    # compute loss
    loss = sum((Y - X %*% LASSO$Coefficients)^2)/2 + 
        lam * sum(abs(ind * LASSO$Coefficients))
    
    returns = list(Call = call, Iterations = LASSO$Iterations, 
        Loss = loss, Coefficients = LASSO$Coefficients)
    return(returns)
    
}




##-----------------------------------------------------------------------------------
