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
#' @param crit criterion for convergence. Criterion \code{obj} will loop until the change in the objective after an iteration over the parameter set is less than \code{tol}. Criterion \code{max} will loop until the maximum change in the estimate after an iteration over the parameter set is less than \code{tol}. Defaults to \code{obj}.

#' @param tol tolerance for algorithm convergence. Defaults to 1e-4
#' @param maxit maximum iterations. Defaults to 1e4
#' @param ind optional matrix specifying which coefficients will be penalized.
#' 
#' @return returns list of returns which includes:
#' \item{Iterations}{number of iterations.}
#' \item{Objective}{value of the objective function.}
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
lasso = function(X, Y, lam = 0.1, crit = "obj", tol = 1e-04, 
    maxit = 10000, ind = NULL) {
    
    # checks
    if (lam <= 0) {
        stop("lam must be positive!")
    }
    
    # initialization
    p = dim(X)[2]
    r = dim(Y)[2]
    X = as.matrix(X)
    Y = as.matrix(Y)
    if (is.null(ind)) {
        ind = matrix(1, nrow = p, ncol = r)
    }
    
    
    # execute lassoc
    LASSO = lassoc(X = X, Y = Y, ind = ind, lam = lam, crit = crit, 
        tol = tol, maxit = maxit)
    
    returns = list(Iterations = LASSO$Iterations, Objective = LASSO$Objective, 
        Coefficients = LASSO$Coefficients)
    return(returns)
    
}




##-----------------------------------------------------------------------------------
