## Matt Galloway


#' @title Penalized precision matrix estimation via ADMM
#' 
#' @description Penalized precision matrix estimation using the graphical lasso (glasso) algorithm.
#' Consider the case where \eqn{X_{1}, ..., X_{n}} are iid \eqn{N_{p}(\mu,
#' \Sigma)} and we are tasked with estimating the precision matrix,
#' denoted \eqn{\Omega \equiv \Sigma^{-1}}. This function solves the
#' following optimization problem:
#' \describe{
#' \item{Objective:}{
#' \eqn{\hat{\Omega}_{\lambda} = \arg\min_{\Omega \in S_{+}^{p}}
#' \left\{ Tr\left(S\Omega\right) - \log \det\left(\Omega \right) +
#' \lambda \left\| \Omega \right\|_{1} \right\}}}
#' }
#' where \eqn{\lambda > 0} and we define
#' \eqn{\left\|A \right\|_{1} = \sum_{i, j} \left| A_{ij} \right|}.
#' 
#' @details For details on the implementation of 'GLASSOO', see the vignette
#' \url{https://mgallow.github.io/GLASSOO/}.
#' 
#' @param X option to provide a nxp data matrix. Each row corresponds to a single observation and each column contains n observations of a single feature/variable.
#' @param S option to provide a pxp sample covariance matrix (denominator n). If argument is \code{NULL} and \code{X} is provided instead then \code{S} will be computed automatically.
#' @param lam tuning parameter for elastic net penalty. Defaults to grid of values \code{10^seq(-5, 5, 0.5)}.
#' @param diagonal option to penalize the diagonal elements of the estimated precision matrix (\eqn{\Omega}). Defaults to \code{FALSE}.
#' @param crit.out criterion for convergence in outer (blockwise) loop. Criterion \code{avg} will loop until the average absolute parameter change is less than \code{tol.out} times tolerance multiple. Criterion \code{max} will loop until the maximum change in the estimated Sigma after an iteration over the parameter set is less than \code{tol.out}. Defaults to \code{avg}.
#' @param crit.in criterion for convergence in inner (lasso) loop. Criterion for convergence. Criterion \code{loss} will loop until the change in the objective for each response after an iteration is less than \code{tol.in}. Criterion \code{avg} will loop until the average absolute change for each response is less than \code{tol.in} times tolerance multiple. Similary, criterion \code{max} will loop until the maximum absolute change is less than \code{tol.in} times tolerance multiple. Defaults to \code{loss}.
#' @param tol.out convergence tolerance for outer (blockwise) loop. Defaults to 1e-4.
#' @param tol.in convergence tolerance for inner (lasso) loop. Defaults to 1e-4.
#' @param maxit.out maximum number of iterations for outer (blockwise) loop. Defaults to 1e4.
#' @param maxit.in maximum number of iterations for inner (lasso) loop. Defaults to 1e4.
#' 
#' @return returns class object \code{ADMMsigma} which includes:
#' \item{Call}{function call.}
#' \item{Iterations}{number of iterations}
#' \item{Lambdas}{grid of lambda values for CV.}
#' \item{maxit.out}{maximum number of iterations for outer (blockwise) loop.}
#' \item{maxit.in}{maximum number of iterations for inner (lasso) loop.}
#' \item{Omega}{estimated penalized precision matrix.}
#' \item{Sigma}{estimated covariance matrix from the penalized precision matrix (inverse of Omega).}
#' \item{Loglik}{penalized log-likelihood for Omega}
#' 
#' @references
#' \itemize{
#' \item 
#' For more information on the graphical lasso algorithm, see: \cr
#' Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. 'Sparse inverse covariance estimation with the graphical lasso.' \emph{Biostatistics} 9.3 (2008): 432-441.\cr
#' \url{http://statweb.stanford.edu/~tibs/ftp/glasso-bio.pdf}
#' }
#' 
#' @author Matt Galloway \email{gall0441@@umn.edu}
#' 
#' @seealso \code{\link{plot.GLASSO}}
#' 
#' @export
#' 
#' @examples
#' # generate data from a tridiagonal precision matrix
#' # first compute covariance matrix (can confirm inverse is tridiagonal)
#' S = matrix(0.7, nrow = 5, ncol = 5)
#' for (i in 1:5){
#'  for (j in 1:5){
#'    S[i, j] = 0.7^abs(i - j)
#'  }
#' }
#'
#' # generate 100 x 5 matrix with rows drawn from iid N_p(0, S)
#' Z = matrix(rnorm(100*5), nrow = 100, ncol = 5)
#' out = eigen(S, symmetric = TRUE)
#' S.sqrt = out$vectors %*% diag(out$values^0.5)
#' S.sqrt = S.sqrt %*% t(out$vectors)
#' X = Z %*% S.sqrt
#'
#' # lasso penalty (lam = 0.1)
#' GLASSO(X, lam = 0.1)
#'
#' # produce CV heat map for GLASSO
#' plot(GLASSO(X))

# we define the GLASSO precision matrix
# estimation function
GLASSO = function(X = NULL, S = NULL, lam = 0.1, 
    diagonal = FALSE, crit.out = c("avg", "max"), 
    crit.in = c("loss", "avg", "max"), tol.out = 1e-04, 
    tol.in = 1e-04, maxit.out = 10000, maxit.in = 10000) {
    
    # checks
    if (is.null(X) && is.null(S)) {
        stop("Must provide entry for X or S!")
    }
    if (!all(lam > 0)) {
        stop("lam must be positive!")
    }
    if (!(all(c(tol.out, tol.in, maxit.out, maxit.in) > 
        0))) {
        stop("Entry must be positive!")
    }
    if (all(c(maxit.out, maxit.in)%%1 != 0)) {
        stop("Entry must be an integer!")
    }
    if (length(lam) > 1) {
        stop("Must provide single value for lam.")
    }
    
    # match values
    crit.out = match.arg(crit.out)
    crit.in = match.arg(crit.in)
    lam = sort(lam)
    call = match.call()
    n = ifelse(is.null(X), nrow(S), nrow(X))
    
    # compute sample covariance matrix, if necessary
    if (is.null(S)) {
        S = (nrow(X) - 1)/nrow(X) * cov(X)
    }
    
    # specify initial estimate for Sigma
    alpha = min(c(lam/max(abs(S - diag(S))), 1))
    init = (1 - alpha) * S
    diag(init) = diag(S)
    
    # compute estimate at specified tuning parameter
    GLASSO = GLASSOc(S = S, initSigma = init, lam = lam, 
        crit_out = crit.out, crit_in = crit.in, tol_out = tol.out, 
        tol_in = tol.in, maxit_out = maxit.out, maxit_in = maxit.in)
    
    
    # option to penalize diagonal
    if (diagonal) {
        C = 1
    } else {
        C = 1 - diag(ncol(S))
    }
    
    # compute loglik
    n = ifelse(is.null(X), nrow(S), nrow(X))
    loglik = (-n/2) * (sum(GLASSO$Omega * S) - determinant(GLASSO$Omega, 
        logarithm = TRUE)$modulus[1] + GLASSO$lam * 
        sum(abs(C * GLASSO$Omega)))
    
    
    # return values
    tuning = matrix(c(log10(lam), lam), ncol = 2)
    colnames(tuning) = c("log10(lam)", "alpha")
    returns = list(Call = call, Iterations = GLASSO$Iterations, 
        Lambdas = GLASSO$lam, maxit.out = maxit.out, 
        maxit.in = maxit.in, Omega = GLASSO$Omega, 
        Sigma = GLASSO$Sigma, Loglik = loglik)
    
    class(returns) = "GLASSO"
    return(returns)
    
}






##-----------------------------------------------------------------------------------



#' @title Print GLASSO object
#' @description Prints GLASSO object and suppresses output if needed.
#' @param x class object GLASSO
#' @param ... additional arguments.
#' @keywords internal
#' @export
print.GLASSO = function(x, ...) {
    
    # print warning if maxit reached
    if (x$maxit.out <= x$Iterations) {
        print("Maximum iterations reached...!")
    }
    
    # print call
    cat("\nCall: ", paste(deparse(x$Call), sep = "\n", 
        collapse = "\n"), "\n", sep = "")
    
    # print iterations
    cat("\nIterations:\n")
    print.default(x$Iterations, quote = FALSE)
    
    # print optimal tuning parameters cat('\nTuning
    # parameters:\n') print.default(round(x$Tuning,
    # 3), print.gap = 2L, quote = FALSE)
    
    # print loglik
    cat("\nLog-likelihood: ", paste(round(x$Loglik, 
        5), sep = "\n", collapse = "\n"), "\n", sep = "")
    
    # print Omega if dim <= 10
    if (nrow(x$Omega) <= 10) {
        cat("\nOmega:\n")
        print.default(round(x$Omega, 5))
    } else {
        cat("\n(...output suppressed due to large dimension!)\n")
    }
    
}
