## Matt Galloway


#' @title Parallel CV (uses CV_GLASSOc)
#' @description Parallel implementation of cross validation.
#'
#' @param X nxp data matrix. Each row corresponds to a single observation and each column contains n observations of a single feature/variable.
#' @param lam positive tuning parameters for elastic net penalty. If a vector of parameters is provided, they should be in increasing order. Defaults to grid of values \code{10^seq(-2, 2, 0.2)}.
#' @param diagonal option to penalize the diagonal elements of the estimated precision matrix (\eqn{\Omega}). Defaults to \code{FALSE}.
#' @param crit.out criterion for convergence in outer (blockwise) loop. Criterion \code{avg} will loop until the average absolute parameter change is less than \code{tol.out} times tolerance multiple. Criterion \code{max} will loop until the maximum change in the estimated Sigma after an iteration over the parameter set is less than \code{tol.out}. Defaults to \code{avg}.
#' @param crit.in criterion for convergence in inner (lasso) loop. Criterion for convergence. Criterion \code{loss} will loop until the relative change in the objective for each response after an iteration is less than \code{tol.in}. Criterion \code{avg} will loop until the average absolute change for each response is less than \code{tol.in} times tolerance multiple. Similary, criterion \code{max} will loop until the maximum absolute change is less than \code{tol.in} times tolerance multiple. Defaults to \code{loss}.
#' @param tol.out convergence tolerance for outer (blockwise) loop. Defaults to 1e-4.
#' @param tol.in convergence tolerance for inner (lasso) loop. Defaults to 1e-4.
#' @param maxit.out maximum number of iterations for outer (blockwise) loop. Defaults to 1e4.
#' @param maxit.in maximum number of iterations for inner (lasso) loop. Defaults to 1e4.
#' @param adjmaxit.out adjusted maximum number of iterations. During cross validation this option allows the user to adjust the maximum number of iterations after the first \code{lam} tuning parameter has converged (for each \code{alpha}). This option is intended to be paired with \code{warm} starts and allows for 'one-step' estimators. Defaults to NULL.
#' @param K specify the number of folds for cross validation.
#' @param crit.cv cross validation criterion (\code{loglik}, \code{AIC}, or \code{BIC}). Defaults to \code{loglik}.
#' @param start specify \code{warm} or \code{cold} start for cross validation. Default is \code{warm}.
#' @param cores option to run CV in parallel. Defaults to \code{cores = 1}.
#' @param trace option to display progress of CV. Choose one of \code{progress} to print a progress bar, \code{print} to print completed tuning parameters, or \code{none}.
#' 
#' @return returns list of returns which includes:
#' \item{lam}{optimal tuning parameter.}
#' \item{min.error}{minimum average cross validation error (cv.crit) for optimal parameters.}
#' \item{avg.error}{average cross validation error (cv.crit) across all folds.}
#' \item{cv.error}{cross validation errors (cv.crit).}
#' 
#' @keywords internal

# we define the CV_GLASSOc function
CVP_GLASSO = function(X = NULL, lam = 10^seq(-2, 2, 0.2), 
    diagonal = FALSE, crit.out = c("avg", "max"), crit.in = c("loss", 
        "avg", "max"), tol.out = 1e-04, tol.in = 1e-04, maxit.out = 10000, 
    maxit.in = 10000, adjmaxit.out = NULL, K = 5, crit.cv = c("loglik", 
        "AIC", "BIC"), start = c("warm", "cold"), cores = 1, 
    trace = c("progress", "print", "none")) {
    
    # match values
    crit.out = match.arg(crit.out)
    crit.in = match.arg(crit.in)
    crit.cv = match.arg(crit.cv)
    start = match.arg(start)
    trace = match.arg(trace)
    lam = sort(lam)
    
    # make cluster and register cluster
    num_cores = detectCores()
    if (cores > num_cores) {
        print(paste("Only detected", num_cores, "cores...", 
            sep = " "))
    }
    if (cores > K) {
        print("Number of cores exceeds K... setting cores = K")
        cores = K
    }
    
    cluster = makeCluster(cores)
    registerDoParallel(cluster)
    
    # use cluster for each fold in CV
    n = nrow(X)
    ind = sample(n)
    k = NULL
    CV = foreach(k = 1:K, .packages = "GLASSOO", .combine = "cbind", 
        .inorder = FALSE) %dopar% {
        
        leave.out = ind[(1 + floor((k - 1) * n/K)):floor(k * 
            n/K)]
        
        # training set
        X.train = X[-leave.out, , drop = FALSE]
        X_bar = apply(X.train, 2, mean)
        X.train = scale(X.train, center = X_bar, scale = FALSE)
        
        # validation set
        X.valid = X[leave.out, , drop = FALSE]
        X.valid = scale(X.valid, center = X_bar, scale = FALSE)
        
        # sample covariances
        S.train = crossprod(X.train)/(dim(X.train)[1])
        S.valid = crossprod(X.valid)/(dim(X.valid)[1])
        
        # run foreach loop on CVP_ADMMc
        CVP_GLASSOc(n = nrow(X.valid), S_train = S.train, 
            S_valid = S.valid, lam = lam, diagonal = diagonal, 
            crit_out = crit.out, crit_in = crit.in, tol_out = tol.out, 
            tol_in = tol.in, maxit_out = maxit.out, maxit_in = maxit.in, 
            adjmaxit_out = adjmaxit.out, crit_cv = crit.cv, 
            start = start, trace = trace)
        
    }
    
    # determine optimal tuning parameters
    AVG = as.matrix(apply(CV, 1, mean))
    best = which(AVG == min(AVG), arr.ind = TRUE)
    error = min(AVG)
    best_lam = lam[best[1]]
    
    # stop cluster
    stopCluster(cluster)
    
    # return best lam and alpha values
    return(list(lam = best_lam, min.error = error, avg.error = AVG, 
        cv.error = CV))
    
}




