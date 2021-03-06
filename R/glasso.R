## Matt Galloway


#' @title Penalized precision matrix estimation
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
#' @param nlam number of \code{lam} tuning parameters for penalty term generated from \code{lam.min.ratio} and \code{lam.max} (automatically generated). Defaults to 10.
#' @param lam.min.ratio smallest \code{lam} value provided as a fraction of \code{lam.max}. The function will automatically generate \code{nlam} tuning parameters from \code{lam.min.ratio*lam.max} to \code{lam.max} in log10 scale. \code{lam.max} is calculated to be the smallest \code{lam} such that all off-diagonal entries in \code{Omega} are equal to zero (\code{alpha} = 1). Defaults to 1e-2.
#' @param lam option to provide positive tuning parameters for penalty term. This will cause \code{nlam} and \code{lam.min.ratio} to be disregarded. If a vector of parameters is provided, they should be in increasing order. Defaults to NULL.
#' @param diagonal option to penalize the diagonal elements of the estimated precision matrix (\eqn{\Omega}). Defaults to \code{FALSE}.
#' @param path option to return the regularization path. This option should be used with extreme care if the dimension is large. If set to TRUE, cores must be set to 1 and errors and optimal tuning parameters will based on the full sample. Defaults to FALSE.
#' @param crit.out criterion for convergence in outer (blockwise) loop. Criterion \code{avg} will loop until the average absolute parameter change is less than \code{tol.out} times tolerance multiple. Criterion \code{max} will loop until the maximum change in the estimated Sigma after an iteration over the parameter set is less than \code{tol.out}. Defaults to \code{avg}.
#' @param crit.in criterion for convergence in inner (lasso) loop. Criterion for convergence. Criterion \code{loss} will loop until the relative change in the objective for each response after an iteration is less than \code{tol.in}. Criterion \code{avg} will loop until the average absolute change for each response is less than \code{tol.in} times tolerance multiple. Similary, criterion \code{max} will loop until the maximum absolute change is less than \code{tol.in} times tolerance multiple. Defaults to \code{loss}.
#' @param tol.out convergence tolerance for outer (blockwise) loop. Defaults to 1e-4.
#' @param tol.in convergence tolerance for inner (lasso) loop. Defaults to 1e-4.
#' @param maxit.out maximum number of iterations for outer (blockwise) loop. Defaults to 1e4.
#' @param maxit.in maximum number of iterations for inner (lasso) loop. Defaults to 1e4.
#' @param adjmaxit.out adjusted maximum number of iterations. During cross validation this option allows the user to adjust the maximum number of iterations after the first \code{lam} tuning parameter has converged. This option is intended to be paired with \code{warm} starts and allows for 'one-step' estimators. Defaults to NULL.
#' @param K specify the number of folds for cross validation.
#' @param crit.cv cross validation criterion (\code{loglik}, \code{AIC}, or \code{BIC}). Defaults to \code{loglik}.
#' @param start specify \code{warm} or \code{cold} start for cross validation. Default is \code{warm}.
#' @param cores option to run CV in parallel. Defaults to \code{cores = 1}.
#' @param trace option to display progress of CV. Choose one of \code{progress} to print a progress bar, \code{print} to print completed tuning parameters, or \code{none}.
#' 
#' @return returns class object \code{GLASSOO} which includes:
#' \item{Call}{function call.}
#' \item{Iterations}{number of iterations/}
#' \item{Tuning}{optimal tuning parameters.}
#' \item{Lambdas}{grid of lambda values for CV.}
#' \item{maxit.out}{maximum number of iterations for outer (blockwise) loop.}
#' \item{maxit.in}{maximum number of iterations for inner (lasso) loop.}
#' \item{Omega}{estimated penalized precision matrix.}
#' \item{Sigma}{estimated covariance matrix from the penalized precision matrix (inverse of Omega).}
#' \item{Path}{array containing the solution path. Solutions will be ordered by ascending lambda values.}
#' \item{MIN.error}{minimum average cross validation error (cv.crit) for optimal parameters.}
#' \item{AVG.error}{average cross validation error (cv.crit) across all folds.}
#' \item{CV.error}{cross validation errors (cv.crit).}
#' 
#' @references
#' \itemize{
#' \item Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. 'Sparse inverse covariance estimation with the graphical lasso.' \emph{Biostatistics} 9.3 (2008): 432-441.
#' \item Banerjee, Onureen, Ghauoui, Laurent El, and d'Aspremont, Alexandre. 2008. 'Model Selection through Sparse Maximum Likelihood Estimation for Multivariate Gaussian or Binary Data.' \emph{Journal of Machine Learning Research} 9: 485-516.
#' \item Tibshirani, Robert. 1996. 'Regression Shrinkage and Selection via the Lasso.' \emph{Journal of the Royal Statistical Society. Series B (Methodological)}. JSTOR: 267-288.
#' \item Meinshausen, Nicolai and Buhlmann, Peter. 2006. 'High-Dimensional Graphs and Variable Selection with the Lasso.' \emph{The Annals of Statistics}. JSTOR: 1436-1462.
#' \item Witten, Daniela M, Friedman, Jerome H, and Simon, Noah. 2011. 'New Insights and Faster computations for the Graphical Lasso.' \emph{Journal of Computation and Graphical Statistics}. Taylor and Francis: 892-900.
#' \item Tibshirani, Robert, Bien, Jacob, Friedman, Jerome, Hastie, Trevor, Simon, Noah, Jonathan, Taylor, and Tibshirani, Ryan J. 'Strong Rules for Discarding Predictors in Lasso-Type Problems.' \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}. Wiley Online Library 74 (2): 245-266.
#' \item Ghaoui, Laurent El, Viallon, Vivian, and Rabbani, Tarek. 2010. 'Safe Feature Elimination for the Lasso and Sparse Supervised Learning Problems.' \emph{arXiv preprint arXiv: 1009.4219}.
#' \item Osborne, Michael R, Presnell, Brett, and Turlach, Berwin A. 'On the Lasso and its Dual.' \emph{Journal of Computational and Graphical Statistics}. Taylor and Francis 9 (2): 319-337.
#' \item Rothman, Adam. 2017. 'STAT 8931 notes on an algorithm to compute the Lasso-penalized Gausssian likelihood precision matrix estimator.' 
#' }
#' 
#' @author Matt Galloway \email{gall0441@@umn.edu}
#' 
#' @seealso \code{\link{plot.GLASSO}}
#' 
#' @export
#' 
#' @examples
#' # generate data from a sparse matrix
#' # first compute covariance matrix
#' S = matrix(0.7, nrow = 5, ncol = 5)
#' for (i in 1:5){
#'  for (j in 1:5){
#'    S[i, j] = S[i, j]^abs(i - j)
#'  }
#'  }
#'
#' # generate 100 x 5 matrix with rows drawn from iid N_p
#' set.seed(123)
#' Z = matrix(rnorm(100*5), nrow = 100, ncol = 5)
#' out = eigen(S, symmetric = TRUE)
#' S.sqrt = out$vectors %*% diag(out$values^0.5)
#' S.sqrt = S.sqrt %*% t(out$vectors)
#' X = Z %*% S.sqrt
#'
#' # lasso penalty CV
#' GLASSO(X)

# we define the GLASSO precision matrix estimation
# function
GLASSO = function(X = NULL, S = NULL, nlam = 10, lam.min.ratio = 0.01, 
    lam = NULL, diagonal = FALSE, path = FALSE, crit.out = c("avg", 
        "max"), crit.in = c("loss", "avg", "max"), tol.out = 1e-04, 
    tol.in = 1e-04, maxit.out = 10000, maxit.in = 10000, 
    adjmaxit.out = NULL, K = 5, crit.cv = c("loglik", "AIC", 
        "BIC"), start = c("warm", "cold"), cores = 1, trace = c("progress", 
        "print", "none")) {
    
    # checks
    if (is.null(X) && is.null(S)) {
        stop("Must provide entry for X or S!")
    }
    if (!all(lam > 0)) {
        stop("lam must be positive!")
    }
    if (!(all(c(tol.out, tol.in, maxit.out, maxit.in, adjmaxit.out, 
        K, cores) > 0))) {
        stop("Entry must be positive!")
    }
    if (!(all(sapply(c(tol.out, tol.in, maxit.out, maxit.in, 
        adjmaxit.out, K, cores, nlam, lam.min.ratio), length) <= 
        1))) {
        stop("Entry must be single value!")
    }
    if (all(c(maxit.out, maxit.in, adjmaxit.out, K, cores)%%1 != 
        0)) {
        stop("Entry must be an integer!")
    }
    if (cores < 1) {
        stop("Number of cores must be positive!")
    }
    if (cores > 1 && path) {
        cat("Parallelization not possible when producing solution path. Setting cores = 1...\n\n")
        cores = 1
    }
    K = ifelse(path, 1, K)
    if (cores > K) {
        cat("Number of cores exceeds K... setting cores = K\n\n")
        cores = K
    }
    if (is.null(adjmaxit.out)) {
        adjmaxit.out = maxit.out
    }
    
    # match values
    crit.out = match.arg(crit.out)
    crit.in = match.arg(crit.in)
    crit.cv = match.arg(crit.cv)
    start = match.arg(start)
    trace = match.arg(trace)
    call = match.call()
    MIN.error = AVG.error = CV.error = NULL
    n = ifelse(is.null(X), nrow(S), nrow(X))
    
    # compute sample covariance matrix, if necessary
    if (is.null(S)) {
        S = (nrow(X) - 1)/nrow(X) * cov(X)
    }
    
    Sminus = S
    diag(Sminus) = 0
    
    # compute grid of lam values, if necessary
    if (is.null(lam)) {
        if (!((lam.min.ratio <= 1) && (lam.min.ratio > 
            0))) {
            cat("lam.min.ratio must be in (0, 1]... setting to 1e-2!\n\n")
            lam.min.ratio = 0.01
        }
        if (!((nlam > 0) && (nlam%%1 == 0))) {
            cat("nlam must be a positive integer... setting to 10!\n\n")
            nlam = 10
        }
        
        # calculate lam.max and lam.min
        lam.max = max(abs(Sminus))
        lam.min = lam.min.ratio * lam.max
        
        # calculate grid of lambda values
        lam = 10^seq(log10(lam.min), log10(lam.max), length = nlam)
        
    } else {
        
        # sort lambda values
        lam = sort(lam)
        
    }
    
    # perform cross validation, if necessary
    if ((length(lam) > 1) & (!is.null(X) || path)) {
        
        # run CV in parallel?
        if (cores > 1) {
            
            # execute CVP_GLASSO
            GLASSO = CVP_GLASSO(X = X, lam = lam, diagonal = diagonal, 
                crit.out = crit.out, crit.in = crit.in, 
                tol.out = tol.out, tol.in = tol.in, maxit.out = maxit.out, 
                maxit.in = maxit.in, adjmaxit.out = adjmaxit.out, 
                K = K, crit.cv = crit.cv, start = start, 
                cores = cores, trace = trace)
            MIN.error = GLASSO$min.error
            AVG.error = GLASSO$avg.error
            CV.error = GLASSO$cv.error
            
        } else {
            
            # execute CV_ADMMc
            if (is.null(X)) {
                X = matrix(0)
            }
            GLASSO = CV_GLASSOc(X = X, S = S, lam = lam, 
                diagonal = diagonal, path = path, crit_out = crit.out, 
                crit_in = crit.in, tol_out = tol.out, tol_in = tol.in, 
                maxit_out = maxit.out, maxit_in = maxit.in, 
                adjmaxit_out = adjmaxit.out, K = K, crit_cv = crit.cv, 
                start = start, trace = trace)
            MIN.error = GLASSO$min.error
            AVG.error = GLASSO$avg.error
            CV.error = GLASSO$cv.error
            Path = GLASSO$path
            
        }
        
        # print warning if lam on boundary
        if ((GLASSO$lam == lam[1]) && (length(lam) != 1) && 
            !path) {
            cat("\nOptimal tuning parameter on boundary... consider providing a smaller lam value or decreasing lam.min.ratio!")
        }
        
        # specify initial estimate for Sigma
        if (diagonal) {
            
            # simply force init to be positive definite final
            # diagonal elements will be increased by lam
            init = S + GLASSO$lam
            
        } else {
            
            # provide estimate that is pd and dual feasible
            alpha = min(c(GLASSO$lam/max(abs(Sminus)), 
                1))
            init = (1 - alpha) * S
            diag(init) = diag(S)
            
        }
        
        # compute final estimate at best tuning parameters
        GLASSO = GLASSOc(S = S, initSigma = init, initOmega = diag(ncol(S)), 
            lam = GLASSO$lam, crit_out = crit.out, crit_in = crit.in, 
            tol_out = tol.out, tol_in = tol.in, maxit_out = maxit.out, 
            maxit_in = maxit.in)
        
        
    } else {
        
        # execute ADMM_sigmac
        if (length(lam) > 1) {
            stop("Must set specify X, set path = TRUE, or provide single value for lam.")
        }
        
        # specify initial estimate for Sigma
        if (diagonal) {
            
            # simply force init to be positive definite final
            # diagonal elements will be increased by lam
            init = S + lam
            
        } else {
            
            # provide estimate that is pd and dual feasible
            alpha = min(c(lam/max(abs(Sminus)), 1))
            init = (1 - alpha) * S
            diag(init) = diag(S)
            
        }
        
        GLASSO = GLASSOc(S = S, initSigma = init, initOmega = diag(ncol(S)), 
            lam = lam, crit_out = crit.out, crit_in = crit.in, 
            tol_out = tol.out, tol_in = tol.in, maxit_out = maxit.out, 
            maxit_in = maxit.in)
        
    }
    
    
    # option to penalize diagonal
    if (diagonal) {
        C = 1
    } else {
        C = 1 - diag(ncol(S))
    }
    
    # compute penalized loglik
    loglik = (-n/2) * (sum(GLASSO$Omega * S) - determinant(GLASSO$Omega, 
        logarithm = TRUE)$modulus[1] + GLASSO$lam * sum(abs(C * 
        GLASSO$Omega)))
    
    
    # return values
    tuning = matrix(c(log10(GLASSO$lam), GLASSO$lam), ncol = 2)
    colnames(tuning) = c("log10(lam)", "lam")
    if (!path) {
        Path = NULL
    }
    
    returns = list(Call = call, Iterations = GLASSO$Iterations, 
        Tuning = tuning, Lambdas = lam, maxit.out = maxit.out, 
        maxit.in = maxit.in, Omega = GLASSO$Omega, Sigma = GLASSO$Sigma, 
        Path = Path, Loglik = loglik, MIN.error = MIN.error, 
        AVG.error = AVG.error, CV.error = CV.error)
    
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
        cat("\nMaximum iterations reached...!")
    }
    
    # print call
    cat("\nCall: ", paste(deparse(x$Call), sep = "\n", 
        collapse = "\n"), "\n", sep = "")
    
    # print iterations
    cat("\nIterations:\n")
    print.default(x$Iterations, quote = FALSE)
    
    # print optimal tuning parameters
    cat("\nTuning parameter:\n")
    print.default(round(x$Tuning, 3), print.gap = 2L, quote = FALSE)
    
    # print loglik
    cat("\nLog-likelihood: ", paste(round(x$Loglik, 5), 
        sep = "\n", collapse = "\n"), "\n", sep = "")
    
    # print Omega if dim <= 10
    if (nrow(x$Omega) <= 10) {
        cat("\nOmega:\n")
        print.default(round(x$Omega, 5))
    } else {
        cat("\n(...output suppressed due to large dimension!)\n")
    }
    
}



##-----------------------------------------------------------------------------------




#' @title Plot GLASSO object
#' @description Produces a plot for the cross validation errors, if available.
#' @param x class object GLASSO
#' @param type produce either 'heatmap' or 'line' graph
#' @param footnote option to print footnote of optimal values. Defaults to TRUE.
#' @param ... additional arguments.
#' @export
#' @examples
#' # generate data from a sparse matrix
#' # first compute covariance matrix
#' S = matrix(0.7, nrow = 5, ncol = 5)
#' for (i in 1:5){
#'  for (j in 1:5){
#'    S[i, j] = S[i, j]^abs(i - j)
#'  }
#'  }
#'
#' # generate 100 x 5 matrix with rows drawn from iid N_p(0, S)
#' set.seed(123)
#' Z = matrix(rnorm(100*5), nrow = 100, ncol = 5)
#' out = eigen(S, symmetric = TRUE)
#' S.sqrt = out$vectors %*% diag(out$values^0.5)
#' S.sqrt = S.sqrt %*% t(out$vectors)
#' X = Z %*% S.sqrt
#' 
#' # produce line graph for GLASSO
#' plot(GLASSO(X))
#' 
#' # produce CV heat map for GLASSO
#' plot(GLASSO(X), type = 'heatmap')

plot.GLASSO = function(x, type = c("line", "heatmap"), 
    footnote = TRUE, ...) {
    
    # check
    type = match.arg(type)
    Means = NULL
    if (is.null(x$CV.error)) {
        stop("No cross validation errors to plot!")
    }
    
    if (type == "line") {
        
        # gather values to plot
        cv = cbind(expand.grid(lam = x$Lambdas, alpha = 0), 
            Errors = as.data.frame.table(x$CV.error)$Freq)
        
        # produce line graph
        graph = ggplot(summarise(group_by(cv, lam), Means = mean(Errors)), 
            aes(log10(lam), Means)) + geom_jitter(width = 0.2, 
            color = "navy blue") + theme_minimal() + geom_line(color = "red") + 
            labs(title = "Cross-Validation Errors", y = "Error") + 
            geom_vline(xintercept = x$Tuning[1], linetype = "dotted")
        
    } else {
        
        # augment values for heat map (helps visually)
        lam = x$Lambdas
        cv = expand.grid(lam = lam, alpha = 0)
        Errors = 1/(c(x$AVG.error) + abs(min(x$AVG.error)) + 
            1)
        cv = cbind(cv, Errors)
        
        # design color palette
        bluetowhite <- c("#000E29", "white")
        
        # produce ggplot heat map
        graph = ggplot(cv, aes(alpha, log10(lam))) + geom_raster(aes(fill = Errors)) + 
            scale_fill_gradientn(colours = colorRampPalette(bluetowhite)(2), 
                guide = "none") + theme_minimal() + labs(title = "Heatmap of Cross-Validation Errors") + 
            theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
                axis.ticks.x = element_blank())
        
    }
    
    if (footnote) {
        
        # produce with footnote
        graph + labs(caption = paste("**Optimal: log10(lam) = ", 
            round(x$Tuning[1], 3), sep = ""))
        
    } else {
        
        # produce without footnote
        graph
        
    }
}
