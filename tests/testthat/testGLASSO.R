

# generate data from a sparse matrix
# first compute covariance matrix
S = matrix(0.7, nrow = 5, ncol = 5)
for (i in 1:5){
  for (j in 1:5){
    S[i, j] = S[i, j]^abs(i - j)
  }
}

# generate 100 x 5 matrix with rows drawn from iid N_p(0, S)
Z = matrix(rnorm(100*5), nrow = 100, ncol = 5)
out = eigen(S, symmetric = TRUE)
S.sqrt = out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)
X = Z %*% S.sqrt

# calculate sample covariance
(nrow(X) - 1)/nrow(X)*cov(X)

# elastic-net type penalty (use CV for optimal lambda and alpha)
expect_error(GLASSO(X), NA)
expect_warning(GLASSO(X), NA)

# lasso penalty (lam = 0.1)
expect_error(GLASSO(X, lam = 0.1), NA)
expect_warning(GLASSO(X, lam = 0.1), NA)

expect_error(GLASSO(S = S, lam = 0.1), NA)
expect_warning(GLASSO(S = S, lam = 0.1), NA)

# parallel CV
expect_error(GLASSO(X, cores = 2), NA)
expect_warning(GLASSO(X, cores = 2), NA)

# adjmaxit
expect_error(GLASSO(X, adjmaxit.out = 2), NA)
expect_warning(GLASSO(X, adjmaxit.out = 2), NA)

# parallel adjmaxit
expect_error(GLASSO(X, adjmaxit.out = 2, cores = 2), NA)
expect_warning(GLASSO(X, adjmaxit.out = 2, cores = 2), NA)

# path
expect_error(GLASSO(X, path = TRUE), NA)
expect_warning(GLASSO(X, path = TRUE), NA)

expect_error(GLASSO(S = S, path = TRUE), NA)
expect_warning(GLASSO(S = S, path = TRUE), NA)
