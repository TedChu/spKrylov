#' @export
krylov_neg_loglik <- function(theta = c(2, 0.2),
                              y,
                              X = NULL,
                              dist_mat,
                              cov.model = "exponential",
                              cov.taper = "wend1",
                              delta = 2,
                              ctrl = list(),
                              gls = FALSE,
                              nu = NULL) {
  n <- length(y)
  ctrl <- do.call("krylov.ctrl", ctrl)
  
  cov.fun <- switch(cov.model,
    "exponential" = cov.exp,
    "matern" = cov.mat,
    "spherical" = cov.sph
  )
  taper.fun <- switch(cov.taper,
    "spherical" = cov.sph,
    "wend1" = cov.wend1,
    "wend2" = cov.wend2
  )
  cg.method <- switch(ctrl$cg.precond,
    "no_precond" = 1,
    "ICHOL0" = 2,
    "ILUT" = 3,
    "Jacobi" = 4,
    "row scaling" = 5,
    "Chow-Patel ICHOL0" = 6,
    "ILU0" = 7,
    "Block-ILU0" = 8,
    "Block-ILUT" = 9
  )
  
  model.theta <- switch(cov.model,
    "exponential" = c(theta[1], 1 - theta[2], theta[2]),
    "matern" = c(theta[1], 1 - theta[2], nu, theta[2]),
    "spherical" = c(theta[1], 1 - theta[2], theta[2])
  )

  psi <- cov.fun(dist_mat, theta = model.theta) *
    taper.fun(dist_mat, theta = delta)
  # calculate log determinant
  if (abs(sum(psi) - sum(spam::diag(psi))) / n^2 < 1e-8) {
    logdet <- sum(log(spam::diag(psi)))
    psi <- as.dgRMatrix.spam(psi)
  } else {
    psi <- as.dgRMatrix.spam(psi)
    logdet <- krylov.logdet(psi, ctrl$logdet.m, ctrl$logdet.nr, rad.seed = ctrl$rad.seed)
  }
  # calculate matrix inverse
  if (!is.null(X)) {
    if (gls) {
      tmp <- vclCG_mrhs(psi, X, cg.method, ctrl$cg.max_iter, ctrl$cg.tol)
      beta <- solve(crossprod(X, tmp), crossprod(tmp, y))
      print(head(tmp))
    } else {
      beta <- solve(crossprod(X), crossprod(X, y))
    }
    y_tilde <- y - X %*% beta
    z <- vclCG_m(psi, y_tilde, cg.method, ctrl$cg.max_iter, ctrl$cg.tol)
    sigma2 <- crossprod(y_tilde, z) / n
  } else {
    z <- vclCG_m(psi, y, cg.method, ctrl$cg.max_iter, ctrl$cg.tol)
    sigma2 <- crossprod(y, z) / n
  }
  print(theta)
  loglik <- n / 2 * log(2 * pi) + n / 2 + 1 / 2 * logdet + n / 2 * log(sigma2)
  print(loglik)
  return(loglik)
}
