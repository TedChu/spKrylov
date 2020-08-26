#' @export
krylov.ctrl <- function(cg.tol = 1e-6, cg.max_iter = 1000,
                        cg.precond = "ichol",
                        logdet.m = 30, logdet.nr = 1,
                        rad.seed = 2019) {
  if (!is.numeric(cg.tol) || cg.tol <= 0) {
    stop("tolerance for conjugate gradient solver must be > 0")
  }
  if (!is.numeric(cg.max_iter) || cg.max_iter <= 0) {
    stop("maximum number of iterations for conjugate gradient must be > 0")
  }
  if (!is.numeric(logdet.m) || logdet.m <= 0) {
    stop("The size of the krylov space for lanczos algorithm must be > 0")
  }
  if (!is.numeric(logdet.nr) || logdet.nr <= 0) {
    stop("Number of Monte Carlo replications must be > 0")
  }
  ctrl <- list(
    cg.tol = cg.tol, cg.max_iter = cg.max_iter,
    cg.precond = cg.precond,
    logdet.m = logdet.m, logdet.nr = logdet.nr,
    rad.seed = rad.seed
  )
  return(ctrl)
}

krylov.rademacher <- function(n, nr, rad.seed = 2019) {

  # Save the old random seed and use the new one, if any
  if (!is.na(rad.seed)) {
    if (exists(".Random.seed")) {
      saved.seed <- .Random.seed
    }
    else {
      saved.seed <- NA
    }
    set.seed(rad.seed)
  }

  logdet.rhs <- matrix(sign(rnorm(n * nr)), n, nr)

  # Restore the old random seed, if any
  if (!is.na(rad.seed) && !is.na(saved.seed)) {
    .Random.seed <- saved.seed
  }

  return(logdet.rhs)
}

krylov.logdet <- function(A, m, nr, rad.seed = 2019) {
  n <- nrow(A)
  r <- krylov.rademacher(n, nr, rad.seed = rad.seed)
  g <- foreach(i = 1:nr, .combine = c) %do% {
    vcl_logdet(A, r[, i], m)
  }
  return(sum(g) * n / nr)
}

vclCG_mrhs <- function(eigen_sparsematrix, eigen_rhs, method, solver_iters, solver_tolerance) {
  if (!is.matrix(eigen_rhs)) {
    result <- vclCG_m(eigen_sparsematrix, eigen_rhs, method, solver_iters, solver_tolerance)
  } else {
    result <- foreach(i = 1:ncol(eigen_rhs), .combine = cbind) %dopar% {
      vclCG_m(eigen_sparsematrix, eigen_rhs[, i], method, solver_iters, solver_tolerance)
    }
  }
}

krylov_estimate <- function(theta = c(2, 0.2),
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
  psi <- as.dgRMatrix.spam(psi)
  if (!is.null(X)) {
    if (gls) {
      tmp <- vclCG_mrhs(psi, X, cg.method, ctrl$cg.max_iter, ctrl$cg.tol)
      beta <- solve(crossprod(X, tmp), crossprod(tmp, y))
    } else {
      beta <- solve(crossprod(X), crossprod(X, y))
    }
    y_tilde <- y - X %*% beta
    z <- vclCG_m(psi, y_tilde, cg.method, ctrl$cg.max_iter, ctrl$cg.tol)
    sigma2 <- crossprod(y_tilde, z) / n
  } else {
    beta <- NULL
    z <- vclCG_m(psi, y, cg.method, ctrl$cg.max_iter, ctrl$cg.tol)
    sigma2 <- crossprod(y, z) / n
  }

  theta.hat <- c(sigma2, theta)
  return(list(beta = beta, theta = theta.hat, z = z))
}
