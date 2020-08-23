#' @export
chol_neg_loglik <- function(theta,
                            y,
                            X = NULL,
                            dist_mat,
                            cov.model = "exponential",
                            nu = NULL) {
  n <- length(y)
  cov.fun <- switch(cov.model,
    exponential = cov.exp,
    matern = cov.mat,
    spherical = cov.sph
  )
  model.theta <- switch(cov.model,
    "exponential" = c(theta[1], 1 - theta[2], theta[2]),
    "matern" = c(theta[1], 1 - theta[2], nu, theta[2]),
    "spherical" = c(theta[1], 1 - theta[2], theta[2])
  )
  
  psi <- cov.fun(dist_mat, theta = model.theta)
  U <- chol(psi)
  U.inv <- backsolve(U, diag(1, nrow = n))
  Gaminv <- tcrossprod(U.inv)
  if (!is.null(X)) {
    tmp <- crossprod(X, Gaminv)
    beta <- solve(tmp %*% X, tmp %*% y)
    resid <- y - X %*% beta
    z <- Gaminv %*% resid
    sigma2 <- drop(crossprod(resid, z)) / n
  } else {
    beta <- NULL
    z <- Gaminv %*% y
    sigma2 <- drop(crossprod(y, z)) / n
  }

  neg.loglik <- n / 2 * log(2 * pi) + n / 2 + sum(log(diag(U))) + n / 2 * log(sigma2)
  return(neg.loglik)
}

chol_est <- function(theta,
                     y,
                     X = NULL,
                     dist_mat,
                     cov.model = "exponential",
                     nu = NULL) {
  n <- length(y)
  cov.fun <- switch(cov.model,
    exponential = cov.exp,
    matern = cov.mat,
    spherical = cov.sph
  )
  model.theta <- switch(cov.model,
    "exponential" = c(theta[1], 1 - theta[2], theta[2]),
    "matern" = c(theta[1], 1 - theta[2], nu, theta[2]),
    "spherical" = c(theta[1], 1 - theta[2], theta[2])
  )
  
  psi <- cov.fun(dist_mat, theta = model.theta)
  U <- chol(psi)
  U.inv <- backsolve(U, diag(1, nrow = n))
  Gaminv <- tcrossprod(U.inv)
  if (!is.null(X)) {
    tmp <- crossprod(X, Gaminv)
    beta <- solve(tmp %*% X, tmp %*% y)
    resid <- y - X %*% beta
    z <- Gaminv %*% resid
    sigma2 <- drop(crossprod(resid, z)) / n
  } else {
    beta <- NULL
    z <- Gaminv %*% y
    sigma2 <- drop(crossprod(y, z)) / n
  }

  theta.hat <- c(sigma2, theta)
  return(list(beta = beta, theta = theta.hat, z = z))
}

chol_mle <- function(y,
                     X = NULL,
                     dist_mat,
                     cov.model = "exponential",
                     theta.init = c(2, 0.2),
                     theta.lower = c(1e-04, 0),
                     theta.upper = c(2, 1),
                     nu = NULL) {
  n <- length(y)
  fit <- nlminb(theta.init,
    chol_neg_loglik,
    lower = theta.lower,
    upper = theta.upper,
    y = y,
    X = X,
    dist_mat = dist_mat,
    cov.model = cov.model,
    nu = nu
  )
  opt <- fit$par
  obj <- fit$objective
  conv <- ifelse(fit$conv == 0, T, F)
  result <- chol_est(
    theta = opt,
    y = y,
    X = X,
    dist_mat = dist_mat,
    cov.model = cov.model,
    nu = nu
  )
  return(list(result = result, conv = conv, neg.loglik = obj))
}
