#' @export
krylov.predict <- function(object,
                           coords = NULL,
                           coords.ho = NULL,
                           X.ho = NULL,
                           method = "krylov",
                           cov.model = "exponential",
                           cov.taper = "wend1",
                           delta = 2,
                           dist.mat = NULL,
                           nu = NULL) {
  beta <- object$beta
  theta <- object$theta
  z <- object$z
  cov.fun <- switch(cov.model,
    exponential = cov.exp,
    matern = cov.mat,
    spherical = cov.sph
  )
  if (is.null(dist.mat) & is.null(coords.ho) & is.null(coords)) {
    stop("error: either dist.mat or coords.ho and coords must be specified")
  }
  if (is.null(dist.mat) & !is.null(coords.ho) & !is.null(coords)) {
    if (method %in% c("krylov", "spam")) {
      dist.mat <- nearest.dist(coords.ho,
        coords,
        miles = FALSE,
        delta = delta
      )
    } else {
      dist.mat <- rdist(coords.ho, coords)
    }
  }
  if (method %in% c("krylov", "spam")) {
    taper.fun <- switch(cov.taper,
      wend1 = cov.wend1,
      wend2 = cov.wend2
    )
    model.theta <- switch(cov.model,
      "exponential" = c(theta[2], 1 - theta[3], theta[3]),
      "matern" = c(theta[2], 1 - theta[3], nu, theta[3]),
      "spherical" = c(theta[2], 1 - theta[3], theta[3])
    )
    psi0 <- cov.fun(dist.mat, theta = model.theta) *
      taper.fun(dist.mat, theta = delta)
  } else {
    model.theta <- switch(cov.model,
      "exponential" = c(theta[2], 1 - theta[3], theta[3]),
      "matern" = c(theta[2], 1 - theta[3], nu, theta[3]),
      "spherical" = c(theta[2], 1 - theta[3], theta[3])
    )

    psi0 <- cov.fun(dist.mat, theta = model.theta)
  }
  if (!is.null(beta)) {
    pred <- X.ho %*% beta + psi0 %*% z
  } else {
    pred <- psi0 %*% z
  }
  return(pred)
}
