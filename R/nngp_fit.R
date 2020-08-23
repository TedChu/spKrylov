#' @export
#-------------------#
#        NNGP
#-------------------#
library(spNNGP)
nngp_score <- function(theta,
                       y,
                       X = NULL,
                       coords,
                       cov.model = "exponential",
                       score.rule = "rmspe",
                       nngp.seed = 2019,
                       n.neighbors = 10,
                       nu = NULL,
                       n.omp.threads = 1) {
  # theta = c(sigma2,c,r)
  # Save the old random seed and use the new one, if any
  if (!is.na(nngp.seed)) {
    if (exists(".Random.seed")) {
      saved.seed <- .Random.seed
    }
    else {
      saved.seed <- NA
    }
    set.seed(nngp.seed)
  }
  
  theta.alpha <- switch(cov.model,
    "exponential" = matrix(c(1 / theta[3], theta[2] / (1 - theta[2])), 1, 2),
    "matern" = matrix(c(1 / theta[3], theta[2] / (1 - theta[2]), nu), 1, 3),
    "spherical" = matrix(c(1 / theta[3], theta[2] / (1 - theta[2])), 1, 2)
  )

  colnames(theta.alpha) <- switch(cov.model,
    "exponential" = c("phi", "alpha"),
    "matern" = c("phi", "alpha", "nu"),
    "spherical" = c("phi", "alpha")
  )
  m.c <- spConjNNGP(y ~ X - 1,
    coords = coords, n.neighbors = n.neighbors,
    k.fold = 5, score.rule = score.rule,
    n.omp.threads = n.omp.threads,
    theta.alpha = theta.alpha,
    sigma.sq.IG = c(2, theta[1] * theta[2]),
    cov.model = cov.model
  )$k.fold.scores
  score <- m.c[colnames(m.c) == score.rule]

  # Restore the old random seed, if any
  if (!is.na(nngp.seed) && !is.na(saved.seed)) {
    .Random.seed <- saved.seed
  }

  return(score)
}


# Test
# system.time(nngp_crps <- nngp_score(theta,
#                                     y=train[,3],
#                                     X=X_train,
#                                     coords=train[,1:2]))

nngp_est <- function(theta,
                     y,
                     X = NULL,
                     coords,
                     cov.model = "exponential",
                     score.rule = "rmspe",
                     coords.ho = NULL,
                     X.ho = NULL,
                     nngp.seed = 2019,
                     n.neighbors = 10,
                     nu = NULL,
                     n.omp.threads = 1) {
  # theta = c(sigma2,c,r)
  # Save the old random seed and use the new one, if any
  if (!is.na(nngp.seed)) {
    if (exists(".Random.seed")) {
      saved.seed <- .Random.seed
    }
    else {
      saved.seed <- NA
    }
    set.seed(nngp.seed)
  }

  theta.alpha <- switch(cov.model,
    "exponential" = matrix(c(1 / theta[3], theta[2] / (1 - theta[2])), 1, 2),
    "matern" = matrix(c(1 / theta[3], theta[2] / (1 - theta[2]), nu), 1, 3),
    "spherical" = matrix(c(1 / theta[3], theta[2] / (1 - theta[2])), 1, 2)
  )

  colnames(theta.alpha) <- switch(cov.model,
    "exponential" = c("phi", "alpha"),
    "matern" = c("phi", "alpha", "nu"),
    "spherical" = c("phi", "alpha")
  )

  if (!is.null(X)) {
    result <- spConjNNGP(y ~ X - 1,
      coords = coords, n.neighbors = n.neighbors,
      X.0 = X.ho, coords.0 = coords.ho,
      k.fold = 5, score.rule = score.rule,
      n.omp.threads = n.omp.threads,
      theta.alpha = theta.alpha,
      sigma.sq.IG = c(2, theta[1] * theta[2]),
      cov.model = cov.model
    )
  } else {
    if (!is.null(coords.ho)) {
      X.ho <- matrix(rep(1, nrow(coords.ho)), nrow(coords.ho), 1)
    }
    result <- spConjNNGP(y ~ 1,
      coords = coords, n.neighbors = n.neighbors,
      X.0 = X.ho, coords.0 = coords.ho,
      k.fold = 5, score.rule = score.rule,
      n.omp.threads = n.omp.threads,
      theta.alpha = theta.alpha,
      sigma.sq.IG = c(2, theta[1] * theta[2]),
      cov.model = cov.model
    )
  }
  # Restore the old random seed, if any
  if (!is.na(nngp.seed) && !is.na(saved.seed)) {
    .Random.seed <- saved.seed
  }

  return(result)
}

# Test
# system.time(nngp_crps <- nngp_est(theta,
#                                   y=train[,3],
#                                   X=X_train,
#                                   coords=train[,1:2],
#                                   cov.model="exponential",
#                                   coords.ho=test[,1:2],
#                                   X.ho=X_test))

nngp_fit <- function(y,
                     X = NULL,
                     coords,
                     cov.model = "exponential",
                     score.rule = "rmspe",
                     theta.init = c(9, 0.2, 2),
                     theta.lower = c(9, 1e-4, 1e-4),
                     theta.upper = c(9, 1, max(dist_mat0)),
                     coords.ho = NULL,
                     X.ho = NULL,
                     nngp.seed = 2019,
                     n.neighbors = 10,
                     nu = nu,
                     n.omp.threads = 1) {
  n <- length(y)
  fit <- nlminb(theta.init,
    nngp_score,
    lower = theta.lower,
    upper = theta.upper,
    y = y,
    X = X,
    coords = coords,
    cov.model = cov.model,
    score.rule = score.rule,
    nngp.seed = nngp.seed,
    n.neighbors = n.neighbors,
    nu = nu,
    n.omp.threads = n.omp.threads
  )
  opt <- fit$par
  obj <- fit$objective
  conv <- ifelse(fit$conv == 0, T, F)
  tmp <- nngp_est(
    theta = opt,
    y = y,
    X = X,
    coords = coords,
    cov.model = cov.model,
    score.rule = score.rule,
    coords.ho = coords.ho,
    X.ho = X.ho,
    nngp.seed = nngp.seed,
    n.neighbors = n.neighbors,
    nu = nu,
    n.omp.threads = n.omp.threads
  )
  result <- list()
  result$beta <- tmp[[1]]
  attr(result$beta, "dimnames") <- NULL
  result$beta.sd <- sqrt(diag(tmp[[2]]))
  attr(result$beta.sd, "names") <- NULL
  nugget <- 1 - 1 / (tmp[[6]][2] + 1)
  result$theta.hat <- c(tmp[[7]] / (1 - nugget), nugget, 1 / tmp[[6]][1])
  result$pred <- tmp[[3]]
  return(list(result = result, conv = conv, neg.loglik = obj))
}
