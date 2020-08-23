#' @export
krylov_mle <- function(y,
                       X = NULL,
                       dist_mat,
                       cov.model = "exponential",
                       cov.taper = "wend1",
                       theta.init = c(2, 0.2),
                       theta.lower = c(0.0001, 0),
                       theta.upper = c(2, 1),
                       delta = 2,
                       ctrl = list(),
                       gls = FALSE,
                       nu = NULL) {
  n <- length(y)
  ctrl <- do.call("krylov.ctrl", ctrl)
  fit <- nlminb(theta.init,
    krylov_neg_loglik,
    lower = theta.lower,
    upper = theta.upper,
    y = y,
    X = X,
    dist_mat = dist_mat,
    cov.model = cov.model,
    cov.taper = cov.taper,
    delta = delta,
    ctrl = ctrl,
    gls = gls,
    nu = nu
  )
  opt <- fit$par
  obj <- fit$objective
  conv <- ifelse(fit$conv == 0, T, F)
  result <- krylov_estimate(
    theta = opt, y = y, X = X,
    dist_mat = dist_mat,
    cov.model = cov.model,
    cov.taper = cov.taper,
    delta = delta, ctrl = ctrl,
    gls = gls, nu = nu
  )
  return(list(result = result, conv = conv, neg.loglik = obj))
}
