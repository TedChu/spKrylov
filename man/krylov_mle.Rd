\name{krylov_mle}
\alias{krylov_mle}
\title{The function \code{krylov_mle} is used to fit spatial regression models.}
\usage{
krylov_mle(y, X = NULL, dist_mat, cov.model = "exponential",
  cov.taper = "wend1", theta.init = c(2, 0.2), theta.lower = c(1e-04,
  0), theta.upper = c(2, 1), delta = 2, ctrl = list(), gls = FALSE,
  nu = NULL)
}
\arguments{
\item{y}{the response variable.}

\item{X}{the model matrix.}

\item{dist_mat}{the distance matrix.}

\item{cov.model}{a quoted keyword that specifies the covariance function used to model the spatial dependence structure among the observations. Supported keywords are: \code{"exponential"}, \code{"matern"} and \code{"spherical"}. }

\item{cov.taper}{cov.taper a quoted keyword that specifies the tapering function. Supported keywords are: \code{"spherical"}, \code{"wend1"} and \code{"wend2"}. }

\item{theta.init}{numeric vector, initial values for the covariance parameters to be optimized.}

\item{theta.lower}{vector of lower limit of covariance parameters.}

\item{theta.upper}{vector of upper limit of covariance parameters.}

\item{delta}{tapering threshold parameter.}

\item{ctrl}{A list of control parameters. See 'Details'.}

\item{gls}{If \code{gls==TURE}, use generalized least square method. Otherwise, use ordinary least square method.}

\item{nu}{the smooth parameter for the Matern covariance function.}
}
\details{ The \code{ctrl} argument is a list that can supply any of the following components:
\describe{
\item{cg.tol}{convergence tolerance for the conjugate gradient algorithm.}
\item{cg.max_iter}{maximum number of iterations for the conjugate gradient algorithm.}
\item{cg.precond}{a quoted keyword that specifies preconditioner. Supported keywords are:
\code{"no_precond"},\code{"ICHOL0"},\code{"ILUT"},\code{"Jacobi"},\code{"row scaling"},
\code{"ILU0"},\code{"Block-ILU0"} and \code{"Block-ILUT"}.}
\item{logdet.m}{the order of the Gaussian quadrature rule.}
\item{logdet.nr}{ the number of Monte Carlo iterations.}
\item{rad.seed}{ random seed.}
}
}
\value{
A list with components:
\describe{
\item{result}{ A list consists of the regression coefficient beta, the covariance parameter theta, and the solution to the linear system \eqn{z=\Gamma^{-1}(y-X\beta)} at the \eqn{k}th iteration.} 
\item{conv}{ An integer code. 0 indicates successful convergence.}
\item{neg.loglik}{ Negative log-likelihood.}
}
}
\description{
The function \code{krylov_mle} is used to fit spatial regression models.
}
\examples{
library(spKrylov)
library(foreach)
library(spam)

set.seed(2019)
n <- 70^2 # sample size
delta <- 6 # 	only distances smaller than delta are recorded.
rep <- 1
data("sim_n_4900rep_1delta_6")
nr <- 1

dist_mat <- nearest.dist(train[, 1:2],
  miles = FALSE,
  upper = NULL,
  delta = delta
)

system.time(
  fit.mle <- krylov_mle(
    y = train[, 3],
    dist_mat = dist_mat,
    cov.model = "exponential",
    cov.taper = "wend1",
    theta.init = c(3, 0.15),
    theta.lower = c(0.0001, 0.0001),
    theta.upper = c(max(dist_mat), 1),
    delta = delta,
    ctrl = list(
      cg.tol = 1e-6, cg.max_iter = 1000,
      cg.precond = "no_precond",
      logdet.m = 30, logdet.nr = nr
    ),
    gls = FALSE
  )
)

fit.mle$result$theta
}
