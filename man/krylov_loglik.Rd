\name{krylov_loglik}
\alias{krylov_loglik}
\title{The function \code{krylov_loglik} is used to calculate the approximate negative log-likelihood function by Krylov methods.}
\usage{
krylov_loglik(theta = c(2, 0.2), y, X = NULL, dist_mat,
  cov.model = "exponential", cov.taper = "wend1", delta = 2,
  ctrl = list(), gls = FALSE, nu = NULL)
}
\arguments{
\item{theta}{the vector of covariance parameters.}

\item{y}{response variable.}

\item{X}{the model matrix.}

\item{dist_mat}{the distance matrix.}

\item{cov.model}{a quoted keyword that specifies the covariance function used to model the spatial dependence structure among the observations. Supported keywords are: \code{"exponential"}, \code{"matern"} and \code{"spherical"}. }

\item{cov.taper}{cov.taper a quoted keyword that specifies the tapering function. Supported keywords are: \code{"spherical"}, \code{"wend1"} and \code{"wend2"}. }

\item{delta}{tapering threshold parameter.}

\item{ctrl}{A list of control parameters. See 'Details'.}

\item{gls}{If \code{gls==TURE}, use generalized least square method. Otherwise, use ordinary least square method.}

\item{nu}{the smooth parameter for the Matern covariance function.}
}
\details{The \code{ctrl} argument is a list that can supply any of the following components:
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
the negative log-likelihood.
}
\description{
The function \code{krylov_loglik} is used to calculate the approximate negative log-likelihood function by Krylov methods.
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
theta <- c(2, 0.2)

dist_mat <- nearest.dist(train[, 1:2],
  miles = FALSE,
  upper = NULL,
  delta = delta
)
#########################################################
### Approximate Negative Log-likelihood by Krylov Method
#########################################################
system.time(
  lik2 <- krylov_neg_loglik(
    theta = theta,
    y = train[, 3],
    dist_mat = dist_mat,
    cov.model = "exponential",
    cov.taper = "wend1",
    delta = delta,
    ctrl = list(
      cg.tol = 1e-6, cg.max_iter = 1000,
      cg.precond = "Jacobi",
      logdet.m = 50, logdet.nr = 1
    ),
    gls = FALSE
  )
)

lik2
}
