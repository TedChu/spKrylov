# spKrylov

Krylov Methods for Large Spatial Datasets

## Installation
This package can be installed directly from R by running the following code:
```{r}
library(devtools)
install_github("liujl93/spKrylov")
```

## Downloads:
Alternatively, you can download the binary file:
- Package Source: [spKrylov_1.0.tar.gz](https://liujl93.github.io/files/spKrylov_1.0.tar.gz)
- Windows binaries: [spKrylov_1.0.zip](https://liujl93.github.io/files/spKrylov_1.0.zip)
- OS X binaries: [spKrylov_1.0.tgz](https://liujl93.github.io/files/spKrylov_1.0.tgz)

```{r}
## MacOS
install.packages("/path/to/spKrylov_1.0.tgz", repos = NULL, type = .Platform$pkgType)
## Windows
install.packages("/path/to/spKrylov_1.0.zip", repos = NULL, type = "win.binary")
```

## Example

This package includes a data file `sim_n_4900rep_1delta_6`. It includes the following objects:
- train: N by 3 matrix, first two columns stores the coordinates and the last column stores the observations
- test: N_ho by 3 matrix, first two columns are coordinates and the last column are observations
- delta: only distances smaller than delta are recorded, used in nearest.dist()
- dist_mat: distance matrix of train data, only distances smaller than delta are recorded, stored in compressed form (class: spam)

To run this example, we will need the following `R` packages.
- spatstat: runifpoint
- RandomFields: RMexp, RMnugget, RFsimulate
- spam/spam64: nearest.dist
- Matrix: sparse Matrix manipulations
- doRNG/foreach/doParallel: parallel computing support

###### Install Packages
```
rm(list=ls(all=TRUE))
library(spKrylov)
library(foreach)
library(spam)
```    

###### Approximate Negative Log-likelihood
```
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
```

###### Find MLE
```
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
      cg.precond = "Jacobi",
      logdet.m = 30, logdet.nr = nr
    ),
    gls = FALSE
  )
)

fit.mle$result$theta
```
