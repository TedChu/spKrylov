# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

vclCG_m <- function(eigen_sparsematrix, eigen_rhs, method, solver_iters, solver_tolerance) {
  .Call(`_spKrylov_vclCG_m`, eigen_sparsematrix, eigen_rhs, method, solver_iters, solver_tolerance)
}

vcl_logdet <- function(eigen_sparsematrix, eigen_r, m) {
  .Call(`_spKrylov_vcl_logdet`, eigen_sparsematrix, eigen_r, m)
}
