// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// vclCG_m
Eigen::VectorXd vclCG_m(const Eigen::SparseMatrix<double, Eigen::RowMajor>& eigen_sparsematrix, const Eigen::VectorXd& eigen_rhs, const int method, const int solver_iters, const double solver_tolerance);
RcppExport SEXP _spKrylov_vclCG_m(SEXP eigen_sparsematrixSEXP, SEXP eigen_rhsSEXP, SEXP methodSEXP, SEXP solver_itersSEXP, SEXP solver_toleranceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::SparseMatrix<double, Eigen::RowMajor>& >::type eigen_sparsematrix(eigen_sparsematrixSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type eigen_rhs(eigen_rhsSEXP);
    Rcpp::traits::input_parameter< const int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< const int >::type solver_iters(solver_itersSEXP);
    Rcpp::traits::input_parameter< const double >::type solver_tolerance(solver_toleranceSEXP);
    rcpp_result_gen = Rcpp::wrap(vclCG_m(eigen_sparsematrix, eigen_rhs, method, solver_iters, solver_tolerance));
    return rcpp_result_gen;
END_RCPP
}
// vcl_logdet
double vcl_logdet(const Eigen::SparseMatrix<double, Eigen::RowMajor>& eigen_sparsematrix, const Eigen::VectorXd& eigen_r, const int m);
RcppExport SEXP _spKrylov_vcl_logdet(SEXP eigen_sparsematrixSEXP, SEXP eigen_rSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::SparseMatrix<double, Eigen::RowMajor>& >::type eigen_sparsematrix(eigen_sparsematrixSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type eigen_r(eigen_rSEXP);
    Rcpp::traits::input_parameter< const int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(vcl_logdet(eigen_sparsematrix, eigen_r, m));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_spKrylov_vclCG_m", (DL_FUNC) &_spKrylov_vclCG_m, 5},
    {"_spKrylov_vcl_logdet", (DL_FUNC) &_spKrylov_vcl_logdet, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_spKrylov(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
