// [[Rcpp::depends(RcppEigen)]]

#define VIENNACL_HAVE_EIGEN 1
#include <RcppEigen.h>
#include <Eigen/Sparse>

//#include <viennacl/ocl/backend.hpp>
#include <viennacl/matrix.hpp>
#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/compressed_matrix.hpp"
//#include "viennacl/context.hpp"
#include "viennacl/linalg/cg.hpp"
#include "viennacl/linalg/ilu.hpp"
#include "viennacl/linalg/ichol.hpp"
#include "logdet.hpp"
#include <iostream>
#include <vector>

#ifdef VIENNACL_WITH_OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
double vcl_logdet(const Eigen::SparseMatrix<double, Eigen::RowMajor> &eigen_sparsematrix,
                      const Eigen::VectorXd &eigen_r, const int m){

  //viennacl::context ctx;
  int size = eigen_r.size();

  // Set up some ViennaCL objects:
  viennacl::compressed_matrix<double> vcl_compressed_matrix(size, size);
  viennacl::vector<double> vcl_r(size);

  // copy from eigen to ViennaCL
  viennacl::copy(eigen_sparsematrix,  vcl_compressed_matrix);
  viennacl::copy(eigen_r,  vcl_r);

  double logdet0 = logdet(vcl_compressed_matrix, vcl_r, m);

  return logdet0;
}
