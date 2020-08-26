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
#include "viennacl/linalg/bicgstab.hpp"
#include "viennacl/linalg/gmres.hpp"
#include "viennacl/linalg/mixed_precision_cg.hpp"

#include "viennacl/linalg/ilu.hpp"
#include "viennacl/linalg/ichol.hpp"
#include "viennacl/linalg/jacobi_precond.hpp"
#include "viennacl/linalg/row_scaling.hpp"

#include <iostream>
#include <vector>

#ifdef VIENNACL_WITH_OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
Eigen::VectorXd vclCG_m(const Eigen::SparseMatrix<double, Eigen::RowMajor> &eigen_sparsematrix,
                      const Eigen::VectorXd &eigen_rhs,
                      const int method,
                      const int solver_iters,
                      const double solver_tolerance){

  //viennacl::context ctx;
  int size = eigen_rhs.size();

  // maximum iterations and tolerance
  //unsigned int solver_iters = 1000;
  //double solver_tolerance = 1e-6;

  // Set up some ViennaCL objects:
  viennacl::compressed_matrix<double> vcl_compressed_matrix(size, size);
  viennacl::vector<double> vcl_rhs(size);
  viennacl::vector<double> vcl_result(size);

  // copy from eigen to ViennaCL
  viennacl::copy(eigen_sparsematrix,  vcl_compressed_matrix);
  viennacl::copy(eigen_rhs,  vcl_rhs);

  ///////////////////////////////////////////////////////////////////////////////
  //////////////////////              CG solver                //////////////////
  ///////////////////////////////////////////////////////////////////////////////

  viennacl::linalg::cg_tag cg_solver(solver_tolerance, solver_iters);

  switch (method) {
    case 1:{
      //std::cout << "------- CG solver (no preconditioner) via ViennaCL, compressed_matrix ----------" << std::endl;
      vcl_result = viennacl::linalg::solve(vcl_compressed_matrix, vcl_rhs, cg_solver, viennacl::linalg::no_precond());
      break;
    }
    case 2:{
      //std::cout << "------- CG solver (ICHOL0 preconditioner) via ViennaCL, compressed_matrix ----------" << std::endl;
      viennacl::linalg::ichol0_precond< viennacl::compressed_matrix<double> > vcl_ichol0(vcl_compressed_matrix, viennacl::linalg::ichol0_tag());
      vcl_result = viennacl::linalg::solve(vcl_compressed_matrix, vcl_rhs, cg_solver, vcl_ichol0);
      break;
    }
    case 3:{
      viennacl::linalg::ilut_precond< viennacl::compressed_matrix<double> > vcl_ilut(vcl_compressed_matrix, viennacl::linalg::ilut_tag());
      //std::cout << "------- CG solver (ILUT preconditioner) via ViennaCL, compressed_matrix ----------" << std::endl;
      vcl_result = viennacl::linalg::solve(vcl_compressed_matrix, vcl_rhs, cg_solver, vcl_ilut);
      break;
    }
    case 4:{
      viennacl::linalg::jacobi_precond< viennacl::compressed_matrix<double> > vcl_jacobi_csr(vcl_compressed_matrix, viennacl::linalg::jacobi_tag());
      //std::cout << "------- CG solver (Jacobi preconditioner) via ViennaCL, compressed_matrix ----------" << std::endl;
      vcl_result = viennacl::linalg::solve(vcl_compressed_matrix, vcl_rhs, cg_solver, vcl_jacobi_csr);
      break;
    }
    case 5:{
      viennacl::linalg::row_scaling< viennacl::compressed_matrix<double> > vcl_row_scaling_csr(vcl_compressed_matrix, viennacl::linalg::row_scaling_tag(1));
      //std::cout << "------- CG solver (row scaling preconditioner) via ViennaCL, compressed_matrix ----------" << std::endl;
      vcl_result = viennacl::linalg::solve(vcl_compressed_matrix, vcl_rhs, cg_solver, vcl_row_scaling_csr);
      break;
    }
    case 6:{
      viennacl::linalg::chow_patel_icc_precond< viennacl::compressed_matrix<double> > vcl_chow_patel_icc(vcl_compressed_matrix, viennacl::linalg::chow_patel_tag());
      //std::cout << "------- CG solver (Chow-Patel ICHOL0 preconditioner) via ViennaCL, compressed_matrix ----------" << std::endl;
      vcl_result = viennacl::linalg::solve(vcl_compressed_matrix, vcl_rhs, cg_solver, vcl_chow_patel_icc);
      break;
    }
    case 7:{
      viennacl::linalg::ilu0_precond< viennacl::compressed_matrix<double> > vcl_ilu0(vcl_compressed_matrix, viennacl::linalg::ilu0_tag());
      //std::cout << "------- CG solver (ILU0 preconditioner) via ViennaCL, compressed_matrix ----------" << std::endl;
      vcl_result = viennacl::linalg::solve(vcl_compressed_matrix, vcl_rhs, cg_solver, vcl_ilu0);
      break;
    }
    case 8:{
      viennacl::linalg::block_ilu_precond< viennacl::compressed_matrix<double>,
                                           viennacl::linalg::ilu0_tag>          vcl_block_ilu0(vcl_compressed_matrix, viennacl::linalg::ilu0_tag());
      //std::cout << "------- CG solver (Block-ILU0 preconditioner) via ViennaCL, compressed_matrix ----------" << std::endl;
      vcl_result = viennacl::linalg::solve(vcl_compressed_matrix, vcl_rhs, cg_solver, vcl_block_ilu0);
      break;
    }
    case 9:{
      viennacl::linalg::block_ilu_precond< viennacl::compressed_matrix<double>,
                                           viennacl::linalg::ilut_tag>          vcl_block_ilut(vcl_compressed_matrix, viennacl::linalg::ilut_tag());
      //std::cout << "------- CG solver (Block-ILUT preconditioner) via ViennaCL, compressed_matrix ----------" << std::endl;
      vcl_result = viennacl::linalg::solve(vcl_compressed_matrix, vcl_rhs, cg_solver, vcl_block_ilut);
      break;
    }
    default:{
      std::cout << "Invalid response" << std::endl;
    }
  }
  viennacl::vector<double> residual(vcl_rhs);
  residual -= viennacl::linalg::prod(vcl_compressed_matrix, vcl_result);
  //std::cout << "Relative residual: " << viennacl::linalg::norm_2(residual) / viennacl::linalg::norm_2(vcl_rhs) << std::endl;
  //std::cout << "Estimated rel. residual: " << cg_solver.error() << std::endl;
  //std::cout << "Iterations: " << cg_solver.iters() << std::endl;

  Eigen::VectorXd eigen_result(size);
  viennacl::copy(vcl_result,  eigen_result);
  return eigen_result;
}
