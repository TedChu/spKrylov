/**
*   @brief Implementation of the Lanczos algorithm
*
*   @param A            The system matrix
*   @param r            Random start vector
*   @param krylov_dim   Size of krylov-space
*   @param tag          The Lanczos tag holding tolerances, etc.
*   @return             Returns the approximated log determinant of A
*/

#include <cmath>
#include <vector>
#include <numeric>      // std::inner_product
#include <valarray>     // std::valarray, std::log(valarray)
#include "viennacl/vector.hpp"
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/io/matrix_market.hpp"
#include "viennacl/tools/random.hpp"

#include <viennacl/linalg/prod.hpp>
#include <viennacl/linalg/inner_prod.hpp>
#include <viennacl/linalg/norm_2.hpp>
#include <viennacl/linalg/sum.hpp>
#include <viennacl/linalg/scalar_operations.hpp>
#include <viennacl/linalg/vector_operations.hpp>
#include <viennacl/linalg/matrix_operations.hpp>
//#include <viennacl/ocl/kernel.hpp>
#include "viennacl/linalg/qr-method.hpp"
#include "viennacl/linalg/qr-method-common.hpp"

using namespace viennacl::linalg;

template< typename MatrixT, typename NumericT>
NumericT logdet(MatrixT const& A, viennacl::vector<NumericT> & r,
                int krylov_dim)
{
  std::vector<NumericT> alphas, betas;
  viennacl::vector<NumericT> Aq(r.size());
  viennacl::matrix<NumericT, viennacl::column_major> Q(r.size(), krylov_dim + 1);  // Krylov basis (each Krylov vector is one column)

  NumericT norm_r = norm_2(r);
  NumericT beta = norm_r;
  r /= norm_r;

  // first Krylov vector:
  viennacl::vector_base<NumericT> q0(Q.handle(), Q.size1(), 0, 1);
  q0 = r;

  //
  // Step 1: Run Lanczos' method to obtain tridiagonal matrix
  //
  for (unsigned int i = 0; i < krylov_dim; i++)
  {
    betas.push_back(beta);
    // last available vector from Krylov basis:
    viennacl::vector_base<NumericT> q_i(Q.handle(), Q.size1(), i * Q.internal_size1(), 1);

    // Lanczos algorithm:
    // - Compute A * q:
    Aq = viennacl::linalg::prod(A, q_i);

    // - Form Aq <- Aq - <Aq, q_i> * q_i - beta * q_{i-1}, where beta is ||q_i|| before normalization in previous iteration
    NumericT alpha = viennacl::linalg::inner_prod(Aq, q_i);
    Aq -= alpha * q_i;

    if (i > 0)
    {
      viennacl::vector_base<NumericT> q_iminus1(Q.handle(), Q.size1(), (i-1) * Q.internal_size1(), 1);
      Aq -= beta * q_iminus1;
    }

    // normalize Aq and add to Krylov basis at column i+1 in Q:
    beta = viennacl::linalg::norm_2(Aq);
    viennacl::vector_base<NumericT> q_iplus1(Q.handle(), Q.size1(), (i+1) * Q.internal_size1(), 1);
    q_iplus1 = Aq / beta;

    alphas.push_back(alpha);
  }

  //
  // Step 2: Compute eigenvalues and eigenvectors of tridiagonal matrix obtained during Lanczos iterations:
  //
  //std::vector<NumericT> eigenvalues = bisect(alphas, betas);
  viennacl::matrix<NumericT> QQ = viennacl::identity_matrix<NumericT>(alphas.size());
  // fill d and e with data here
  viennacl::linalg::tql2(QQ, alphas, betas);
  //
  // Step 3: get sum(eigenvectors[1,]^2*log(eigenvalues))
  //
  viennacl::vector<NumericT> qq1 = viennacl::row(QQ, 0);
  NumericT logdet = 0;
  for(unsigned int i = 0; i < alphas.size(); i++){
    logdet+=qq1[i]*qq1[i]*log(alphas[i]);
  }

  return logdet;
}


//#include <iomanip> // setprecision
//int main()
//{
//  // If you GPU does not support double precision, use `float` instead of `double`:
//  typedef double     ScalarType;
//  ScalarType logdet0;
//  /**
//  *  Create the sparse matrix and read data from a Matrix-Market file:
//  **/
//  std::vector< std::map<unsigned int, ScalarType> > host_A;
//  std::cout << "Reading Matrix file..." << std::endl;
//  if (!viennacl::io::read_matrix_market_file(host_A, "/Users/lolofter/Documents/Projects/krylov_ref/ViennaCL-1.7.1/examples/testdata/mat65k.mtx"))
//  {
//    std::cout << "Error reading Matrix file" << std::endl;
//    return EXIT_FAILURE;
//  }
//  std::cout << "End Reading Matrix file..." << std::endl;
//  viennacl::compressed_matrix<ScalarType> A;
//  viennacl::copy(host_A, A);
//  std::cout << "End copying Matrix file..." << std::endl;
//
//  viennacl::vector<ScalarType> r(A.size1());
//  viennacl::tools::uniform_random_numbers<ScalarType> randomNumber;
//  for (unsigned int i = 0; i < A.size1(); ++i)
//    r[i] = (randomNumber() > 0.5 ? 1 : 0);
//  std::cout << "End assigning vector..." << std::endl;
//  unsigned int krylov_dim = 20;
//
//    /**
//    *  Run the Lanczos method for computing log deterimant.
//    **/
//  std::cout << "Running Lanczos algorithm. This might take a while..." << std::endl;
//  logdet0 = logdet(A, r, krylov_dim);
//    /**
//    *  Print the computed logdet and exit:
//    **/
//  std::cout << "logdet" << std::setprecision(7) << logdet0;
//    // test approximated eigenvector by computing A*v:
//
//  return EXIT_SUCCESS;
//}
