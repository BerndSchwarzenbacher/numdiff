#ifndef FILE_NGBLA
#define FILE_NGBLA

#include "ngstd.hpp"

/// namespace for basic linear algebra
namespace ngbla
{
  using namespace std;
  using namespace ngstd;

  using ngstd::CArray;

  typedef std::complex<double> Complex;
  inline double fabs (Complex v) { return std::abs (v); }

  inline bool IsComplex(double v) { return false; }
  inline bool IsComplex(Complex v) { return true; }
}


#include "expr.hpp"
#include "vector.hpp"
#include "matrix.hpp"

#include "exception.hpp"
#include "localheap.hpp"
#include "array.hpp"
#include "autodiff.hpp"


//#include "cholesky.hpp"
//#include "symmetricmatrix.hpp"
//#include "bandmatrix.hpp"
//#include "tensor.hpp"

//#include "cuda_bla.hpp"

namespace ngbla
{

  /// Computes eigenvalues and vectors of the symmetric matrix mat.
  extern NGS_DLL_HEADER void CalcEigenSystem (const FlatMatrix<double> & mat,
			       FlatVector<double> & lami,
			       FlatMatrix<double> & eigenvecs);



  /// Computes the Schur Complement.
  extern NGS_DLL_HEADER void CalcSchurComplement (const FlatMatrix<double> a, 
				   FlatMatrix<double> s,
				   const BitArray & used,
				   LocalHeap & lh);
}

#include "ng_lapack.hpp"

#endif
