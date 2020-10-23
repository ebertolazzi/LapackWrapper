/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#ifdef __GNUC__ 
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wpadded"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wweak-template-vtables"
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wpadded"
#endif

#include <iomanip>
#include <vector>

#include "lapack_wrapper++.hh"

#include "code/wrapper.cxx"
#include "code/sparse.cxx"
#include "code/rank_estimate.cxx"

#include "code++/band.cxx"
#include "code++/block_trid.cxx"
#include "code++/eig.cxx"
#include "code++/ls.cxx"
#include "code++/lsc.cxx"
#include "code++/lu.cxx"
#include "code++/pinv.cxx"
#include "code++/qn.cxx"
#include "code++/qr.cxx"
#include "code++/svd.cxx"
#include "code++/trid.cxx"

namespace lapack_wrapper {

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template integer rankEstimate( integer   M,
                                 integer   N,
                                 real      A[],
                                 integer   LDA,
                                 real      RCOND,
                                 real      SVAL[3] );

  template integer rankEstimate( integer    M,
                                 integer    N,
                                 doublereal A[],
                                 integer    LDA,
                                 doublereal RCOND,
                                 doublereal SVAL[3] );

  template class LU_no_alloc<real>;
  template class LU_no_alloc<doublereal>;
  template class LU<real>;
  template class LU<doublereal>;

  template class LUPQ_no_alloc<real>;
  template class LUPQ_no_alloc<doublereal>;
  template class LUPQ<real>;
  template class LUPQ<doublereal>;

  template class QR_no_alloc<real>;
  template class QR_no_alloc<doublereal>;
  template class QR<real>;
  template class QR<doublereal>;

  template class QRP_no_alloc<real>;
  template class QRP_no_alloc<doublereal>;
  template class QRP<real>;
  template class QRP<doublereal>;

  template class SVD_no_alloc<real>;
  template class SVD_no_alloc<doublereal>;
  template class SVD<real>;
  template class SVD<doublereal>;

  template class BandedLU<real>;
  template class BandedLU<doublereal>;

  template class BandedSPD<real>;
  template class BandedSPD<doublereal>;

  template class BlockTridiagonalSymmetic<real>;
  template class BlockTridiagonalSymmetic<doublereal>;

  template class Eigenvalues<real>;
  template class Eigenvalues<doublereal>;

  template class Eigenvectors<real>;
  template class Eigenvectors<doublereal>;

  template class GeneralizedEigenvalues<real>;
  template class GeneralizedEigenvalues<doublereal>;

  template class GeneralizedEigenvectors<real>;
  template class GeneralizedEigenvectors<doublereal>;

  template class GeneralizedSVD<real>;
  template class GeneralizedSVD<doublereal>;

  template class LSS_no_alloc<real>;
  template class LSS_no_alloc<doublereal>;
  template class LSS<real>;
  template class LSS<doublereal>;

  template class LSY_no_alloc<real>;
  template class LSY_no_alloc<doublereal>;
  template class LSY<real>;
  template class LSY<doublereal>;

  template class LSC<real>;
  template class LSC<doublereal>;

  template class PINV_no_alloc<real>;
  template class PINV_no_alloc<doublereal>;
  template class PINV<real>;
  template class PINV<doublereal>;

  template class QN<real>;
  template class QN<doublereal>;

  template class BFGS<real>;
  template class BFGS<doublereal>;

  template class DFP<real>;
  template class DFP<doublereal>;

  template class TridiagonalSPD<real>;
  template class TridiagonalSPD<doublereal>;

  template class TridiagonalLU<real>;
  template class TridiagonalLU<doublereal>;

  template class TridiagonalQR<real>;
  template class TridiagonalQR<doublereal>;

} // end namespace lapack_wrapper

///
/// eof: lapack_wrapper++.cc
///
