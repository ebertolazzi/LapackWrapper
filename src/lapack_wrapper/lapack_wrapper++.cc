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

namespace lapack_wrapper {

  /*\
   *
   *   RELEASE 3.0, WGS COPYRIGHT 1997.
   *
   *   PURPOSE
   *
   *   To compute (optionally) a rank-revealing QR factorization of a
   *   real general M-by-N matrix  A,  which may be rank-deficient,
   *   and estimate its effective rank using incremental condition
   *   estimation.
   *
   *   The routine uses a QR factorization with column pivoting:
   *      A * P = Q * R,  where  R = [ R11 R12 ],
   *                                 [  0  R22 ]
   *   with R11 defined as the largest leading submatrix whose estimated
   *   condition number is less than 1/RCOND.  The order of R11, RANK,
   *   is the effective rank of A.
   *
   *   MB03OD  does not perform any scaling of the matrix A.
   *
   *   ARGUMENTS
   *
   *   Mode Parameters
   *
   *   JOBQR   CHARACTER*1
   *           = 'Q':  Perform a QR factorization with column pivoting;
   *           = 'N':  Do not perform the QR factorization (but assume
   *                   that it has been done outside).
   *
   *   Input/Output Parameters
   *
   *   M       (input) INTEGER
   *           The number of rows of the matrix A.  M >= 0.
   *
   *   N       (input) INTEGER
   *           The number of columns of the matrix A.  N >= 0.
   *
   *   A       (input/output) DOUBLE PRECISION array, dimension
   *           ( LDA, N )
   *           On entry with JOBQR = 'Q', the leading M by N part of this
   *           array must contain the given matrix A.
   *           On exit with JOBQR = 'Q', the leading min(M,N) by N upper
   *           triangular part of A contains the triangular factor R,
   *           and the elements below the diagonal, with the array TAU,
   *           represent the orthogonal matrix Q as a product of
   *           min(M,N) elementary reflectors.
   *           On entry and on exit with JOBQR = 'N', the leading
   *           min(M,N) by N upper triangular part of A contains the
   *           triangular factor R, as determined by the QR factorization
   *           with pivoting.  The elements below the diagonal of A are
   *           not referenced.
   *
   *   LDA     INTEGER
   *           The leading dimension of the array A.  LDA >= max(1,M).
   *
   *   RCOND   (input) DOUBLE PRECISION
   *           RCOND is used to determine the effective rank of A, which
   *           is defined as the order of the largest leading triangular
   *           submatrix R11 in the QR factorization with pivoting of A,
   *           whose estimated condition number is less than 1/RCOND.
   *           RCOND >= 0.
   *           NOTE that when SVLMAX > 0, the estimated rank could be
   *           less than that defined above (see SVLMAX).
   *
   *   TAU     (output) DOUBLE PRECISION array, dimension ( MIN( M, N ) )
   *           On exit with JOBQR = 'Q', the leading min(M,N) elements of
   *           TAU contain the scalar factors of the elementary
   *           reflectors.
   *           Array TAU is not referenced when JOBQR = 'N'.
   *
   *   RANK    (output) INTEGER
   *           The effective (estimated) rank of A, i.e. the order of
   *           the submatrix R11.
   *
   *   SVAL    (output) DOUBLE PRECISION array, dimension ( 3 )
   *           The estimates of some of the singular values of the
   *           triangular factor R:
   *           SVAL(1): largest singular value of R(1:RANK,1:RANK);
   *           SVAL(2): smallest singular value of R(1:RANK,1:RANK);
   *           SVAL(3): smallest singular value of R(1:RANK+1,1:RANK+1),
   *                    if RANK < MIN( M, N ), or of R(1:RANK,1:RANK),
   *                    otherwise.
   *           If the triangular factorization is a rank-revealing one
   *           (which will be the case if the leading columns were well-
   *           conditioned), then SVAL(1) will also be an estimate for
   *           the largest singular value of A, and SVAL(2) and SVAL(3)
   *           will be estimates for the RANK-th and (RANK+1)-st singular
   *           values of A, respectively.
   *           By examining these values, one can confirm that the rank
   *           is well defined with respect to the chosen value of RCOND.
   *           The ratio SVAL(1)/SVAL(2) is an estimate of the condition
   *           number of R(1:RANK,1:RANK).
   *
   *   Workspace
   *
   *   DWORK   DOUBLE PRECISION array, dimension ( LDWORK )
   *           where LDWORK = max( 1, 2*min( M, N ) )
   *
   *   Error Indicator
   *
   *   INFO    INTEGER
   *           = 0:  successful exit
   *           < 0:  if INFO = -i, the i-th argument had an illegal
   *                 value.
   *
   *   METHOD
   *
   *   The routine computes or uses a QR factorization with column
   *   pivoting of A,  A * P = Q * R,  with  R  defined above, and then
   *   finds the largest leading submatrix whose estimated condition
   *   number is less than 1/RCOND, taking the possible positive value of
   *   SVLMAX into account.  This is performed using the LAPACK
   *   incremental condition estimation scheme and a slightly modified
   *   rank decision test.
   *
   *   CONTRIBUTOR
   *
   *   V. Sima, Katholieke Univ. Leuven, Belgium, Nov. 1996.
   *
   *  ******************************************************************
  \*/

  template <typename T>
  integer
  rankEstimate(
    integer   M,
    integer   N,
    T         A[],
    integer   LDA,
    T         RCOND,
    T         SVAL[3]
  ) {

    integer MN = std::min(M, N);
    std::vector<T> Wmin( MN ), Wmax( MN );

    // Test the input scalar arguments.
    LAPACK_WRAPPER_ASSERT(
      M >= 0 && N >= 0,
      "rankEstimate, bad size matrix " << M << " x " << N
    );
    LAPACK_WRAPPER_ASSERT(
      LDA >= max_index(1,M),
      "rankEstimate, bad leading dimension ldA = " << LDA
    );
    LAPACK_WRAPPER_ASSERT(
      RCOND >= 0,
      "rankEstimate, bad condision number rcond = " << RCOND
    );

    // Quick return if possible
    SVAL[0] = 0;
    SVAL[1] = 0;
    SVAL[2] = 0;
    if ( MN == 0 ) return 0;

    // Determine RANK using incremental condition estimation
    integer RANK = 0;
    T SMAX = std::abs( A[0] );
    if ( SMAX > 0 ) {
      T SMIN   = SMAX;
      T SMINPR = SMIN;
      Wmin[0] = Wmax[0] = 1;
      while ( ++RANK < MN ) {
        T SMAXPR, S1, C1, S2, C2;
        T * pA0r = A + RANK * LDA;
        T & Arr  = pA0r[RANK];
        laic1( 2, RANK, &Wmin.front(), SMIN, pA0r, Arr, SMINPR, S1, C1 );
        laic1( 1, RANK, &Wmax.front(), SMAX, pA0r, Arr, SMAXPR, S2, C2 );

        if ( SMAXPR*RCOND > SMINPR ) break;

        for ( integer i=0; i < RANK; ++i )
          { Wmin[i] *= S1; Wmax[i] *= S2; }

        Wmin[RANK] = C1;
        Wmax[RANK] = C2;
        SMIN = SMINPR;
        SMAX = SMAXPR;
      }
      SVAL[0] = SMAX;
      SVAL[1] = SMIN;
      SVAL[2] = SMINPR;
    }
    return RANK;
  }

}

#include "code++/lu.cxx"
#include "code++/qr.cxx"
#include "code++/svd.cxx"
#include "code++/ls.cxx"
#include "code++/trid.cxx"
#include "code++/band.cxx"
#include "code++/qn.cxx"
#include "code++/eig.cxx"
#include "code++/block_trid.cxx"

namespace lapack_wrapper {

  /*\
  :|:   _     ____   ___   ____
  :|:  | |   / ___| / _ \ / ___|
  :|:  | |   \___ \| | | | |
  :|:  | |___ ___) | |_| | |___
  :|:  |_____|____/ \__\_\\____|
  \*/

  #if 0

  template <typename T>
  void
  LSQC<T>::factorize(
    integer         NR,
    integer         NC,
    valueType const M[],
    integer         ldM,
    bool const      row_select[] // row selected for minimization
  ) {
    integer NN = NR+NC;
    lu.allocate( NN, NN );
    lu.zero();
    NR1 = NR2 = 0;
    for ( integer i = 0; i < NR; ++i ) {
      if ( row_select[i] ) to_row1[NR1++] = i;
      else                 to_row2[NR2++] = i;
    }

    allocReals.allocate( size_t(NN) );
    rhs = allocReals( size_t(NN) );

    allocIntegers.allocate( size_t(NR) );
    to_row1 = allocIntegers( size_t(NR1) );
    to_row2 = allocIntegers( size_t(NR2) );

    /*\
    :|:  / 0   A^T  B^T \
    :|:  |              |
    :|:  | A   -I    0  |
    :|:  |              |
    :|:  \ B    0    0  /
    \*/
    // A
    for ( integer i = 0; i < NR1; ++i ) {
      integer irow = to_row1[i];
      integer NCi  = NC+i;
      lapack_wrapper::copy(
        NC,
        M + irow, ldM,
        lu.Apointer()+NCi, lu.numRow()
      );
      lapack_wrapper::copy(
        NC,
        M + irow, ldM,
        lu.Apointer()+NCi*lu.numRow(), 1
      );
      *(lu.Apointer()+(NCi+1)*lu.numRow()) = -1;
    }
    // B
    for ( integer i = 0; i < NR2; ++i ) {
      integer irow = to_row2[i];
      integer NCi  = NC+NR1+i;
      lapack_wrapper::copy(
        NC,
        M + irow, ldM,
        lu.Apointer()+NCi, lu.numRow()
      );
      lapack_wrapper::copy(
        NC,
        M + irow, ldM,
        lu.Apointer()+NCi*lu.numRow(), 1
      );
    }

    lu.factorize( "LSQC" );
  }

  template <typename T>
  void
  LSQC<T>::solve( valueType xb[] ) const {
    integer NR = lu.numRow();
    lapack_wrapper::zero( NR, rhs, 1 );
    // A^T b
    for ( integer i = 0; i < NR1; ++i ) {
      integer ii = to_row1[i];
      rhs[i] = lapack_wrapper::dot( NR, xb, 1, this->M+ii*this->ldM, 1);
    }
    // c
    for ( integer i = 0; i < NR2; ++i ) rhs[NR1+i] = xb[to_row2[i]];
  }

  template <typename T>
  void
  LSQC<T>::solve(
    integer   nrhs,
    valueType B[],
    integer   ldB
  ) const {

  }

  #endif

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

  template class LU<real>;
  template class LU<doublereal>;

  template class LUPQ<real>;
  template class LUPQ<doublereal>;

  template class QR<real>;
  template class QR<doublereal>;

  template class QRP<real>;
  template class QRP<doublereal>;

  template class SVD<real>;
  template class SVD<doublereal>;

  template class LSS<real>;
  template class LSS<doublereal>;

  template class LSY<real>;
  template class LSY<doublereal>;

  template class TridiagonalSPD<real>;
  template class TridiagonalSPD<doublereal>;

  template class TridiagonalLU<real>;
  template class TridiagonalLU<doublereal>;

  template class TridiagonalQR<real>;
  template class TridiagonalQR<doublereal>;

  template class BlockTridiagonalSymmetic<real>;
  template class BlockTridiagonalSymmetic<doublereal>;

  template class BandedLU<real>;
  template class BandedLU<doublereal>;

  template class BandedSPD<real>;
  template class BandedSPD<doublereal>;

#if 0
  template class LSQC<real>;
  template class LSQC<doublereal>;
#endif

  template class QN<real>;
  template class QN<doublereal>;

  template class BFGS<real>;
  template class BFGS<doublereal>;

  template class DFP<real>;
  template class DFP<doublereal>;

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

} // end namespace lapack_wrapper

///
/// eof: lapack_wrapper++.cc
///
