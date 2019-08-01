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

  /*\
   |   _    _   _
   |  | |  | | | |
   |  | |  | | | |
   |  | |__| |_| |
   |  |_____\___/
   |
  \*/

  template <typename T>
  LU<T>::LU()
  : Factorization<T>()
  , allocReals("allocReals")
  , allocIntegers("allocIntegers")
  {}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  LU<T>::~LU() {
    allocReals.free();
    allocIntegers.free();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU<T>::allocate( integer NR, integer NC ) {
    if ( nRow != NR || nCol != NC ) {
      nRow = NR;
      nCol = NC;
      allocReals.allocate( size_t(nRow*nCol+2*(nRow+nCol)) );
      allocIntegers.allocate( size_t(2*nRow) );
      Amat    = allocReals( size_t(nRow*nCol) );
      Work    = allocReals( size_t(2*(nRow+nCol)) );
      i_pivot = allocIntegers( size_t(nRow) );
      Iwork   = allocIntegers( size_t(nRow) );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU<T>::factorize( char const who[] ) {
    integer info = getrf( nRow, nCol, Amat, nRow, i_pivot );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "LU::factorize[" << who << "] getrf INFO = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU<T>::factorize(
    char const      who[],
    integer         NR,
    integer         NC,
    valueType const A[],
    integer LDA
  ) {
    allocate( NR, NC );
    integer info = gecopy( nRow, nCol, A, LDA, Amat, nRow );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "LU::factorize[" << who << "] gecopy INFO = " << info
    );
    factorize( who );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU<T>::check_ls( char const who[] ) const {
    LAPACK_WRAPPER_ASSERT(
      nRow == nCol,
      "LU<T>::" << who << ", rectangular matrix " <<
      nRow << " x " << nCol
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU<T>::solve( valueType xb[] ) const {
    check_ls("solve");
    integer info = getrs(
      NO_TRANSPOSE, nRow, 1, Amat, nRow, i_pivot, xb, nRow
    );
    LAPACK_WRAPPER_ASSERT( info == 0, "LU::solve getrs INFO = " << info );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU<T>::t_solve( valueType xb[] ) const {
    check_ls("t_solve");
    integer info = getrs(
      TRANSPOSE, nRow, 1, Amat, nRow, i_pivot, xb, nRow
    );
    LAPACK_WRAPPER_ASSERT( info == 0, "LU::t_solve getrs INFO = " << info );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU<T>::solve( integer nrhs, valueType B[], integer ldB ) const {
    check_ls("solve");
    integer info = getrs(
      NO_TRANSPOSE, nRow, nrhs, Amat, nRow, i_pivot, B, ldB
    );
    LAPACK_WRAPPER_ASSERT( info == 0, "LU::solve getrs INFO = " << info );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LU<T>::t_solve( integer nrhs, valueType B[], integer ldB ) const {
    check_ls("t_solve");
    integer info = getrs(
      TRANSPOSE, nRow, nrhs, Amat, nRow, i_pivot, B, ldB
    );
    LAPACK_WRAPPER_ASSERT( info >= 0, "LU::t_solve getrs INFO = " << info );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  typename LU<T>::valueType
  LU<T>::cond1( valueType norm1 ) const {
    valueType rcond;
    integer info = gecon1(
      nRow, Amat, nRow, norm1, rcond, Work, Iwork
    );
    LAPACK_WRAPPER_ASSERT( info == 0, "LU::cond1, gecon1 return info = " << info );
    return rcond;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  typename LU<T>::valueType
  LU<T>::condInf( valueType normInf ) const {
    valueType rcond;
    integer info = geconInf(
      nRow, Amat, nRow, normInf, rcond, Work, Iwork
    );
    LAPACK_WRAPPER_ASSERT( info == 0, "LU::condInf, geconInf return info = " << info );
    return rcond;
  }

  /*\
   |   _    _   _ ____   ___
   |  | |  | | | |  _ \ / _ \
   |  | |  | | | | |_) | | | |
   |  | |__| |_| |  __/| |_| |
   |  |_____\___/|_|    \__\_\
   |
  \*/

  template <typename T>
  LUPQ<T>::LUPQ()
  : Factorization<T>()
  , allocReals("LUPQ-allocReals")
  , allocIntegers("LUPQ-allocIntegers")
  {}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  LUPQ<T>::~LUPQ() {
    allocReals.free();
    allocIntegers.free();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LUPQ<T>::allocate( integer NR, integer NC ) {
    LAPACK_WRAPPER_ASSERT(
      NR == NC,
      "LUPQ<T>::allocate, cannot allocate rectangular matrix " <<
      NR << " x " << NC
    );
    if ( nRow != NR || nCol != NC ) {
      nRow = NR;
      nCol = NC;
      allocReals.allocate( size_t(nRow*nCol) );
      Amat = allocReals( size_t(nRow*nCol) );
      allocIntegers.allocate( size_t(2*nRow) );
      ipiv = allocIntegers( size_t(nRow) );
      jpiv = allocIntegers( size_t(nRow) );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LUPQ<T>::check_ls( char const who[] ) const {
    LAPACK_WRAPPER_ASSERT(
      nRow == nCol,
      "LUPQ<T>::" << who << ", rectangular matrix " <<
      nRow << " x " << nCol
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LUPQ<T>::factorize( char const who[] ) {
    integer info = getc2( nRow, Amat, nRow, ipiv, jpiv );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "LUPQ::factorize[" << who << "] getrf INFO = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LUPQ<T>::factorize(
    char const      who[],
    integer         NR,
    integer         NC,
    valueType const A[],
    integer         LDA
  ) {
    LAPACK_WRAPPER_ASSERT(
      NR == NC,
      "LUPQ<T>::factorize[" << who << 
      "], cannot factorize rectangular matrix " <<
      NR << " x " << NC
    );
    allocate( NR, NC );
    integer info = gecopy( nRow, nCol, A, LDA, Amat, nRow );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "LUPQ::factorize[" << who << "] gecopy INFO = " << info
    );
    factorize( who );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LUPQ<T>::solve( valueType xb[] ) const {
    // Apply permutations IPIV to RHS
    swaps( 1, xb, nRow, 0, nRow-2, ipiv, 1 );

    // Solve for L part
    trsv( LOWER, NO_TRANSPOSE, UNIT, nRow, Amat, nRow, xb, 1 );

    // Solve for U part
    trsv( UPPER, NO_TRANSPOSE, NON_UNIT, nRow, Amat, nRow, xb, 1 );

    // Apply permutations JPIV to the solution (RHS)
    swaps( 1, xb, nRow, 0, nRow-2, jpiv, -1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LUPQ<T>::solve( integer nrhs, valueType B[], integer ldB ) const {
    // Apply permutations IPIV to RHS
    swaps( nrhs, B, ldB, 0, nRow-2, ipiv, 1 );

    // Solve for L part
    trsm( LEFT, LOWER, NO_TRANSPOSE, UNIT, nRow, nrhs, 1.0, Amat, nRow, B, ldB );

    // Solve for U part
    trsm( LEFT, UPPER, NO_TRANSPOSE, NON_UNIT, nRow, nrhs, 1.0, Amat, nRow, B, ldB );

    // Apply permutations JPIV to the solution (RHS)
    swaps( nrhs, B, ldB, 0, nRow-2, jpiv, -1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LUPQ<T>::t_solve( valueType xb[] ) const {
    // Apply permutations JPIV to the solution (RHS)
    swaps( 1, xb, nRow, 0, nRow-2, jpiv, 1 );

    // Solve for U part
    trsv( UPPER, TRANSPOSE, NON_UNIT, nRow, Amat, nRow, xb, 1 );

    // Solve for L part
    trsv( LOWER, TRANSPOSE, UNIT, nRow, Amat, nRow, xb, 1 );

    // Apply permutations IPIV to RHS
    swaps( 1, xb, nRow, 0, nRow-2, ipiv, -1 );

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LUPQ<T>::t_solve( integer nrhs, valueType B[], integer ldB ) const {

    // Apply permutations JPIV to the solution (RHS)
    swaps( nrhs, B, ldB, 0, nRow-2, jpiv, 1 );

    // Solve for U part
    trsm( LEFT, UPPER, TRANSPOSE, NON_UNIT, nRow, nrhs, 1.0, Amat, nRow, B, ldB );

    // Solve for L part
    trsm( LEFT, LOWER, TRANSPOSE, UNIT, nRow, nrhs, 1.0, Amat, nRow, B, ldB );

    // Apply permutations IPIV to RHS
    swaps( nrhs, B, ldB, 0, nRow-2, ipiv, -1 );
  }

  /*\
   |    ___  ____
   |   / _ \|  _ \
   |  | | | | |_) |
   |  | |_| |  _ <
   |   \__\_\_| \_\
   |
  \*/

  template <typename T>
  void
  QR<T>::setMaxNrhs( integer mnrhs ) {
    LAPACK_WRAPPER_ASSERT(
      mnrhs > 0,
      "lapack_wrapper::QR::setMaxNrhs, maxNrhs = " << mnrhs
    );
    maxNrhs = mnrhs;
  }

  template <typename T>
  void
  QR<T>::allocate( integer NR, integer NC, integer Lwrk ) {
    nRow       = NR;
    nCol       = NC;
    nReflector = std::min(nRow,nCol);
    Lwork      = Lwrk;
    allocReals.allocate( size_t(nRow*nCol+Lwork+nReflector) );
    Amat = allocReals( size_t(nRow*nCol) );
    Work = allocReals( size_t(Lwork) );
    Tau  = allocReals( size_t(nReflector) );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR<T>::allocate( integer NR, integer NC ) {
    if ( nRow != NR || nCol != NC || Lwork < maxNrhs ) {
      valueType tmp; // get optimal allocation
      integer info = geqrf( NR, NC, nullptr, NR, nullptr, &tmp, -1 );
      LAPACK_WRAPPER_ASSERT(
        info == 0,
        "QR::allocate call lapack_wrapper::geqrf return info = " << info
      );
      integer L = std::max( integer(tmp), maxNrhs );
      if ( L < NR ) L = NR;
      if ( L < NC ) L = NC;
      allocate( NR, NC, L );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR<T>::applyQ(
    SideMultiply  SIDE,
    Transposition TRANS,
    integer       nRefl,
    integer       NR,
    integer       NC,
    valueType     C[],
    integer       ldC
  ) const {
    LAPACK_WRAPPER_ASSERT(
      (SIDE == lapack_wrapper::LEFT  && NR == nRow) ||
      (SIDE == lapack_wrapper::RIGHT && NC == nRow),
      "QR::applyQ NR = " << NR << " NC = " << NC << " nRow = " << nRow
    );
    integer info = ormqr(
      SIDE, TRANS,
      NR, NC,
      nRefl,  // numero riflettori usati nel prodotto Q
      Amat, nRow /*ldA*/,
      Tau,
      C, ldC,
      Work, Lwork
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "QR::applyQ call lapack_wrapper::ormqr return info = " << info <<
      " Lwork = " << Lwork
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR<T>::getR( valueType R[], integer ldR ) const {
    integer minRC = std::min( nRow, nCol );
    gezero( minRC, minRC, R, ldR );
    for ( integer i = 0; i < minRC; ++i )
      for ( integer j = i; j < minRC; ++j )
        R[i+j*ldR] = Amat[ i+j*nRow];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR<T>::solve( valueType xb[] ) const {
    LAPACK_WRAPPER_ASSERT(
      nRow == nCol,
      "in QR::solve, factored matrix must be square"
    );
    Qt_mul(xb);
    invR_mul(xb);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR<T>::t_solve( valueType xb[] ) const {
    LAPACK_WRAPPER_ASSERT(
      nRow == nCol,
      "in QR::solve_t, factored matrix must be square"
    );
    invRt_mul(xb);
    Q_mul(xb);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR<T>::solve( integer nrhs, valueType XB[], integer ldXB ) const {
    LAPACK_WRAPPER_ASSERT(
      nRow == nCol,
      "in QR::solve, factored matrix must be square"
    );
    Qt_mul( nRow, nrhs, XB, ldXB );
    invR_mul( nRow, nrhs, XB, ldXB );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QR<T>::t_solve( integer nrhs, valueType XB[], integer ldXB ) const {
    LAPACK_WRAPPER_ASSERT(
      nRow == nCol,
      "in QR::solve_t, factored matrix must be square"
    );
    invRt_mul( nRow, nrhs, XB, ldXB );
    Q_mul( nRow, nrhs, XB, ldXB );
  }

  /*\
   |    ___  ____  ____
   |   / _ \|  _ \|  _ \
   |  | | | | |_) | |_) |
   |  | |_| |  _ <|  __/
   |   \__\_\_| \_\_|
   |
  \*/

  template <typename T>
  void
  QRP<T>::allocate( integer NR, integer NC ) {
    if ( nRow != NR || nCol != NC || Lwork < maxNrhs ) {
      valueType tmp; // get optimal allocation
      integer info = geqp3( NR, NC, nullptr, NR, nullptr, nullptr, &tmp, -1 );
      LAPACK_WRAPPER_ASSERT(
        info == 0,
        "QRP::allocate call lapack_wrapper::geqp3 return info = " << info
      );
      integer L = std::max( integer(tmp), maxNrhs );
      if ( L < NR ) L = NR;
      if ( L < NC ) L = NC;
      QR<T>::allocate( NR, NC, L );
    }
    allocIntegers.allocate(size_t(NC));
    JPVT = allocIntegers(size_t(NC));
  }

  template <typename T>
  void
  QRP<T>::permute( valueType x[] ) const {
    // applico permutazione
    for ( integer i = 0; i < nCol; ++i ) Work[JPVT[i]-1] = x[i];
    copy( nCol, Work, 1, x, 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QRP<T>::inv_permute( valueType x[] ) const {
    // applico permutazione
    for ( integer i = 0; i < nCol; ++i ) Work[i] = x[JPVT[i]-1];
    copy( nCol, Work, 1, x, 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QRP<T>::solve( valueType xb[] ) const {
    LAPACK_WRAPPER_ASSERT(
      nRow == nCol,
      "in QRP::solve, factored matrix must be square"
    );
    Qt_mul(xb);
    invR_mul(xb);
    permute(xb); // da aggiungere!
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QRP<T>::t_solve( valueType xb[] ) const {
    LAPACK_WRAPPER_ASSERT(
      nRow == nCol,
      "in QRP::solve_t, factored matrix must be square"
    );
    inv_permute(xb); // da aggiungere!
    invRt_mul(xb);
    Q_mul(xb);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QRP<T>::solve( integer nrhs, valueType XB[], integer ldXB ) const {
    LAPACK_WRAPPER_ASSERT(
      nRow == nCol,
      "in QRP::solve, factored matrix must be square"
    );
    Qt_mul( nRow, nrhs, XB, ldXB );
    invR_mul( nRow, nrhs, XB, ldXB );
    permute_rows( nRow, nrhs, XB, ldXB ); // da aggiungere!
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QRP<T>::t_solve( integer nrhs, valueType XB[], integer ldXB ) const {
    LAPACK_WRAPPER_ASSERT(
      nRow == nCol,
      "in QRP::solve_t, factored matrix must be square"
    );
    inv_permute_rows( nRow, nrhs, XB, ldXB ); // da aggiungere!
    invRt_mul( nRow, nrhs, XB, ldXB );
    Q_mul( nRow, nrhs, XB, ldXB );
  }

  /*\
   |   ______     ______
   |  / ___\ \   / /  _ \
   |  \___ \\ \ / /| | | |
   |   ___) |\ V / | |_| |
   |  |____/  \_/  |____/
   |
  \*/

  template <typename T>
  void
  SVD<T>::allocate( integer NR, integer NC ) {

    if ( nRow != NR || nCol != NC ) {
      nRow  = NR;
      nCol  = NC;
      minRC = std::min(NR,NC);
      valueType tmp;
      integer info = gesvd(
        REDUCED, REDUCED,
        NR, NC,
        nullptr, NR,
        nullptr,
        nullptr, NR,
        nullptr, minRC,
        &tmp, -1
      );
      LAPACK_WRAPPER_ASSERT(
        info == 0, "SVD::allocate, in gesvd info = " << info
      );
      Lwork = integer(tmp);
      info = gesdd(
        REDUCED,
        NR, NC,
        nullptr, NR,
        nullptr,
        nullptr, NR,
        nullptr, minRC,
        &tmp, -1, nullptr
      );
      if ( integer(tmp) > Lwork ) Lwork = integer(tmp);
      allocReals.allocate( size_t(nRow*nCol+minRC*(nRow+nCol+1)+Lwork) );
      Amat = allocReals( size_t(nRow*nCol));
      Svec  = allocReals( size_t(minRC) );
      Umat  = allocReals( size_t(minRC*nRow) );
      VTmat = allocReals( size_t(minRC*nCol) );
      Work  = allocReals( size_t(Lwork) );
      allocIntegers.allocate( size_t(8*minRC) );
      IWork = allocIntegers( size_t(8*minRC) );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  SVD<T>::factorize( char const who[] ) {
    integer info;
    switch ( svd_used ) {
    case USE_GESVD:
      info = gesvd(
        REDUCED,
        REDUCED,
        nRow, nCol, Amat, nRow,
        Svec,
        Umat, nRow,
        VTmat, minRC,
        Work, Lwork
      );
      LAPACK_WRAPPER_ASSERT(
        info == 0,
        "SVD::factorize[" << who << 
        "] call lapack_wrapper::gesvd return info = " << info
      );
      break;
    case USE_GESDD:
      info = gesdd(
        REDUCED,
        nRow, nCol, Amat, nRow,
        Svec,
        Umat, nRow,
        VTmat, minRC,
        Work, Lwork, IWork
      );
      LAPACK_WRAPPER_ASSERT(
        info == 0,
        "SVD::factorize[" << who << 
        "] call lapack_wrapper::gesdd return info = " << info
      );
      break;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  SVD<T>::solve( valueType xb[] ) const {
    // A = U*S*VT
    // U*S*VT*x=b --> VT^T S^+ U^T b
    // U  nRow x minRC
    // VT minRC x nCol
    valueType smin = rcond*Svec[0];
    Ut_mul( 1.0, xb, 1, 0.0, Work, 1 );
    for ( integer i = 0; i < minRC; ++i ) Work[i] /= std::max(Svec[i],smin);
    V_mul( 1.0, Work, 1, 0.0, xb, 1 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  SVD<T>::t_solve( valueType xb[] ) const {
    // A = U*S*VT
    // U*S*VT*x=b --> VT^T S^+ U^T b
    // U  nRow x minRC
    // VT minRC x nCol
    valueType smin = rcond*Svec[0];
    Vt_mul( 1.0, xb, 1, 0.0, Work, 1 );
    for ( integer i = 0; i < minRC; ++i ) Work[i] /= std::max(Svec[i],smin);
    U_mul( 1.0, Work, 1, 0.0, xb, 1 );
  }

  /*\
   |   _     ____ ____
   |  | |   / ___/ ___|
   |  | |   \___ \___ \
   |  | |___ ___) |__) |
   |  |_____|____/____/
   |
  \*/

  template <typename T>
  void
  LSS<T>::setMaxNrhs( integer mnrhs ) {
    LAPACK_WRAPPER_ASSERT(
      mnrhs > 0,
      "LSS::setMaxNrhs, maxNrhs = " << mnrhs
    );
    maxNrhs         = mnrhs;
    maxNrhs_changed = true;
  }

  template <typename T>
  void
  LSS<T>::allocate( integer NR, integer NC ) {

    if ( nRow != NR || nCol != NC || maxNrhs_changed ) {
      nRow = NR;
      nCol = NC;
      valueType tmp;
      integer info = gelss(
        NR, NC, maxNrhs,
        nullptr, NR, nullptr, NR, nullptr,
        rcond, rank, &tmp, -1
      );
      LAPACK_WRAPPER_ASSERT(
        info == 0, "LSS::allocate, in gelss info = " << info
      );

      Lwork = integer(tmp);
      if ( NR != NC ) {
        info = gelss(
          NC, NR, maxNrhs,
          nullptr, NC, nullptr, NC, nullptr,
          rcond, rank, &tmp, -1
        );
        LAPACK_WRAPPER_ASSERT(
          info == 0, "LSS::allocate, in gelss info = " << info
        );
        if ( Lwork < integer(tmp) ) Lwork = integer(tmp);
      }

      integer minRC = std::min(NR,NC);

      allocReals.allocate( size_t(2*NR*NC+Lwork+minRC) );
      Amat     = allocReals( size_t(2*nRow*nCol) );
      Work     = allocReals( size_t(Lwork) );
      sigma    = allocReals( size_t(minRC) );
      AmatWork = Amat+nRow*nCol;

      maxNrhs_changed = false;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSS<T>::solve( valueType xb[] ) const {
    // save matrix
    copy( nRow*nCol, Amat, 1, AmatWork, 1);
    integer info = gelss(
      nRow, nCol, 1,
      AmatWork, nRow,
      xb, nRow,
      sigma, rcond, rank,
      Work, Lwork
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0, "LSS::solve (rhs=1), in gelss info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSS<T>::t_solve( valueType xb[] ) const {
    // save matrix
    for ( integer i = 0; i < nCol; ++i )
      copy( nRow, Amat+i*nRow, 1, AmatWork+i, nCol );
    integer info = gelss(
      nCol, nRow, 1,
      AmatWork, nCol,
      xb, nCol,
      sigma, rcond, rank,
      Work, Lwork
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0, "LSS::t_solve (rhs=1), in gelss info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSS<T>::solve(
    integer   nrhs,
    valueType B[],
    integer   ldB
  ) const {
    // save matrix
    copy( nRow*nCol, Amat, 1, AmatWork, 1 );
    integer info = gelss(
      nRow, nCol, nrhs,
      AmatWork, nRow,
      B, ldB,
      sigma, rcond, rank,
      Work, Lwork
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "LSS::solve (rhs=" << nrhs << "), in gelss info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSS<T>::t_solve(
    integer   nrhs,
    valueType B[],
    integer   ldB
  ) const {
    // save matrix
    for ( integer i = 0; i < nCol; ++i )
      copy( nRow, Amat+i*nRow, 1, AmatWork+i, nCol );
    integer info = gelss(
      nCol, nRow, nrhs,
      AmatWork, nCol,
      B, ldB,
      sigma, rcond, rank,
      Work, Lwork
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "LSS::t_solve (rhs=" << nrhs << "), in gelss info = " << info
    );
  }

  /*\
   |  _     ______   __
   | | |   / ___\ \ / /
   | | |   \___ \\ V /
   | | |___ ___) || |
   | |_____|____/ |_|
   |
  \*/

  template <typename T>
  void
  LSY<T>::setMaxNrhs( integer mnrhs ) {
    LAPACK_WRAPPER_ASSERT(
      mnrhs > 0, "LSY::setMaxNrhs, maxNrhs = " << mnrhs
    );
    maxNrhs         = mnrhs;
    maxNrhs_changed = true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSY<T>::allocate( integer NR, integer NC ) {

    if ( nRow != NR || nCol != NC || maxNrhs_changed ) {
      nRow = NR;
      nCol = NC;
      valueType tmp;
      integer info = gelsy(
        NR, NC, maxNrhs, nullptr, NR, nullptr, NR, nullptr,
        rcond, rank, &tmp, -1
      );
      LAPACK_WRAPPER_ASSERT(
        info == 0, "LSY::allocate, in gelss info = " << info
      );

      Lwork = integer(tmp);
      if ( NR != NC ) {
        info = gelsy(
          NC, NR, maxNrhs, nullptr, NC, nullptr, NC, nullptr,
          rcond, rank, &tmp, -1
        );
        LAPACK_WRAPPER_ASSERT(
          info == 0, "LSY::allocate, in gelss info = " << info
        );
        if ( Lwork < integer(tmp) ) Lwork = integer(tmp);
      }

      allocReals.allocate( size_t(2*NR*NC+Lwork) );
      Amat     = allocReals( size_t(2*nRow*nCol) );
      Work     = allocReals( size_t(Lwork) );
      AmatWork = Amat+nRow*nCol;
      allocInts.allocate( size_t(NC) );
      jpvt = allocInts( size_t(NC) );

      maxNrhs_changed = false;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSY<T>::solve( valueType xb[] ) const {
    // save matrix
    copy( nRow*nCol, Amat, 1, AmatWork, 1);
    integer info = gelsy(
      nRow, nCol, 1,
      AmatWork, nRow,
      xb, nRow, jpvt,
      rcond, rank,
      Work, Lwork
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0, "LSS::solve (rhs=1), in gelss info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSY<T>::t_solve( valueType xb[] ) const {
    // save matrix
    for ( integer i = 0; i < nCol; ++i )
      copy( nRow, Amat+i*nRow, 1, AmatWork+i, nCol );
    integer info = gelsy(
      nCol, nRow, 1,
      AmatWork, nCol,
      xb, nCol, jpvt,
      rcond, rank,
      Work, Lwork
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0, "LSS::t_solve (rhs=1), in gelss info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSY<T>::solve(
    integer   nrhs,
    valueType B[],
    integer   ldB
  ) const {
    // save matrix
    copy( nRow*nCol, Amat, 1, AmatWork, 1 );
    integer info = gelsy(
      nRow, nCol, nrhs,
      AmatWork, nRow,
      B, ldB, jpvt,
      rcond, rank,
      Work, Lwork
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "LSD::solve (rhs=" << nrhs << "), in gelsd info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  LSY<T>::t_solve(
    integer   nrhs,
    valueType B[],
    integer   ldB
  ) const {
    // save matrix
    for ( integer i = 0; i < nCol; ++i )
      copy( nRow, Amat+i*nRow, 1, AmatWork+i, nCol );
    integer info = gelsy(
      nCol, nRow, nrhs,
      AmatWork, nCol,
      B, ldB, jpvt,
      rcond, rank,
      Work, Lwork
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "LSD::t_solve (rhs=" << nrhs << "), in gelsd info = " << info
    );
  }

  //============================================================================

  template <typename valueType>
  inline
  void
  tridiag_axpy(
    integer         N,
    valueType       alpha,
    valueType const L[],
    valueType const D[],
    valueType const U[],
    valueType const x[],
    valueType       beta,
    valueType       y[]
  ) {
    if ( isZero(beta) ) {
      y[0] = alpha*(D[0]*x[0] + U[0] * x[1]);
      for ( integer i = 1; i < N-1; ++i )
        y[i] = alpha*(D[i]*x[i] + U[i] * x[i+1] + L[i-1] * x[i-1]);
      y[N-1] = alpha*(D[N-1]*x[N-1] + L[N-2] * x[N-2]);
    } else {
      y[0] = beta*y[0] + alpha*(D[0]*x[0] + U[0] * x[1]);
      for ( integer i = 1; i < N-1; ++i )
        y[i] = beta*y[i] + alpha*(D[i]*x[i] + U[i] * x[i+1] + L[i-1] * x[i-1]);
      y[N-1] = beta*y[N-1] + alpha*(D[N-1]*x[N-1] + L[N-2] * x[N-2]);
    }
  }

  //============================================================================
  /*\
   |   _____     _     _ _                               _ ____  ____  ____
   |  |_   _| __(_) __| (_) __ _  __ _  ___  _ __   __ _| / ___||  _ \|  _ \
   |    | || '__| |/ _` | |/ _` |/ _` |/ _ \| '_ \ / _` | \___ \| |_) | | | |
   |    | || |  | | (_| | | (_| | (_| | (_) | | | | (_| | |___) |  __/| |_| |
   |    |_||_|  |_|\__,_|_|\__,_|\__, |\___/|_| |_|\__,_|_|____/|_|   |____/
   |                             |___/
  \*/
  template <typename T>
  void
  TridiagonalSPD<T>::factorize(
    char const      who[],
    integer         N,
    valueType const _L[],
    valueType const _D[]
  ) {
    if ( nRC != N ) {
      nRC = N;
      allocReals.allocate(3*N);
      L    = allocReals(N);
      D    = allocReals(N);
      WORK = allocReals(N);
    }
    copy( N, _L, 1, L, 1 );
    copy( N, _D, 1, D, 1 );
    integer info = pttrf( N, L, D );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "TridiagonalSPD::factorize[" << who << 
      "], return info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  T
  TridiagonalSPD<T>::cond1( valueType norm1 ) const {
    valueType rcond;
    integer info = ptcon1( nRC, D, L, norm1, rcond, WORK );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "TridiagonalSPD::cond1, return info = " << info
    );
    return rcond;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalSPD<T>::solve( valueType xb[] ) const {
    integer info = pttrs( nRC, 1, D, L, xb, nRC );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "TridiagonalSPD::solve, return info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalSPD<T>::t_solve( valueType xb[] ) const {
    integer info = pttrs( nRC, 1, D, L, xb, nRC );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "TridiagonalSPD::solve, return info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalSPD<T>::solve( integer nrhs, valueType xb[], integer ldXB ) const {
    integer info = pttrs( nRC, nrhs, D, L, xb, ldXB );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "TridiagonalSPD::solve, return info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalSPD<T>::t_solve( integer nrhs, valueType xb[], integer ldXB ) const {
    integer info = pttrs( nRC, nrhs, D, L, xb, ldXB );
    LAPACK_WRAPPER_ASSERT(
      info == 0, 
      "TridiagonalSPD::solve, return info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalSPD<T>::axpy(
    integer         N,
    valueType       alpha,
    valueType const _L[],
    valueType const _D[],
    valueType const x[],
    valueType       beta,
    valueType       y[]
  ) const {
    tridiag_axpy( N, alpha, _L, _D, _L, x, beta, y );
  }

  //============================================================================
  /*\
   |   _____     _     _ _                               _ _    _   _
   |  |_   _| __(_) __| (_) __ _  __ _  ___  _ __   __ _| | |  | | | |
   |    | || '__| |/ _` | |/ _` |/ _` |/ _ \| '_ \ / _` | | |  | | | |
   |    | || |  | | (_| | | (_| | (_| | (_) | | | | (_| | | |__| |_| |
   |    |_||_|  |_|\__,_|_|\__,_|\__, |\___/|_| |_|\__,_|_|_____\___/
   |                             |___/
  \*/
  template <typename T>
  void
  TridiagonalLU<T>::factorize(
    char const      who[],
    integer         N,
    valueType const _L[],
    valueType const _D[],
    valueType const _U[]
  ) {
    if ( nRC != N ) {
      nRC = N;
      allocReals.allocate(6*N);
      allocIntegers.allocate(2*N);
      L     = allocReals(N);
      D     = allocReals(N);
      U     = allocReals(N);
      U2    = allocReals(N);
      WORK  = allocReals(2*N);
      IPIV  = allocIntegers(N);
      IWORK = allocIntegers(N);
    }
    copy( N, _L, 1, L, 1 );
    copy( N, _D, 1, D, 1 );
    copy( N, _U, 1, U, 1 );
    integer info = gttrf( N, L, D, U, U2, IPIV );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "TridiagonalLU::factorize[" << who << 
      "], return info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  T
  TridiagonalLU<T>::cond1( valueType norm1 ) const {
    valueType rcond;
    integer info = gtcon1( nRC, L, D, U, U2, IPIV, norm1, rcond, WORK, IWORK );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "TridiagonalLU::cond1, return info = " << info
    );
    return rcond;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  T
  TridiagonalLU<T>::condInf( valueType normInf ) const {
    valueType rcond;
    integer info = gtconInf( nRC, L, D, U, U2, IPIV, normInf, rcond, WORK, IWORK );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "TridiagonalLU::cond1, return info = " << info
    );
    return rcond;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalLU<T>::solve( valueType xb[] ) const {
    integer info = gttrs( NO_TRANSPOSE, nRC, 1, L, D, U, U2, IPIV, xb, nRC );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "TridiagonalLU::solve, return info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalLU<T>::t_solve( valueType xb[] ) const {
    integer info = gttrs( TRANSPOSE, nRC, 1, L, D, U, U2, IPIV, xb, nRC );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "TridiagonalLU::solve, return info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalLU<T>::solve( integer nrhs, valueType xb[], integer ldXB ) const {
    integer info = gttrs( NO_TRANSPOSE, nRC, nrhs, L, D, U, U2, IPIV, xb, ldXB );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "TridiagonalLU::solve, return info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalLU<T>::t_solve( integer nrhs, valueType xb[], integer ldXB ) const {
    integer info = gttrs( TRANSPOSE, nRC, nrhs, L, D, U, U2, IPIV, xb, ldXB );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "TridiagonalLU::solve, return info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalLU<T>::axpy(
    integer         N,
    valueType       alpha,
    valueType const _L[],
    valueType const _D[],
    valueType const _U[],
    valueType const x[],
    valueType       beta,
    valueType       y[]
  ) const {
    tridiag_axpy( N, alpha, _L, _D, _U, x, beta, y );
  }

  //============================================================================
  /*\
   |   _____     _     _ _                               _  ___  ____
   |  |_   _| __(_) __| (_) __ _  __ _  ___  _ __   __ _| |/ _ \|  _ \
   |    | || '__| |/ _` | |/ _` |/ _` |/ _ \| '_ \ / _` | | | | | |_) |
   |    | || |  | | (_| | | (_| | (_| | (_) | | | | (_| | | |_| |  _ <
   |    |_||_|  |_|\__,_|_|\__,_|\__, |\___/|_| |_|\__,_|_|\__\_\_| \_\
   |                             |___/
  \*/

  template <typename T>
  void
  TridiagonalQR<T>::factorize(
    char const      /* who */[],
    integer         N,
    valueType const L[],
    valueType const D[],
    valueType const U[]
  ) {
    allocReals.allocate(size_t(5*(N-1)));
    nRC = N;
    C   = allocReals(size_t(N-1));
    S   = allocReals(size_t(N-1));
    BD  = allocReals(size_t(N));
    BU  = allocReals(size_t(N-1));
    BU2 = allocReals(size_t(N-2));

    /*\
      | d u       | d u @     | d u @     | d u @     | d u @     |
      | l d u     | 0 d u     | 0 d u @   | 0 d u @   | 0 d u @   |
      |   l d u   |   l d u   |   0 d u   |   0 d u @ |   0 d u @ |
      |     l d u |     l d u |     l d u |     0 d u |     0 d u |
      |       l d |       l d |       l d |       l d |       0 d |
    \*/

    lapack_wrapper::copy( N,   D, 1, BD, 1 );
    lapack_wrapper::copy( N-1, U, 1, BU, 1 );
    lapack_wrapper::zero( N-2, BU2, 1 );

    normInfA = 0;
    integer i = 0;
    for (; i < N-2; ++i ) {
      valueType Li = L[i];
      rotg( BD[i], Li, C[i], S[i] );
      rot( 1, &BU[i],  1, &BD[i+1], 1, C[i], S[i] );
      rot( 1, &BU2[i], 1, &BU[i+1], 1, C[i], S[i] );
      valueType sum = std::abs(BD[i]) + std::abs(BU[i]) + std::abs(BU2[i]);
      if ( sum > normInfA ) normInfA = sum;
    }
    valueType Li = L[i];
    rotg( BD[i], Li, C[i], S[i] );
    rot( 1, &BU[i], 1, &BD[i+1], 1, C[i], S[i] );

    valueType sum = std::abs(BD[i]) + std::abs(BU[i]);
    if ( sum > normInfA ) normInfA = sum;
    sum = std::abs(BD[i+1]);
    if ( sum > normInfA ) normInfA = sum;

    // Q A = R
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalQR<T>::Rsolve( valueType xb[] ) const {
    xb[nRC-1] /= BD[nRC-1];
    xb[nRC-2] = (xb[nRC-2]-BU[nRC-2]*xb[nRC-1])/BD[nRC-2];
    for ( integer i = nRC-3; i >= 0; --i )
      xb[i] = (xb[i]-BU[i]*xb[i+1]-BU2[i]*xb[i+2])/BD[i];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalQR<T>::RsolveTransposed( valueType xb[] ) const {
    xb[0] /= BD[0];
    xb[1] = (xb[1]-BU[0]*xb[0])/BD[1];
    for ( integer i = 2; i < nRC; ++i )
      xb[i] = (xb[i]-BU[i]*xb[i-1]-BU2[i]*xb[i-2])/BD[i];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalQR<T>::solve( valueType xb[] ) const {
    // A x = b --> Q A x = Q b --> R x = Q b
    // applico Q b
    for ( integer i = 0; i < nRC-1; ++i )
      rot( 1, &xb[i], 1, &xb[i+1], 1, C[i], S[i] );
    Rsolve( xb );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalQR<T>::t_solve( valueType xb[] ) const {
    // A^T x = b --> A^T Q^T Q x = b --> R^T Q x = b --> R^T y = b  x = Q^T y
    RsolveTransposed( xb );
    // applico Q^T b
    for ( integer i = nRC-2; i >= 0; --i )
      rot( 1, &xb[i], 1, &xb[i+1], 1, C[i], -S[i] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalQR<T>::solve( integer nrhs, valueType xb[], integer ldXB ) const {
    // A x = b --> Q A x = Q b --> R x = Q b
    // applico Q b
    for ( integer i = 0; i < nRC-1; ++i )
      rot( nrhs, &xb[i], ldXB, &xb[i+1], ldXB, C[i], S[i] );
    for ( integer i = 0; i < nrhs; ++i )
      Rsolve( xb+i*ldXB );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalQR<T>::t_solve( integer nrhs, valueType xb[], integer ldXB ) const {
    // A^T x = b --> A^T Q^T Q x = b --> R^T Q x = b --> R^T y = b  x = Q^T y
    for ( integer i = 0; i < nrhs; ++i )
      RsolveTransposed(xb+i*ldXB);
    for ( integer i = nRC-2; i >= 0; --i )
      rot( nrhs, &xb[i], ldXB, &xb[i+1], ldXB, C[i], -S[i] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  TridiagonalQR<T>::axpy(
    integer         N,
    valueType       alpha,
    valueType const L[],
    valueType const D[],
    valueType const U[],
    valueType const x[],
    valueType       beta,
    valueType       y[]
  ) const {
    tridiag_axpy( N, alpha, L, D, U, x, beta, y );
  }

  /*\
   *
   *  Solve
   *
   *  min || T x - b ||^2 + lambda ||x||^2
   *
   *
   *  / * * *       \
   *  |   * * *     |
   *  |     * * *   |
   *  |       * * * |
   *  |         * * |
   *  |           * |
   *  | x - - - - - |
   *  |   x         |
   *  |     x       |
   *  |       x     |
   *  |         x   |
   *  \           x /
   *
  \*/

  template <typename T>
  void
  TridiagonalQR<T>::lsq(
    integer nrhs,
    T       RHS[],
    integer ldRHS,
    T       lambda_in
  ) const {

    valueType lambda = normInfA * lambda_in;
    std::vector<T> D(nRC),
                   U(nRC-1),
                   U2(nRC-2),
                   tmp(nrhs);
    T CC, SS;

    for ( integer i = 0; i < nRC-1; ++i )
      rot( nrhs, RHS+i, ldRHS, RHS+i+1, ldRHS, C[i], S[i] );

    copy( nRC,   BD,  1, &D.front(),  1 );
    copy( nRC-1, BU,  1, &U.front(),  1 );
    copy( nRC-2, BU2, 1, &U2.front(), 1 );
    T line[3];
    integer i = 0;
    while ( i < nRC-1 ) {
      line[0] = line[2] = 0; line[1] = lambda;
      std::fill( tmp.begin(), tmp.end(), T(0) );
      integer j = i;
      while ( j < nRC-2 ) {
        line[0] = line[1];
        line[1] = line[2];
        line[2] = 0;
        rotg( D[j], line[0], CC, SS );
        rot( 1, &U[j],  1, &line[1], 1, CC, SS );
        rot( 1, &U2[j], 1, &line[2], 1, CC, SS );
        rot( nrhs, RHS+j, ldRHS, &tmp.front(), 1, CC, SS );
        ++j;
      }
      // penultima
      rotg( D[j], line[1], CC, SS );
      rot( 1, &U[j],  1, &line[2], 1, CC, SS );
      rot( nrhs, RHS+j, ldRHS, &tmp.front(), 1, CC, SS );
      ++j;

      rotg( D[j], line[2], CC, SS );
      rot( nrhs, RHS+j, ldRHS, &tmp.front(), 1, CC, SS );

      // ultima
      line[2] = lambda;
      rotg( D[j], line[2], CC, SS );
      rot( nrhs, RHS+j, ldRHS, &tmp.front(), 1, CC, SS );

      ++i; // next lambda
    }
    if ( nRC > 0 ) {
      integer j = nRC-1;
      line[0] = lambda;
      std::fill( tmp.begin(), tmp.end(), T(0) );
      rotg( D[j], line[0], CC, SS );
      rot( nrhs, RHS+j, ldRHS, &tmp.front(), 1, CC, SS );
    }

    for ( integer j = 0; j < nrhs; ++j ) {
      T * xb = RHS + j*ldRHS;
      xb[nRC-1] /= D[nRC-1];
      xb[nRC-2] = (xb[nRC-2]-U[nRC-2]*xb[nRC-1])/D[nRC-2];
      for ( integer k = nRC-3; k >= 0; --k )
        xb[k] = (xb[k]-U[k]*xb[k+1]-U2[k]*xb[k+2])/D[k];
    }
  }

  //============================================================================
  /*\
   |  ___ _         _     _____    _    _ _                         _
   | | _ ) |___  __| |__ |_   _| _(_)__| (_)__ _ __ _ ___ _ _  __ _| |
   | | _ \ / _ \/ _| / /   | || '_| / _` | / _` / _` / _ \ ' \/ _` | |
   | |___/_\___/\__|_\_\   |_||_| |_\__,_|_\__,_\__, \___/_||_\__,_|_|
   |                                            |___/
   |  ___                _       _
   | / __|_  _ _ __  ___| |_ _ _(_)__
   | \__ \ || | '  \/ -_)  _| '_| / _|
   | |___/\_, |_|_|_\___|\__|_| |_\__|
   |      |__/
  \*/

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::setup(
    integer       nblks,
    integer const rBlocks[]
  ) {
    integer N = rBlocks[nblks], nrmax = 0;
    allocIntegers.allocate( N + nblks+1 );
    allocRpointers.allocate( 2*nblks-1 );
    allocIpointers.allocate( nblks );
    this->row_blocks    = allocIntegers( nblks+1 );
    this->D_blocks      = allocRpointers( nblks );
    this->L_blocks      = allocRpointers( nblks-1 );
    this->B_permutation = allocIpointers( nblks );
    // evalute the memry usage for the L and D blocks
    integer nr0 = rBlocks[1] - rBlocks[0];
    nrmax = nr0;
    this->nnz = nr0*nr0;

    LAPACK_WRAPPER_ASSERT(
      nr0 >= 0,
      "BlockTridiagonalSymmetic::setup, nr0 = " << nr0
    );

    for ( integer i = 1; i < nblks; ++i ) {
      integer nr = rBlocks[i+1] - rBlocks[i];
      LAPACK_WRAPPER_ASSERT(
        nr >= 0,
        "BlockTridiagonalSymmetic::setup, nr = " << nr
      );
      this->nnz += nr*(nr+nr0);
      if ( nr > nrmax ) nrmax = nr;
      nr0 = nr;
    }
    allocReals.allocate( this->nnz + nrmax * nrmax );
    nr0 = rBlocks[1] - rBlocks[0];
    this->D_blocks[0] = allocReals( nr0 * nr0 );
    B_permutation[0] = allocIntegers( nr0 );
    for ( integer i = 1; i < nblks; ++i ) {
      integer nr = rBlocks[i+1] - rBlocks[i];
      this->D_blocks[i]   = allocReals( nr * nr );
      this->L_blocks[i-1] = allocReals( nr * nr0 );
      nr0 = nr;
      this->B_permutation[i] = allocIntegers( nr );
    }

    Work = allocReals( nrmax * nrmax );

    this->zero();
    this->nBlocks = nblks;
    std::copy( rBlocks, rBlocks+nblks+1, this->row_blocks );
    is_factorized = false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::setup(
    integer nblks, integer const block_size
  ) {
    // use a temporary vector
    std::vector<integer> rBlocks;
    rBlocks.reserve(nblks+1);
    integer n = 0;
    for ( integer k = 0; k <= nblks; ++k, n += block_size )
      rBlocks.push_back(n);
    this->setup( nblks, &rBlocks.front() );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::zero() { // fill to 0 all the blocks
    std::fill( D_blocks[0], D_blocks[0]+this->nnz, 0 );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  T const &
  BlockTridiagonalSymmetic<T>::operator () ( integer ii, integer jj ) const {
    integer const * ib = this->row_blocks;
    integer const * ie = this->row_blocks+this->nBlocks+1;
    integer iBlock = integer(std::upper_bound( ib, ie, ii ) - ib) - 1;
    integer jBlock = integer(std::upper_bound( ib, ie, jj ) - ib) - 1;

    LAPACK_WRAPPER_ASSERT(
      row_blocks[iBlock] <= ii && ii < row_blocks[iBlock+1],
      "bad iBlock"
    );
    LAPACK_WRAPPER_ASSERT(
      row_blocks[jBlock] <= jj && jj < row_blocks[jBlock+1],
      "bad iBlock"
    );

    integer nr = row_blocks[iBlock+1] - row_blocks[iBlock];
    integer i  = ii - row_blocks[iBlock];
    integer j  = jj - row_blocks[jBlock];
    LAPACK_WRAPPER_ASSERT(
      jBlock == iBlock || jBlock+1 == iBlock,
      "BlockTridiagonalSymmetic:find( " << ii << " , " << jj <<
      " ) --> ( iBlock = " << iBlock << ", jBlock = " << jBlock <<
      " ) --> ( i = " << i << ", j = " << j << " ) out of range"
    );

    integer ij = i+j*nr;
    if ( iBlock == jBlock ) {
      return this->D_blocks[jBlock][ij];
    } else {
      return this->L_blocks[jBlock][ij];
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  T &
  BlockTridiagonalSymmetic<T>::operator () ( integer ii, integer jj ) {
    integer const * ib = this->row_blocks;
    integer const * ie = this->row_blocks+this->nBlocks+1;
    integer iBlock = integer(std::upper_bound( ib, ie, ii ) - ib) - 1;
    integer jBlock = integer(std::upper_bound( ib, ie, jj ) - ib) - 1;

    LAPACK_WRAPPER_ASSERT(
      row_blocks[iBlock] <= ii && ii < row_blocks[iBlock+1],
      "bad iBlock"
    );
    LAPACK_WRAPPER_ASSERT(
      row_blocks[jBlock] <= jj && jj < row_blocks[jBlock+1],
      "bad iBlock"
    );

    integer nr = row_blocks[iBlock+1] - row_blocks[iBlock];
    integer i  = ii - row_blocks[iBlock];
    integer j  = jj - row_blocks[jBlock];
    LAPACK_WRAPPER_ASSERT(
      jBlock == iBlock || jBlock+1 == iBlock,
      "BlockTridiagonalSymmetic:find( " << ii << " , " << jj <<
      " ) --> ( iBlock = " << iBlock << ", jBlock = " << jBlock <<
      " ) --> ( i = " << i << ", j = " << j << " ) out of range"
    );

    integer ij = i+j*nr;
    if ( iBlock == jBlock ) {
      return this->D_blocks[jBlock][ij];
    } else {
      return this->L_blocks[jBlock][ij];
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::insert(
    integer   ii,
    integer   jj,
    valueType v,
    bool      sym
  ) {
    integer const * ib = this->row_blocks;
    integer const * ie = this->row_blocks+this->nBlocks+1;
    integer jBlock = integer(std::upper_bound( ib, ie, jj ) - ib) - 1;
    integer jmin   = row_blocks[jBlock];
    integer jmax   = row_blocks[jBlock+1];
    integer j      = jj - jmin;
    LAPACK_WRAPPER_ASSERT(
      jmin <= jj && jj < jmax && jmin <= ii,
      "bad iBlock ii = " << ii << " jj = " << jj
    );
    if ( ii < jmax ) {
      integer i  = ii - jmin;
      integer nr = jmax-jmin;
      valueType * D = this->D_blocks[jBlock];
      D[i+j*nr] = v;
      if ( sym ) D[j+i*nr] = v;
    } else {
      integer imin = jmax;
      integer imax = row_blocks[jBlock+2];
      integer nr   = imax-imin;
      integer i    = ii - imin;
      this->L_blocks[jBlock][i+j*nr] = v;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::setD(
    integer         n,
    valueType const data[],
    integer         ldData,
    bool            transposed
  ) {
    integer nr = this->DnumRows(n);
    if ( transposed ) {
      getranspose( nr, nr, data, ldData, this->D_blocks[n], nr );
    } else {
      integer ierr = gecopy( nr, nr, data, ldData, this->D_blocks[n], nr );
      LAPACK_WRAPPER_ASSERT(
        ierr == 0,
        "BlockTridiagonalSymmetic::setD, gecopy return ierr = " << ierr
      );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::setL(
    integer         n,
    valueType const data[],
    integer         ldData,
    bool            transposed
  ) {
    integer nr = this->LnumRows(n);
    integer nc = this->LnumCols(n);
    if ( transposed ) {
      getranspose( nr, nc, data, ldData, this->L_blocks[n], nr );
    } else {
      integer ierr = gecopy( nr, nc, data, ldData, this->L_blocks[n], nr );
      LAPACK_WRAPPER_ASSERT(
        ierr == 0,
        "BlockTridiagonalSymmetic::setL, gecopy return ierr = " << ierr
      );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::setD(
    integer           n,
    valueType const * data,
    integer           ldData,
    integer           beginRow,
    integer           beginCol,
    integer           nrow,
    integer           ncol,
    bool              transposed
  ) {
    integer ldD = this->DnumRows(n);
    valueType * D = D_blocks[n]+beginRow+beginCol*ldD;

    if ( transposed ) {
      getranspose( nrow, ncol, data, ldData, D, ldD );
    } else {
      integer ierr = gecopy( nrow, ncol, data, ldData, D, ldD );
      LAPACK_WRAPPER_ASSERT(
        ierr == 0,
        "BlockTridiagonalSymmetic::setD (block), gecopy return ierr = " << ierr
      );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::setL(
    integer           n,
    valueType const * data,
    integer           ldData,
    integer           beginRow,
    integer           beginCol,
    integer           nrow,
    integer           ncol,
    bool              transposed
  ) {
    integer ldL = this->LnumRows(n);
    valueType * L = L_blocks[n]+beginRow+beginCol*ldL;
    if ( transposed ) {
      getranspose( nrow, ncol, data, ldData, L, ldL );
    } else {
      integer ierr = gecopy( nrow, ncol, data, ldData, L, ldL );
      LAPACK_WRAPPER_ASSERT(
        ierr == 0,
        "BlockTridiagonalSymmetic::setL (block), gecopy return ierr = " << ierr
      );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::factorize( char const who[] ) {
    LAPACK_WRAPPER_ASSERT(
      !is_factorized,
      "BlockTridiagonalSymmetic::factorize[" << who << 
      "], already factored"
    );

    integer info = lapack_wrapper::getrf(
      this->DnumRows(0), this->DnumCols(0),
      this->D_blocks[0], this->DnumRows(0),
      B_permutation[0]
    );

    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "BlockTridiagonalSymmetic::factorize[" << who << 
      "] getrf INFO = " << info
    );

    for ( integer k=1; k < this->nBlocks; ++k ) {
      integer nr0 = this->DnumRows(k-1);
      integer nr1 = this->DnumRows(k);
      valueType * L0 = this->L_blocks[k-1];
      valueType * D0 = this->D_blocks[k-1];
      valueType * D1 = this->D_blocks[k];
      // Work = L^T nr0 x nr1
      getranspose( nr1, nr0, L0, nr1, Work, nr0 );
      // solve
      info = getrs( TRANSPOSE, nr0, nr1, D0, nr0, B_permutation[k-1], Work, nr0 );
      LAPACK_WRAPPER_ASSERT(
        info == 0,
        "BlockTridiagonalSymmetic::factorize[" << who << 
        "] getrs INFO = " << info
      );

      // DD{k}   = DD{k} - (LL{k-1}*DD{k-1}) *LL{k-1}.';

      gemm(
        TRANSPOSE,
        TRANSPOSE,
        nr1, nr1, nr0,
        -1.0, Work, nr0,
        L0, nr1,
        1.0, D1, nr1
      );

      // Work --> L nr1 x nr0
      getranspose( nr0, nr1, Work, nr0, L0, nr1 );

      info = getrf( nr1, nr1, D1, nr1, B_permutation[k] );

      LAPACK_WRAPPER_ASSERT(
        info == 0,
        "BlockTridiagonalSymmetic::factorize[" << who << 
        "] getrf INFO = " << info
      );

    }
    is_factorized = true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::solve( valueType xb[] ) const {
    LAPACK_WRAPPER_ASSERT(
      is_factorized,
      "BlockTridiagonalSymmetic::solve, matrix not factored"
    );

    // RR{k} = RR{k}-LL{k-1}*RR{k-1};
    integer k = 0;
    integer nr0 = this->DnumRows(0), nr1;
    valueType * xkm1 = xb, *xk;
    while ( ++k < this->nBlocks ) {
      nr1 = this->DnumRows(k);
      valueType const * L0 = this->L_blocks[k-1];
      xk = xkm1 + nr0;
      gemv(
        NO_TRANSPOSE, nr1, nr0,
        -1.0, L0, nr1,
        xkm1, 1,
        1.0,
        xk, 1
      );
      xkm1 = xk;
      nr0  = nr1;
    }
    // RR{k} = DD{k}\RR{k};
    xk = xb;
    for ( k = 0; k < this->nBlocks; ++k ) {
      nr1 = this->DnumRows(k);
      // solve
      valueType const * D1 = this->D_blocks[k];
      integer info = getrs( NO_TRANSPOSE, nr1, 1, D1, nr1, B_permutation[k], xk, nr1 );
      LAPACK_WRAPPER_ASSERT(
        info == 0,
        "BlockTridiagonalSymmetic::solve getrs INFO = " << info
      );
      xk += nr1;
    }
    // RR{k} = RR{k}-LL{k}.'*RR{k+1};
    nr1 = this->DnumRows(k-1);
    xk -= nr1;
    while ( --k > 0 ) {
      nr0 = this->DnumRows(k-1);
      valueType const * L0 = this->L_blocks[k-1];
      xkm1 = xk - nr0;
      gemv(
        TRANSPOSE, nr1, nr0,
        -1.0, L0, nr1,
        xk, 1,
        1.0,
        xkm1, 1
      );
      xk  = xkm1;
      nr1 = nr0;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::solve(
    integer   nrhs,
    valueType B[],
    integer   ldB
  ) const {
    LAPACK_WRAPPER_ASSERT(
      is_factorized,
      "BlockTridiagonalSymmetic::solve, matrix not factored"
    );

    // RR{k} = RR{k}-LL{k-1}*RR{k-1};
    integer k = 0;
    integer nr0 = this->DnumRows(0), nr1;
    valueType * Bkm1 = B, *Bk;
    while ( ++k < this->nBlocks ) {
      nr1 = this->DnumRows(k);
      valueType const * L0 = this->L_blocks[k-1];
      Bk = Bkm1 + nr0;
      gemm(
        NO_TRANSPOSE, NO_TRANSPOSE, nr1, nrhs, nr0,
        -1.0, L0, nr1,
        Bkm1, ldB,
        1.0,
        Bk, ldB
      );
      Bkm1 = Bk;
      nr0  = nr1;
    }
    // RR{k} = DD{k}\RR{k};
    Bk = B;
    for ( k = 0; k < this->nBlocks; ++k ) {
      nr1 = this->DnumRows(k);
      // solve
      valueType const * D1 = this->D_blocks[k];
      integer info = getrs(
        NO_TRANSPOSE, nr1, nrhs,
        D1, nr1, B_permutation[k],
        Bk, ldB
      );
      LAPACK_WRAPPER_ASSERT(
        info == 0,
        "BlockTridiagonalSymmetic::solve getrs INFO = " << info
      );
      Bk += nr1;
    }
    // RR{k} = RR{k}-LL{k}.'*RR{k+1};
    nr1 = this->DnumRows(k-1);
    Bk -= nr1;
    while ( --k > 0 ) {
      nr0 = this->DnumRows(k-1);
      valueType const * L0 = this->L_blocks[k-1];
      Bkm1 = Bk - nr0;
      gemm(
        TRANSPOSE, NO_TRANSPOSE, nr0, nrhs, nr1,
        -1.0, L0, nr1,
        Bk, ldB,
        1.0,
        Bkm1, ldB
      );
      Bk  = Bkm1;
      nr1 = nr0;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::t_solve( valueType xb[] ) const
  { BlockTridiagonalSymmetic<T>::solve( xb ); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BlockTridiagonalSymmetic<T>::t_solve(
    integer   nrhs,
    valueType xb[],
    integer   ldXB
  ) const {
    BlockTridiagonalSymmetic<T>::solve( nrhs, xb, ldXB );
  }

  /*\
   |   ____                  _          _ __  __       _        _
   |  | __ )  __ _ _ __   __| | ___  __| |  \/  | __ _| |_ _ __(_)_  __
   |  |  _ \ / _` | '_ \ / _` |/ _ \/ _` | |\/| |/ _` | __| '__| \ \/ /
   |  | |_) | (_| | | | | (_| |  __/ (_| | |  | | (_| | |_| |  | |>  <
   |  |____/ \__,_|_| |_|\__,_|\___|\__,_|_|  |_|\__,_|\__|_|  |_/_/\_\
  \*/

  template <typename T>
  BandedLU<T>::BandedLU()
  : allocReals("_BandedLU_reals")
  , allocIntegers("_BandedLU_integers")
  , m(0)
  , n(0)
  , nL(0)
  , nU(0)
  , ldAB(0)
  , is_factorized(false)
  {}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  BandedLU<T>::~BandedLU()
  {}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::setup(
    integer _m,
    integer _n,
    integer _nL,
    integer _nU
  ) {
    m    = _m;
    n    = _n;
    nL   = _nL;
    nU   = _nU;
    ldAB = 2*nL+nU+1;
    integer nnz = n*ldAB;
    allocReals.allocate( nnz );
    allocIntegers.allocate(m);
    AB   = allocReals( nnz );
    ipiv = allocIntegers( m );
    is_factorized = false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::solve( valueType xb[] ) const {
    LAPACK_WRAPPER_ASSERT( is_factorized, "BandedLU::solve, matrix not yet factorized" );
    LAPACK_WRAPPER_ASSERT( m == n, "BandedLU::solve, matrix must be square" );
    integer info = gbtrs( NO_TRANSPOSE, m, nL, nU, 1, AB, ldAB, ipiv, xb, m );
    LAPACK_WRAPPER_ASSERT( info == 0, "BandedLU::solve, info = " << info );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::t_solve( valueType xb[] ) const {
    LAPACK_WRAPPER_ASSERT( is_factorized, "BandedLU::solve, matrix not yet factorized" );
    LAPACK_WRAPPER_ASSERT( m == n, "BandedLU::solve, matrix must be square" );
    integer info = gbtrs( TRANSPOSE, m, nL, nU, 1, AB, ldAB, ipiv, xb, m );
    LAPACK_WRAPPER_ASSERT( info == 0, "BandedLU::t_solve, info = " << info );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::solve( integer nrhs, valueType B[], integer ldB ) const {
    LAPACK_WRAPPER_ASSERT( is_factorized, "BandedLU::solve, matrix not yet factorized" );
    LAPACK_WRAPPER_ASSERT( m == n, "BandedLU::solve, matrix must be square" );
    integer info = gbtrs( NO_TRANSPOSE, m, nL, nU, nrhs, AB, ldAB, ipiv, B, ldB );
    LAPACK_WRAPPER_ASSERT( info == 0, "BandedLU::solve, info = " << info );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::t_solve( integer nrhs, valueType B[], integer ldB ) const {
    LAPACK_WRAPPER_ASSERT( is_factorized, "BandedLU::solve, matrix not yet factorized" );
    LAPACK_WRAPPER_ASSERT( m == n, "BandedLU::solve, matrix must be square" );
    integer info = gbtrs( TRANSPOSE, m, nL, nU, nrhs, AB, ldAB, ipiv, B, ldB );
    LAPACK_WRAPPER_ASSERT( info == 0, "BandedLU::t_solve, info = " << info );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::factorize( char const      who[] ) {
    LAPACK_WRAPPER_ASSERT(
      !is_factorized,
      "BandedLU::factorize[" << who << "], matrix yet factorized"
    );
    LAPACK_WRAPPER_ASSERT(
      m == n,
      "BandedLU::factorize[" << who << "], matrix must be square"
    );
    integer info = gbtrf( m, n, nL, nU, AB, ldAB, ipiv );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "BandedLU::factorize[" << who << "], info = " << info
    );
    is_factorized = true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::zero() {
    integer nnz = m*(2*nL+nU+1);
    lapack_wrapper::zero( nnz, AB, 1 );
    is_factorized = false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::check( integer i, integer j ) const {
    LAPACK_WRAPPER_ASSERT(
      i >= 0 && i < m && j >= 0 && j < n,
      "BandedLU::check( " << i << " , " << j << " ) out of range"
    );
    LAPACK_WRAPPER_ASSERT(
      j >= i-nL && j <= i+nU,
      "BandedLU::check( " << i << " , " << j << " ) out of band"
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::insert(
    integer   i,
    integer   j,
    valueType v,
    bool      sym
  ) {
    LAPACK_WRAPPER_ASSERT(
      i >= 0 && i < m && j >= 0 && j < n,
      "BandedLU::insert( " << i << " , " << j << " ) out of range"
    );
    LAPACK_WRAPPER_ASSERT(
      i-nL <= j && j <= i+nU,
      "BandedLU::insert( " << i << " , " << j << " ) out of band"
    );
    (*this)(i,j) = v;
    if ( sym && i != j ) {
      LAPACK_WRAPPER_ASSERT(
        j-nL <= i && i <= j+nU,
        "BandedLU::insert( " << i << " , " << j << " ) out of band"
      );
      (*this)(j,i) = v;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::load_block(
    integer         nr,
    integer         nc,
    valueType const B[],
    integer         ldB,
    integer         irow,
    integer         icol
  ) {

    LAPACK_WRAPPER_ASSERT(
      !is_factorized,
      "BandedLU::load_block, matrix is factorized"
    );

    #if 1
    for ( integer r = 0; r < nr; ++r )
      for ( integer c = 0; c < nc; ++c )
        AB[iaddr( irow+r, icol+c )] = B[ r + c * ldB ];
    #else
    // must be checked
    check( irow,      icol      );
    check( irow+nr-1, icol+nc-1 );
    check( irow,      icol+nc-1 );
    check( irow+nr-1, icol      );

    // copy by diagonal
    for ( integer r = 0; r < nr; ++r ) {
      integer ia = iaddr( irow+r, icol );
      copy( std::min(nr-r,nc), B+r, ldB+1, AB+ia, ldAB );
    }
    for ( integer c = 1; c < nc; ++c ) {
      integer ia = iaddr( irow, icol+c );
      copy( std::min(nc-c,nr), B+c*ldB, ldB+1, AB+ia, ldAB );
    }
    #endif
  }

  // y <- beta*y + alpha*A*x
  /*
    +---------+
    | \       |
    |  \      |
    +---+-----+
  */
  template <typename T>
  void
  BandedLU<T>::aAxpy(
    valueType       alpha,
    valueType const x[],
    valueType       y[]
  ) const {

    valueType const * col = AB + nL;
    for ( integer j = 0; j < n; ++j, col += ldAB ) {
      integer imin  = j-nU;
      integer imax  = std::min(j+nL,m-1);
      integer imin0 = imin > 0 ? imin : 0;
      lapack_wrapper::axpy(
        imax-imin0+1,
        alpha*x[j],
        col+imin0-imin, 1,
        y+imin0,        1
      );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedLU<T>::dump( ostream_type & stream ) const {
    for ( integer i = 0; i <= nL+nU; ++i ) {
      valueType const * col = AB + nL + i;
      for ( integer j = 0; j < n; ++j, col += ldAB )
        stream << std::setw(10) << col[0] << ' ';
      stream << '\n';
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  BandedSPD<T>::BandedSPD()
  : allocReals("_BandedSPD_reals")
  , n(0)
  , nD(0)
  , ldAB(0)
  , is_factorized(false)
  {}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  BandedSPD<T>::~BandedSPD()
  {}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedSPD<T>::setup(
    ULselect _UPLO,
    integer  _N,
    integer  _nD
  ) {
    UPLO = _UPLO;
    n    = _N;
    nD   = _nD;
    ldAB = nD+1;
    integer nnz = n*ldAB;
    allocReals.allocate( nnz );
    AB   = allocReals( nnz );
    is_factorized = false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedSPD<T>::solve( valueType xb[] ) const {
    LAPACK_WRAPPER_ASSERT(
      is_factorized,
      "BandedSPD::solve, matrix not yet factorized"
    );
    integer info = pbtrs( UPLO, n, nD, 1, AB, ldAB, xb, n );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "BandedSPD::solve, info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedSPD<T>::t_solve( valueType xb[] ) const {
    LAPACK_WRAPPER_ASSERT(
      is_factorized,
      "BandedSPD::solve, matrix not yet factorized"
    );
    integer info = pbtrs( UPLO, n, nD, 1, AB, ldAB, xb, n );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "BandedSPD::t_solve, info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedSPD<T>::solve( integer nrhs, valueType B[], integer ldB ) const {
    LAPACK_WRAPPER_ASSERT(
      is_factorized,
      "BandedSPD::solve, matrix not yet factorized"
    );
    integer info = pbtrs( UPLO, n, nD, nrhs, AB, ldAB, B, ldB );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "BandedSPD::solve, info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedSPD<T>::t_solve( integer nrhs, valueType B[], integer ldB ) const {
    LAPACK_WRAPPER_ASSERT(
      is_factorized,
      "BandedSPD::solve, matrix not yet factorized"
    );
    integer info = pbtrs( UPLO, n, nD, nrhs, AB, ldAB, B, ldB );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "BandedSPD::t_solve, info = " << info
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedSPD<T>::factorize( char const who[] ) {
    LAPACK_WRAPPER_ASSERT(
      !is_factorized,
      "BandedSPD::factorize[" << who << "], matrix yet factorized"
    );
    integer info = pbtrf( UPLO, n, nD, AB, ldAB );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "BandedSPD::factorize[" << who << "], info = " << info
    );
    is_factorized = true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedSPD<T>::zero() {
    lapack_wrapper::zero( n*ldAB, AB, 1 );
    is_factorized = false;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  BandedSPD<T>::insert(
    integer   i,
    integer   j,
    valueType v,
    bool      sym
  ) {
    (*this)(i,j) = v;
    if ( sym && i != j ) (*this)(j,i) = v;
  }

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

  /*\
   |    ___                  _ _   _               _
   |   / _ \ _   _  __ _ ___(_) \ | | _____      _| |_ ___  _ __
   |  | | | | | | |/ _` / __| |  \| |/ _ \ \ /\ / / __/ _ \| '_ \
   |  | |_| | |_| | (_| \__ \ | |\  |  __/\ V  V /| || (_) | | | |
   |   \__\_\\__,_|\__,_|___/_|_| \_|\___| \_/\_/  \__\___/|_| |_|
  \*/

  template <typename T>
  void
  QN<T>::allocate( integer N ) {
    LAPACK_WRAPPER_ASSERT(
      N > 0 && N <= 1000,
      "QN<T>::allocate, N = " << N << " must be > 0 and <= 1000"
    );
    n = N;
    allocReals.allocate( size_t(n*(n+3)) );
    H = allocReals( size_t(n*n) );
    s = allocReals( size_t(n) );
    y = allocReals( size_t(n) );
    z = allocReals( size_t(n) );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  QN<T>::print( ostream_type & stream ) const {
    for ( integer i = 0; i < this->n; ++i ) {
      for ( integer j = 0; j < this->n; ++j )
        stream
          << std::setw(14)
          << this->H[i+j*this->n]
          << ' ';
      stream << '\n';
    }
  }

  /*\
   |  ____  _____ ____ ____
   | | __ )|  ___/ ___/ ___|
   | |  _ \| |_ | |  _\___ \
   | | |_) |  _|| |_| |___) |
   | |____/|_|   \____|____/
  \*/

  template <typename T>
  void
  BFGS<T>::update(
    valueType const y[],
    valueType const s[]
  ) {
    valueType sy = dot( n, s, 1, y, 1 );
    if ( sy > 0 ) {
      mult( y, z );
      valueType yHy = dot( n, z, 1, y, 1 );
      syr( LOWER, n, (1+yHy/sy)/sy, s, 1, H, n );
      syr2( LOWER, n, -1/sy, s, 1, z, 1, H, n );
    }
  }

  /*\
   |   ____  _____ ____
   |  |  _ \|  ___|  _ \
   |  | | | | |_  | |_) |
   |  | |_| |  _| |  __/
   |  |____/|_|   |_|
  \*/

  template <typename T>
  void
  DFP<T>::update(
    valueType const y[],
    valueType const s[]
  ) {
    valueType sy = dot( n, s, 1, y, 1 );
    if ( sy > 0 ) {
      mult( y, z );
      valueType yHy = dot( n, z, 1, y, 1 );
      syr( LOWER, n, -1/yHy, z, 1, H, n );
      syr( LOWER, n,   1/sy, s, 1, H, n );
    }
  }

  /*\
  :|:   _____ _                            _
  :|:  | ____(_) __ _  ___ _ ____   ____ _| |_   _  ___  ___
  :|:  |  _| | |/ _` |/ _ \ '_ \ \ / / _` | | | | |/ _ \/ __|
  :|:  | |___| | (_| |  __/ | | \ V / (_| | | |_| |  __/\__ \
  :|:  |_____|_|\__, |\___|_| |_|\_/ \__,_|_|\__,_|\___||___/
  :|:           |___/
  \*/

  template <typename T>
  Eigenvalues<T>::Eigenvalues()
  : mem_real("Eigenvalues::mem_real")
  , N(0)
  , Re(nullptr)
  , Im(nullptr)
  , Work(nullptr)
  , A_saved(nullptr)
  {}

  template <typename T>
  void
  Eigenvalues<T>::allocate( integer Nin ) {
    this->N = Nin;
    // calcolo memoria ottimale
    valueType Lworkdummy;
    integer info = geev(
      false, false, Nin, nullptr, Nin,
      nullptr, nullptr, nullptr, Nin, nullptr, Nin, &Lworkdummy, -1
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "Eigenvalues<T>::allocate, call geev return info = " << info
    );
    this->Lwork = integer(Lworkdummy);
    this->mem_real.allocate( size_t( this->Lwork + (2+this->N) * this->N) ) ;
    this->Re      = this->mem_real( size_t(this->N) );
    this->Im      = this->mem_real( size_t(this->N) );
    this->Work    = this->mem_real( size_t(this->Lwork) );
    this->A_saved = this->mem_real( size_t(this->N*this->N) );
  }

  template <typename T>
  void
  Eigenvalues<T>::compute( ) {
    integer info = geev(
      false, false, this->N, this->A_saved, this->N,
      this->Re, this->Im, nullptr, this->N, nullptr, this->N,
      this->Work, this->Lwork
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "Eigenvalues<T>::compute, call geev return info = " << info
    );
  }

  template <typename T>
  Eigenvalues<T>::Eigenvalues(
    integer         NRC,
    valueType const data[],
    integer         ldData
  )
  : mem_real("Eigenvalues::mem_real")
  , N(0)
  , Re(nullptr)
  , Im(nullptr)
  , Work(nullptr)
  , A_saved(nullptr)
  {
    this->setup( NRC, data, ldData );
  }

  template <typename T>
  Eigenvalues<T>::Eigenvalues( MatW const & M )
  : mem_real("Eigenvalues::mem_real")
  , N(0)
  , Re(nullptr)
  , Im(nullptr)
  , Work(nullptr)
  , A_saved(nullptr)
  {
    this->setup( M );
  }

  template <typename T>
  Eigenvalues<T>::Eigenvalues(
    integer         NRC,
    integer         nnz,
    valueType const values[],
    integer   const row[],
    integer   const col[]
  )
  : mem_real("Eigenvalues::mem_real")
  , N(0)
  , Re(nullptr)
  , Im(nullptr)
  , Work(nullptr)
  , A_saved(nullptr)
  {
    this->setup( NRC, nnz, values, row, col );
  }

  template <typename T>
  void
  Eigenvalues<T>::setup(
    integer         NRC,
    valueType const data[],
    integer         ldData
  ) {
    this->allocate( NRC );
    integer info = gecopy( NRC, NRC, data, ldData, A_saved, NRC );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "Eigenvalues<T>::setup, call gecopy return info = " << info
    );
    this->compute();
  }

  template <typename T>
  void
  Eigenvalues<T>::setup( MatW const & M ) {
    this->allocate( M.numRows() );
    integer info = gecopy( this->N, this->N, M.get_data(), M.lDim(),  this->A_saved, this->N );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "Eigenvalues<T>::setup, call gecopy return info = " << info
    );
    this->compute();
  }

  template <typename T>
  void
  Eigenvalues<T>::setup(
    integer         NRC,
    integer         nnz,
    valueType const values[],
    integer   const row[],
    integer   const col[]
  ) {
    this->allocate( NRC );
    std::fill(  this->A_saved, this->A_saved + NRC*NRC, 0 );
    for ( integer i = 0; i < nnz; ++i )  this->A_saved[row[i]+col[i]*NRC] += values[i];
    this->compute();
  }

  template <typename T>
  void
  Eigenvalues<T>::getEigenvalue(
    integer n, valueType & re, valueType & im
  ) const {
    re = Re[n];
    im = Im[n];
  }

  template <typename T>
  void
  Eigenvalues<T>::getEigenvalue(
    integer n, std::complex<valueType> & eig
  ) const {
    eig = std::complex<valueType>(Re[n],Im[n]);
  }

  template <typename T>
  void
  Eigenvalues<T>::getEigenvalues(
    std::vector<valueType> & re,
    std::vector<valueType> & im
  ) const {
    re.clear(); re.reserve( this->N );
    im.clear(); im.reserve( this->N );
    for ( int i = 0;i < this->N; ++i ) {
      re.push_back( Re[i] );
      im.push_back( Im[i] );
    }
  }

  template <typename T>
  void
  Eigenvalues<T>::getEigenvalues(
    std::vector<std::complex<valueType> > & eigs
  ) const {
    eigs.clear(); eigs.reserve( this->N );
    for ( int i = 0;i < this->N; ++i )
      eigs.push_back( std::complex<valueType>(Re[i],Im[i]) );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  Eigenvectors<T>::Eigenvectors()
  : mem_real("Eigenvectors::mem_real")
  , N(0)
  , Re(nullptr)
  , Im(nullptr)
  , A_saved(nullptr)
  , VL(nullptr)
  , VR(nullptr)
  , Work(nullptr)
  {}

  template <typename T>
  void
  Eigenvectors<T>::allocate( integer Nin ) {
    this->N = Nin;
    // calcolo memoria ottimale
    integer doLwork    = -1;
    T       Lworkdummy = 1;
    integer info = geev(
      true, true,
      this->N,
      nullptr, this->N,
      nullptr, nullptr,
      this->VL, this->N,
      this->VR, this->N,
      &Lworkdummy, doLwork
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "Eigenvectors::allocate, call geev return info = " << info
    );
    this->Lwork = integer(Lworkdummy);
    this->mem_real.allocate( size_t( this->Lwork + (2+3*this->N) * this->N) ) ;
    this->Re      = this->mem_real( size_t(this->N) );
    this->Im      = this->mem_real( size_t(this->N) );
    this->A_saved = this->mem_real( size_t(this->N*this->N) );
    this->VL      = this->mem_real( size_t(this->N*this->N) );
    this->VR      = this->mem_real( size_t(this->N*this->N) );
    this->Work    = this->mem_real( size_t(this->Lwork) );
  }

  template <typename T>
  void
  Eigenvectors<T>::compute( ) {
    integer info = geev(
      this->VL != nullptr,
      this->VR != nullptr,
      this->N,
      this->A_saved, this->N,
      this->Re,      this->Im,
      this->VL,      this->N,
      this->VR,      this->N,
      this->Work,    this->Lwork
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "GeneralizedEigenvectors::compute, call ggevx return info = " << info
    );
  }

  template <typename T>
  Eigenvectors<T>::Eigenvectors(
    integer         NRC,
    valueType const A_data[],
    integer         ldA
  )
  : mem_real("Eigenvectors::mem_real")
  , N(0)
  , Re(nullptr)
  , Im(nullptr)
  , A_saved(nullptr)
  , VL(nullptr)
  , VR(nullptr)
  , Work(nullptr)
  {
    this->setup( NRC, A_data, ldA );
  }

  template <typename T>
  Eigenvectors<T>::Eigenvectors( MatW const & A )
  : mem_real("Eigenvectors::mem_real")
  , N(0)
  , Re(nullptr)
  , Im(nullptr)
  , A_saved(nullptr)
  , VL(nullptr)
  , VR(nullptr)
  , Work(nullptr)
  {
    this->setup( A );
  }

  template <typename T>
  Eigenvectors<T>::Eigenvectors(
    integer         NRC,
    integer         A_nnz,
    valueType const A_values[],
    integer   const A_row[],
    integer   const A_col[]
  )
  : mem_real("Eigenvectors::mem_real")
  , N(0)
  , Re(nullptr)
  , Im(nullptr)
  , A_saved(nullptr)
  , VL(nullptr)
  , VR(nullptr)
  , Work(nullptr)
  {
    this->setup( NRC, A_nnz, A_values, A_row, A_col );
  }

  template <typename T>
  void
  Eigenvectors<T>::setup(
    integer         NRC,
    valueType const A_data[],
    integer         ldA
  ) {
    this->allocate( NRC );
    integer info = gecopy( NRC, NRC, A_data, ldA, A_saved, NRC );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "Eigenvectors::setup, call gecopy return info = " << info
    );
    this->compute();
  }

  template <typename T>
  void
  Eigenvectors<T>::setup( MatW const & A ) {
    this->allocate( A.numRows() );
    integer info = gecopy( this->N, this->N, A.get_data(), A.lDim(),  this->A_saved, this->N );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "Eigenvectors::setup, call gecopy return info = " << info
    );
    this->compute();
  }

  template <typename T>
  void
  Eigenvectors<T>::setup(
    integer         NRC,
    integer         A_nnz,
    valueType const A_values[],
    integer   const A_row[],
    integer   const A_col[]
  ) {
    this->allocate( NRC );
    std::fill(  this->A_saved, this->A_saved + NRC*NRC, 0 );
    for ( integer i = 0; i < A_nnz; ++i )
      this->A_saved[A_row[i]+A_col[i]*NRC] += A_values[i];
    this->compute();
  }

  template <typename T>
  void
  Eigenvectors<T>::getEigenvalue(
    integer n, valueType & re, valueType & im
  ) const {
    re = Re[n];
    im = Im[n];
  }

  template <typename T>
  void
  Eigenvectors<T>::getEigenvalue(
    integer n, std::complex<valueType> & eig
  ) const {
    eig = std::complex<valueType>(Re[n],Im[n]);
  }

  template <typename T>
  void
  Eigenvectors<T>::getEigenvalues(
    std::vector<valueType> & re,
    std::vector<valueType> & im
  ) const {
    re.clear(); re.reserve( this->N );
    im.clear(); im.reserve( this->N );
    for ( int i = 0;i < this->N; ++i ) {
      re.push_back( Re[i] );
      im.push_back( Im[i] );
    }
  }

  template <typename T>
  void
  Eigenvectors<T>::getEigenvalues(
    std::vector<std::complex<valueType> > & eigs
  ) const {
    eigs.clear(); eigs.reserve( this->N );
    for ( int i = 0;i < this->N; ++i )
      eigs.push_back( std::complex<valueType>(Re[i],Im[i]) );
  }

  template <typename T>
  void
  Eigenvectors<T>::getLeftEigenvector(
    std::vector<std::vector<complexType> > & vecs
  ) const {
    vecs.resize( size_t(this->N) );
    for ( integer n = 0; n < this->N; ++n ) {
      std::vector<complexType> & v = vecs[n];
      v.clear(); v.reserve( this->N );
      T const * vr = this->VL + n * this->N;
      if ( Im[n] > 0 ) {
        std::vector<complexType> & v1 = vecs[++n];
        v1.clear(); v1.reserve( this->N );
        T const * vi = vr + this->N;
        for ( integer j = 0; j < this->N; ++j ) {
          // salvo vettore "gia" coniugato
          v.push_back( complexType( vr[j], -vi[j] ) );
          v1.push_back( complexType( vr[j], vi[j] ) );
        }
      } else {
        for ( integer j = 0; j < this->N; ++j )
          v.push_back( complexType( vr[j], 0 ) );
      }
    }
  }

  template <typename T>
  void
  Eigenvectors<T>::getRightEigenvector(
    std::vector<std::vector<complexType> > & vecs
  ) const {
    vecs.resize( size_t(this->N) );
    for ( integer n = 0; n < this->N; ++n ) {
      std::vector<complexType> & v = vecs[n];
      v.clear(); v.reserve( this->N );
      T const * vr = this->VR + n * this->N;
      if ( Im[n] > 0 ) {
        std::vector<complexType> & v1 = vecs[++n];
        v1.clear(); v1.reserve( this->N );
        T const * vi = vr + this->N;
        for ( integer j = 0; j < this->N; ++j ) {
          v.push_back( complexType( vr[j], vi[j] ) );
          v1.push_back( complexType( vr[j], -vi[j] ) );
        }
      } else {
        for ( integer j = 0; j < this->N; ++j )
          v.push_back( complexType( vr[j], 0 ) );
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  GeneralizedEigenvalues<T>::GeneralizedEigenvalues()
  : mem_real("GeneralizedEigenvalues::mem_real")
  , N(0)
  , alphaRe(nullptr)
  , alphaIm(nullptr)
  , beta(nullptr)
  , Work(nullptr)
  , A_saved(nullptr)
  , B_saved(nullptr)
  {}

  template <typename T>
  void
  GeneralizedEigenvalues<T>::allocate( integer Nin ) {
    this->N = Nin;
    // calcolo memoria ottimale
    valueType Lworkdummy;
    integer info = ggev(
      false, false, Nin, nullptr, Nin, nullptr, Nin,
      nullptr, nullptr, nullptr,
      nullptr, Nin, nullptr, Nin, &Lworkdummy, -1
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "GeneralizedEigenvalues::allocate, call geev return info = " << info
    );
    this->Lwork = integer(Lworkdummy);
    this->mem_real.allocate( size_t( this->Lwork + (3+2*this->N) * this->N) ) ;
    this->alphaRe = this->mem_real( size_t(this->N) );
    this->alphaIm = this->mem_real( size_t(this->N) );
    this->beta    = this->mem_real( size_t(this->N) );
    this->Work    = this->mem_real( size_t(this->Lwork) );
    this->A_saved = this->mem_real( size_t(this->N*this->N) );
    this->B_saved = this->mem_real( size_t(this->N*this->N) );
  }

  template <typename T>
  void
  GeneralizedEigenvalues<T>::compute( ) {
    integer info = ggev(
      false, false,
      this->N,
      this->A_saved, this->N,
      this->B_saved, this->N,
      this->alphaRe, this->alphaIm, this->beta,
      nullptr, this->N, nullptr, this->N,
      this->Work, this->Lwork
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "GeneralizedEigenvalues::compute, call geev return info = " << info
    );
  }

  template <typename T>
  GeneralizedEigenvalues<T>::GeneralizedEigenvalues(
    integer         NRC,
    valueType const A_data[],
    integer         ldA,
    valueType const B_data[],
    integer         ldB
  )
  : mem_real("GeneralizedEigenvalues::mem_real")
  , N(0)
  , alphaRe(nullptr)
  , alphaIm(nullptr)
  , beta(nullptr)
  , Work(nullptr)
  , A_saved(nullptr)
  , B_saved(nullptr)
  {
    this->setup( NRC, A_data, ldA, B_data, ldB );
  }

  template <typename T>
  GeneralizedEigenvalues<T>::GeneralizedEigenvalues(
    MatW const & A, MatW const & B
  )
  : mem_real("GeneralizedEigenvalues::mem_real")
  , N(0)
  , alphaRe(nullptr)
  , alphaIm(nullptr)
  , beta(nullptr)
  , Work(nullptr)
  , A_saved(nullptr)
  , B_saved(nullptr)
  {
    this->setup( A, B );
  }

  template <typename T>
  GeneralizedEigenvalues<T>::GeneralizedEigenvalues(
    integer         NRC,
    integer         A_nnz,
    valueType const A_values[],
    integer   const A_row[],
    integer   const A_col[],
    integer         B_nnz,
    valueType const B_values[],
    integer   const B_row[],
    integer   const B_col[]
  )
  : mem_real("GeneralizedEigenvalues::mem_real")
  , N(0)
  , alphaRe(nullptr)
  , alphaIm(nullptr)
  , beta(nullptr)
  , Work(nullptr)
  , A_saved(nullptr)
  , B_saved(nullptr)
  {
    this->setup(
      NRC,
      A_nnz, A_values, A_row, A_col,
      B_nnz, B_values, B_row, B_col
    );
  }

  template <typename T>
  void
  GeneralizedEigenvalues<T>::setup(
    integer         NRC,
    valueType const A_data[],
    integer         ldA,
    valueType const B_data[],
    integer         ldB
  ) {
    this->allocate( NRC );
    integer info1 = gecopy( NRC, NRC, A_data, ldA, A_saved, NRC );
    integer info2 = gecopy( NRC, NRC, B_data, ldB, B_saved, NRC );
    LAPACK_WRAPPER_ASSERT(
      info1 == 0 && info2 == 0,
      "GeneralizedEigenvalues::setup, call gecopy return info1 = " << info1 <<
      ", info2 = " << info2
    );
    this->compute();
  }

  template <typename T>
  void
  GeneralizedEigenvalues<T>::setup( MatW const & A, MatW const & B ) {
    this->allocate( A.numRows() );
    integer info1 = gecopy( this->N, this->N, A.get_data(), A.lDim(),  this->A_saved, this->N );
    integer info2 = gecopy( this->N, this->N, B.get_data(), B.lDim(),  this->B_saved, this->N );
    LAPACK_WRAPPER_ASSERT(
      info1 == 0 && info2 == 0,
      "GeneralizedEigenvalues::setup, call gecopy return info1 = " << info1 <<
      ", info2 = " << info2
    );
    this->compute();
  }

  template <typename T>
  void
  GeneralizedEigenvalues<T>::setup(
    integer         NRC,
    integer         A_nnz,
    valueType const A_values[],
    integer   const A_row[],
    integer   const A_col[],
    integer         B_nnz,
    valueType const B_values[],
    integer   const B_row[],
    integer   const B_col[]
  ) {
    this->allocate( NRC );
    std::fill(  this->A_saved, this->A_saved + NRC*NRC, 0 );
    std::fill(  this->B_saved, this->B_saved + NRC*NRC, 0 );
    for ( integer i = 0; i < A_nnz; ++i )
      this->A_saved[A_row[i]+A_col[i]*NRC] += A_values[i];
    for ( integer i = 0; i < B_nnz; ++i )
      this->B_saved[B_row[i]+B_col[i]*NRC] += B_values[i];
    this->compute();
  }

  template <typename T>
  void
  GeneralizedEigenvalues<T>::getEigenvalue(
    integer n, valueType & re, valueType & im
  ) const {
    re = alphaRe[n]/beta[n];
    im = alphaIm[n]/beta[n];
  }

  template <typename T>
  void
  GeneralizedEigenvalues<T>::getEigenvalue(
    integer n, std::complex<valueType> & eig
  ) const {
    eig = std::complex<valueType>(alphaRe[n],alphaIm[n])/beta[n];
  }

  template <typename T>
  void
  GeneralizedEigenvalues<T>::getEigenvalues(
    std::vector<valueType> & re,
    std::vector<valueType> & im
  ) const {
    re.clear(); re.reserve( this->N );
    im.clear(); im.reserve( this->N );
    for ( int i = 0;i < this->N; ++i ) {
      re.push_back( alphaRe[i]/beta[i] );
      im.push_back( alphaIm[i]/beta[i] );
    }
  }

  template <typename T>
  void
  GeneralizedEigenvalues<T>::getEigenvalues(
    std::vector<std::complex<valueType> > & eigs
  ) const {
    eigs.clear(); eigs.reserve( this->N );
    for ( int i = 0;i < this->N; ++i )
      eigs.push_back( std::complex<valueType>(alphaRe[i],alphaIm[i])/beta[i] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  GeneralizedEigenvectors<T>::GeneralizedEigenvectors()
  : mem_real("GeneralizedEigenvectors::mem_real")
  , mem_int("GeneralizedEigenvectors::mem_int")
  , N(0)
  , alphaRe(nullptr)
  , alphaIm(nullptr)
  , beta(nullptr)
  , A_saved(nullptr)
  , B_saved(nullptr)
  , VL(nullptr)
  , VR(nullptr)
  , lscale(nullptr)
  , rscale(nullptr)
  , rconde(nullptr)
  , rcondv(nullptr)
  , Work(nullptr)
  , iWork(nullptr)
  , bWork(nullptr)
  {}

  template <typename T>
  void
  GeneralizedEigenvectors<T>::allocate( integer Nin ) {
    this->N = Nin;
    // calcolo memoria ottimale
    integer doLwork    = -1;
    T       Lworkdummy = 1;
    integer info = ggevx(
      lapack_wrapper::PERMUTE_AND_SCALE,
      false, false,
      lapack_wrapper::EIGENVALUES_AND_EIGENVECTORS,
      this->N, nullptr, this->N, nullptr, this->N,
      nullptr, nullptr, nullptr, this->VL, this->N, this->VR, this->N,
      ilo, ihi,
      nullptr, nullptr,
      abnorm,  bbnorm,
      nullptr, nullptr,
      &Lworkdummy, doLwork,
      nullptr, nullptr
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "GeneralizedEigenvectors::allocate, call geev return info = " << info
    );
    this->Lwork = integer(Lworkdummy);
    this->mem_real.allocate( size_t( this->Lwork + (7+4*this->N) * this->N) ) ;
    this->mem_int.allocate( size_t( 2*this->N + 6 ) ) ;
    this->alphaRe = this->mem_real( size_t(this->N) );
    this->alphaIm = this->mem_real( size_t(this->N) );
    this->beta    = this->mem_real( size_t(this->N) );
    this->A_saved = this->mem_real( size_t(this->N*this->N) );
    this->B_saved = this->mem_real( size_t(this->N*this->N) );
    this->VL      = this->mem_real( size_t(this->N*this->N) );
    this->VR      = this->mem_real( size_t(this->N*this->N) );
    this->lscale  = this->mem_real( size_t(this->N) );
    this->rscale  = this->mem_real( size_t(this->N) );
    this->rconde  = this->mem_real( size_t(this->N) );
    this->rcondv  = this->mem_real( size_t(this->N) );
    this->Work    = this->mem_real( size_t(this->Lwork) );
    this->iWork   = this->mem_int( size_t(this->N+6) );
    this->bWork   = this->mem_int( size_t(this->N) );
  }

  template <typename T>
  void
  GeneralizedEigenvectors<T>::compute( ) {
    integer info = ggevx(
      lapack_wrapper::PERMUTE_ONLY, // "B", // ATTENZIONE LA SCALATURA NON FUNZIONA
      this->VL != nullptr,
      this->VR != nullptr,
      lapack_wrapper::EIGENVALUES_AND_EIGENVECTORS,
      this->N,
      this->A_saved, this->N,
      this->B_saved, this->N,
      this->alphaRe, this->alphaIm, this->beta,
      this->VL, this->N,
      this->VR, this->N,
      ilo, ihi,
      this->lscale, this->rscale,
      abnorm,       bbnorm,
      this->rconde, this->rcondv,
      this->Work,   this->Lwork,
      this->iWork,  this->bWork
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "GeneralizedEigenvectors::compute, call ggevx return info = " << info
    );
  }

  template <typename T>
  GeneralizedEigenvectors<T>::GeneralizedEigenvectors(
    integer         NRC,
    valueType const A_data[],
    integer         ldA,
    valueType const B_data[],
    integer         ldB
  )
  : mem_real("GeneralizedEigenvectors::mem_real")
  , mem_int("GeneralizedEigenvectors::mem_int")
  , N(0)
  , alphaRe(nullptr)
  , alphaIm(nullptr)
  , beta(nullptr)
  , A_saved(nullptr)
  , B_saved(nullptr)
  , VL(nullptr)
  , VR(nullptr)
  , lscale(nullptr)
  , rscale(nullptr)
  , rconde(nullptr)
  , rcondv(nullptr)
  , Work(nullptr)
  , iWork(nullptr)
  , bWork(nullptr)
  {
    this->setup( NRC, A_data, ldA, B_data, ldB );
  }

  template <typename T>
  GeneralizedEigenvectors<T>::GeneralizedEigenvectors(
    MatW const & A, MatW const & B
  )
  : mem_real("GeneralizedEigenvectors::mem_real")
  , mem_int("GeneralizedEigenvectors::mem_int")
  , N(0)
  , alphaRe(nullptr)
  , alphaIm(nullptr)
  , beta(nullptr)
  , A_saved(nullptr)
  , B_saved(nullptr)
  , VL(nullptr)
  , VR(nullptr)
  , lscale(nullptr)
  , rscale(nullptr)
  , rconde(nullptr)
  , rcondv(nullptr)
  , Work(nullptr)
  , iWork(nullptr)
  , bWork(nullptr)
  {
    this->setup( A, B );
  }

  template <typename T>
  GeneralizedEigenvectors<T>::GeneralizedEigenvectors(
    integer         NRC,
    integer         A_nnz,
    valueType const A_values[],
    integer   const A_row[],
    integer   const A_col[],
    integer         B_nnz,
    valueType const B_values[],
    integer   const B_row[],
    integer   const B_col[]
  )
  : mem_real("GeneralizedEigenvectors::mem_real")
  , mem_int("GeneralizedEigenvectors::mem_int")
  , N(0)
  , alphaRe(nullptr)
  , alphaIm(nullptr)
  , beta(nullptr)
  , A_saved(nullptr)
  , B_saved(nullptr)
  , VL(nullptr)
  , VR(nullptr)
  , lscale(nullptr)
  , rscale(nullptr)
  , rconde(nullptr)
  , rcondv(nullptr)
  , Work(nullptr)
  , iWork(nullptr)
  , bWork(nullptr)
  {
    this->setup(
      NRC,
      A_nnz, A_values, A_row, A_col,
      B_nnz, B_values, B_row, B_col
    );
  }

  template <typename T>
  void
  GeneralizedEigenvectors<T>::setup(
    integer         NRC,
    valueType const A_data[],
    integer         ldA,
    valueType const B_data[],
    integer         ldB
  ) {
    this->allocate( NRC );
    integer info1 = gecopy( NRC, NRC, A_data, ldA, A_saved, NRC );
    integer info2 = gecopy( NRC, NRC, B_data, ldB, B_saved, NRC );
    LAPACK_WRAPPER_ASSERT(
      info1 == 0 && info2 == 0,
      "GeneralizedEigenvectors::setup, call gecopy return info1 = " << info1 <<
      ", info2 = " << info2
    );
    this->compute();
  }

  template <typename T>
  void
  GeneralizedEigenvectors<T>::setup( MatW const & A, MatW const & B ) {
    this->allocate( A.numRows() );
    integer info1 = gecopy( this->N, this->N, A.get_data(), A.lDim(),  this->A_saved, this->N );
    integer info2 = gecopy( this->N, this->N, B.get_data(), B.lDim(),  this->B_saved, this->N );
    LAPACK_WRAPPER_ASSERT(
      info1 == 0 && info2 == 0,
      "GeneralizedEigenvectors::setup, call gecopy return info1 = " << info1 <<
      ", info2 = " << info2
    );
    this->compute();
  }

  template <typename T>
  void
  GeneralizedEigenvectors<T>::setup(
    integer         NRC,
    integer         A_nnz,
    valueType const A_values[],
    integer   const A_row[],
    integer   const A_col[],
    integer         B_nnz,
    valueType const B_values[],
    integer   const B_row[],
    integer   const B_col[]
  ) {
    this->allocate( NRC );
    std::fill(  this->A_saved, this->A_saved + NRC*NRC, 0 );
    std::fill(  this->B_saved, this->B_saved + NRC*NRC, 0 );
    for ( integer i = 0; i < A_nnz; ++i )
      this->A_saved[A_row[i]+A_col[i]*NRC] += A_values[i];
    for ( integer i = 0; i < B_nnz; ++i )
      this->B_saved[B_row[i]+B_col[i]*NRC] += B_values[i];
    this->compute();
  }

  template <typename T>
  void
  GeneralizedEigenvectors<T>::getEigenvalue(
    integer n, valueType & re, valueType & im
  ) const {
    re = alphaRe[n]/beta[n];
    im = alphaIm[n]/beta[n];
  }

  template <typename T>
  void
  GeneralizedEigenvectors<T>::getEigenvalue(
    integer n, std::complex<valueType> & eig
  ) const {
    eig = std::complex<valueType>(alphaRe[n],alphaIm[n])/beta[n];
  }

  template <typename T>
  void
  GeneralizedEigenvectors<T>::getEigenvalues(
    std::vector<valueType> & re,
    std::vector<valueType> & im
  ) const {
    re.clear(); re.reserve( this->N );
    im.clear(); im.reserve( this->N );
    for ( int i = 0;i < this->N; ++i ) {
      re.push_back( alphaRe[i]/beta[i] );
      im.push_back( alphaIm[i]/beta[i] );
    }
  }

  template <typename T>
  void
  GeneralizedEigenvectors<T>::getEigenvalues(
    std::vector<std::complex<valueType> > & eigs
  ) const {
    eigs.clear(); eigs.reserve( this->N );
    for ( int i = 0;i < this->N; ++i )
      eigs.push_back( std::complex<valueType>(alphaRe[i],alphaIm[i])/beta[i] );
  }

  template <typename T>
  void
  GeneralizedEigenvectors<T>::getLeftEigenvector(
    std::vector<std::vector<complexType> > & vecs
  ) const {
    vecs.resize( size_t(this->N) );
    for ( integer n = 0; n < this->N; ++n ) {
      std::vector<complexType> & v = vecs[n];
      v.clear(); v.reserve( this->N );
      T const * vr = this->VL + n * this->N;
      if ( alphaIm[n] > 0 ) {
        std::vector<complexType> & v1 = vecs[++n];
        v1.clear(); v1.reserve( this->N );
        T const * vi = vr + this->N;
        for ( integer j = 0; j < this->N; ++j ) {
          // salvo vettore "gia" coniugato
          v.push_back( complexType( vr[j], -vi[j] ) );
          v1.push_back( complexType( vr[j], vi[j] ) );
        }
      } else {
        for ( integer j = 0; j < this->N; ++j )
          v.push_back( complexType( vr[j], 0 ) );
      }
    }
  }

  template <typename T>
  void
  GeneralizedEigenvectors<T>::getRightEigenvector(
    std::vector<std::vector<complexType> > & vecs
  ) const {
    vecs.resize( size_t(this->N) );
    for ( integer n = 0; n < this->N; ++n ) {
      std::vector<complexType> & v = vecs[n];
      v.clear(); v.reserve( this->N );
      T const * vr = this->VR + n * this->N;
      if ( alphaIm[n] > 0 ) {
        std::vector<complexType> & v1 = vecs[++n];
        v1.clear(); v1.reserve( this->N );
        T const * vi = vr + this->N;
        for ( integer j = 0; j < this->N; ++j ) {
          v.push_back( complexType( vr[j], vi[j] ) );
          v1.push_back( complexType( vr[j], -vi[j] ) );
        }
      } else {
        for ( integer j = 0; j < this->N; ++j )
          v.push_back( complexType( vr[j], 0 ) );
      }
    }
  }

  /*\
  :|:
  :|:    ____                           _ _             _ ______     ______
  :|:   / ___| ___ _ __   ___ _ __ __ _| (_)_______  __| / ___\ \   / /  _ \
  :|:  | |  _ / _ \ '_ \ / _ \ '__/ _` | | |_  / _ \/ _` \___ \\ \ / /| | | |
  :|:  | |_| |  __/ | | |  __/ | | (_| | | |/ /  __/ (_| |___) |\ V / | |_| |
  :|:   \____|\___|_| |_|\___|_|  \__,_|_|_/___\___|\__,_|____/  \_/  |____/
  :|:
  \*/

  template <typename T>
  GeneralizedSVD<T>::GeneralizedSVD()
  : mem_real("GeneralizedSVD(real)")
  , mem_int("GeneralizedSVD(int)")
  , M(0)
  , N(0)
  , P(0)
  , K(0)
  , L(0)
  , Lwork(0)
  , Work(nullptr)
  , IWork(nullptr)
  , alpha_saved(nullptr)
  , beta_saved(nullptr)
  , A_saved(nullptr)
  , B_saved(nullptr)
  , U_saved(nullptr)
  , V_saved(nullptr)
  , Q_saved(nullptr)
  {}

  template <typename T>
  GeneralizedSVD<T>::GeneralizedSVD(
    integer         m,
    integer         n,
    integer         p,
    valueType const A[], integer ldA_in,
    valueType const B[], integer ldB_in
  )
  : mem_real("GeneralizedSVD(real)")
  , mem_int("GeneralizedSVD(int)")
  , M(0)
  , N(0)
  , P(0)
  , K(0)
  , L(0)
  , Lwork(0)
  , Work(nullptr)
  , IWork(nullptr)
  , alpha_saved(nullptr)
  , beta_saved(nullptr)
  , A_saved(nullptr)
  , B_saved(nullptr)
  , U_saved(nullptr)
  , V_saved(nullptr)
  , Q_saved(nullptr)
  { this->setup( m, n, p, A, ldA_in, B, ldB_in ); }

  template <typename T>
  GeneralizedSVD<T>::GeneralizedSVD( MatW const & A, MatW const & B )
  : mem_real("GeneralizedSVD(real)")
  , mem_int("GeneralizedSVD(int)")
  , M(0)
  , N(0)
  , P(0)
  , K(0)
  , L(0)
  , Lwork(0)
  , Work(nullptr)
  , IWork(nullptr)
  , alpha_saved(nullptr)
  , beta_saved(nullptr)
  , A_saved(nullptr)
  , B_saved(nullptr)
  , U_saved(nullptr)
  , V_saved(nullptr)
  , Q_saved(nullptr)
  { this->setup( A, B ); }

  template <typename T>
  GeneralizedSVD<T>::GeneralizedSVD(
    integer         m,
    integer         n,
    integer         p,
    integer         A_nnz,
    valueType const A_values[],
    integer   const A_row[],
    integer   const A_col[],
    integer         B_nnz,
    valueType const B_values[],
    integer   const B_row[],
    integer   const B_col[]
  )
  : mem_real("GeneralizedSVD(real)")
  , mem_int("GeneralizedSVD(int)")
  , M(0)
  , N(0)
  , P(0)
  , K(0)
  , L(0)
  , Lwork(0)
  , Work(nullptr)
  , IWork(nullptr)
  , alpha_saved(nullptr)
  , beta_saved(nullptr)
  , A_saved(nullptr)
  , B_saved(nullptr)
  , U_saved(nullptr)
  , V_saved(nullptr)
  , Q_saved(nullptr)
  { this->setup(
      m, n, p,
      A_nnz, A_values, A_row, A_col,
      B_nnz, B_values, B_row, B_col
    );
  }

  template <typename T>
  void
  GeneralizedSVD<T>::allocate( integer m, integer n, integer p ) {
    integer k, l;
    real    wL;
    integer info = ggsvd(
      true, true, true, m, n, p, k, l,
      nullptr, m,
      nullptr, p,
      nullptr,
      nullptr,
      nullptr, m,
      nullptr, p,
      nullptr, n,
      &wL,
      -1,
      nullptr
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "GeneralizedSVD<T>::allocate(m=" << m << ",n" << n << ",p=" << p <<
      ") failed, info = " << info
    );
    this->M     = m;
    this->N     = n;
    this->P     = p;
    this->Lwork = integer(wL);

    this->mem_int.allocate( n );
    this->IWork = mem_int( n );

    this->mem_real.allocate( Lwork + (m+p+2)*n + m*m + p*p + n*n );
    this->Work        = mem_real( Lwork );
    this->alpha_saved = mem_real( n );
    this->beta_saved  = mem_real( n );
    this->A_saved     = mem_real( m*n );
    this->B_saved     = mem_real( p*n );
    this->U_saved     = mem_real( m*m );
    this->V_saved     = mem_real( p*p );
    this->Q_saved     = mem_real( n*n );

    this->U.setup( this->U_saved, m, m, m );
    this->V.setup( this->V_saved, p, p, p );
    this->Q.setup( this->Q_saved, n, n, n );
    this->Dbeta.setup( this->beta_saved, n );
    this->Dalpha.setup( this->alpha_saved, n );
  }

  template <typename T>
  void
  GeneralizedSVD<T>::compute() {
    integer info = ggsvd(
      true, true, true,
      this->M, this->N, this->P, this->K, this->L,
      this->A_saved, this->M,
      this->B_saved, this->P,
      this->alpha_saved,
      this->beta_saved,
      this->U_saved, this->M,
      this->V_saved, this->P,
      this->Q_saved, this->N,
      Work,
      Lwork,
      IWork
    );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "GeneralizedSVD<T>::compute() failed, info = " << info
    );
    // R is stored in A(1:K+L,N-K-L+1:N) on exit.
    this->R.setup( this->A_saved + this->M * ( this->N-this->K-this->L ),
                   this->N, this->K+this->L, this->M );
  }

  template <typename T>
  void
  GeneralizedSVD<T>::setup(
    integer         m,
    integer         n,
    integer         p,
    valueType const A[], integer ldA,
    valueType const B[], integer ldB
  ) {
    this->allocate( m, n, p );
    integer info = gecopy( m, n, A, ldA, A_saved, m );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "GeneralizedSVD<T>::setup(...) failed to copy A, info = " << info
    );
    info = gecopy( p, n, B, ldB, B_saved, p );
    LAPACK_WRAPPER_ASSERT(
      info == 0,
      "GeneralizedSVD<T>::setup(...) failed to copy B, info = " << info
    );
    compute();
  }

  template <typename T>
  void
  GeneralizedSVD<T>::setup( MatW const & A, MatW const & B ) {
    integer m = A.numRows();
    integer n = A.numCols();
    integer p = B.numRows();
    LAPACK_WRAPPER_ASSERT(
      n == B.numCols(),
      "GeneralizedSVD<T>::setup( A, B ) incompatible matrices\n" <<
      "A is " << A.numRows() << " x " << A.numCols() << "\n"
      "B is " << B.numRows() << " x " << B.numCols()
    );
    this->setup( m, n, p, A.get_data(), A.lDim(), B.get_data(), B.lDim() );
  }

  template <typename T>
  void
  GeneralizedSVD<T>::setup(
    integer         m,
    integer         n,
    integer         p,
    integer         A_nnz,
    valueType const A_values[],
    integer   const A_row[],
    integer   const A_col[],
    integer         B_nnz,
    valueType const B_values[],
    integer   const B_row[],
    integer   const B_col[]
  ) {
    this->allocate( m, n, p );
    lapack_wrapper::zero( this->N*this->M, this->A_saved, 1 );
    lapack_wrapper::zero( this->P*this->M, this->B_saved, 1 );
    for ( integer k = 0; k < A_nnz; ++k ) A(A_row[k],A_col[k]) = A_values[k];
    for ( integer k = 0; k < B_nnz; ++k ) B(B_row[k],B_col[k]) = B_values[k];
    compute();
  }

  template <typename T>
  void
  GeneralizedSVD<T>::info( ostream_type & stream, valueType eps ) const {
    stream
      << "A = " << this->M << " x " << this->N << '\n'
      << "B = " << this->P << " x " << this->N << '\n';
    for ( integer i = 0; i < this->N; ++i ) {
      T a =  this->alpha_saved[i];
      T b =  this->beta_saved[i];
      stream
        << "alpha[" << i << "]=" << std::setw(14) << a
        << ", beta[" << i << "]=" << std::setw(14) << b
        << ", alpha^2+beta^2 = " << a*a+b*b << '\n';
    }
    stream << "U\n"; this->U.print0( stream, eps );
    stream << "V\n"; this->V.print0( stream, eps );
    stream << "Q\n"; this->Q.print0( stream, eps );
    stream << "R\n"; this->R.print0( stream, eps );
    stream << "C\n"; this->Dalpha.print( stream );
    stream << "S\n"; this->Dbeta.print( stream );
    stream << '\n';
  }

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
